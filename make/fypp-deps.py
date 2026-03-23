from collections.abc import Iterable
from io import TextIOBase
from pathlib import Path
from dataclasses import dataclass
from textwrap import indent

import re
import argparse
import sys


INCLUDE_EX = re.compile(r"\s*#:include\s+(?P<quote>['\"])(?P<filename>[^'\"]+)(?P=quote)")


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="fypp-deps",
        description="""Scans a FYPP file (Fortran Preprocessor Language) for
        #:include statements that fortdepend doesn't register. Outputs Makefile
        compatible dependency statements that can be directly included by Make.

        The output (if any) is written to a file with the same basename as the
        target file, but the suffix changed to `.d`.
        """,
        epilog="""
        Conditional Includes: Fypp processes #:include statements in a first pass,
        resulting in a larger "virtual" code, in which other preprocessing statements
        are evaluated. There is no such thing as conditional inclusion.

        Any #:include statement is unconditional, even if the evaluation of the
        resulting expanded code may be conditional, so an #:include always
        leads to a dependency even when the included content goes unused.
        """)
    parser.add_argument("input")
    parser.add_argument("target")

    return parser.parse_args()


@dataclass
class CircularDependencyError(Exception):
    visited: list[Path]
    file: Path

    def __str__(self) -> str:
        cycle = "\n".join(f"includes `{f.relative_to(self.visited[0].parent, walk_up=True)}`" for f in self.visited[1:] + [self.file])
        return f"Circular dependency detected in `{self.file}`: starting from `{self.visited[0]}`:\n" + indent(cycle, "  - ")


def scan_deps(input_path: Path, visited: list[Path] | None = None) -> set[Path]:
    visited = visited or list()
    if input_path in visited:
        raise CircularDependencyError(visited, input_path)

    visited.append(input_path)
    deps = set()

    with input_path.open("r") as fid:
        for line in fid:
            if m := INCLUDE_EX.match(line):
                d = (input_path.parent / m["filename"]).resolve()
                deps.add(d)
                deps.update(scan_deps(d, visited))

    visited.pop()

    return deps


if __name__ == "__main__":
    args = parse_arguments()
    input_path = Path(args.input)

    try:
        deps = scan_deps(input_path)
    except (CircularDependencyError, FileNotFoundError) as e:
        print(f"Error scanning `{input_path}`:\n   ", e)
        sys.exit(-1)

    if deps:
        target_path = Path(args.target)
        target_path.parent.mkdir(parents=True, exist_ok=True)
        with target_path.with_suffix(".d").open("w") as fid:
            print(args.target + ":", " ".join(map(str, deps)), file=fid)
