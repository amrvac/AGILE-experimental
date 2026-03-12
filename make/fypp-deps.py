from collections.abc import Iterable
from io import TextIOBase
from pathlib import Path

import re
import argparse


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


if __name__ == "__main__":
    args = parse_arguments()
    deps = []
    input_path = Path(args.input)
    with input_path.open("r") as fid:
        for line in fid:
            if m := INCLUDE_EX.match(line):
                deps.append(input_path.parent / m["filename"])
    if deps:
        target_path = Path(args.target)
        target_path.parent.mkdir(parents=True, exist_ok=True)
        with target_path.with_suffix(".d").open("w") as fid:
            print(args.target + ":", " ".join(map(str, deps)), file=fid)
