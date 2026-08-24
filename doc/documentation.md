# About the AGILE documentation

[TOC]

# Generating the documentation {#doc-gen}

API documentation for AGILE is generated from the source code using
[FORD](https://forddocs.readthedocs.io/) (the Fortran Documenter). Unlike
upstream MPI-AMRVAC, this fork no longer uses Doxygen.

To generate it locally:

    uv sync --extra docs
    source .venv/bin/activate
    doc/ford/build_docs.sh

This first expands `src/*.fpp` with `fypp` (FORD cannot parse fypp directives
directly) using the compile-time flags in `doc/ford/docgen.par`, then runs
FORD on the result. Open the output with

    firefox doc/ford/html/index.html

Because AGILE's physics modules (`hd`, `mhd`, `ffhd`, `srhd`) are
mutually-exclusive fypp branches spliced into the same files
(`mod_physics`, `mod_finite_volume`, `mod_dt`, `mod_source`, ...) at
compile time, one `fypp`/FORD pass can only show one physics module's
variant of those files. `docgen.par` currently selects `mhd` with most of
its optional features enabled, to maximize what's documented; see
`doc/ford/project.md` for details.

# How to write documentation {#doc-howto}

FORD understands the same `!>`/`!<` Doxygen-style comment syntax already
used throughout AGILE, so existing source comments don't need to change.
For the full set of features (special commands like `@note`, `@warning`,
`@todo`, cross-references, LaTeX math, etc.) see the
[FORD user guide](https://forddocs.readthedocs.io/en/latest/user_guide/writing_documentation.html).

## Documenting source code {#doc-src}

You can write documentation comments almost in the same way as regular
comments, using the following syntax:

    ! The number of iterations (normal comment, ignored by FORD)
    integer :: x

    !> The number of iterations (variant 1)
    integer :: bum_its

    integer :: x !< The number of iterations (variant 2)

Note that `!>` describes the next statement, whereas `!<` describes the previous statement.
Multi-line comments can be formed in the following way:

    !> a long line
    !> of text
    !> it is really long

    !> a long line
    !! of text
    !! it is really long

You can document variables, functions, subroutines, modules, types and arguments.
Here are some examples to demonstrate the syntax:

    !> Compute the square of x
    subroutine square(x, x2)
        real, intent(in)  :: x  !< The number we will square
        real, intent(out) :: x2 !< The square of x

        ! This comment will not appear in the documentation
        x2 = x**2
    end subroutine square

    !> A module that contains nothing
    !>
    !> A longer description for the module that does nothing,
    !> although it seems hard to be very elaborate.
    module meaning_of_life
    end module meaning_of_life

    !> A 2D point coordinate
    type coordinate
        real :: x !< The x-coordinate
        real :: y !< The y-coordinate
    end type coordinate

## Documentation in markdown files {#doc-md}

FORD can also fold separate markdown pages (like this one) into the
generated site, via its `page_dir` setting -- see the
[FORD page documentation](https://forddocs.readthedocs.io/en/latest/user_guide/writing_documentation.html#writing-pages).
This is **not currently wired up**: the `doc/*.md` files here still use
Doxygen-specific syntax (`{#label}` section anchors, `@ref label` links,
bare `[TOC]`) that would need to be converted to FORD's page-linking
conventions first. Until that migration happens, these markdown files are
only readable as plain files (e.g. on GitHub), not as part of the generated
FORD site.
