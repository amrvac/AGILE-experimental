# MPI-AMRVAC data file format
> this doc describes the latest version of the data file format, which is numbered "6"

All data files consist of a single snapshot, and they can be used for restart
and/or further conversion to other data formats directly usable for
visualization. Note that restart is possible on a differing number of CPUs,
and may suddenly allow more refinement levels. Also, note that the individual
snapshots will typically have different lengths, as the number of grid blocks
will vary dynamically. The data is saved in binary format (double precision
up to version 5; version 6 uses single precision, see below).
You can find the exact implementation in `src/mod_input_output.t`, more specifically in `snapshot_write_header()` and `snapshot_read_header()`.

The normal output (`<base_filename>NNNN.dat`, file type 2, `dtsave_dat`)
is always written in the double precision version 5 format and is the file
to restart from. In addition, single precision **version 6 snapshots**
(`<base_filename>_snapNNNN.dat`, file type 6) can be written for analysis on
their own schedule (`dtsave_snap`/`ditsave_snap` in the `savelist`
namelist); they are roughly half the size.

A snapshot (`.dat` file) contains a header, mesh tree information, and block
data, in the following order:

# Header

```{fortran}
integer           :: Version number
integer           :: Byte offset where tree information starts
integer           :: Byte offset where block data starts
integer           :: nw
integer           :: ndir
integer           :: ndim
integer           :: levmax
integer           :: nleafs
integer           :: nparents
integer           :: it
double precision  :: global_time
double precision  :: xprobmin(ndim)
double precision  :: xprobmax(ndim)
integer           :: domain_nx(ndim)
integer           :: block_nx(ndim)
logical           :: periodic(ndim)
character(len=16) :: geometry
logical           :: staggered
integer           :: size_real         ! from version 6: bytes per real in block data (4 or 8)
character(len=16) :: w_names(nw)
character(len=16) :: physics_type
! The physics parameters, such as gamma
integer           :: n_params
double precision  :: parameters(n_params)
character(len=16) :: parameter_names(n_params)
! Indexes for file output (for restarting)
integer           :: snapshotnext
integer           :: slicenext
integer           :: collapsenext
```

# Tree information

```{fortran}
logical :: leaf(nleafs+nparents)
integer :: refinement_level(nleafs)
integer :: spatial_index(ndim, nleafs)
integer(kind=MPI_OFFSET_KIND) :: offset_block(nleafs)
```

# Block 1 to nleafs

```{fortran}
integer :: n_ghost_lo(ndim) ! number of ghost cells on lower boundaries
integer :: n_ghost_hi(ndim) ! number of ghost cells on upper boundaries
! block_shape = 1-n_ghost_lo:block_nx+n_ghost_hi
! double precision up to version 5, single precision from version 6
double precision :: w(block_shape, nw)
```

# Version history

## Version 1

Version 1 contained the following information

    1. block data (nw variables)
    2. leaf/parent logical array
    3. header:
        nx^D
        domain_nx^D
        xprobmin^D
        xprobmax^D
        nleafs
        levmax
        ndim
        ndir
        nw
        it
        global_time

The idea is that you can reconstruct the full grid when you know the Morton
order used for the leaf/parent logical array.

## Version 2

Version 2 had the same information as version 1, but changes were made to the
Morton order on the coarse grid, causing incompatibility.

## Version 5

Version 5 introduced the `geometry` parameter (e.g. "polar_2D",
"cartesian_1.75D"... ), as well as a `periodic` and a `staggered` flags.
The periodic flag a `ndir`-long 1D boolean array. For each direction,
it defaults to `.false.` but is set to `.true.` if at least one
quantity is using a periodic boundary condition in the corresponding
direction. This change was motivated by improving the (upcoming)
compatibility with `yt`.  The `staggered` flag is added to support
staggered grid for constrained transport MHD.

## Version 6 (current)

Version 6 has the same header and layout as version 5, but stores the block
data (including any staggered-field data) in single precision, roughly halving
the file size. A single `size_real` integer (bytes per real, = 4) is added to
the header after the `staggered` flag, so that readers decode the block data
from this field rather than mapping the version number to a precision. It is
written by the snapshot output stream (file type 6,
`<base_filename>_snapNNNN.dat`) and is intended for analysis; use the normal
output (always v5, double precision) for restarts.
