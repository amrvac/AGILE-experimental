"""
Module containing reading and processing methods for an MPI-AMRVAC .datfiles file.

@author: Jannis Theunissen (original)
         Clément Robert (extensions, modifications)
         Adrian Kelly (v6)
         Olaf Willocx (v6)
Last edit: 1 August 2026
"""

import struct
import numpy as np

# Size of basic types (in bytes)
SIZE_LOGICAL = 4
SIZE_INT = 4
SIZE_DOUBLE = 8
NAME_LEN = 16

# For un-aligned data, use '=' (for aligned data set to '')
ALIGN = '='


def get_header(istream):
    """Read header from an MPI-AMRVAC 2.1 snapshot. This is compatible with versions down to 2.0.
    :param: istream     open datfile buffer with 'rb' mode
    :return: h          header information contained in a dictionary
    """
    istream.seek(0)
    h = {}

    fmt = ALIGN + 'i'
    [h['datfile_version']] = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))

    if h['datfile_version'] < 3:
        raise IOError("Unsupported AMRVAC .datfiles file version: %d", h['datfile_version'])

    # Read scalar data at beginning of file
    fmt = ALIGN + 9 * 'i' + 'd'
    hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    [h['offset_tree'], h['offset_blocks'], h['nw'],
     h['ndir'], h['ndim'], h['levmax'], h['nleafs'], h['nparents'],
     h['it'], h['time']] = hdr

    # Read min/max coordinates
    fmt = ALIGN + h['ndim'] * 'd'
    h['xmin'] = np.array(
        struct.unpack(fmt, istream.read(struct.calcsize(fmt))))
    h['xmax'] = np.array(
        struct.unpack(fmt, istream.read(struct.calcsize(fmt))))

    # Read domain and block size (in number of cells)
    fmt = ALIGN + h['ndim'] * 'i'
    h['domain_nx'] = np.array(
        struct.unpack(fmt, istream.read(struct.calcsize(fmt))))
    h['block_nx'] = np.array(
        struct.unpack(fmt, istream.read(struct.calcsize(fmt))))

    if h['datfile_version'] < 5:
        h['size_real'] = SIZE_DOUBLE

    if h['datfile_version'] >= 5:
        # Read periodicity
        fmt = ALIGN + h['ndim'] * 'i' # Fortran logical is 4 byte int
        h['periodic'] = np.array(
            struct.unpack(fmt, istream.read(struct.calcsize(fmt))), dtype=bool)

        # Read geometry name
        fmt = ALIGN + NAME_LEN * 'c'
        hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
        h['geometry'] = b''.join(hdr).strip().decode()

        # Read staggered flag
        fmt = ALIGN + 'i' # Fortran logical is 4 byte int
        h['staggered'] = bool(
            struct.unpack(fmt, istream.read(struct.calcsize(fmt)))[0])

        if h['datfile_version'] >= 6:
            # Bytes per real in the block data (4 for single precision)
            fmt = ALIGN + 'i'
            [h['size_real']] = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))

    # Read w_names
    w_names = []
    for i in range(h['nw']):
        fmt = ALIGN + NAME_LEN * 'c'
        hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
        w_names.append(b''.join(hdr).strip().decode())
    h['w_names'] = w_names

    # Read physics type
    fmt = ALIGN + NAME_LEN * 'c'
    hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    h['physics_type'] = b''.join(hdr).strip().decode()

    # Read number of physics-defined parameters
    fmt = ALIGN + 'i'
    [n_pars] = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))

    # First physics-parameter values are given, then their names
    fmt = ALIGN + n_pars * 'd'
    vals = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))

    fmt = ALIGN + n_pars * NAME_LEN * 'c'
    names = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    # Split and join the name strings (from one character array)
    names = [b''.join(names[i:i+NAME_LEN]).strip().decode()
             for i in range(0, len(names), NAME_LEN)]

    # Store the values corresponding to the names
    for val, name in zip(vals, names):
        h[name] = val
    return h


def get_tree_info(istream):
    """
    Read levels, morton-curve indices, and byte offsets for each block as stored in the datfile
    :param: istream         open datfile buffer with 'rb' mode
    :return: block_lvls     numpy array with block levels
             block_ixs      numpy array with block morton-curve indices
             block_offsets  numpy array with block (data) offset in the datfile
             block_nghost   numpy array (nleafs, 2*ndim) of per-block ghost-cell
                            counts: [n_ghost_lo(ndim), n_ghost_hi(ndim)]. Zero for
                            files written without a ghost halo.
    """
    istream.seek(0)
    header = get_header(istream)
    nleafs = header['nleafs']
    nparents = header['nparents']
    ndim = header['ndim']

    # Read tree info. Skip 'leaf' array
    istream.seek(header['offset_tree'] + (nleafs+nparents) * SIZE_LOGICAL)

    # Read block levels
    fmt = ALIGN + nleafs * 'i'
    block_lvls = np.array(struct.unpack(fmt, istream.read(struct.calcsize(fmt))))

    # Read block indices
    fmt = ALIGN + nleafs * header['ndim'] * 'i'
    block_ixs = np.reshape(struct.unpack(fmt, istream.read(struct.calcsize(fmt))),
                           [nleafs, header['ndim']])

    gfmt = ALIGN + (2 * ndim) * 'i'
    bcsize = struct.calcsize(gfmt)

    if header['datfile_version'] >= 6:
        # From version 6, the per-block ghost-cell counts are stored in the tree
        # as Fortran n_ghost(ndim, 2, nleafs).
        fmt = ALIGN + nleafs * 2 * ndim * 'i'
        block_nghost = np.reshape(
            struct.unpack(fmt, istream.read(struct.calcsize(fmt))),
            [nleafs, 2 * ndim])

        # Offsets point directly at the block data.
        fmt = ALIGN + nleafs * 'q'
        block_offsets = np.array(struct.unpack(fmt, istream.read(struct.calcsize(fmt))))
    else:
        # Up to version 5, each stored offset points at that block's ghost-count,
        # `2*ndim` ints.
        fmt = ALIGN + nleafs * 'q'
        raw_offsets = np.array(struct.unpack(fmt, istream.read(struct.calcsize(fmt))))
        block_offsets = raw_offsets + bcsize          # start of block data

        # Per-block ghost-cell counts, read from each block's header.
        block_nghost = np.zeros((nleafs, 2 * ndim), dtype=int)
        for i, raw in enumerate(raw_offsets):
            istream.seek(int(raw))
            block_nghost[i] = struct.unpack(gfmt, istream.read(bcsize))
    return block_lvls, block_ixs, block_offsets, block_nghost


def block_shape_with_ghost(header, nghost):
    """Per-block array shape including its ghost halo: (block_nx + lo + hi, ..., nw)."""
    ndim = header['ndim']
    lo = np.asarray(nghost[:ndim])
    hi = np.asarray(nghost[ndim:])
    return np.append(np.asarray(header['block_nx']) + lo + hi, header['nw'])


def strip_ghost(block, header, nghost):
    """Return the interior (block_nx) part of a ghosted block, dropping the halo."""
    ndim = header['ndim']
    lo = np.asarray(nghost[:ndim])
    hi = np.asarray(nghost[ndim:])
    sl = tuple(slice(int(lo[k]), block.shape[k] - int(hi[k])) for k in range(ndim))
    return block[sl + (slice(None),)]


def get_single_block_data(istream, byte_offset, block_shape, size_real=SIZE_DOUBLE):
    """"
    Retrieve a specific block from the datfile
    :param: istream       open datfile buffer in 'rb' mode
    :param: byte_offset   offset of the given block in the datfile
    :param: block_shape   the shape of the block (list containing dimensions + number of variables)
    :param: size_real     bytes per real in the block data (header['size_real'], 4 or 8)
    :return: block_data   numpy array containing the block data, with dimensions equal to block_shape
    """
    istream.seek(byte_offset)
    # Read actual data
    fmt = ALIGN + int(np.prod(block_shape)) * ('f' if size_real == 4 else 'd')
    d = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    # Fortran ordering
    block_data = np.reshape(d, block_shape, order='F')
    return block_data


def read_block(dataset, ileaf, keep_ghost=False):
    """Read leaf block `ileaf` at its true (ghosted) shape, stripping the halo
    unless keep_ghost=True. Handles both ghosted and ghost-free files."""
    nghost = dataset.block_nghost[ileaf]
    shape = block_shape_with_ghost(dataset.header, nghost)
    size_real = dataset.header.get('size_real', SIZE_DOUBLE)
    block = get_single_block_data(dataset.file, dataset.block_offsets[ileaf], shape, size_real)
    if not keep_ghost:
        block = strip_ghost(block, dataset.header, nghost)
    return block


def get_blocks(dataset, keep_ghost=False):
    """
    Reads all block data from an MPI-AMRVAC 2.0 snapshot.
    :param dataset      instance of 'amrvac_reader.load_file' class
    :param keep_ghost   if True, return each block with its ghost halo (for
                        per-block gradients / periodic-edge stencils); if False
                        (default) strip the halo so blocks tile the domain.
    :return list containing block data as dictionaries with level, morton index,
            data ('w') and ghost counts ('nghost').
    """
    size_real = dataset.header.get('size_real', SIZE_DOUBLE)
    hdr = dataset.header
    blocks = []
    for ileaf, offset in enumerate(dataset.block_offsets):
        nghost = dataset.block_nghost[ileaf]
        shape = block_shape_with_ghost(hdr, nghost)
        block = get_single_block_data(dataset.file, offset, shape, size_real)
        if not keep_ghost:
            block = strip_ghost(block, hdr, nghost)
        lvl = dataset.block_lvls[ileaf]
        ix = dataset.block_ixs[ileaf]
        b = {'lvl': lvl, 'ix': ix, 'w': block, 'nghost': nghost}
        blocks.append(b)

    return blocks



def get_uniform_data(dataset):
    """
    Retrieves the data for a uniform data set.
    :param dataset: instance of 'amrvac_reader.load_file' class
    :return The raw data as a NumPy array.
    """
    blocks = get_blocks(dataset)
    hdr = dataset.header

    refined_nx = 2 ** (hdr['levmax'] - 1) * hdr['domain_nx']
    domain_shape = np.append(refined_nx, hdr['nw'])
    d = np.zeros(domain_shape, order='F')

    for b in blocks:
        i0 = (b['ix'] - 1) * hdr['block_nx']
        i1 = i0 + hdr['block_nx']
        if hdr['ndim'] == 1:
            d[i0[0]:i1[0], :] = b['w']
        elif hdr['ndim'] == 2:
            d[i0[0]:i1[0], i0[1]:i1[1], :] = b['w']
        elif hdr['ndim'] == 3:
            d[i0[0]:i1[0], i0[1]:i1[1], i0[2]:i1[2], :] = b['w']
    return d
