module mod_zfp

#:if defined('COMPRESS_ZFP')
  use, intrinsic :: iso_c_binding
  use zfp
  use mod_comm_lib, only: mpistop
  implicit none
  private

  public :: zfp_padded_dim
  public :: zfp_max_field_bytes
  public :: zfp_compress_field_dp
  public :: zfp_compress_field_sp

  !> Staging arrays, dynamically reallocated.
  double precision, allocatable, target, save :: pad_dp(:,:,:)
  real(kind=4),     allocatable, target, save :: pad_sp(:,:,:)

  !> Corrected interface for zfp_write_header: the zFORp binding of ZFP is
  !> missing the `value` attribute on `mask`, causes garbage writes.
  interface
    function zfp_write_header_c(stream, field, mask) result(num_bits)&
        bind(c, name="zfp_write_header")
      import :: c_ptr, c_int, c_size_t
      type(c_ptr), value    :: stream, field
      integer(c_int), value :: mask
      integer(c_size_t)     :: num_bits
    end function zfp_write_header_c
  end interface

contains

  !> Need to pad, ZFP uses 4^3 blocks.
  elemental function zfp_padded_dim(n) result(np)
    integer, intent(in) :: n
    integer             :: np
    np = ((n + 3) / 4) * 4
  end function zfp_padded_dim

  !> Worst case (compression ratio 1) number of bytes of a compressed field,
  !> with header. (nx, ny, nz) are unpadded.
  function zfp_max_field_bytes(nx, ny, nz, write_sp, prec) result(nbytes)
    integer, intent(in) :: nx, ny, nz, prec
    logical, intent(in) :: write_sp
    integer(kind=8)     :: nbytes

    integer(kind=1), target :: scratch(64)
    type(zFORp_bitstream)   :: bs
    type(zFORp_stream)      :: stream
    type(zFORp_field)       :: field
    integer                 :: iprec

    field = zFORp_field_3d(c_null_ptr,&
        merge(zFORp_type_float, zFORp_type_double, write_sp),&
        zfp_padded_dim(nx), zfp_padded_dim(ny), zfp_padded_dim(nz))
    bs = zFORp_bitstream_stream_open(c_loc(scratch), int(size(scratch), 8))
    stream = zFORp_stream_open(bs)
    iprec = zFORp_stream_set_precision(stream, prec)

    ! zFORp_stream_maximum_size describes the field data. Add the header bits,
    ! but round to the word size. (x+63)/64 is ceiling.
    nbytes = zFORp_stream_maximum_size(stream, field) + &
        ((zFORp_header_max_bits + 63) / 64) * 8

    call zFORp_stream_close(stream)
    call zFORp_bitstream_stream_close(bs)
    call zFORp_field_free(field)
  end function zfp_max_field_bytes

  !> Compress one double precision field into buf(buf_start:), returning the
  !> byte count of the stream.
  function zfp_compress_field_dp(w, nx, ny, nz, prec, buf, buf_start, buf_size)&
      result(nbytes)
    integer, intent(in)          :: nx, ny, nz, prec
    double precision, intent(in) :: w(nx, ny, nz)
    integer(kind=8), intent(in)  :: buf_start, buf_size
    integer(kind=1), target, intent(inout) :: buf(buf_size)
    integer(kind=8)              :: nbytes

    integer :: i, j, k, px, py, pz

    px = zfp_padded_dim(nx)
    py = zfp_padded_dim(ny)
    pz = zfp_padded_dim(nz)

    if (allocated(pad_dp)) then
      if (any(shape(pad_dp) /= [px, py, pz])) deallocate(pad_dp)
    end if
    if (.not. allocated(pad_dp)) allocate(pad_dp(px, py, pz))

    ! Pad by edge replication on the right. It's important to keep the exponent
    ! bits uniform over the data, it improves compressability.
    do k = 1, pz
    do j = 1, py
    do i = 1, px
      pad_dp(i, j, k) = w(min(i, nx), min(j, ny), min(k, nz))
    end do
    end do
    end do

    nbytes = compress_padded(c_loc(pad_dp), zFORp_type_double, px, py, pz,&
        prec, buf, buf_start, buf_size)
  end function zfp_compress_field_dp

  !> Compress one single precision field into buf(buf_start:), returning the
  !> byte count of the stream.
  function zfp_compress_field_sp(w, nx, ny, nz, prec, buf, buf_start, buf_size)&
      result(nbytes)
    integer, intent(in)         :: nx, ny, nz, prec
    real(kind=4), intent(in)    :: w(nx, ny, nz)
    integer(kind=8), intent(in) :: buf_start, buf_size
    integer(kind=1), target, intent(inout) :: buf(buf_size)
    integer(kind=8)             :: nbytes

    integer :: i, j, k, px, py, pz

    px = zfp_padded_dim(nx)
    py = zfp_padded_dim(ny)
    pz = zfp_padded_dim(nz)

    if (allocated(pad_sp)) then
      if (any(shape(pad_sp) /= [px, py, pz])) deallocate(pad_sp)
    end if
    if (.not. allocated(pad_sp)) allocate(pad_sp(px, py, pz))

    ! Pad by edge replication on the right. It's important to keep the exponent
    ! bits uniform over the data, it improves compressability.
    do k = 1, pz
    do j = 1, py
    do i = 1, px
      pad_sp(i, j, k) = w(min(i, nx), min(j, ny), min(k, nz))
    end do
    end do
    end do

    nbytes = compress_padded(c_loc(pad_sp), zFORp_type_float, px, py, pz,&
        prec, buf, buf_start, buf_size)
  end function zfp_compress_field_sp

  !> Compress an already padded field. The header goes first.
  function compress_padded(field_ptr, scalar_type, px, py, pz, prec, buf,&
      buf_start, buf_size) result(nbytes)
    type(c_ptr), intent(in)     :: field_ptr
    integer, intent(in)         :: scalar_type, px, py, pz, prec
    integer(kind=8), intent(in) :: buf_start, buf_size
    integer(kind=1), target, intent(inout) :: buf(buf_size)
    integer(kind=8)             :: nbytes

    type(zFORp_bitstream) :: bs
    type(zFORp_stream)    :: stream
    type(zFORp_field)     :: field
    integer               :: iprec
    integer(kind=8)       :: header_bits

    field = zFORp_field_3d(field_ptr, scalar_type, px, py, pz)
    bs = zFORp_bitstream_stream_open(&
      c_loc(buf(buf_start)), buf_size - buf_start + 1)
    stream = zFORp_stream_open(bs)
    iprec = zFORp_stream_set_precision(stream, prec)

    header_bits = zfp_write_header_c(&
        transfer(stream, c_null_ptr),&
        transfer(field, c_null_ptr),&
        int(zFORp_header_full, c_int))
    if (header_bits == 0) call mpistop("mod_zfp: writing zfp header failed")

    nbytes = zFORp_compress(stream, field)
    if (nbytes == 0) call mpistop("mod_zfp: zfp compression failed")

    call zFORp_stream_close(stream)
    call zFORp_bitstream_stream_close(bs)
    call zFORp_field_free(field)
  end function compress_padded

#:endif

end module mod_zfp
