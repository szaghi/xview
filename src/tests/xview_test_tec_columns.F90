!< xview test: every Tecplot zone carries exactly the variables declared in the file header.
program xview_test_tec_columns
!< xview test: every Tecplot zone carries exactly the variables declared in the file header.
!<
!< Load a grd/icc/rst tuple of the xnavis-level-set case with many auxiliary fields, save it as a cell-centered ascii
!< Tecplot file, read it back and check that the header declares `file_tec%nvar` names and that every zone holds
!< 3 nodal (x, y, z) plus `nvar-3` cell-centered variables. Run from the repository root.
use penf
use xview_file_grd_object
use xview_file_icc_object
use xview_file_rst_object
use xview_file_tec_object

implicit none
character(*), parameter :: ipath='src/tests/xnavis-level-set/input/' !< Path to input files.
character(*), parameter :: opath='exe/'                              !< Path to output files.
character(*), parameter :: oname='xview_test_tec_columns'            !< Output file base name.
type(file_grd_object)   :: file_grd                                  !< File grd.
type(file_icc_object)   :: file_icc                                  !< File icc.
type(file_rst_object)   :: file_rst                                  !< File rst.
type(file_tec_object)   :: file_tec                                  !< File Tecplot.
character(1024)         :: line                                      !< Line buffer.
integer(I4P)            :: file_unit                                 !< File unit.
integer(I4P)            :: iostat                                    !< IO status.
integer(I4P)            :: nvar                                      !< Number of variables in the header.
integer(I4P)            :: zones_number                              !< Number of zones read.
integer(I4P)            :: values_number                             !< Number of values read in the current zone.
integer(I4P)            :: values_expected                           !< Number of values expected in the current zone.
integer(I4P)            :: ni, nj, nk                                !< Zone nodes dimensions.
integer(I4P)            :: b                                         !< Counter.

call file_grd%load_file(filename=ipath//'CC/cc.02.grd.p002', compute_metrics=.true.) ! as the CLI does for gradient fields
call file_icc%load_file(filename=ipath//'CC/cc.02.p002', file_grd=file_grd, correct_metrics=.true.,  &
                        is_cell_centered=.true.)
call file_icc%compute_patches_extents(file_grd=file_grd)
call file_rst%load_file(file_grd=file_grd, file_icc=file_icc, filename=ipath//'RST/boa_00.02.p002',               &
                        is_level_set=.true., is_oneeq=.true., is_cell_centered=.true.,                             &
                        compute_lambda2=.true., compute_qfactor=.true., compute_liutex=.true.,                     &
                        compute_helicity=.true., compute_vorticity=.true., compute_grad_p=.true., compute_LT=.true., &
                        compute_div2LT=.true.)
if (.not.file_rst%is_loaded) error stop 'error: restart file not loaded'

call file_tec%initialize(path=opath, file_name=oname, is_binary=.false., is_cell_centered=.true., &
                         icc=file_icc%blocks(1), sol=file_rst%blocks(1))
do b=1, file_grd%blocks_number
   call file_tec%save_block_file_tec(is_binary=.false., blocks_map=[0_I4P, b], grd=file_grd%blocks(b), &
                                     is_cell_centered=.true., icc=file_icc%blocks(b), sol=file_rst%blocks(b))
enddo
call file_tec%finalize(is_binary=.false.)

open(newunit=file_unit, file=opath//oname//'.dat', action='read')
read(file_unit, '(A)') line
nvar = count_char(trim(line), '"') / 2
if (nvar/=file_tec%nvar) error stop 'error: header names number differs from file_tec%nvar'
if (nvar<=3) error stop 'error: no variable other than the coordinates saved'
zones_number = 0
values_number = 0
values_expected = 0
do
   read(file_unit, '(A)', iostat=iostat) line
   if (is_iostat_end(iostat)) exit
   if (index(line, 'ZONE')>0) then
      call check_zone
      ni = zone_dimension('I=')
      nj = zone_dimension('J=')
      nk = zone_dimension('K=')
      values_expected = 3 * ni * nj * nk + (nvar - 3) * (ni - 1) * (nj - 1) * (nk - 1)
      values_number = 0
      zones_number = zones_number + 1
   else
      values_number = values_number + count_tokens(trim(line))
   endif
enddo
call check_zone
close(file_unit, status='delete')
if (zones_number/=file_grd%blocks_number) error stop 'error: zones number differs from blocks number'
print '(A)', trim(str(zones_number, no_sign=.true.))//' zones checked against '//trim(str(nvar, no_sign=.true.))//' variables'
contains
   subroutine check_zone()
   !< Check the number of values of the zone just read.
   if (values_number/=values_expected) then
      print '(A)', 'zone '//trim(str(zones_number, no_sign=.true.))//': values '//trim(str(values_number, no_sign=.true.))// &
                   ', expected '//trim(str(values_expected, no_sign=.true.))
      error stop 'error: zone values number differs from the declared variables'
   endif
   endsubroutine check_zone

   pure function count_char(string, char) result(n)
   !< Return the number of occurrences of a character into a string.
   character(*), intent(in) :: string !< String.
   character(1), intent(in) :: char   !< Character.
   integer(I4P)             :: n      !< Number of occurrences.
   integer(I4P)             :: c      !< Counter.

   n = 0
   do c=1, len(string)
      if (string(c:c)==char) n = n + 1
   enddo
   endfunction count_char

   pure function count_tokens(string) result(n)
   !< Return the number of blank-separated tokens of a string.
   character(*), intent(in) :: string !< String.
   integer(I4P)             :: n      !< Number of tokens.
   integer(I4P)             :: c      !< Counter.

   n = 0
   do c=1, len(string)
      if (string(c:c)/=' '.and.(c==1.or.string(max(c-1,1):max(c-1,1))==' ')) n = n + 1
   enddo
   endfunction count_tokens

   function zone_dimension(key) result(n)
   !< Return the zone dimension following the given key (e.g. 'I=') in the zone header line.
   character(*), intent(in) :: key !< Dimension key.
   integer(I4P)             :: n   !< Dimension.
   integer(I4P)             :: c   !< Key position.

   c = index(line, key) + len(key)
   read(line(c:index(line(c:), ',')+c-2), *) n
   endfunction zone_dimension
endprogram xview_test_tec_columns
