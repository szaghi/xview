!< xview test: the Lighthill tensor is the outer product of the momentum, LT_ij = u_i * u_j.
program xview_test_lighthill_tensor
!< xview test: the Lighthill tensor is the outer product of the momentum, LT_ij = u_i * u_j.
!<
!< Load a grd/icc/rst tuple of the xnavis-level-set case with the Lighthill tensor computed at cell centers (no
!< interpolation at nodes), then check every cell of every block: LT must be allocated with the momentum bounds, equal
!< bitwise the outer product of the momentum and therefore be symmetric. Run from the repository root.
use penf
use vecfor
use xview_file_grd_object
use xview_file_icc_object
use xview_file_rst_object

implicit none
character(*), parameter :: ipath='src/tests/xnavis-level-set/input/' !< Path to input files.
type(file_grd_object)   :: file_grd                                  !< File grd.
type(file_icc_object)   :: file_icc                                  !< File icc.
type(file_rst_object)   :: file_rst                                  !< File rst.
real(R8P)               :: u(1:3)                                    !< Momentum components.
integer(I4P)            :: cells_number                              !< Number of checked cells.
integer(I4P)            :: b, i, j, k, m, n                          !< Counter.

call file_grd%load_file(filename=ipath//'CC/cc.02.grd.p002')
call file_icc%load_file(filename=ipath//'CC/cc.02.p002', file_grd=file_grd, correct_metrics=.false., &
                        is_cell_centered=.true.)
call file_rst%load_file(file_grd=file_grd, file_icc=file_icc, filename=ipath//'RST/boa_00.02.p002', &
                        is_level_set=.true., is_oneeq=.true., is_cell_centered=.true., compute_LT=.true.)
if (.not.file_rst%is_loaded) error stop 'error: restart file not loaded'

cells_number = 0
do b=1, file_rst%blocks_number
   if (.not.allocated(file_rst%blocks(b)%LT)) error stop 'error: Lighthill tensor not allocated'
   associate(momentum=>file_rst%blocks(b)%momentum, LT=>file_rst%blocks(b)%LT)
   if (any(lbound(LT)/=[1, 1, lbound(momentum)]).or.any(ubound(LT)/=[3, 3, ubound(momentum)])) &
      error stop 'error: Lighthill tensor bounds differ from momentum ones'
   do k=lbound(momentum, dim=3), ubound(momentum, dim=3)
      do j=lbound(momentum, dim=2), ubound(momentum, dim=2)
         do i=lbound(momentum, dim=1), ubound(momentum, dim=1)
            u = [momentum(i,j,k)%x, momentum(i,j,k)%y, momentum(i,j,k)%z]
            do n=1, 3
               do m=1, 3
                  if (LT(m,n,i,j,k)/=u(m)*u(n)) error stop 'error: LT_ij /= u_i * u_j'
                  if (LT(m,n,i,j,k)/=LT(n,m,i,j,k)) error stop 'error: Lighthill tensor not symmetric'
               enddo
            enddo
            cells_number = cells_number + 1
         enddo
      enddo
   enddo
   endassociate
enddo
if (cells_number==0) error stop 'error: no cell checked'
print '(A)', 'Lighthill tensor checked on '//trim(str(cells_number, no_sign=.true.))//' cells'
endprogram xview_test_lighthill_tensor
