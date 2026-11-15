!< xview test: cell-to-node interpolation is exact for a linear field on every node, the last node planes included.
program xview_test_node_interpolation
!< xview test: cell-to-node interpolation is exact for a linear field on every node, the last node planes included.
!<
!< A block with ghost cells (gc = 2, arrays from 1-gc) holds a pressure linear in the cell indices; node i lies between
!< cells i and i+1, so the interpolated value at node (i, j, k), i = 0..Ni, must be the field at (i+1/2, j+1/2, k+1/2).
!< An assumed-shape (0:) dummy in the interpolation rebases the arrays and leaves node planes Ni, Nj, Nk wrong.
!< The momentum components are checked too: they reach the interpolation as non-contiguous actual arguments.
use penf
use xview_block_rst_object

implicit none
integer(I4P), parameter :: Ni=4, Nj=3, Nk=2 !< Block dimensions (cells).
type(block_rst_object)  :: block            !< Block.
real(R8P)               :: expected         !< Expected node value.
real(R8P)               :: err_max          !< Maximum error over the nodes.
integer(I4P)            :: i, j, k          !< Counters.

call block%init(Ni=Ni, Nj=Nj, Nk=Nk, is_oneeq=.false.)
do k=lbound(block%pressure, dim=3), ubound(block%pressure, dim=3)
   do j=lbound(block%pressure, dim=2), ubound(block%pressure, dim=2)
      do i=lbound(block%pressure, dim=1), ubound(block%pressure, dim=1)
         block%pressure(i,j,k) = field(real(i, R8P), real(j, R8P), real(k, R8P))
         block%momentum(i,j,k)%x =  2._R8P * block%pressure(i,j,k)
         block%momentum(i,j,k)%y = -3._R8P * block%pressure(i,j,k)
         block%momentum(i,j,k)%z =  0.5_R8P * block%pressure(i,j,k)
      enddo
   enddo
enddo
call block%interpolate_at_nodes
err_max = 0._R8P
do k=0, Nk
   do j=0, Nj
      do i=0, Ni
         expected = field(i + 0.5_R8P, j + 0.5_R8P, k + 0.5_R8P)
         err_max = max(err_max, abs(block%pressure(i,j,k) - expected))
         err_max = max(err_max, abs(block%momentum(i,j,k)%x -  2._R8P  * expected))
         err_max = max(err_max, abs(block%momentum(i,j,k)%y - (-3._R8P * expected)))
         err_max = max(err_max, abs(block%momentum(i,j,k)%z -  0.5_R8P * expected))
      enddo
   enddo
enddo
print '(A,ES10.3)', 'max |node - exact| (p, u, v, w) over nodes 0..N: ', err_max
if (err_max > 1.e-12_R8P) error stop 'error: cell-to-node interpolation is not exact for a linear field'

contains
   pure function field(x, y, z) result(f)
   !< Linear test field.
   real(R8P), intent(in) :: x, y, z !< Cell (or node) index coordinates.
   real(R8P)             :: f       !< Field value.

   f = 1._R8P + 2._R8P * x + 30._R8P * y + 500._R8P * z
   endfunction field
endprogram xview_test_node_interpolation
