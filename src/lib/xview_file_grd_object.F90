!< xview, file grd class definition.
module xview_file_grd_object
!< xview, file grd class definition.
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf
use stringifor

use xview_block_grd_object
use xview_file_object

implicit none
private
public :: file_grd_object

type, extends(file_object) :: file_grd_object
   !< **File grd** class, parse, manipulate and emit Xnavis grd file.
   integer(I4P)                        :: blocks_number=0 !< Number of blocks contained into the file.
   type(block_grd_object), allocatable :: blocks(:)       !< Blocks contained into the grd file.
   contains
      ! plublic methods
      procedure, pass(self) :: destroy   !< Destroy dynamic memory.
      procedure, pass(self) :: alloc     !< Allocate dynamic memory.
      procedure, pass(self) :: load_file !< Load file.
      procedure, pass(self) :: save_file !< Save file.
endtype file_grd_object

contains
   ! plublic methods
   elemental subroutine destroy(self)
   !< Destroy dynamic memory.
   class(file_grd_object), intent(inout) :: self !< File data.
   integer(I4P)                          :: b    !< Counter.

   call self%file_object%destroy
   self%blocks_number = 0
   if (allocated(self%blocks)) then
      do b=lbound(self%blocks, dim=1), ubound(self%blocks, dim=1)
         call self%blocks(b)%destroy
      enddo
      deallocate(self%blocks)
   endif
   endsubroutine destroy

   elemental subroutine alloc(self)
   !< Allocate dynamic memory.
   class(file_grd_object), intent(inout) :: self !< File data.

   allocate(self%blocks(1:self%blocks_number))
   endsubroutine alloc

   subroutine load_file(self, filename, is_centers_to_compute, is_extents_to_compute, verbose)
   !< Load file.
   class(file_grd_object), intent(inout)        :: self                  !< File data.
   character(*),           intent(in)           :: filename              !< File name.
   logical,                intent(in), optional :: is_centers_to_compute !< Flag to activate the computation of cell centers.
   logical,                intent(in), optional :: is_extents_to_compute !< Flag to activate the computation of extents.
   logical,                intent(in), optional :: verbose               !< Activate verbose mode.
   logical                                      :: verbose_              !< Activate verbose mode, local variable.
   integer(I4P)                                 :: file_unit             !< Logical file unit.
   integer(I4P)                                 :: b                     !< Counter.

   call self%destroy
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   self%filename = trim(adjustl(filename))
   if (self%is_file_present()) then
      open(newunit=file_unit, file=trim(adjustl(filename)), form='unformatted', action='read')
      read(file_unit, end=10, err=10) self%blocks_number
      call self%alloc
      do b=1, self%blocks_number
         call self%blocks(b)%load_dimensions(file_unit=file_unit,                          &
                                             is_centers_to_allocate=is_centers_to_compute, &
                                             is_extents_to_allocate=is_extents_to_compute)
      enddo
      do b=1, self%blocks_number
         call self%blocks(b)%load_nodes(file_unit=file_unit)
      enddo
      10 close(file_unit)
      self%is_loaded = .true.
   else
      if (verbose_) write(stderr, "(A)") 'error: file "'//trim(adjustl(filename))//'" not found!'
      self%is_loaded = .false.
   endif
   endsubroutine load_file

   subroutine save_file(self, filename)
   !< Save file.
   class(file_grd_object), intent(in)           :: self      !< File data.
   character(*),           intent(in), optional :: filename  !< File name.
   integer(I4P)                                 :: file_unit !< Logical file unit.
   integer(I4P)                                 :: b         !< Counter.

   if (present(filename)) then
      open(newunit=file_unit, file=trim(adjustl(filename)), form='unformatted', action='write')
   elseif (self%filename%is_allocated()) then
      open(newunit=file_unit, file=trim(adjustl(self%filename%chars())), form='unformatted', action='write')
   else
      error stop 'error: nor "filename" neither "self%filename" have been specified for "file_grd_object%save_file" method!'
   endif
   write(file_unit) self%blocks_number
   do b=1, self%blocks_number
      call self%blocks(b)%save_dimensions(file_unit=file_unit)
   enddo
   do b=1, self%blocks_number
      call self%blocks(b)%save_nodes(file_unit=file_unit)
   enddo
   close(file_unit)
   endsubroutine save_file
endmodule xview_file_grd_object
