!< xview, file icc class definition.
module xview_file_icc_object
!< xview, file icc class definition.
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf

use xview_block_icc_object
use xview_file_object
use xview_file_grd_object

implicit none
private
public :: file_icc_object

type, extends(file_object) :: file_icc_object
   !< File icc class definition.
   integer(I4P)                        :: blocks_number=0      !< Number of blocks contained into the file.
   type(block_icc_object), allocatable :: blocks(:)            !< Blocks contained into the icc file.
   integer(I4P)                        :: unstruct_dimension=0 !< Dimension of unstructured array of rcc.
   real(R4P), allocatable              :: rcc(:)               !< rcc unstructured array.
   contains
      ! public methods
      procedure, pass(self) :: destroy                 !< Destroy dynamic memory.
      procedure, pass(self) :: alloc                   !< Allocate dynamic memory.
      procedure, pass(self) :: compute_patches_extents !< Compute patches extents for given bc or whole block extents.
      procedure, pass(self) :: load_file               !< Load file.
      procedure, pass(self) :: save_file               !< Save file.
endtype file_icc_object

contains
   elemental subroutine destroy(self)
   !< Destroy dynamic memory.
   class(file_icc_object), intent(inout) :: self !< File data.
   integer(I4P)                          :: b    !< Counter.

   call self%file_object%destroy
   self%blocks_number = 0
   if (allocated(self%blocks)) then
      do b=lbound(self%blocks, dim=1), ubound(self%blocks, dim=1)
         call self%blocks(b)%destroy
      enddo
      deallocate(self%blocks)
   endif
   self%unstruct_dimension = 0
   if (allocated(self%rcc)) deallocate(self%rcc)
   endsubroutine destroy

   elemental subroutine alloc(self, blocks, rcc)
   !< Allocate dynamic memory.
   class(file_icc_object), intent(inout)        :: self   !< File data.
   logical,                intent(in), optional :: blocks !< Flag to allocate blocks.
   logical,                intent(in), optional :: rcc    !< Flag to allocate rcc unstructured array.

   if (present(blocks)) then
      if (blocks) allocate(self%blocks(1:self%blocks_number))
   endif
   if (present(rcc)) then
      if (rcc)  allocate(self%rcc(1:self%unstruct_dimension))
   endif
   endsubroutine alloc

   pure subroutine compute_patches_extents(self, file_grd, patch, offset)
   !< Compute the patches extents of given patch boundary conditions or the use the whole block.
   class(file_icc_object), intent(in)           :: self     !< File icc data.
   type(file_grd_object),  intent(inout)        :: file_grd !< File grd data.
   integer(I4P),           intent(in), optional :: patch    !< Patch bc to be found.
   integer(I4P),           intent(in), optional :: offset   !< Offset from patch.
   integer(I4P)                                 :: b        !< Counter.

   do b=1, file_grd%blocks_number
      if (self%is_loaded.and.present(patch)) then
         ! specific BC has been queried
         call file_grd%blocks(b)%compute_patches_extents(tcc=self%blocks(b)%tcc, patch=patch)
      else
         ! whole blocks extents have been queried
         call file_grd%blocks(b)%compute_patches_extents
      endif
   enddo
   endsubroutine compute_patches_extents

   subroutine load_file(self, filename, file_grd, correct_metrics, is_cell_centered, verbose)
   !< Load file.
   class(file_icc_object), intent(inout)        :: self             !< File data.
   character(*),           intent(in)           :: filename         !< File name.
   type(file_grd_object),  intent(inout)        :: file_grd         !< File grd data.
   logical,                intent(in)           :: correct_metrics  !< Flag to activate the correction of metrics.
   logical,                intent(in), optional :: is_cell_centered !< Define variables at cell centers or nodes.
   logical,                intent(in), optional :: verbose          !< Activate verbose mode.
   logical                                      :: verbose_         !< Activate verbose mode, local variable.
   integer(I4P)                                 :: file_unit        !< Logical file unit.
   integer(I4P)                                 :: b                !< Counter.
   integer(I4P)                                 :: i                !< Counter.

   call self%destroy
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   self%filename = trim(adjustl(filename))
   if (self%is_file_present()) then
      open(newunit=file_unit, file=trim(adjustl(filename)), form='unformatted', action='read')
      read(file_unit) self%blocks_number
      call self%alloc(blocks=.true.)
      do b=1, self%blocks_number
         call self%blocks(b)%load_dimensions(file_unit=file_unit)
      enddo
      do b=1, self%blocks_number
         call self%blocks(b)%load_icc(file_unit=file_unit)
      enddo
      read(file_unit) self%unstruct_dimension
      call self%alloc(rcc=.true.)
      read(file_unit) (self%rcc(i), i=1, self%unstruct_dimension)
      10 close(file_unit)
      self%is_loaded = .true.
   else
      if (verbose_) write(stderr, "(A)") 'error: file "'//trim(adjustl(filename))//'" not found!'
      self%is_loaded = .false.
   endif
   do b=1, self%blocks_number
      call self%blocks(b)%compute_cc(rcc=self%rcc, is_cell_centered=is_cell_centered)
   enddo
   if (correct_metrics) then
      do b=1, file_grd%blocks_number
         call file_grd%blocks(b)%correct_metrics_bc(tcc=self%blocks(b)%tcc)
      enddo
   endif
   endsubroutine load_file

   subroutine save_file(self, filename)
   !< Save file.
   class(file_icc_object), intent(in)           :: self      !< File data.
   character(*),           intent(in), optional :: filename  !< File name.
   integer(I4P)                                 :: file_unit !< Logical file unit.
   integer(I4P)                                 :: b         !< Counter.
   integer(I4P)                                 :: i         !< Counter.

   if (present(filename)) then
      open(newunit=file_unit, file=trim(adjustl(filename)), form='unformatted', action='write')
   elseif (self%filename%is_allocated()) then
      open(newunit=file_unit, file=trim(adjustl(self%filename%chars())), form='unformatted', action='write')
   else
      error stop 'error: nor "filename" neither "self%filename" have been specified for "file_icc_object%save_file" method!'
   endif
   write(file_unit) self%blocks_number
   do b=1, self%blocks_number
      call self%blocks(b)%save_dimensions(file_unit=file_unit)
   enddo
   do b=1, self%blocks_number
      call self%blocks(b)%save_icc(file_unit=file_unit)
   enddo
   write(file_unit) self%unstruct_dimension
   do i=1, self%unstruct_dimension
      write(file_unit) self%rcc(i)
   enddo
   close(file_unit)
   endsubroutine save_file
endmodule xview_file_icc_object
