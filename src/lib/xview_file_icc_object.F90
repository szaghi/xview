!< xview, file icc class definition.
module xview_file_icc_object
!< xview, file icc class definition.
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf

use xview_block_icc_object

implicit none
private
public :: file_icc_object

type :: file_icc_object
   !< File icc class definition.
   character(len=:), allocatable       :: file_name            !< File name.
   integer(I4P)                        :: blocks_number=0      !< Number of blocks contained into the file.
   logical                             :: is_loaded=.false.    !< Flag for checking if the file has been correctly loaded.
   type(block_icc_object), allocatable :: blocks(:)            !< Blocks contained into the icc file.
   integer(I4P)                        :: unstruct_dimension=0 !< Dimension of unstructured array of rcc.
   real(R4P), allocatable              :: rcc(:)               !< rcc unstructured array.
   contains
      ! public methods
      procedure, pass(self) :: destroy         !< Destroy dynamic memory.
      procedure, pass(self) :: alloc           !< Allocate dynamic memory.
      procedure, pass(self) :: is_file_present !< Inquire if the file path is valid.
      procedure, pass(self) :: load_file       !< Load file.
      procedure, pass(self) :: save_file       !< Save file.
endtype file_icc_object

contains
   elemental subroutine destroy(self)
   !< Destroy dynamic memory.
   class(file_icc_object), intent(inout) :: self !< File data.
   integer(I4P)                          :: b    !< Counter.

   if (allocated(self%file_name)) deallocate(self%file_name)
   self%blocks_number = 0
   self%is_loaded = .false.
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

   function is_file_present(self) result(is_present)
   !< Inquire if the file path is valid.
   class(file_icc_object), intent(in) :: self       !< File data.
   logical                            :: is_present !< Inquiring result.

   is_present = allocated(self%file_name)
   if (is_present) inquire(file=trim(adjustl(self%file_name)), exist=is_present)
   endfunction is_file_present

   subroutine load_file(self, file_name, is_cell_centered)
   !< Load file.
   class(file_icc_object), intent(inout)        :: self             !< File data.
   character(*),           intent(in)           :: file_name        !< File name.
   logical,                intent(in), optional :: is_cell_centered !< Define variables at cell centers or nodes.
   integer(I4P)                                 :: file_unit        !< Logical file unit.
   integer(I4P)                                 :: b                !< Counter.
   integer(I4P)                                 :: i                !< Counter.

   call self%destroy
   self%file_name = trim(adjustl(file_name))
   if (self%is_file_present()) then
      open(newunit=file_unit, file=trim(adjustl(file_name)), form='unformatted', action='read')
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
      write(stderr, "(A)")'error: file "'//trim(adjustl(file_name))//'" not found!'
      self%is_loaded = .false.
   endif
   do b=1, self%blocks_number
      call self%blocks(b)%compute_cc(rcc=self%rcc, is_cell_centered=is_cell_centered)
   enddo
   endsubroutine load_file

   subroutine save_file(self, file_name)
   !< Save file.
   class(file_icc_object), intent(in)           :: self      !< File data.
   character(*),           intent(in), optional :: file_name !< File name.
   integer(I4P)                                 :: file_unit !< Logical file unit.
   integer(I4P)                                 :: b         !< Counter.
   integer(I4P)                                 :: i         !< Counter.

   if (present(file_name)) then
      open(newunit=file_unit, file=trim(adjustl(file_name)), form='unformatted', action='write')
   elseif (allocated(self%file_name)) then
      open(newunit=file_unit, file=trim(adjustl(self%file_name)), form='unformatted', action='write')
   else
      error stop 'error: nor "file_name" neither "self%file_name" have been specified for "file_icc_object%save_file" method!'
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
