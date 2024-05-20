!< xview, file rst class definition.
module xview_file_rst_object
!< xview, file rst class definition.
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf

use xview_block_rst_object

implicit none
private
public :: file_rst_object

type :: file_rst_object
   !< File rst class definition.
   character(len=:), allocatable       :: file_name         !< File name of sol file.
   integer(I4P)                        :: blocks_number=0   !< Number of blocks contained into the file.
   logical                             :: is_loaded=.false. !< Flag for checking if the file has been correctly loaded.
   type(block_rst_object), allocatable :: blocks(:)         !< Blocks contained into the sol file.
   real(R8P)                           :: time=0._R8P       !< Time of solution.
   contains
      ! plublic methods
      procedure, pass(self) :: destroy         !< Destroy dynamic memory.
      procedure, pass(self) :: alloc           !< Allocate dynamic memory.
      procedure, pass(self) :: is_file_present !< Inquire if the file path is valid.
      procedure, pass(self) :: load_file       !< Load file.
      procedure, pass(self) :: save_file       !< Save file.
endtype file_rst_object

contains
   ! plublic methods
   elemental subroutine destroy(self)
   !< Destroy dynamic memory.
   class(file_rst_object), intent(inout) :: self !< File data.
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
   self%time = 0._R8P
   endsubroutine destroy

   elemental subroutine alloc(self)
   !< Allocate dynamic memory.
   class(file_rst_object), intent(inout) :: self !< File data.

   allocate(self%blocks(1:self%blocks_number))
   endsubroutine alloc

   function is_file_present(self) result(is_present)
   !< Inquire if the file path is valid.
   class(file_rst_object), intent(in) :: self     !< File data.
   logical                            :: is_present !< Inquiring result.

   is_present = allocated(self%file_name)
   if (is_present) inquire(file=trim(adjustl(self%file_name)), exist=is_present)
   endfunction is_file_present

   subroutine load_file(self, blocks_number, file_name, is_level_set, turbulence_model, is_cell_centered)
   !< Load file.
   class(file_rst_object), intent(inout)        :: self             !< File data.
   integer(I4P),           intent(in)           :: blocks_number    !< Number of blocks contained into the file.
   character(*),           intent(in)           :: file_name        !< File name of geo file.
   logical,                intent(in), optional :: is_level_set     !< Flag for level set function presence.
   character(*),           intent(in), optional :: turbulence_model !< Turbulence model: 'zero_eq', 'one_eq', 'two_eq'.
   logical,                intent(in), optional :: is_cell_centered !< Define variables at cell centers or nodes.
   integer(I4P)                                 :: file_unit        !< Logical file unit of geo file.
   integer(I4P)                                 :: b                !< Counter.

   call self%destroy
   self%file_name = trim(adjustl(file_name))
   if (self%is_file_present()) then
      self%blocks_number = blocks_number
      call self%alloc
      open(newunit=file_unit, file=trim(adjustl(file_name)), form='unformatted', action='read')
      read(file_unit) self%time
      do b=1, self%blocks_number
         call self%blocks(b)%load_dimensions(file_unit=file_unit, is_level_set=is_level_set, turbulence_model=turbulence_model)
      enddo
      do b=1, self%blocks_number
         call self%blocks(b)%load_solution(file_unit=file_unit, is_cell_centered=is_cell_centered)
      enddo
      close(file_unit)
      self%is_loaded = .true.
   else
      write(stderr, "(A)")'error: file "'//trim(adjustl(file_name))//'" not found!'
      self%is_loaded = .false.
   endif
   endsubroutine load_file

   subroutine save_file(self, file_name)
   !< Save file.
   class(file_rst_object), intent(in)           :: self      !< File data.
   character(*),           intent(in), optional :: file_name !< File name.
   integer(I4P)                                 :: file_unit !< Logical file unit.
   integer(I4P)                                 :: b         !< Counter.

   if (present(file_name)) then
      open(newunit=file_unit, file=trim(adjustl(file_name)), form='unformatted', action='write')
   elseif (allocated(self%file_name)) then
      open(newunit=file_unit, file=trim(adjustl(self%file_name)), form='unformatted', action='write')
   else
      error stop 'error: nor "file_name" neither "self%file_name" have been specified for "file_rst_object%save_file" method!'
   endif
   write(file_unit) self%time
   do b=1, self%blocks_number
      call self%blocks(b)%save_dimensions(file_unit=file_unit)
   enddo
   do b=1, self%blocks_number
      call self%blocks(b)%save_solution(file_unit=file_unit)
   enddo
   close(file_unit)
   endsubroutine save_file
endmodule xview_file_rst_object
