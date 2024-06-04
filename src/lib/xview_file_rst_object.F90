!< xview, file rst class definition.
module xview_file_rst_object
!< xview, file rst class definition.
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf

use xview_block_rst_object
use xview_file_object

implicit none
private
public :: file_rst_object

type, extends(file_object) :: file_rst_object
   !< File rst class definition.
   integer(I4P)                        :: blocks_number=0 !< Number of blocks contained into the file.
   type(block_rst_object), allocatable :: blocks(:)       !< Blocks contained into the sol file.
   real(R8P)                           :: time=0._R8P     !< Time of solution.
   contains
      ! plublic methods
      procedure, pass(self) :: destroy   !< Destroy dynamic memory.
      procedure, pass(self) :: alloc     !< Allocate dynamic memory.
      procedure, pass(self) :: load_file !< Load file.
      procedure, pass(self) :: save_file !< Save file.
endtype file_rst_object

contains
   ! plublic methods
   elemental subroutine destroy(self)
   !< Destroy dynamic memory.
   class(file_rst_object), intent(inout) :: self !< File data.
   integer(I4P)                          :: b    !< Counter.

   call self%file_object%destroy
   self%blocks_number = 0
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

   subroutine load_file(self, blocks_number, filename, is_level_set, is_zeroeq, is_oneeq, is_twoeq, is_cell_centered, verbose)
   !< Load file.
   class(file_rst_object), intent(inout)        :: self             !< File data.
   integer(I4P),           intent(in)           :: blocks_number    !< Number of blocks contained into the file.
   character(*),           intent(in)           :: filename         !< File name of geo file.
   logical,                intent(in), optional :: is_level_set     !< Flag for level set function presence.
   logical,                intent(in), optional :: is_zeroeq        !< Use *zero* equations turbulence model.
   logical,                intent(in), optional :: is_oneeq         !< Use *one* equations turbulence model.
   logical,                intent(in), optional :: is_twoeq         !< Use *two* equations turbulence model.
   logical,                intent(in), optional :: is_cell_centered !< Define variables at cell centers or nodes.
   logical,                intent(in), optional :: verbose          !< Activate verbose mode.
   logical                                      :: verbose_         !< Activate verbose mode, local variable.
   integer(I4P)                                 :: file_unit        !< Logical file unit of geo file.
   integer(I4P)                                 :: b                !< Counter.

   call self%destroy
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   self%filename = trim(adjustl(filename))
   if (self%is_file_present()) then
      self%blocks_number = blocks_number
      call self%alloc
      open(newunit=file_unit, file=trim(adjustl(filename)), form='unformatted', action='read')
      read(file_unit) self%time
      do b=1, self%blocks_number
         call self%blocks(b)%load_dimensions(file_unit=file_unit, &
                                             is_level_set=is_level_set, is_zeroeq=is_zeroeq, is_oneeq=is_oneeq, is_twoeq=is_twoeq)
      enddo
      do b=1, self%blocks_number
         call self%blocks(b)%load_solution(file_unit=file_unit, is_cell_centered=is_cell_centered)
      enddo
      close(file_unit)
      self%is_loaded = .true.
   else
      if (verbose_) write(stderr, "(A)") 'warning: file "'//trim(adjustl(filename))//'" not found!'
      self%is_loaded = .false.
   endif
   endsubroutine load_file

   subroutine save_file(self, filename)
   !< Save file.
   class(file_rst_object), intent(in)           :: self      !< File data.
   character(*),           intent(in), optional :: filename  !< File name.
   integer(I4P)                                 :: file_unit !< Logical file unit.
   integer(I4P)                                 :: b         !< Counter.

   if (present(filename)) then
      open(newunit=file_unit, file=trim(adjustl(filename)), form='unformatted', action='write')
   elseif (self%filename%is_allocated()) then
      open(newunit=file_unit, file=trim(adjustl(self%filename%chars())), form='unformatted', action='write')
   else
      error stop 'error: nor "filename" neither "self%filename" have been specified for "file_rst_object%save_file" method!'
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
