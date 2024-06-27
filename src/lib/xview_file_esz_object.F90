!< xview, file extracted subzone class definition.
module xview_file_esz_object
!< xview, file extracted subzone class definition.
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf

use xview_block_esz_object
use xview_file_object

implicit none
private
public :: file_esz_object

type, extends(file_object) :: file_esz_object
   !< File rst class definition.
   type(block_esz_object) :: block_esz !< Extracted block subzone.
   contains
      ! plublic methods
      procedure, pass(self) :: destroy   !< Destroy dynamic memory.
      procedure, pass(self) :: load_file !< Load file.
endtype file_esz_object

contains
   ! plublic methods
   elemental subroutine destroy(self)
   !< Destroy dynamic memory.
   class(file_esz_object), intent(inout) :: self !< File data.

   call self%file_object%destroy
   call self%block_esz%destroy
   endsubroutine destroy

   subroutine load_file(self, filename, is_level_set, is_zeroeq, is_oneeq, is_twoeq, verbose)
   !< Load file.
   class(file_esz_object), intent(inout)        :: self         !< File data.
   character(*),           intent(in)           :: filename     !< File name of geo file.
   logical,                intent(in), optional :: is_level_set !< Flag for level set function presence.
   logical,                intent(in), optional :: is_zeroeq    !< Use *zero* equations turbulence model.
   logical,                intent(in), optional :: is_oneeq     !< Use *one* equations turbulence model.
   logical,                intent(in), optional :: is_twoeq     !< Use *two* equations turbulence model.
   logical,                intent(in), optional :: verbose      !< Activate verbose mode.
   logical                                      :: verbose_     !< Activate verbose mode, local variable.
   integer(I4P)                                 :: file_unit    !< Logical file unit.

   call self%destroy
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   self%filename = trim(adjustl(filename))
   if (self%is_file_present()) then
      open(newunit=file_unit, file=trim(adjustl(filename)), form='unformatted', action='read')
      call self%block_esz%load_solution(file_unit=file_unit, &
                                        is_level_set=is_level_set, is_zeroeq=is_zeroeq, is_oneeq=is_oneeq, is_twoeq=is_twoeq)
      close(file_unit)
      self%is_loaded = .true.
   else
      if (verbose_) write(stderr, "(A)") 'warning: file "'//trim(adjustl(filename))//'" not found!'
      self%is_loaded = .false.
   endif
   endsubroutine load_file
endmodule xview_file_esz_object
