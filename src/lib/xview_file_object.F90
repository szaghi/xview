!< xview, file (base) class definition.
module xview_file_object
!< xview, file (base) class definition.
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf
use stringifor

implicit none
private
public :: file_object

type :: file_object
   !< **File base** class: base object to be extended by other file classes.
   type(string) :: filename          !< File name.
   logical      :: is_loaded=.false. !< Flag for checking if the file has been correctly loaded.
   contains
      ! plublic methods
      procedure, pass(self) :: destroy         !< Destroy dynamic memory.
      procedure, pass(self) :: is_file_present !< Inquire if the file path is valid.
endtype file_object

contains
   ! plublic methods
   elemental subroutine destroy(self)
   !< Destroy dynamic memory.
   class(file_object), intent(inout) :: self !< File data.

   call self%filename%free
   self%is_loaded = .false.
   endsubroutine destroy

   function is_file_present(self) result(is_present)
   !< Inquire if the file path is valid.
   class(file_object), intent(in) :: self       !< File data.
   logical                        :: is_present !< Inquiring result.

   is_present = .false.
   if (self%filename%is_allocated()) inquire(file=trim(adjustl(self%filename%chars())), exist=is_present)
   endfunction is_file_present
endmodule xview_file_object
