!< xview, mb.par class definition.
module xview_file_mbpar_object
!< xview, mb.par class definition.
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf
use stringifor

use xview_file_object

implicit none
private
public :: file_mbpar_object

type, extends(file_object) :: file_mbpar_object
   !< mb.par class definition.
   type(string) :: basename_grd     !< Basename of grid files.
   type(string) :: basename_icc     !< Basename of topological files.
   type(string) :: basename_rst     !< Basename of solution restart files.
   type(string) :: basename_ini     !< Basename of initial conditions files.
   type(string) :: turbulence_model !< Turbulence model.
   real(R8P)    :: RE=0._R8P        !< Reynolds number.
   real(R8P)    :: FR=0._R8P        !< Froude number.
   real(R8P)    :: WE=0._R8P        !< Weber number.
   contains
      ! public methods
      procedure, pass(self) :: description !< Return pretty-printed object description.
      procedure, pass(self) :: destroy     !< Destroy dynamic memory.
      procedure, pass(self) :: load_file   !< Load mb.par file data.
endtype file_mbpar_object

contains
   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(file_mbpar_object), intent(in) :: self             !< File data.
   character(len=:), allocatable        :: desc             !< Description.
   character(len=1), parameter          :: NL=new_line('a') !< New line character.

   if (self%is_loaded) then
      desc =       self%filename//' parsed data:'//NL
      desc = desc//'  Basename of grid files:               '//self%basename_grd    //NL
      desc = desc//'  Basename of topological files:        '//self%basename_icc    //NL
      desc = desc//'  Basename of solution restart files:   '//self%basename_rst    //NL
      desc = desc//'  Basename of initial conditions files: '//self%basename_ini    //NL
      desc = desc//'  Turbulence model:                     '//self%turbulence_model//NL
      desc = desc//'  Reynolds number:                      '//str(self%RE)         //NL
      desc = desc//'  Froude number:                        '//str(self%FR)         //NL
      desc = desc//'  Weber number:                         '//str(self%WE)
   else
      desc = 'warning: file "'//self%filename//'" not loaded!'
   endif
   endfunction description

   elemental subroutine destroy(self)
   !< Destroy dynamic memory.
   class(file_mbpar_object), intent(inout) :: self !< File data.

   call self%file_object%destroy
   call self%basename_grd%free
   call self%basename_icc%free
   call self%basename_rst%free
   call self%basename_ini%free
   call self%turbulence_model%free
   self%RE=0._R8P
   self%FR=0._R8P
   self%WE=0._R8P
   endsubroutine destroy

   subroutine load_file(self, path, filename, verbose)
   !< Load mb.par file data.
   class(file_mbpar_object), intent(inout)        :: self          !< File data.
   character(*),             intent(in), optional :: path          !< Path to mb.par.
   character(*),             intent(in), optional :: filename      !< File name of mb.par, optionally customizable.
   logical,                  intent(in), optional :: verbose       !< Activate verbose mode.
   character(:), allocatable                      :: path_         !< Path to mb.par, local variable.
   character(:), allocatable                      :: filename_     !< File name of mb.par, local varibale
   logical                                        :: verbose_      !< Activate verbose mode, local variable.
   type(string)                                   :: mbpar_str     !< mb.par parsing string.
   type(string), allocatable                      :: mbpar_rows(:) !< mb.par rows.
   integer(I4P)                                   :: gln           !< Grid levels number.

   path_     = ''       ; if (present(path    )) path_     = trim(adjustl(path))
   filename_ = 'mb.par' ; if (present(filename)) filename_ = trim(adjustl(filename))
   verbose_  = .false.  ; if (present(verbose))  verbose_  = verbose
   self%filename = path_//filename_

   self%is_loaded = .false.
   if (self%is_file_present()) then
      call mbpar_str%read_file(file=path_//filename)             ! read mb.par file as a single stream
      call mbpar_str%split(tokens=mbpar_rows, sep=new_line('a')) ! get mb.par file rows
      ! parse mb.par rows
      call get_token(row=mbpar_rows(1       ), sep=' ', n=1, var_s=self%basename_grd    ) ! grid files basename
      call get_token(row=mbpar_rows(2       ), sep=' ', n=1, var_s=self%basename_ini    ) ! initial conditions files basename
      call get_token(row=mbpar_rows(3       ), sep=' ', n=1, var_s=self%basename_icc    ) ! topological files basename
      call get_token(row=mbpar_rows(6       ), sep=' ', n=1, var_s=self%basename_rst    ) ! solution restart files basename
      call get_token(row=mbpar_rows(9       ), sep=' ', n=1, var_i=gln                  ) ! grid levels number
      call get_token(row=mbpar_rows(9+gln+13), sep=' ', n=1, var_s=self%turbulence_model) ! turbulence model
      call get_token(row=mbpar_rows(9+gln+10), sep=' ', n=1, var_r=self%RE              ) ! Reynolds number
      call get_token(row=mbpar_rows(9+gln+11), sep=' ', n=1, var_r=self%FR              ) ! Froude number
      call get_token(row=mbpar_rows(9+gln+11), sep=' ', n=2, var_r=self%WE              ) ! Weber number
      ! forward mb.par relative path
      self%basename_grd = path_//self%basename_grd
      self%basename_ini = path_//self%basename_ini
      self%basename_icc = path_//self%basename_icc
      self%basename_rst = path_//self%basename_rst
      ! transform turbulence model string to upper case
      self%turbulence_model = self%turbulence_model%upper()
      self%is_loaded = .true.
   else
      if (verbose_) write(stderr, "(A)") 'warning: file "'//self%filename//'" not found!'
   endif
   endsubroutine load_file

   ! non TBP methods
   subroutine get_token(row, sep, n, var_s, var_i, var_r)
   !< Get token in given row by column (n) and columns separator.
   type(string), intent(in)            :: row        !< Row to be parsed.
   character(*), intent(in)            :: sep        !< Columns seperator.
   integer(I4P), intent(in)            :: n          !< Column number.
   type(string), intent(out), optional :: var_s      !< Outuput variable, string type.
   integer(I4P), intent(out), optional :: var_i      !< Outuput variable, integer type.
   real(R8P),    intent(out), optional :: var_r      !< Outuput variable, real type.
   type(string), allocatable           :: columns(:) !< Columns.

   call row%split(tokens=columns, sep=sep)
   if     (present(var_s)) then
      var_s = columns(n)
   elseif (present(var_i)) then
      var_i = columns(n)%to_number(kind=1_I4P)
   elseif (present(var_r)) then
      var_r = columns(n)%to_number(kind=1._R8P)
   endif
   endsubroutine get_token
endmodule xview_file_mbpar_object
