!< xview, UI class definition.
module xview_ui_object
!< xview, UI class definition.

use flap
use penf

implicit none
private
public :: ui_object

integer(I4P), parameter :: MAX_CHAR_LENGTH=999 !< Maximum length of strings.

type :: ui_object
   !< UI class definition.
   type(command_line_interface) :: cli                !< Command line interface handler.
   character(MAX_CHAR_LENGTH)   :: ipath              !< Path to mb.par file.
   character(MAX_CHAR_LENGTH)   :: procinput_filename !< proc.input file name.
   character(MAX_CHAR_LENGTH)   :: mbpar_filename     !< mb.par file name.
   contains
      ! public methods
      procedure, pass(self) :: parse_cli !< Parse command line interface.
      ! private methods
      procedure, pass(self), private :: set_cli !< Set command line interface.
endtype ui_object

contains
   ! public methods
   subroutine parse_cli(self)
   !< Parse command line interface.
   class(ui_object), intent(inout) :: self  !< User inteface.
   integer(I4P)                    :: error !< Error trapping flag.

   call self%set_cli
   call self%cli%parse(error=error)
   if (error/=0) stop

   call self%cli%get(switch='--ipath',     val=self%ipath,              error=error) ; if (error/=0) stop
   call self%cli%get(switch='--procinput', val=self%procinput_filename, error=error) ; if (error/=0) stop
   call self%cli%get(switch='--mbpar',     val=self%mbpar_filename,     error=error) ; if (error/=0) stop
   endsubroutine parse_cli

   ! private methods
   subroutine set_cli(self)
   !< Set command line interface.
   class(ui_object), intent(inout) :: self  !< User inteface.
   integer(I4P)                    :: error !< Error trapping flag.

   call self%cli%init(progname    = 'xview',                                                               &
                      version     = 'v0.0.1',                                                              &
                      authors     = 'Stefano Zaghi',                                                       &
                      license     = 'GPL v3',                                                              &
                      help        = 'Usage: ',                                                             &
                      description = 'xview, post-processing and analysis tool for Xall/Xnavis CFD solver', &
                      examples    = ["xview --mbpar mb.par"],                                              &
                      epilog      = new_line('a')//"all done")

   call self%cli%add(switch='--ipath',           &
                     switch_ab='-ipath',         &
                     help='path to input files', &
                     required=.false.,           &
                     act='store',                &
                     def='./',                   &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--procinput',          &
                     switch_ab='-procinput',        &
                     help='Xnavis input file name', &
                     required=.false.,              &
                     act='store',                   &
                     def='proc.input',              &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--mbpar',              &
                     switch_ab='-mbpar',            &
                     help='Xnavis input file name', &
                     required=.false.,              &
                     act='store',                   &
                     def='mb.par',                  &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--grd',                &
                     switch_ab='-g',                &
                     help='basename of grid files', &
                     required=.false.,              &
                     act='store',                   &
                     def='cc',                      &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--icc',                    &
                     switch_ab='-i',                    &
                     help='basename of topology files', &
                     required=.false.,                  &
                     act='store',                       &
                     def='cc',                          &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--out',                  &
                     switch_ab='-o',                  &
                     help='basename of output files', &
                     required=.false.,                &
                     act='store',                     &
                     def='output',                    &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--ascii',                &
                     help='write ascii output files', &
                     required=.false.,                &
                     act='store_true',                &
                     def='.false.',                   &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--tec',                    &
                     help='write output Tecplot files', &
                     required=.false.,                  &
                     act='store',                       &
                     def='no',                          &
                     choices='yes,no',                  &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--vtk',                &
                     help='write output VTK files', &
                     required=.false.,              &
                     act='store',                   &
                     def='yes',                     &
                     choices='yes,no',              &
                     error=error)
   if (error/=0) stop
   endsubroutine set_cli
endmodule xview_ui_object
