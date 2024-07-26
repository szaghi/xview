!< xview, UI class definition.
module xview_ui_object
!< xview, UI class definition.
use flap
use penf

use xview_file_mbpar_object
use xview_file_procinput_object
use xview_file_blksmap_object

implicit none
private
public :: ui_object
public :: OPT_UNSET

integer(I4P), parameter :: MAX_CHAR_LENGTH=999 !< Maximum length of strings.
character(6), parameter :: OPT_UNSET='unset!'  !< Option unset flag.

type :: ui_object
   !< UI class definition.
   type(command_line_interface) :: cli                      !< Command line interface handler.
   type(file_procinput_object)  :: file_procinput           !< proc.input.
   type(file_blksmap_object)    :: file_blksmap             !< File block groups map.
   integer(I4P)                 :: myrank=-1                !< Input files process number in proc.input data.
   type(file_mbpar_object)      :: file_mbpar               !< mb.par.
   character(MAX_CHAR_LENGTH)   :: mbpar_file_name          !< File name of mb.par.
   character(MAX_CHAR_LENGTH)   :: procinput_file_name      !< File name of proc.input.
   character(MAX_CHAR_LENGTH)   :: blksmap_file_name        !< File name of blocks groups map ini file.
   character(MAX_CHAR_LENGTH)   :: ipath                    !< Path to input files.
   character(MAX_CHAR_LENGTH)   :: filename_grd             !< File name of file grd.
   character(MAX_CHAR_LENGTH)   :: filename_icc             !< File name of file icc.
   character(MAX_CHAR_LENGTH)   :: filename_rst             !< File name of file sol.
   character(MAX_CHAR_LENGTH)   :: opath                    !< Path to output files.
   character(MAX_CHAR_LENGTH)   :: basename_out             !< Base name of output files.
   character(2)                 :: grid_level               !< Grid level to be postprocessed.
   real(R8P)                    :: RE=0._R8P                !< Reynolds number.
   real(R8P)                    :: FR=-1._R8P               !< Froude number.
   real(R8P)                    :: rFR2=0._R8P              !< 1/(Froude number)^2.
   real(R8P)                    :: zfs=0._R8P               !< Z quote of free surface.
   logical                      :: is_cell_centered=.false. !< Data are cell-centered, otherwise node-centered.
   logical                      :: is_level_set=.false.     !< Solution has level set variable.
   logical                      :: is_dns=.false.           !< Solution has not turbulent module (DNS).
   logical                      :: is_zeroeq=.false.        !< Solution has zero equations turbulent variables.
   logical                      :: is_oneeq=.true.          !< Solution has one  equations turbulent variables.
   logical                      :: is_twoeq=.false.         !< Solution has two  equations turbulent variables.
   logical                      :: is_ascii=.false.         !< Output is ascii formatted.
   logical                      :: is_tec=.false.           !< Output is Tecplot formatted.
   logical                      :: is_vtk=.true.            !< Output is VTK formatted.
   logical                      :: is_patch=.false.         !< Extract patch instead of whole volume.
   integer(I4P)                 :: patch                    !< Patch boundary conditions.
   logical                      :: is_extsubzone=.false.    !< Input files contain extracted subzones instead regular files.
   logical                      :: compute_metrics=.false.  !< Compute metrics.
   logical                      :: compute_lambda2=.false.  !< Compute lamda2 field.
   logical                      :: compute_qfactor=.false.  !< Compute qfactor field.
   logical                      :: compute_helicity=.false. !< Compute helicity field.
   logical                      :: compute_vorticity=.false.!< Compute vorticity field.
   logical                      :: compute_div2LT=.false.   !< Compute double divergence of Lighthill tensor.
   logical                      :: compute_k_ratio=.false.  !< Compute kinetic energy ratio.
   logical                      :: compute_yplus=.false.    !< Compute y+ field.
   logical                      :: compute_tau=.false.      !< Compute tau field.
   logical                      :: compute_div_tau=.false.  !< Compute divergence of tau field.
   logical                      :: compute_loads=.false.    !< Compute loads (forces and torques).
   logical                      :: do_glob=.false.          !< Do glob files search, input file names become base name.
   logical                      :: do_postprocess=.false.   !< Do postprocess for Paraview/Tecplo visualization.
   logical                      :: do_interpolate=.false.   !< Do interpolate on new grids.
   logical                      :: do_aeroacustic=.false.   !< Do aeroacustic propagation analysis.
   logical                      :: verbose=.false.          !< Activate verbose mode.
   contains
      ! public methods
      procedure, pass(self) :: get_options !< Get (user) options.
      ! private methods
      procedure, pass(self), private :: parse_cli !< Parse command line interface.
      procedure, pass(self), private :: set_cli   !< Set command line interface.
endtype ui_object

contains
   ! public methods
   subroutine get_options(self)
   !< Get (user) options.
   class(ui_object), intent(inout) :: self !< User inteface.

   call self%parse_cli
   call self%file_mbpar%load_file(path=self%ipath, filename=self%mbpar_file_name, verbose=self%verbose)
   if (self%file_mbpar%is_loaded) then
      ! override cli with mb.par data
      self%RE  = self%file_mbpar%RE
      self%FR  = self%file_mbpar%FR
      self%zfs = self%file_mbpar%zfs
      self%is_level_set = self%file_mbpar%FR > 0._R8P
      select case(self%file_mbpar%turbulence_model%slice(1,3))
      case('BAL')
         self%is_zeroeq = .true.
      case('SGS')
         self%is_zeroeq = .true.
      case('SPA')
         self%is_oneeq = .true.
      case('LAM')
         self%is_twoeq = .true.
      case('CHA')
         self%is_twoeq = .true.
      case('DES')
         self%is_oneeq = .true.
      case('DDE')
         self%is_oneeq = .true.
      endselect
      if (self%verbose) print '(A)', self%file_mbpar%description()
   endif
   if (self%FR> 0._R8P) self%rFR2 = 1._R8P/(self%FR**2)
   call self%file_procinput%load_file(path=self%ipath, filename=self%procinput_file_name, verbose=self%verbose)
   if (self%blksmap_file_name/=OPT_unset) &
      call self%file_blksmap%load_file(path=self%ipath, filename=self%blksmap_file_name, verbose=self%verbose)
   if (self%file_blksmap%is_loaded) then
      if (self%verbose) print '(A)', self%file_blksmap%description()
   endif
   ! create output dirs if necessary, not portable (only for *nix system)
   if (trim(adjustl(self%opath))/='') call execute_command_line('mkdir -p '//trim(adjustl(self%opath)))
   if (self%is_vtk) call execute_command_line('mkdir -p '//trim(adjustl(self%opath))//trim(adjustl(self%basename_out))//'-vts')
   endsubroutine get_options

   ! private methods
   subroutine parse_cli(self)
   !< Parse command line interface.
   class(ui_object), intent(inout) :: self         !< User inteface.
   integer(I4P)                    :: error        !< Error trapping flag.
   integer(I4P)                    :: turbulent_eq !< Number equations of Turbulent model.
   character(3)                    :: is_tec       !< Flag for setting Tecplot output files.
   character(3)                    :: is_vtk       !< Flag for setting VTK output files.

   call self%set_cli
   call self%cli%parse(error=error) ; if (error/=0) stop

   call self%cli%get(switch='--verbose',            val=self%verbose,             error=error) ; if (error/=0) stop
   call self%cli%get(switch='--ipath',              val=self%ipath,               error=error) ; if (error/=0) stop
   call self%cli%get(switch='--procinput',          val=self%procinput_file_name, error=error) ; if (error/=0) stop
   call self%cli%get(switch='--blocks-map',         val=self%blksmap_file_name,   error=error) ; if (error/=0) stop
   call self%cli%get(switch='--myrank',             val=self%myrank,              error=error) ; if (error/=0) stop
   call self%cli%get(switch='--mbpar',              val=self%mbpar_file_name,     error=error) ; if (error/=0) stop
   call self%cli%get(switch='--grid-level',         val=self%grid_level,          error=error) ; if (error/=0) stop
   call self%cli%get(switch='--glob',               val=self%do_glob,             error=error) ; if (error/=0) stop
   call self%cli%get(switch='--grd',                val=self%filename_grd,        error=error) ; if (error/=0) stop
   call self%cli%get(switch='--icc',                val=self%filename_icc,        error=error) ; if (error/=0) stop
   call self%cli%get(switch='--rst',                val=self%filename_rst,        error=error) ; if (error/=0) stop
   call self%cli%get(switch='--opath',              val=self%opath,               error=error) ; if (error/=0) stop
   call self%cli%get(switch='--out',                val=self%basename_out,        error=error) ; if (error/=0) stop
   call self%cli%get(switch='--RE',                 val=self%RE,                  error=error) ; if (error/=0) stop
   call self%cli%get(switch='--FR',                 val=self%FR,                  error=error) ; if (error/=0) stop
   call self%cli%get(switch='--zfs',                val=self%zfs,                 error=error) ; if (error/=0) stop
   call self%cli%get(switch='--level-set',          val=self%is_level_set,        error=error) ; if (error/=0) stop
   call self%cli%get(switch='--no-turbulent-model', val=self%is_dns,              error=error) ; if (error/=0) stop
   call self%cli%get(switch='--turbulent-eq',       val=turbulent_eq,             error=error) ; if (error/=0) stop
   if (self%cli%run_command('postprocess')) then
      call self%cli%get(group='postprocess',switch='--ascii',    val=self%is_ascii,          error=error) ; if (error/=0) stop
      call self%cli%get(group='postprocess',switch='--tec',      val=is_tec,                 error=error) ; if (error/=0) stop
      call self%cli%get(group='postprocess',switch='--vtk',      val=is_vtk,                 error=error) ; if (error/=0) stop
      call self%cli%get(group='postprocess',switch='--cell',     val=self%is_cell_centered,  error=error) ; if (error/=0) stop
      call self%cli%get(group='postprocess',switch='--metrics',  val=self%compute_metrics,   error=error) ; if (error/=0) stop
      call self%cli%get(group='postprocess',switch='--lambda2',  val=self%compute_lambda2,   error=error) ; if (error/=0) stop
      call self%cli%get(group='postprocess',switch='--qfactor',  val=self%compute_qfactor,   error=error) ; if (error/=0) stop
      call self%cli%get(group='postprocess',switch='--helicity', val=self%compute_helicity,  error=error) ; if (error/=0) stop
      call self%cli%get(group='postprocess',switch='--vorticity',val=self%compute_vorticity, error=error) ; if (error/=0) stop
      call self%cli%get(group='postprocess',switch='--div2LT',   val=self%compute_div2LT,    error=error) ; if (error/=0) stop
      call self%cli%get(group='postprocess',switch='--k-ratio',  val=self%compute_k_ratio,   error=error) ; if (error/=0) stop
      call self%cli%get(group='postprocess',switch='--yplus',    val=self%compute_yplus,     error=error) ; if (error/=0) stop
      call self%cli%get(group='postprocess',switch='--tau',      val=self%compute_tau,       error=error) ; if (error/=0) stop
      call self%cli%get(group='postprocess',switch='--div-tau',  val=self%compute_div_tau,   error=error) ; if (error/=0) stop
      call self%cli%get(group='postprocess',switch='--loads',    val=self%compute_loads,     error=error) ; if (error/=0) stop
      call self%cli%get(group='postprocess',switch='--patch',      val=self%patch,           error=error) ; if (error/=0) stop
      call self%cli%get(group='postprocess',switch='--ext-subzone',val=self%is_extsubzone,   error=error) ; if (error/=0) stop
      if (self%cli%is_passed(group='postprocess', switch='--patch')) self%is_patch=.true.
   endif
   if (self%is_dns) then
      self%is_zeroeq = .false.
      self%is_oneeq = .false.
      self%is_twoeq = .false.
   else
      select case(turbulent_eq)
      case(0)
         self%is_zeroeq = .true.
         self%is_oneeq = .false.
         self%is_twoeq = .false.
      case(1)
         self%is_zeroeq = .false.
         self%is_oneeq = .true.
         self%is_twoeq = .false.
      case(2)
         self%is_zeroeq = .false.
         self%is_oneeq = .false.
         self%is_twoeq = .true.
      endselect
   endif
   self%is_tec = (adjustl(trim(is_tec))=='YES'.or.adjustl(trim(is_tec))=='yes')
   self%is_vtk = (adjustl(trim(is_vtk))=='YES'.or.adjustl(trim(is_vtk))=='yes')
   endsubroutine parse_cli

   subroutine set_cli(self)
   !< Set command line interface.
   class(ui_object), intent(inout) :: self  !< User inteface.
   integer(I4P)                    :: error !< Error trapping flag.

   call self%cli%init(progname    = 'xview',                                                                                   &
                      version     = 'v0.0.1',                                                                                  &
                      authors     = 'Stefano Zaghi',                                                                           &
                      license     = 'GPL v3',                                                                                  &
                      help        = 'Usage: ',                                                                                 &
                      description = 'xview, post-processing and analysis tool for Xall/Xnavis CFD solver',                     &
                      examples    = ["xview postprocess -g grid -i topology -r srestart -o sol_tec --cell --tec yes --vtk no", &
                                     "xview interpolate -g grid -i topology -r srestart -o interpolated --igrd new_grid     ", &
                                     "xview aeroacustic -g grid -i topology -r srestart -o aero_analysis --aero fwh         "],&
                      epilog      = new_line('a')//"all done")

   call self%cli%add(switch='--verbose',           &
                     help='activate verbose mode', &
                     required=.false.,             &
                     act='store_true',             &
                     def='.false.',                &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--ipath',           &
                     switch_ab='-ip',            &
                     help='path to input files', &
                     required=.false.,           &
                     act='store',                &
                     def='',                     &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--procinput',    &
                     help='proc.input input', &
                     required=.false.,        &
                     act='store',             &
                     def='proc.input',        &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--blocks-map',                 &
                     help='blocks groups map ini filename', &
                     required=.false.,                      &
                     act='store',                           &
                     def=OPT_UNSET,                         &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--myrank',                                   &
                     help='input file process number in proc.input data', &
                     required=.false.,                                    &
                     act='store',                                         &
                     def='-1',                                            &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--mbpar',           &
                     help='mb.par Xnavis input', &
                     required=.false.,           &
                     act='store',                &
                     def='mb.par',               &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--grid-level',        &
                     switch_ab='-gl',              &
                     help='grid refinement level', &
                     required=.false.,             &
                     act='store',                  &
                     def='*',                      &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--glob',                                                 &
                     help='do glob files search, input file names become base names', &
                     required=.false.,                                                &
                     act='store_true',                                                &
                     def='.false.',                                                   &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--grd',                &
                     switch_ab='-g',                &
                     help='filename of grid files', &
                     required=.false.,              &
                     act='store',                   &
                     def=OPT_UNSET,                 &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--icc',                    &
                     switch_ab='-i',                    &
                     help='filename of topology files', &
                     required=.false.,                  &
                     act='store',                       &
                     def=OPT_UNSET,                     &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--rst',                            &
                     switch_ab='-r',                            &
                     help='filename of solution restart files', &
                     required=.false.,                          &
                     act='store',                               &
                     def=OPT_UNSET,                             &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--opath',            &
                     switch_ab='-op',             &
                     help='path to output files', &
                     required=.false.,            &
                     act='store',                 &
                     def='',                      &
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

   call self%cli%add(switch='--RE',          &
                     help='Reynolds number', &
                     required=.false.,       &
                     act='store',            &
                     def='1.0',              &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--FR',        &
                     help='Froude number', &
                     required=.false.,     &
                     act='store',          &
                     def='-1.0',           &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--zfs',                 &
                     help='Z quote of free surface', &
                     required=.false.,               &
                     act='store',                    &
                     def='0.0',                      &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--level-set',           &
                     help='solution with level set', &
                     required=.false.,               &
                     act='store_true',               &
                     def='.false.',                  &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--no-turbulent-model',  &
                     help='no turbulent model used', &
                     required=.false.,               &
                     act='store_true',               &
                     def='.false.',                  &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--turbulent-eq',                 &
                     help='number equations turbulent model', &
                     required=.false.,                        &
                     act='store',                             &
                     def='1',                                 &
                     choices='0,1,2',                         &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--compute-aux',                                  &
                     help='compute auxiliary variables (metrics, forces...)', &
                     required=.false.,                                        &
                     act='store_true',                                        &
                     def='.false.',                                           &
                     error=error)
   if (error/=0) stop

   ! commands

   ! postprocess
   call self%cli%add_group(help="usage: ",                                                           &
                           description="postprocess input files for ParaView/TecPlot visualization", &
                           group="postprocess")

   call self%cli%add(switch='--cell',                                            &
                     help='variables other than nodes coord. are cell centered', &
                     required=.false.,                                           &
                     act='store_true',                                           &
                     def='.false.',                                              &
                     group='postprocess',                                        &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--ascii',                &
                     help='write ascii output files', &
                     required=.false.,                &
                     act='store_true',                &
                     def='.false.',                   &
                     group='postprocess',             &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--tec',                    &
                     help='write output Tecplot files', &
                     required=.false.,                  &
                     act='store',                       &
                     def='no',                          &
                     choices='yes,no',                  &
                     group='postprocess',               &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--vtk',                &
                     help='write output VTK files', &
                     required=.false.,              &
                     act='store',                   &
                     def='yes',                     &
                     choices='yes,no',              &
                     group='postprocess',           &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--metrics',          &
                     help='compute mesh metrics', &
                     required=.false.,            &
                     act='store_true',            &
                     def='.false.',               &
                     group='postprocess',         &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--lambda2',           &
                     help='compute lambda2 field', &
                     required=.false.,             &
                     act='store_true',             &
                     def='.false.',                &
                     group='postprocess',          &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--qfactor',           &
                     help='compute qfactor field', &
                     required=.false.,             &
                     act='store_true',             &
                     def='.false.',                &
                     group='postprocess',          &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--helicity',           &
                     help='compute helicity field', &
                     required=.false.,              &
                     act='store_true',              &
                     def='.false.',                 &
                     group='postprocess',           &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--vorticity',           &
                     help='compute vorticity field', &
                     required=.false.,               &
                     act='store_true',               &
                     def='.false.',                  &
                     group='postprocess',            &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--div2LT',                                          &
                     help='compute double divergence of Lighthill Tensor field', &
                     required=.false.,                                           &
                     act='store_true',                                           &
                     def='.false.',                                              &
                     group='postprocess',                                        &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--k-ratio',                        &
                     help='compute kinetic energy ratio field', &
                     required=.false.,                          &
                     act='store_true',                          &
                     def='.false.',                             &
                     group='postprocess',                       &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--patch',               &
                     switch_ab='-p',                 &
                     help='extract specified patch', &
                     required=.false.,               &
                     act='store',                    &
                     def='-999',                     &
                     group='postprocess',            &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--yplus',           &
                     help='compute yplus field', &
                     required=.false.,           &
                     act='store_true',           &
                     def='.false.',              &
                     group='postprocess',        &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--tau',           &
                     help='compute tau field', &
                     required=.false.,         &
                     act='store_true',         &
                     def='.false.',            &
                     group='postprocess',      &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--div-tau',           &
                     help='compute div-tau field', &
                     required=.false.,             &
                     act='store_true',             &
                     def='.false.',                &
                     group='postprocess',          &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--loads',           &
                     help='compute loads field', &
                     required=.false.,           &
                     act='store_true',           &
                     def='.false.',              &
                     group='postprocess',        &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--ext-subzone',                                         &
                     help='postprocess extracted subzones instead of regular files', &
                     required=.false.,                                               &
                     act='store_true',                                               &
                     def='.false.',                                                  &
                     group='postprocess',                                            &
                     error=error)
   if (error/=0) stop

   ! interpolate
   call self%cli%add_group(help="usage: ",                                       &
                           description="interpolate input files on given grids", &
                           group="interpolate")

   call self%cli%add(switch='--igrd',                            &
                     switch_ab='-ig',                            &
                     help='basename of interpolated grid files', &
                     required=.false.,                           &
                     act='store',                                &
                     def='cc',                                   &
                     group='interpolate',                        &
                     error=error)
   if (error/=0) stop

   ! aeroacustic analysis
   call self%cli%add_group(help="usage: ",                     &
                           description="aeroacustic analysis", &
                           group="aeroacustic")

   call self%cli%add(switch='--aero',                      &
                     help='aeroacustic propagation model', &
                     required=.false.,                     &
                     act='store',                          &
                     def='fwh',                            &
                     choices='fwh,fwh-2',                  &
                     group='aeroacustic',                  &
                     error=error)
   if (error/=0) stop
   endsubroutine set_cli
endmodule xview_ui_object
