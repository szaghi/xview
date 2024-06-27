!< xview.
program xview
!< xview.
USE, intrinsic :: iso_fortran_env, only : stderr => error_unit

use penf
use stringifor

use xview_file_grd_object
use xview_file_icc_object
use xview_file_rst_object
use xview_file_esz_object
use xview_file_tec_object
use xview_file_vtk_object
use xview_ui_object

implicit none

type(ui_object) :: ui !< User Interface.

call ui%get_options
if (ui%cli%run_command('postprocess')) then
   if (ui%is_patch) then
      call postprocess(ui=ui, patch=ui%patch)
   else
      call postprocess(ui=ui)
   endif
endif

contains
   subroutine get_filenames(ui, filenames_grd, filenames_icc, filenames_rst, filenames_esz)
   type(ui_object),           intent(in)    :: ui                 !< User Interface.
   type(string), allocatable, intent(inout) :: filenames_grd(:)   !< List of grd files.
   type(string), allocatable, intent(inout) :: filenames_icc(:)   !< List of icc files.
   type(string), allocatable, intent(inout) :: filenames_rst(:)   !< List of rst files.
   type(string), allocatable, intent(inout) :: filenames_esz(:,:) !< List of rst files.
   type(string)                             :: buffer             !< Buffer string.
   integer(I4P)                             :: files_number       !< Number of file.
   integer(I4P)                             :: blocks_number      !< Number of blocks.
   type(string), allocatable                :: tokens(:)          !< String tokens buffer.
   integer(I4P)                             :: Nt                 !< Time steps number.
   integer(I4P)                             :: f, t, t_old, b     !< Counter.

   if (ui%is_extsubzone) then
      if (trim(ui%filename_rst)/=OPT_UNSET)  then
         if (ui%do_glob) then
            ! populate filenames_esz
            ! call buffer%glob(pattern=trim(ui%ipath)//trim(ui%filename_rst)//'*I*', list=filenames_rst)
            call buffer%glob(pattern=trim(ui%ipath)//trim(ui%filename_rst)//'*', list=filenames_rst)
            if (allocated(filenames_rst)) then
               files_number = size(filenames_rst, dim=1)
               ! remove ipath from file names
               do f=1, files_number
                  filenames_rst(f) = filenames_rst(f)%replace(old=trim(ui%ipath), new='')
               enddo
               ! count time steps
               Nt = 1
               call filenames_rst(1)%split(tokens=tokens,sep='_') ; buffer=tokens(2) ; call buffer%split(tokens=tokens,sep='_')
               t_old = tokens(1)%to_number(kind=1_I4P)
               do f=2, files_number
                  call filenames_rst(f)%split(tokens=tokens,sep='_') ; buffer=tokens(2) ; call buffer%split(tokens=tokens,sep='_')
                  t = tokens(1)%to_number(kind=1_I4P)
                  if (t/=t_old) then
                     t_old = t
                     Nt = Nt + 1_I4P
                  endif
               enddo
               blocks_number = files_number / Nt
               ! store file names as blocks_number * Nt matrix
               allocate(filenames_esz(0:blocks_number,1:Nt))
               f = 0_I4P
               do t=1, Nt
                  do b=1, blocks_number
                     f = f + 1_I4P
                     if (b==1) then ! store time step number
                        call filenames_rst(f)%split(tokens=tokens,sep='_');buffer=tokens(2);call buffer%split(tokens=tokens,sep='_')
                        filenames_esz(0,t) = tokens(1)
                     endif
                     filenames_esz(b,t) = filenames_rst(f)
                  enddo
               enddo
               deallocate(filenames_rst)
            endif
         else
            ! fill filenames_esz
            allocate(filenames_esz(0:1,1:1))
            filenames_esz(1,1) = trim(adjustl(ui%filename_rst))
            call filenames_esz(1,1)%split(tokens=tokens,sep='_');buffer=tokens(2);call buffer%split(tokens=tokens,sep='_')
            filenames_esz(0,1) = tokens(1)
         endif
      else
         write(stderr, '(A)') 'error: extracted subzone file (base) name(s) not provided'
      endif
   else
      if (trim(ui%filename_grd)/=OPT_UNSET)  then
         if (ui%do_glob) then
            ! populate filenames_grd, filenames_icc, filenames_rst
            call buffer%glob(pattern=trim(ui%ipath)//trim(ui%filename_grd)//'.*'//trim(adjustl(ui%grid_level))//'.grd.p???', &
                             list=filenames_grd)
            if (allocated(filenames_grd)) then
               files_number = size(filenames_grd, dim=1)
               do f=1, files_number
                  filenames_grd(f) = filenames_grd(f)%replace(old=trim(ui%ipath), new='')
               enddo
               if (trim(ui%filename_icc)/=OPT_UNSET)  then
                  call buffer%glob(pattern=trim(ui%ipath)//trim(ui%filename_icc)//'.*'//trim(adjustl(ui%grid_level))//&
                                   '.p??? | grep -vi .grd', list=filenames_icc)
                  if (allocated(filenames_icc)) then
                     if (size(filenames_icc, dim=1) /= files_number) then
                        write(stderr, '(A)') 'error: number of icc files different from grd files one!'
                        write(stderr, '(A)') '  number of icc files: '//trim(str(size(filenames_icc, dim=1)))
                        write(stderr, '(A)') '  number of grd files: '//trim(str(files_number))
                        stop
                     endif
                     do f=1, files_number
                        filenames_icc(f) = filenames_icc(f)%replace(old=trim(ui%ipath), new='')
                     enddo
                  endif
               else
                  allocate(filenames_icc(1:files_number))
                  do f=1, files_number
                     filenames_icc(f) = OPT_UNSET
                  enddo
               endif
               if (trim(ui%filename_rst)/=OPT_UNSET)  then
                  call buffer%glob(pattern=trim(ui%ipath)//trim(ui%filename_rst)//'_*'//trim(adjustl(ui%grid_level))//'.p???',&
                                   list=filenames_rst)
                  if (allocated(filenames_rst)) then
                     if (size(filenames_rst, dim=1) /= files_number) then
                        write(stderr, '(A)') 'error: number of rst files different from grd files one!'
                        write(stderr, '(A)') '  number of rst files: '//trim(str(size(filenames_rst, dim=1)))
                        write(stderr, '(A)') '  number of grd files: '//trim(str(files_number))
                        stop
                     endif
                     do f=1, files_number
                        filenames_rst(f) = filenames_rst(f)%replace(old=trim(ui%ipath), new='')
                     enddo
                  endif
               else
                  allocate(filenames_rst(1:files_number))
                  do f=1, files_number
                     filenames_rst(f) = OPT_UNSET
                  enddo
               endif
            endif
         else
            ! fill filenames_grd, filenames_icc, filenames_rst
            allocate(filenames_grd(1))
            allocate(filenames_icc(1))
            allocate(filenames_rst(1))
            filenames_grd(1) = trim(adjustl(ui%filename_grd))
            filenames_icc(1) = trim(adjustl(ui%filename_icc))
            filenames_rst(1) = trim(adjustl(ui%filename_rst))
         endif
      else
         write(stderr, '(A)') 'error: at leaset grid file (base) name(s) must be provided'
      endif
   endif
   endsubroutine

   subroutine postprocess(ui, patch)
   !< Post process single or glob tuples of Grd/Icc/rSt (gis) input files.
   type(ui_object), intent(in)           :: ui                 !< User Interface.
   integer(I4P),    intent(in), optional :: patch              !< Patch boundary conditions.
   type(string), allocatable             :: filenames_grd(:)   !< List of grd files.
   type(string), allocatable             :: filenames_icc(:)   !< List of icc files.
   type(string), allocatable             :: filenames_rst(:)   !< List of rst files.
   type(string), allocatable             :: filenames_esz(:,:) !< List of esz files.
   type(string), allocatable             :: files_blocks(:,:)  !< Files blocks grouping.
   type(file_vtk_object)                 :: file_vtk           !< File VTK.
   integer(I4P)                          :: blocks_number      !< Number of blocks.
   integer(I4P)                          :: Nt                 !< Time steps number.
   integer(I4P)                          :: files_number       !< Number of file.
   integer(I4P)                          :: myrank             !< Current process ID.
   integer(I4P)                          :: f, t, b            !< Counter.

   call get_filenames(ui=ui, filenames_grd=filenames_grd, filenames_icc=filenames_icc, filenames_rst=filenames_rst,&
                      filenames_esz=filenames_esz)
   if (ui%is_extsubzone) then
      blocks_number = ubound(filenames_esz, dim=1)
      Nt            = ubound(filenames_esz, dim=2)
      allocate(files_blocks(3,blocks_number))
      do t=1, Nt
         do b=1, blocks_number
            if (ui%verbose) print '(A)', ' Postprocess files "'//trim(ui%ipath)//filenames_esz(b,t)//'"'
            call postprocess_esz(ui           = ui,                                                              &
                                 myrank       = myrank,                                                          &
                                 ipath        = trim(adjustl(ui%ipath)),                                         &
                                 filename_esz = filenames_esz(b,t)%chars(),                                      &
                                 opath        = trim(adjustl(ui%opath)),                                         &
                                 basename_out = trim(adjustl(ui%basename_out)), &
                                 files_blocks = files_blocks(:,b))
         enddo
         if (ui%is_vtk) then
            call file_vtk%save_blocks_file_vtm(file_name=trim(adjustl(ui%opath))//                                               &
                                                         trim(adjustl(ui%basename_out))//'-'//filenames_esz(0,t)%chars()//'.vtm',&
                                               files_blocks=files_blocks)
         endif
      enddo
   else
      myrank = ui%myrank
      files_number = size(filenames_grd, dim=1)
      allocate(files_blocks(3,files_number))
      do f=1, files_number
         if (ui%do_glob) myrank = f-1
         if (ui%verbose) print '(A)', ' Postprocess files "'//trim(ui%ipath)//filenames_grd(f)//', '// &
                                                              trim(ui%ipath)//filenames_icc(f)//', '// &
                                                              trim(ui%ipath)//filenames_rst(f)//'"'
         call postprocess_gis(ui           = ui,                             &
                              myrank       = myrank,                         &
                              ipath        = trim(adjustl(ui%ipath)),        &
                              filename_grd = filenames_grd(f)%chars(),       &
                              filename_icc = filenames_icc(f)%chars(),       &
                              filename_rst = filenames_rst(f)%chars(),       &
                              opath        = trim(adjustl(ui%opath)),        &
                              basename_out = trim(adjustl(ui%basename_out)), &
                              files_blocks = files_blocks(:,f),              &
                              patch=patch)
      enddo
      if (ui%is_vtk) then
         if (ui%file_procinput%is_loaded.and.ui%file_blksmap%is_loaded) &
            call ui%file_blksmap%get_files_blocks(files_blocks=files_blocks)
         call file_vtk%save_blocks_file_vtm(file_name=trim(adjustl(ui%opath))//trim(adjustl(ui%basename_out))//'.vtm',&
                                            files_blocks=files_blocks)
      endif
   endif
   endsubroutine postprocess

   subroutine postprocess_esz(ui,myrank,ipath,filename_esz,opath,basename_out,files_blocks)
   !< Post process single Extracted SubZone (esz) input file.
   type(ui_object), intent(in)    :: ui              !< User Interface.
   integer(I4P),    intent(in)    :: myrank          !< Input files process rank in proc.input.
   character(*),    intent(in)    :: ipath           !< Path to input files.
   character(*),    intent(in)    :: filename_esz    !< File name of file sol (subzone solution).
   character(*),    intent(in)    :: opath           !< Path to output files.
   character(*),    intent(in)    :: basename_out    !< Base name of output files.
   type(string),    intent(inout) :: files_blocks(3) !< Files blocks grouping.
   integer(I4P)                   :: blocks_map(2)   !< Processors/Blocks maps.
   type(string)                   :: buffer          !< String buffer.
   type(string), allocatable      :: tokens(:)       !< String tokens buffer.
   type(file_esz_object)          :: file_esz        !< File esz.
   type(file_vtk_object)          :: file_vtk        !< File VTK.

   if (filename_esz/=OPT_UNSET) call file_esz%load_file(filename=ipath//filename_esz, verbose=ui%verbose, &
                                                        is_level_set=ui%is_level_set,                     &
                                                        is_zeroeq=ui%is_zeroeq,                           &
                                                        is_oneeq=ui%is_oneeq,                             &
                                                        is_twoeq=ui%is_twoeq)
   if (ui%is_vtk) then
      call file_vtk%initialize(path=trim(adjustl(ui%opath)), file_name=trim(adjustl(basename_out)))
      files_blocks(1) = ''
      files_blocks(2) = ''
      files_blocks(3) = file_esz%filename%basename()
      buffer = filename_esz
      call buffer%split(tokens=tokens, sep='_I')
      blocks_map(1) = 0_I4P
      blocks_map(2) = tokens(2)%to_number(kind=1_I4P)
      if (file_esz%is_loaded) &
         call file_vtk%save_block_file_vtk(is_binary=.not.ui%is_ascii, &
                                           save_aux=.true.,            &
                                           blocks_map=blocks_map,      &
                                           grd=file_esz%block_esz%grd, &
                                           is_cell_centered=.false.,   &
                                           icc=file_esz%block_esz%icc, &
                                           sol=file_esz%block_esz%sol, &
                                           files_blocks=files_blocks)
   endif
   endsubroutine postprocess_esz

   subroutine postprocess_gis(ui,myrank,ipath,filename_grd,filename_icc,filename_rst,opath,basename_out,files_blocks,patch)
   !< Post process single tuple of Grd/Icc/rSt (gis) input files.
   type(ui_object), intent(in)           :: ui              !< User Interface.
   integer(I4P),    intent(in)           :: myrank          !< Input files process rank in proc.input.
   character(*),    intent(in)           :: ipath           !< Path to input files.
   character(*),    intent(in)           :: filename_grd    !< File name of file grd.
   character(*),    intent(in)           :: filename_icc    !< File name of file icc.
   character(*),    intent(in)           :: filename_rst    !< File name of file sol (restart solution).
   character(*),    intent(in)           :: opath           !< Path to output files.
   character(*),    intent(in)           :: basename_out    !< Base name of output files.
   type(string),    intent(inout)        :: files_blocks(3) !< Files blocks grouping.
   integer(I4P),    intent(in), optional :: patch           !< Patch boundary conditions.
   type(file_grd_object)                 :: file_grd        !< File grd.
   type(file_icc_object)                 :: file_icc        !< File icc.
   type(file_rst_object)                 :: file_rst        !< File icc.
   integer(I4P), allocatable             :: blocks_map(:,:) !< Processors/Blocks maps.
   type(file_tec_object)                 :: file_tec        !< File Tecplot.
   type(file_vtk_object)                 :: file_vtk        !< File VTK.
   type(string)                          :: basename        !< Current file basename.
   integer(I4P)                          :: b               !< Counter.

   if (filename_grd/=OPT_UNSET) call file_grd%load_file(filename=ipath//filename_grd, &
                                                        is_metrics_to_compute=ui%do_compute_aux, verbose=ui%verbose)
   if (file_grd%is_loaded) then
      if (filename_icc/=OPT_UNSET) call file_icc%load_file(filename=ipath//filename_icc, verbose=ui%verbose)
      ! compute patches extents if patch is queried or whole blocks extensts
      do b=1, file_grd%blocks_number
         if (file_icc%is_loaded.and.present(patch)) then
            ! specific BC has been queried
            call file_grd%blocks(b)%compute_patches_extents(tcc=file_icc%blocks(b)%tcc, patch=patch)
         else
            ! whole blocks extents have been queried
            call file_grd%blocks(b)%compute_patches_extents
         endif
      enddo
      if (file_icc%is_loaded) then
         if (ui%do_compute_aux) then
            do b=1, file_grd%blocks_number
               call file_grd%blocks(b)%correct_metrics_bc(tcc=file_icc%blocks(b)%tcc)
            enddo
          endif
      endif
      if (filename_rst/=OPT_UNSET) call file_rst%load_file(file_grd=file_grd,                                &
                                                           filename=ipath//filename_rst, verbose=ui%verbose, &
                                                           is_level_set=ui%is_level_set,                     &
                                                           is_zeroeq=ui%is_zeroeq,                           &
                                                           is_oneeq=ui%is_oneeq,                             &
                                                           is_twoeq=ui%is_twoeq,                             &
                                                           is_cell_centered=ui%is_cell_centered,             &
                                                           is_aux_to_compute=ui%do_compute_aux)
      if (file_icc%is_loaded.and.file_rst%is_loaded.and.ui%do_compute_aux.and.present(patch)) then
         do b=1, file_grd%blocks_number
            call file_rst%blocks(b)%compute_loads(grd=file_grd%blocks(b), icc=file_icc%blocks(b), patch=patch, RE=ui%RE, &
                                                  rFR2=ui%rFR2, zfs=ui%zfs)
            call file_rst%blocks(b)%compute_yplus(grd=file_grd%blocks(b), icc=file_icc%blocks(b), patch=patch, RE=ui%RE)
            call file_rst%blocks(b)%compute_tau(  grd=file_grd%blocks(b), icc=file_icc%blocks(b), patch=patch, RE=ui%RE)
         enddo
      endif
      if (ui%is_vtk) then
         call file_vtk%initialize(path=trim(adjustl(ui%opath)), file_name=trim(adjustl(basename_out)))
         files_blocks(1) = ''
         files_blocks(2) = ''
         basename = file_grd%filename%basename()
         basename = basename%replace(old='cc.', new='')
         basename = basename%replace(old='grd.', new='')
         files_blocks(3) = basename
         blocks_map = ui%file_procinput%processor_map(processor=myrank, blocks_number=file_grd%blocks_number)
         if (file_icc%is_loaded) then
            if (file_rst%is_loaded) then
               ! save all grd, icc and sol
               do b=1, file_grd%blocks_number
                  call file_vtk%save_block_file_vtk(is_binary=.not.ui%is_ascii,           &
                                                    save_aux=ui%do_compute_aux,           &
                                                    blocks_map=blocks_map(:,b),           &
                                                    grd=file_grd%blocks(b),               &
                                                    is_cell_centered=ui%is_cell_centered, &
                                                    icc=file_icc%blocks(b),               &
                                                    sol=file_rst%blocks(b),               &
                                                    patch=patch,                          &
                                                    files_blocks=files_blocks)
               enddo
            else
               ! save grd and icc
               do b=1, file_grd%blocks_number
                  call file_vtk%save_block_file_vtk(is_binary=.not.ui%is_ascii,           &
                                                    save_aux=ui%do_compute_aux,           &
                                                    blocks_map=blocks_map(:,b),           &
                                                    grd=file_grd%blocks(b),               &
                                                    is_cell_centered=ui%is_cell_centered, &
                                                    icc=file_icc%blocks(b),               &
                                                    patch=patch,                          &
                                                    files_blocks=files_blocks)
               enddo
            endif
         else
            ! save only grd without icc and sol
            do b=1, file_grd%blocks_number
               call file_vtk%save_block_file_vtk(is_binary=.not.ui%is_ascii,           &
                                                 save_aux=ui%do_compute_aux,           &
                                                 blocks_map=blocks_map(:,b),           &
                                                 grd=file_grd%blocks(b),               &
                                                 is_cell_centered=ui%is_cell_centered, &
                                                 files_blocks=files_blocks)
            enddo
         endif
      endif
   else
      write(stderr, '(A)') 'error: at leaset grid files must be loaded'
   endif
   endsubroutine postprocess_gis
endprogram xview
