!< xview.
program xview
!< xview.
USE, intrinsic :: iso_fortran_env, only : stderr => error_unit

use penf
use stringifor

use xview_file_grd_object
use xview_file_icc_object
use xview_file_rst_object
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
   subroutine postprocess(ui, patch)
   !< Post process single or glob tuples of Grd/Icc/rSt (gis) input files.
   type(ui_object), intent(in)           :: ui                !< User Interface.
   integer(I4P),    intent(in), optional :: patch             !< Patch boundary conditions.
   type(string), allocatable             :: filenames_grd(:)  !< List of grd files.
   type(string), allocatable             :: filenames_icc(:)  !< List of icc files.
   type(string), allocatable             :: filenames_rst(:)  !< List of rst files.
   type(string), allocatable             :: files_blocks(:,:) !< Files blocks grouping.
   type(file_vtk_object)                 :: file_vtk          !< File VTK.
   type(string)                          :: buffer            !< Buffer string.
   integer(I4P)                          :: files_number      !< Number of file.
   integer(I4P)                          :: f                 !< Counter.

   if (ui%do_glob) then
      if (trim(ui%filename_grd)/=OPT_UNSET)  then
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
            allocate(files_blocks(3,files_number))
            do f=1, files_number
               if (ui%verbose) print '(A)', ' Postprocess files "'//trim(ui%ipath)//filenames_grd(f)//', '// &
                                                                    trim(ui%ipath)//filenames_icc(f)//', '// &
                                                                    trim(ui%ipath)//filenames_rst(f)//'"'
               call postprocess_gis(ui           = ui,                             &
                                    myrank       = f-1,                            &
                                    ipath        = trim(adjustl(ui%ipath)),        &
                                    filename_grd = filenames_grd(f)%chars(),       &
                                    filename_icc = filenames_icc(f)%chars(),       &
                                    filename_rst = filenames_rst(f)%chars(),       &
                                    opath        = trim(adjustl(ui%opath)),        &
                                    basename_out = trim(adjustl(ui%basename_out)), &
                                    files_blocks = files_blocks(:,f),              &
                                    patch=patch)
            enddo
         else
            write(stderr, '(A)') 'error: at leaset grid files must be loaded'
         endif
      endif
   else
      allocate(files_blocks(3,1))
      call postprocess_gis(ui           = ui,                             &
                           myrank       = ui%myrank,                      &
                           ipath        = trim(adjustl(ui%ipath)),        &
                           filename_grd = trim(adjustl(ui%filename_grd)), &
                           filename_icc = trim(adjustl(ui%filename_icc)), &
                           filename_rst = trim(adjustl(ui%filename_rst)), &
                           opath        = trim(adjustl(ui%opath)),        &
                           basename_out = trim(adjustl(ui%basename_out)), &
                           files_blocks = files_blocks(:,1),              &
                           patch=patch)
   endif
   if (ui%is_vtk) then
      if (ui%file_procinput%is_loaded.and.ui%file_blksmap%is_loaded) &
         call ui%file_blksmap%get_files_blocks(files_blocks=files_blocks)
      call file_vtk%save_blocks_file_vtm(file_name=trim(adjustl(ui%opath))//trim(adjustl(ui%basename_out))//'.vtm',&
                                         files_blocks=files_blocks)
   endif
   endsubroutine postprocess

   subroutine postprocess_gis(ui,myrank,ipath,filename_grd,filename_icc,filename_rst,opath,basename_out,files_blocks,patch)
   !< Post process single tuple of Grd/Icc/rSt (gis) input files.
   type(ui_object),       intent(in)           :: ui              !< User Interface.
   integer(I4P),          intent(in)           :: myrank          !< Input files process rank in proc.input.
   character(*),          intent(in)           :: ipath           !< Path to input files.
   character(*),          intent(in)           :: filename_grd    !< File name of file grd.
   character(*),          intent(in)           :: filename_icc    !< File name of file icc.
   character(*),          intent(in)           :: filename_rst    !< File name of file sol.
   character(*),          intent(in)           :: opath           !< Path to output files.
   character(*),          intent(in)           :: basename_out    !< Base name of output files.
   type(string),          intent(inout)        :: files_blocks(3) !< Files blocks grouping.
   integer(I4P),          intent(in), optional :: patch           !< Patch boundary conditions.
   type(file_grd_object)                       :: file_grd        !< File grd.
   type(file_icc_object)                       :: file_icc        !< File icc.
   type(file_rst_object)                       :: file_rst        !< File icc.
   integer(I4P), allocatable                   :: blocks_map(:,:) !< Processors/Blocks maps.
   type(file_tec_object)                       :: file_tec        !< File Tecplot.
   type(file_vtk_object)                       :: file_vtk        !< File VTK.
   type(string)                                :: basename        !< Current file basename.
   integer(I4P)                                :: b               !< Counter.

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
