!< xview, file VTK class definition.
module xview_file_vtk_object
!< xview, file VTK class definition.
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf
use vecfor
use stringifor
use vtk_fortran

use xview_block_grd_object
use xview_block_icc_object
use xview_block_rst_object

implicit none
private
public :: file_vtk_object

type :: file_vtk_object
   !< File VTK class definition.
   character(len=:), allocatable :: path      !< Output path.
   character(len=:), allocatable :: file_name !< File name.
   contains
      ! public methods
      procedure, pass(self) :: destroy              !< Destroy dynamic memory.
      procedure, pass(self) :: initialize           !< Initialize file.
      procedure, pass(self) :: is_file_present      !< Inquire if the file path is valid.
      procedure, pass(self) :: save_block_file_vtk  !< Save one Xnavis-block-data into a VTK file.
      procedure, nopass     :: save_blocks_file_vtm !< Save multi Xnavis-block-data into a VTM file.
endtype file_vtk_object

contains
   elemental subroutine destroy(self)
   !< Destroy dynamic memory.
   class(file_vtk_object), intent(inout) :: self !< File data.

   if (allocated(self%file_name)) deallocate(self%file_name)
   endsubroutine destroy

   elemental subroutine initialize(self, path, file_name)
   !< Initialize file.
   class(file_vtk_object), intent(inout) :: self      !< File data.
   character(*),           intent(in)    :: path      !< Output path.
   character(*),           intent(in)    :: file_name !< File name.

   call self%destroy
   self%path      = trim(adjustl(path))
   self%file_name = trim(adjustl(file_name))
   endsubroutine initialize

   function is_file_present(self) result(is_present)
   !< Inquire if the file path is valid.
   class(file_vtk_object), intent(in) :: self       !< File data.
   logical                            :: is_present !< Inquiring result.

   if (is_present) inquire(file=self%path//self%file_name, exist=is_present)
   endfunction is_file_present

   subroutine save_block_file_vtk(self,is_binary,save_aux,blocks_map,grd,is_cell_centered,files_blocks,&
                                  icc,sol,patch,file_name)
   !< Save one Xnavis-block-data into a VTK file.
   class(file_vtk_object), intent(inout)        :: self              !< File data.
   logical,                intent(in)           :: is_binary         !< Define binary or ascii output.
   logical,                intent(in)           :: save_aux          !< Save auxiliary variables (metrics, forces...).
   integer(I4P),           intent(in)           :: blocks_map(1:)    !< Blocks (processors) map.
   type(block_grd_object), intent(in)           :: grd               !< Grid of block.
   logical,                intent(in)           :: is_cell_centered  !< Define variables at cell centers or nodes.
   type(string),           intent(inout)        :: files_blocks(1:)  !< Files blocks group.
   type(block_icc_object), intent(in), optional :: icc               !< Boundary condititions of block.
   type(block_rst_object), intent(in), optional :: sol               !< Solution of block.
   integer(I4P),           intent(in), optional :: patch             !< Patch boundary conditions.
   character(*),           intent(in), optional :: file_name         !< File name.
   character(len=:), allocatable                :: base_name         !< Base file name.
   character(len=:), allocatable                :: output_format     !< Output format.
   character(len=:), allocatable                :: file_name_        !< File name, local variable.
   type(vtk_file)                               :: output            !< Output VTK file.
   integer(I4P)                                 :: error             !< Error status.
   integer(I4P)                                 :: nodes_number      !< Number of nodes.
   integer(I4P)                                 :: i1,i2,j1,j2,k1,k2 !< Actual extents.
   integer(I4P)                                 :: p                 !< Counter.

   if (.not.grd%is_loaded) error stop 'error: cannot save VTK block file because passed grd block has not been loaded'
   if ((.not.present(file_name)).and.(.not.allocated(self%file_name))) then
      error stop 'error: nor "file_name" neither "self%file_name" have been specified for "file_vtk_object%save_file" method!'
   else
      if (present(file_name)) then
         base_name = trim(adjustl(file_name))
      else
         base_name = trim(adjustl(self%file_name))
      endif
   endif

   if (.not.allocated(grd%patches_extents)) return ! queried a patch bc that is not present into this block

   output_format = 'ascii' ; if (is_binary) output_format = 'raw'

   do p=1, size(grd%patches_extents, dim=1) ! loop over pathces or whole block
      associate(ci1=>grd%patches_extents(p, 1 ), ci2=>grd%patches_extents(p, 2 ), &
                cj1=>grd%patches_extents(p, 3 ), cj2=>grd%patches_extents(p, 4 ), &
                ck1=>grd%patches_extents(p, 5 ), ck2=>grd%patches_extents(p, 6 ), &
                ni1=>grd%patches_extents(p, 7 ), ni2=>grd%patches_extents(p, 8 ), &
                nj1=>grd%patches_extents(p, 9 ), nj2=>grd%patches_extents(p, 10), &
                nk1=>grd%patches_extents(p, 11), nk2=>grd%patches_extents(p, 12))
         if (present(patch)) then
            file_name_ = base_name//'-'//files_blocks(3)//'-patch_'//trim(strz(patch, nz_pad=3))//&
                         '-face_'//trim(str(grd%patches_extents(p, 0), no_sign=.true.))//         &
                         '-blk_'//trim(strz(blocks_map(2), nz_pad=4))//'-grp_'//trim(strz(blocks_map(1), nz_pad=3))
         else
            file_name_ = base_name//'-'//files_blocks(3)//&
                         '-blk_'//trim(strz(blocks_map(2), nz_pad=4))//'-grp_'//trim(strz(blocks_map(1), nz_pad=3))
         endif
         files_blocks(1) = files_blocks(1)//' '//base_name//'-vts/'//file_name_//'.vts'
         files_blocks(2) = files_blocks(2)//' '//file_name_
         nodes_number = (ni2 - ni1 + 1) * (nj2 - nj1 + 1) * (nk2 - nk1 + 1)
         error = output%initialize(format=output_format,                                       &
                                   filename=self%path//base_name//'-vts/'//file_name_//'.vts', &
                                   mesh_topology='StructuredGrid',                             &
                                   nx1=ni1, nx2=ni2, ny1=nj1, ny2=nj2, nz1=nk1, nz2=nk2)
         error = output%xml_writer%write_piece(nx1=ni1, nx2=ni2, ny1=nj1, ny2=nj2, nz1=nk1, nz2=nk2)
         error = output%xml_writer%write_geo(n=nodes_number,                         &
                                             x=grd%nodes(ni1:ni2,nj1:nj2,nk1:nk2)%x, &
                                             y=grd%nodes(ni1:ni2,nj1:nj2,nk1:nk2)%y, &
                                             z=grd%nodes(ni1:ni2,nj1:nj2,nk1:nk2)%z)
         if (is_cell_centered) then
            i1 = ci1 ; i2 = ci2
            j1 = cj1 ; j2 = cj2
            k1 = ck1 ; k2 = ck2
            error = output%xml_writer%write_dataarray(location='cell', action='open')
         else
            i1 = ni1 ; i2 = ni2
            j1 = nj1 ; j2 = nj2
            k1 = nk1 ; k2 = nk2
            error = output%xml_writer%write_dataarray(location='node', action='open')
         endif
         if (present(icc)) then
            if (icc%is_loaded) then
               error = output%xml_writer%write_dataarray(data_name='icc', x=icc%rcc(i1:i2,j1:j2,k1:k2), one_component=.true.)
            endif
         endif
         if (present(sol)) then
            if (sol%is_loaded) then
               error = output%xml_writer%write_dataarray(data_name='velocity', x=sol%momentum(i1:i2,j1:j2,k1:k2)%x, &
                                                                               y=sol%momentum(i1:i2,j1:j2,k1:k2)%y, &
                                                                               z=sol%momentum(i1:i2,j1:j2,k1:k2)%z)
               error = output%xml_writer%write_dataarray(data_name='pressure', x=sol%pressure(i1:i2,j1:j2,k1:k2), &
                                                         one_component=.true.)
               if (sol%has_viscosity()) then
                  error = output%xml_writer%write_dataarray(data_name='viscosity', x=sol%viscosity(i1:i2,j1:j2,k1:k2), &
                                                            one_component=.true.)
               endif
               if (sol%has_turb_viscosity()) then
                  error = output%xml_writer%write_dataarray(data_name='turbulent viscosity', &
                                                            x=sol%turbulent_viscosity(i1:i2,j1:j2,k1:k2), one_component=.true.)
               endif
               if (sol%has_turb_kinetic_energy()) then
                  error = output%xml_writer%write_dataarray(data_name='turbulent kinetic energy',              &
                                                            x=sol%turbulent_kinetic_energy(i1:i2,j1:j2,k1:k2), &
                                                            one_component=.true.)
               endif
               if (sol%has_turb_kinetic_energy_diss()) then
                  error = output%xml_writer%write_dataarray(data_name='turbulent kinetic energy dissipation',              &
                                                            x=sol%turbulent_kinetic_energy_dissipation(i1:i2,j1:j2,k1:k2), &
                                                            one_component=.true.)
               endif
               if (sol%has_level_set()) then
                  error = output%xml_writer%write_dataarray(data_name='level set', x=sol%level_set(i1:i2,j1:j2,k1:k2), &
                                                            one_component=.true.)
                  error = output%xml_writer%write_dataarray(data_name='level set zero', &
                                                            x=sol%level_set_zero(i1:i2,j1:j2,k1:k2), one_component=.true.)
               endif
            endif
         endif
         if (save_aux) then
            ! if (allocated(grd%volume)) &
            ! error=output%xml_writer%write_dataarray(data_name='volume',x=grd%volume(i1:i2,j1:j2,k1:k2),one_component=.true.)
            if (allocated(sol%div2T)) &
            error=output%xml_writer%write_dataarray(data_name='div2T',x=sol%div2T(i1:i2,j1:j2,k1:k2),one_component=.true.)
            if (allocated(sol%lambda2)) &
            error=output%xml_writer%write_dataarray(data_name='lambda2',x=sol%lambda2(i1:i2,j1:j2,k1:k2),one_component=.true.)
            if (allocated(sol%qfactor)) &
            error=output%xml_writer%write_dataarray(data_name='qfactor',x=sol%qfactor(i1:i2,j1:j2,k1:k2),one_component=.true.)
            if (allocated(sol%helicity)) &
            error=output%xml_writer%write_dataarray(data_name='helicity',x=sol%helicity(i1:i2,j1:j2,k1:k2),one_component=.true.)
            if (allocated(sol%vorticity)) &
            error=output%xml_writer%write_dataarray(data_name='vorticity', x=sol%vorticity(i1:i2,j1:j2,k1:k2)%x, &
                                                                           y=sol%vorticity(i1:i2,j1:j2,k1:k2)%y, &
                                                                           z=sol%vorticity(i1:i2,j1:j2,k1:k2)%z)
            if (allocated(sol%force_hydrostatic).and.present(patch)) &
            error=output%xml_writer%write_dataarray(data_name='force_hydrostatic',x=sol%force_hydrostatic(i1:i2,j1:j2,k1:k2)%x, &
                                                                                  y=sol%force_hydrostatic(i1:i2,j1:j2,k1:k2)%y, &
                                                                                  z=sol%force_hydrostatic(i1:i2,j1:j2,k1:k2)%z)
            if (allocated(sol%force_pressure).and.present(patch)) &
            error=output%xml_writer%write_dataarray(data_name='force_pressure',x=sol%force_pressure(i1:i2,j1:j2,k1:k2)%x, &
                                                                               y=sol%force_pressure(i1:i2,j1:j2,k1:k2)%y, &
                                                                               z=sol%force_pressure(i1:i2,j1:j2,k1:k2)%z)
            if (allocated(sol%force_viscous).and.present(patch)) &
            error=output%xml_writer%write_dataarray(data_name='force_viscous',x=sol%force_viscous(i1:i2,j1:j2,k1:k2)%x, &
                                                                              y=sol%force_viscous(i1:i2,j1:j2,k1:k2)%y, &
                                                                              z=sol%force_viscous(i1:i2,j1:j2,k1:k2)%z)
            if (allocated(sol%torque_hydrostatic).and.present(patch)) &
            error=output%xml_writer%write_dataarray(data_name='torque_hydrostatic',x=sol%torque_hydrostatic(i1:i2,j1:j2,k1:k2)%x, &
                                                                                   y=sol%torque_hydrostatic(i1:i2,j1:j2,k1:k2)%y, &
                                                                                   z=sol%torque_hydrostatic(i1:i2,j1:j2,k1:k2)%z)
            if (allocated(sol%torque_pressure).and.present(patch)) &
            error=output%xml_writer%write_dataarray(data_name='torque_pressure',x=sol%torque_pressure(i1:i2,j1:j2,k1:k2)%x, &
                                                                                y=sol%torque_pressure(i1:i2,j1:j2,k1:k2)%y, &
                                                                                z=sol%torque_pressure(i1:i2,j1:j2,k1:k2)%z)
            if (allocated(sol%torque_viscous).and.present(patch)) &
            error=output%xml_writer%write_dataarray(data_name='torque_viscous',x=sol%torque_viscous(i1:i2,j1:j2,k1:k2)%x, &
                                                                               y=sol%torque_viscous(i1:i2,j1:j2,k1:k2)%y, &
                                                                               z=sol%torque_viscous(i1:i2,j1:j2,k1:k2)%z)
            if (allocated(sol%yplus).and.present(patch)) &
            error=output%xml_writer%write_dataarray(data_name='yplus',x=sol%yplus(i1:i2,j1:j2,k1:k2),one_component=.true.)
            if (allocated(sol%tau).and.present(patch)) &
            error=output%xml_writer%write_dataarray(data_name='tau',x=sol%tau(i1:i2,j1:j2,k1:k2)%x, &
                                                                    y=sol%tau(i1:i2,j1:j2,k1:k2)%y, &
                                                                    z=sol%tau(i1:i2,j1:j2,k1:k2)%z)
            if (allocated(sol%div_tau).and.present(patch)) &
            error=output%xml_writer%write_dataarray(data_name='div_tau',x=sol%div_tau(i1:i2,j1:j2,k1:k2),one_component=.true.)
         endif
         if (is_cell_centered) then
            error = output%xml_writer%write_dataarray(location='cell', action='close')
         else
            error = output%xml_writer%write_dataarray(location='node', action='close')
         endif
         error = output%xml_writer%write_piece()
         error = output%finalize()
      endassociate
   enddo
   endsubroutine save_block_file_vtk

   subroutine save_blocks_file_vtm(file_name, files_blocks)
   !< Save multi Xnavis-block-data into a VTM file.
   character(*), intent(in)  :: file_name           !< File name.
   type(string), intent(in)  :: files_blocks(1:,1:) !< Blocks of VTK files names.
   type(vtm_file)            :: file_vtm            !< VTM file handler.
   integer(I4P)              :: error               !< Error status.
   integer(I4P)              :: b                   !< Counter.

   error = file_vtm%initialize(filename=trim(adjustl(file_name)))
   do b=1, size(files_blocks, dim=2)
      error = file_vtm%write_block(filenames=files_blocks(1,b)//'', names=files_blocks(2,b)//'', name=files_blocks(3,b)//'')
   enddo
   error = file_vtm%finalize()
   endsubroutine save_blocks_file_vtm
endmodule xview_file_vtk_object
