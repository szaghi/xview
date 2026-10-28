!< xview, file TEC class definition.
module xview_file_tec_object
!< xview, file TEC class definition.
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf
use vecfor
use stringifor

use xview_block_grd_object
use xview_block_icc_object
use xview_block_rst_object

implicit none
private
public :: file_tec_object

type :: file_tec_object
   !< File TEC class definition.
   character(len=:), allocatable :: path         !< Output path.
   character(len=:), allocatable :: file_name    !< File name.
   character(len=:), allocatable :: tecvarname   !< Variables name for tecplot header file.
   integer(I4P),     allocatable :: tecvarloc(:) !< Tecplot array of variables location.
   character(len=:), allocatable :: tecvarlocstr !< Tecplot string of variables location.
   integer(I4P),     allocatable :: tecnull(:)   !< Tecplot null array.
   integer(I4P)                  :: nvar=3       !< Number of variables saved.
   integer(I4P)                  :: file_unit=0  !< File unit.
   contains
      ! public methods
      procedure, pass(self) :: destroy             !< Destroy dynamic memory.
      procedure, pass(self) :: finalize            !< Finalize file.
      procedure, pass(self) :: initialize          !< Initialize file.
      procedure, pass(self) :: is_file_present     !< Inquire if the file path is valid.
      procedure, pass(self) :: save_block_file_tec !< Save one Xnavis-block-data into a Tecplot file.
endtype file_tec_object

character(1), parameter :: tecendrec=char(0) !< End-character for binary-record end.
integer(I4P), external ::  tecini112,    &
                           tecauxstr112, &
                           teczne112,    &
                           tecdat112,    &
                           tecend112         !< Tecplot *tecio* functions.

contains
   elemental subroutine destroy(self)
   !< Destroy dynamic memory.
   class(file_tec_object), intent(inout) :: self !< File data.

   if (allocated(self%path        )) deallocate(self%path        )
   if (allocated(self%file_name   )) deallocate(self%file_name   )
   if (allocated(self%tecvarname  )) deallocate(self%tecvarname  )
   if (allocated(self%tecvarloc   )) deallocate(self%tecvarloc   )
   if (allocated(self%tecvarlocstr)) deallocate(self%tecvarlocstr)
   if (allocated(self%tecnull     )) deallocate(self%tecnull     )
   self%nvar = 3
   self%file_unit = 0
   endsubroutine destroy

   subroutine finalize(self, is_binary)
   !< Finalize file.
   class(file_tec_object), intent(inout) :: self      !< File data.
   logical,                intent(in)    :: is_binary !< Define binary or ascii output.
   integer(I4P)                          :: error     !< Error status.

   if (is_binary) then
      error = tecend112()
   else
      close(self%file_unit)
   endif
   endsubroutine finalize

   subroutine initialize(self, path, file_name, is_binary, is_cell_centered, icc, sol)
   !< Initialize file.
   !<
   !< @note Variables names are built from the same `allocated` checks of `save_block_file_tec`, in the same order: any field
   !< added there must be added here too.
   class(file_tec_object), intent(inout)        :: self             !< File data.
   character(*),           intent(in)           :: path             !< Output path.
   character(*),           intent(in)           :: file_name        !< File name.
   logical,                intent(in)           :: is_binary        !< Define binary or ascii output.
   logical,                intent(in)           :: is_cell_centered !< Define variables at cell centers or nodes.
   type(block_icc_object), intent(in), optional :: icc              !< Representative icc block.
   type(block_rst_object), intent(in), optional :: sol              !< Representative sol block.
   character(1)                                 :: vd               !< Variables delimiter, '"' or ''.
   integer(I4P)                                 :: error            !< Error status.

   call self%destroy
   self%path      = trim(adjustl(path))
   self%file_name = trim(adjustl(file_name))
   if (is_binary) then
      vd = ' '
      self%tecvarname = 'x y z'
   else
      vd = '"'
      self%tecvarname = ' VARIABLES ="x" "y" "z"'
   endif
   if (present(icc)) then
      if (icc%is_loaded) call add_var('icc')
   endif
   if (present(sol)) then
      if (sol%is_loaded) then
         call add_var('u') ; call add_var('v') ; call add_var('w') ; call add_var('p')
         if (allocated(sol%level_set                           )) call add_var('f')
         if (allocated(sol%level_set_zero                      )) call add_var('f0')
         if (allocated(sol%viscosity                           )) call add_var('visc')
         if (allocated(sol%turbulent_viscosity                 )) call add_var('vitl')
         if (allocated(sol%turbulent_kinetic_energy            )) call add_var('ken')
         if (allocated(sol%turbulent_kinetic_energy_dissipation)) call add_var('eps')
         if (allocated(sol%lambda2                             )) call add_var('lambda2')
         if (allocated(sol%qfactor                             )) call add_var('qfactor')
         if (allocated(sol%liutex                              )) call add_var('liutex')
         if (allocated(sol%helicity                            )) call add_var('helicity')
         if (allocated(sol%vorticity)) then
            call add_var('vorticity_x') ; call add_var('vorticity_y') ; call add_var('vorticity_z')
         endif
         if (allocated(sol%grad_p)) then
            call add_var('grad_p_x') ; call add_var('grad_p_y') ; call add_var('grad_p_z')
         endif
         if (allocated(sol%LT)) then
            call add_var('LT_11') ; call add_var('LT_12') ; call add_var('LT_13')
            call add_var('LT_21') ; call add_var('LT_22') ; call add_var('LT_23')
            call add_var('LT_31') ; call add_var('LT_32') ; call add_var('LT_33')
         endif
         if (allocated(sol%div2LT                              )) call add_var('div2LT')
      endif
   endif
   allocate(self%tecnull(  1:self%nvar))
   allocate(self%tecvarloc(1:self%nvar))
   if (is_binary) then
      self%tecnull = 0
      self%tecvarloc = 1
      if (self%nvar > 3.and.is_cell_centered) self%tecvarloc(4:self%nvar) = 0
   endif
   if (is_binary) then
      error = tecini112(tecendrec,self%tecvarname//tecendrec,self%path//self%file_name//".plt"//tecendrec,'.'//tecendrec,0,0,1)
   else
     open(newunit=self%file_unit, file=self%path//self%file_name//".dat")
     write(self%file_unit, '(A)', iostat=error)self%tecvarname
   endif
   contains
      subroutine add_var(name)
      !< Add a variable to the header names list.
      character(*), intent(in) :: name !< Variable name.

      self%tecvarname = self%tecvarname//' '//trim(vd)//name//trim(vd)
      self%nvar = self%nvar + 1
      endsubroutine add_var
   endsubroutine initialize

   function is_file_present(self) result(is_present)
   !< Inquire if the file path is valid.
   class(file_tec_object), intent(in) :: self       !< File data.
   logical                            :: is_present !< Inquiring result.

   is_present = allocated(self%file_name)
   if (is_present) inquire(file=trim(adjustl(self%file_name)), exist=is_present)
   endfunction is_file_present

   subroutine save_block_file_tec(self, is_binary, blocks_map, grd, is_cell_centered, icc, sol, patch)
   !< Save one Xnavis-block-data into a Tecplot file.
   class(file_tec_object), intent(inout)        :: self              !< File data.
   logical,                intent(in)           :: is_binary         !< Define binary or ascii output.
   integer(I4P),           intent(in)           :: blocks_map(1:)    !< Blocks (processors) map.
   type(block_grd_object), intent(in)           :: grd               !< Grid of block.
   logical,                intent(in)           :: is_cell_centered  !< Define variables at cell centers or nodes.
   type(block_icc_object), intent(in), optional :: icc               !< Boundary condititions of block.
   type(block_rst_object), intent(in), optional :: sol               !< Solution of block.
   integer(I4P),           intent(in), optional :: patch             !< Patch boundary conditions.
   character(len=:), allocatable                :: zone_name         !< Tecplot zone name.
   integer(I4P)                                 :: error             !< Error status.
   integer(I4P)                                 :: nodes_number      !< Number of nodes.
   integer(I4P)                                 :: cells_number      !< Number of cells.
   integer(I4P)                                 :: ncn               !< Number of cells or nodes.
   integer(I4P)                                 :: i1,i2,j1,j2,k1,k2 !< Actual extents.
   integer(I4P)                                 :: p                 !< Counter.

   if (.not.grd%is_loaded) error stop 'error: cannot save TEC block file because passed grd block has not been loaded'
   if (.not.allocated(grd%patches_extents)) return ! queried a patch bc that is not present into this block

   do p=1, size(grd%patches_extents, dim=1) ! loop over pathces or whole block
      associate(ci1=>grd%patches_extents(p, 1 ), ci2=>grd%patches_extents(p, 2 ), &
                cj1=>grd%patches_extents(p, 3 ), cj2=>grd%patches_extents(p, 4 ), &
                ck1=>grd%patches_extents(p, 5 ), ck2=>grd%patches_extents(p, 6 ), &
                ni1=>grd%patches_extents(p, 7 ), ni2=>grd%patches_extents(p, 8 ), &
                nj1=>grd%patches_extents(p, 9 ), nj2=>grd%patches_extents(p, 10), &
                nk1=>grd%patches_extents(p, 11), nk2=>grd%patches_extents(p, 12), &
                file_unit=>self%file_unit, nvar=>self%nvar, face=>grd%patches_extents(p,0))
         if (present(patch)) then
            zone_name = 'patch_'//trim(strz(patch, nz_pad=3))//                      &
                        '-face_'//trim(str(face, no_sign=.true.))// &
                        '-blk_'//trim(strz(blocks_map(2), nz_pad=4))//'-grp_'//trim(strz(blocks_map(1), nz_pad=3))
         else
            zone_name = 'blk_'//trim(strz(blocks_map(2), nz_pad=4))//'-grp_'//trim(strz(blocks_map(1), nz_pad=3))
         endif
         nodes_number = (ni2 - ni1 + 1) * (nj2 - nj1 + 1) * (nk2 - nk1 + 1)
         cells_number = (ci2 - ci1 + 1) * (cj2 - cj1 + 1) * (ck2 - ck1 + 1)
         if (is_cell_centered) then
            i1 = ci1 ; i2 = ci2
            j1 = cj1 ; j2 = cj2
            k1 = ck1 ; k2 = ck2
            ncn = cells_number
         else
            i1 = ni1 ; i2 = ni2
            j1 = nj1 ; j2 = nj2
            k1 = nk1 ; k2 = nk2
            ncn = nodes_number
         endif
         if (is_binary) then
            error = teczne112(zone_name//tecendrec, &
                              0,                    &
                              (ni2 - ni1 + 1),      &
                              (nj2 - nj1 + 1),      &
                              (nk2 - nk1 + 1),      &
                              0,                    &
                              0,                    &
                              0,                    &
                              0.0,                  &
                              0,                    &
                              0,                    &
                              1,                    &
                              0,                    &
                              0,                    &
                              0,                    &
                              0,                    &
                              0,                    &
                              self%tecnull,         &
                              self%tecvarloc,       &
                              self%tecnull,         &
                              0)
         else
            if (nvar > 3.and.is_cell_centered) then
               write(file_unit,'(A)')' ZONE  T="'//zone_name//'"'//                        &
                                     ', I='//trim(str(no_sign=.true.,n=(ni2 - ni1 + 1)))// &
                                     ', J='//trim(str(no_sign=.true.,n=(nj2 - nj1 + 1)))// &
                                     ', K='//trim(str(no_sign=.true.,n=(nk2 - nk1 + 1)))// &
                                     ', DATAPACKING=BLOCK'//                               &
                                     ', VARLOCATION=([1-3]=NODAL,[4-'//trim(str(nvar, .true.))//']=CELLCENTERED)'
            else
               write(file_unit,'(A)')' ZONE  T="'//zone_name//'"'//                        &
                                     ', I='//trim(str(no_sign=.true.,n=(ni2 - ni1 + 1)))// &
                                     ', J='//trim(str(no_sign=.true.,n=(nj2 - nj1 + 1)))// &
                                     ', K='//trim(str(no_sign=.true.,n=(nk2 - nk1 + 1)))// &
                                     ', DATAPACKING=BLOCK'//                               &
                                     ', VARLOCATION=([1-'//trim(str(nvar,.true.))//']=NODAL)'
            endif
         endif
         if (grd%is_loaded) then
            error = tec_var(n=nodes_number, var=grd%nodes(ni1:ni2,nj1:nj2,nk1:nk2)%x, d=1_I4P)
            error = tec_var(n=nodes_number, var=grd%nodes(ni1:ni2,nj1:nj2,nk1:nk2)%y, d=1_I4P)
            error = tec_var(n=nodes_number, var=grd%nodes(ni1:ni2,nj1:nj2,nk1:nk2)%z, d=1_I4P)
         else
            error stop 'error: cannot save TEC block file because passed grd block has not been loaded'
         endif
         if (present(icc)) then
            if (icc%is_loaded) then
               error = tec_var(n=ncn, var=real(icc%rcc(i1:i2,j1:j2,k1:k2), R8P), d=1_I4P)
            endif
         endif
         if (present(sol)) then
            if (sol%is_loaded) then
               error = tec_var(n=ncn, var=sol%momentum(i1:i2,j1:j2,k1:k2)%x, d=1_I4P)
               error = tec_var(n=ncn, var=sol%momentum(i1:i2,j1:j2,k1:k2)%y, d=1_I4P)
               error = tec_var(n=ncn, var=sol%momentum(i1:i2,j1:j2,k1:k2)%z, d=1_I4P)
               error = tec_var(n=ncn, var=sol%pressure(i1:i2,j1:j2,k1:k2),   d=1_I4P)
               if (allocated(sol%level_set)) &
                  error = tec_var(n=ncn, var=sol%level_set(i1:i2,j1:j2,k1:k2), d=1_I4P)
               if (allocated(sol%level_set_zero)) &
                  error = tec_var(n=ncn, var=sol%level_set_zero(i1:i2,j1:j2,k1:k2), d=1_I4P)
               if (allocated(sol%viscosity)) &
                  error=tec_var(n=ncn, var=sol%viscosity(i1:i2,j1:j2,k1:k2), d=1_I4P)
               if (allocated(sol%turbulent_viscosity)) &
                  error=tec_var(n=ncn, var=sol%turbulent_viscosity(i1:i2,j1:j2,k1:k2), d=1_I4P)
               if (allocated(sol%turbulent_kinetic_energy)) &
                  error=tec_var(n=ncn, var=sol%turbulent_kinetic_energy(i1:i2,j1:j2,k1:k2), d=1_I4P)
               if (allocated(sol%turbulent_kinetic_energy_dissipation)) &
                  error=tec_var(n=ncn, var=sol%turbulent_kinetic_energy_dissipation(i1:i2,j1:j2,k1:k2), d=1_I4P)
               if (allocated(sol%lambda2)) &
                  error = tec_var(n=ncn, var=sol%lambda2(i1:i2,j1:j2,k1:k2), d=1_I4P)
               if (allocated(sol%qfactor)) &
                  error = tec_var(n=ncn, var=sol%qfactor(i1:i2,j1:j2,k1:k2), d=1_I4P)
               if (allocated(sol%liutex)) &
                  error = tec_var(n=ncn, var=sol%liutex(i1:i2,j1:j2,k1:k2), d=1_I4P)
               if (allocated(sol%helicity)) &
                  error = tec_var(n=ncn, var=sol%helicity(i1:i2,j1:j2,k1:k2), d=1_I4P)
               if (allocated(sol%vorticity)) then
                  error = tec_var(n=ncn, var=sol%vorticity(i1:i2,j1:j2,k1:k2)%x, d=1_I4P)
                  error = tec_var(n=ncn, var=sol%vorticity(i1:i2,j1:j2,k1:k2)%y, d=1_I4P)
                  error = tec_var(n=ncn, var=sol%vorticity(i1:i2,j1:j2,k1:k2)%z, d=1_I4P)
               endif
               if (allocated(sol%grad_p)) then
                  error = tec_var(n=ncn, var=sol%grad_p(i1:i2,j1:j2,k1:k2)%x, d=1_I4P)
                  error = tec_var(n=ncn, var=sol%grad_p(i1:i2,j1:j2,k1:k2)%y, d=1_I4P)
                  error = tec_var(n=ncn, var=sol%grad_p(i1:i2,j1:j2,k1:k2)%z, d=1_I4P)
               endif
               if (allocated(sol%LT)) then
                  error = tec_var(n=ncn, var=sol%LT(1,1,i1:i2,j1:j2,k1:k2), d=1_I4P)
                  error = tec_var(n=ncn, var=sol%LT(1,2,i1:i2,j1:j2,k1:k2), d=1_I4P)
                  error = tec_var(n=ncn, var=sol%LT(1,3,i1:i2,j1:j2,k1:k2), d=1_I4P)
                  error = tec_var(n=ncn, var=sol%LT(2,1,i1:i2,j1:j2,k1:k2), d=1_I4P)
                  error = tec_var(n=ncn, var=sol%LT(2,2,i1:i2,j1:j2,k1:k2), d=1_I4P)
                  error = tec_var(n=ncn, var=sol%LT(2,3,i1:i2,j1:j2,k1:k2), d=1_I4P)
                  error = tec_var(n=ncn, var=sol%LT(3,1,i1:i2,j1:j2,k1:k2), d=1_I4P)
                  error = tec_var(n=ncn, var=sol%LT(3,2,i1:i2,j1:j2,k1:k2), d=1_I4P)
                  error = tec_var(n=ncn, var=sol%LT(3,3,i1:i2,j1:j2,k1:k2), d=1_I4P)
               endif
               if (allocated(sol%div2LT)) &
                  error = tec_var(n=ncn, var=sol%div2LT(i1:i2,j1:j2,k1:k2), d=1_I4P)
            endif
         endif
      endassociate
   enddo
   contains
      function tec_var(n, var, d) result(err)
      !< Save variables into Tecplot file.
      integer(I4P), intent(in) :: n        !< Number of var elements.
      real(R8P),    intent(in) :: var(1:n) !< Variable to be saved.
      integer(I4P), intent(in) :: d        !< Double precision output (1 yes, 0 no).
      integer(I4P)             :: err      !< Error traping flag: 0 no errors, >0 error occours.
      integer(I4P)             :: e        !< Element counter.

      if (is_binary) then
         err = tecdat112(n, var, d)
      else
         write(self%file_unit, FR8P, iostat=err) (var(e), e=1,n)
      endif
      endfunction tec_var
   endsubroutine save_block_file_tec
endmodule xview_file_tec_object
