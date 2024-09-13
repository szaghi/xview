!< xview, block (rst) class definition.
module xview_block_esz_object
!< xview, block (rst) class definition.

use penf
use vecfor

use xview_block_grd_object
use xview_block_icc_object
use xview_block_rst_object

implicit none
private
public :: block_esz_object

type :: block_esz_object
   !< **Block** class for *file sol*.
   integer(I4P)           :: Ni=0                     !< Number of cells in i direction, into block containing subzone.
   integer(I4P)           :: Nj=0                     !< Number of cells in j direction, into block containing subzone.
   integer(I4P)           :: Nk=0                     !< Number of cells in k direction, into block containing subzone.
   integer(I4P)           :: gc(6)=[2, 2, 2, 2, 2, 2] !< Number of ghost cells, into block containing subzone.
   integer(I4P)           :: subzone_extents(7)=0_I4P !< Extracted subzone extents (blk-global,i1,i2,j1,j2,k1,k2).
   type(block_grd_object) :: grd                      !< Extracted subzone grid.
   type(block_icc_object) :: icc                      !< Extracted subzone icc.
   type(block_rst_object) :: sol                      !< Extracted subzone solution.
   logical                :: is_loaded=.false.        !< Flag for checking if the block is loaded.
   contains
      ! public methods
      procedure, pass(self) :: destroy       !< Destroy dynamic memory.
      procedure, pass(self) :: load_solution !< Load subzone solution from file.
   endtype block_esz_object

contains
   elemental subroutine destroy(self)
   !< Destroy dynamic memory.
   class(block_esz_object), intent(inout) :: self !< Block data.

   self%Ni = 0_I4P
   self%Nj = 0_I4P
   self%Nk = 0_I4P
   self%gc = [2, 2, 2, 2, 2, 2]
   self%subzone_extents = 0_I4P
   call self%grd%destroy
   call self%icc%destroy
   call self%sol%destroy
   self%is_loaded = .false.
   endsubroutine destroy

   subroutine load_solution(self,file_unit,is_cell_centered,patch,RE,rFR2,zfs,is_level_set,is_zeroeq,is_oneeq,is_twoeq, &
                            compute_lambda2,compute_qfactor,compute_helicity,compute_vorticity,                         &
                            compute_div2LT,compute_grad_p,compute_k_ratio,compute_yplus,compute_tau,compute_div_tau,compute_loads)
   !< Load subzone solution from file.
   class(block_esz_object), intent(inout)        :: self              !< Block data.
   integer(I4P),            intent(in)           :: file_unit         !< Logical file unit.
   logical,                 intent(in), optional :: is_cell_centered  !< Define variables at cell centers or nodes.
   integer(I4P),            intent(in), optional :: patch             !< Patch boundary conditions.
   real(R8P),               intent(in), optional :: RE                !< Reynolds number.
   real(R8P),               intent(in), optional :: rFR2              !< 1/(Froude number)^2.
   real(R8P),               intent(in), optional :: zfs               !< Z quote of free surface.
   logical,                 intent(in), optional :: is_level_set      !< Flag for level set function presence.
   logical,                 intent(in), optional :: is_zeroeq         !< Use *zero* equations turbulence model.
   logical,                 intent(in), optional :: is_oneeq          !< Use *one* equations turbulence model.
   logical,                 intent(in), optional :: is_twoeq          !< Use *two* equations turbulence model.
   logical,                 intent(in), optional :: compute_lambda2   !< Compute lamda2 field.
   logical,                 intent(in), optional :: compute_qfactor   !< Compute qfactor field.
   logical,                 intent(in), optional :: compute_helicity  !< Compute helicity field.
   logical,                 intent(in), optional :: compute_vorticity !< Compute vorticity field.
   logical,                 intent(in), optional :: compute_div2LT    !< Compute double divergence of Lighthill tensor.
   logical,                 intent(in), optional :: compute_grad_p    !< Compute pressure gradient.
   logical,                 intent(in), optional :: compute_k_ratio   !< Compute kinetic energy ratio.
   logical,                 intent(in), optional :: compute_yplus     !< Compute y+ field.
   logical,                 intent(in), optional :: compute_tau       !< Compute tau field.
   logical,                 intent(in), optional :: compute_div_tau   !< Compute divergence of tau field.
   logical,                 intent(in), optional :: compute_loads     !< Compute loads (forces and torques).
   real(R8P), allocatable                        :: rcc_buf(:,:,:)    !< Buffer for reading rcc.
   integer(kind=I4P)                             :: i1,j1,k1,i2,j2,k2 !< Subzone extents.
   integer(I4P)                                  :: i,j,k             !< Counter.

   call self%destroy
   self%sol%is_level_set = .false. ; if (present(is_level_set)) self%sol%is_level_set = is_level_set
   self%sol%is_zeroeq    = .false. ; if (present(is_zeroeq   )) self%sol%is_zeroeq    = is_zeroeq
   self%sol%is_oneeq     = .false. ; if (present(is_oneeq    )) self%sol%is_oneeq     = is_oneeq
   self%sol%is_twoeq     = .false. ; if (present(is_twoeq    )) self%sol%is_twoeq     = is_twoeq
   associate(Ni=>self%Ni, Nj=>self%Nj, Nk=>self%Nk, gc=>self%gc, subzone_extents=>self%subzone_extents, &
             grd=>self%grd, icc=>self%icc, sol=>self%sol)
      ! load subzone extents
      read(file_unit) subzone_extents
      i1 = subzone_extents(2)
      i2 = subzone_extents(3)
      j1 = subzone_extents(4)
      j2 = subzone_extents(5)
      k1 = subzone_extents(6)
      k2 = subzone_extents(7)
      Ni = i2-i1+1 ; Nj = j2-j1+1 ; Nk = k2-k1+1
      ! load subzone grid nodes
      grd%Ni = Ni ; grd%Nj = Nj ; grd%Nk = Nk
      call grd%alloc(is_metrics_to_allocate=.true.)
      read(file_unit)(((grd%nodes(i,j,k)%x,i=i1,i2),j=j1,j2),k=k1,k2)
      read(file_unit)(((grd%nodes(i,j,k)%y,i=i1,i2),j=j1,j2),k=k1,k2)
      read(file_unit)(((grd%nodes(i,j,k)%z,i=i1,i2),j=j1,j2),k=k1,k2)
      allocate(grd%patches_extents(1:1, 0:12))
      grd%patches_extents(1, 0 ) = -999 ! whole block subzone, no patches
      grd%patches_extents(1, 1 ) = i1+1 ; grd%patches_extents(1, 2 ) = i2
      grd%patches_extents(1, 3 ) = j1+1 ; grd%patches_extents(1, 4 ) = j2
      grd%patches_extents(1, 5 ) = k1+1 ; grd%patches_extents(1, 6 ) = k2
      grd%patches_extents(1, 7 ) = i1   ; grd%patches_extents(1, 8 ) = i2
      grd%patches_extents(1, 9 ) = j1   ; grd%patches_extents(1, 10) = j2
      grd%patches_extents(1, 11) = k1   ; grd%patches_extents(1, 12) = k2
      call grd%compute_metrics
      grd%is_loaded = .true.
      ! load subzone solution
      call sol%init(Ni=Ni,Nj=Nj,Nk=Nk,                                                                     &
                    is_level_set=is_level_set,is_zeroeq=is_zeroeq,is_oneeq=is_oneeq,is_twoeq=is_twoeq,     &
                    has_lambda2=compute_lambda2,has_qfactor=compute_qfactor,has_helicity=compute_helicity, &
                    has_vorticity=compute_vorticity,has_div2LT=compute_div2LT,has_grad_p=compute_grad_p,   &
                    has_k_ratio=compute_k_ratio,                                                           &
                    has_yplus=compute_yplus,has_tau=compute_tau,has_div_tau=compute_div_tau,has_loads=compute_loads)
      read(file_unit)(((sol%momentum(i,j,k)%x,i=i1,i2),j=j1,j2),k=k1,k2)
      read(file_unit)(((sol%momentum(i,j,k)%y,i=i1,i2),j=j1,j2),k=k1,k2)
      read(file_unit)(((sol%momentum(i,j,k)%z,i=i1,i2),j=j1,j2),k=k1,k2)
      read(file_unit)(((sol%pressure(i,j,k)  ,i=i1,i2),j=j1,j2),k=k1,k2)
      if (sol%is_level_set) then
         read(file_unit)(((sol%level_set(     i,j,k),i=i1,i2),j=j1,j2),k=k1,k2)
         read(file_unit)(((sol%level_set_zero(i,j,k),i=i1,i2),j=j1,j2),k=k1,k2)
      endif
      if (self%sol%is_zeroeq.or.self%sol%is_oneeq.or.self%sol%is_twoeq) then
         read(file_unit)(((sol%viscosity(i,j,k),i=i1,i2),j=j1,j2),k=k1,k2)
      endif
      if (self%sol%is_oneeq) read(file_unit)(((sol%turbulent_viscosity(i,j,k),i=i1,i2),j=j1,j2),k=k1,k2)
      if (self%sol%is_twoeq) then
         read(file_unit)(((sol%turbulent_kinetic_energy(i,j,k),i=i1,i2),j=j1,j2),k=k1,k2)
         read(file_unit)(((sol%turbulent_kinetic_energy_dissipation(i,j,k),i=i1,i2),j=j1,j2),k=k1,k2)
      endif
      sol%is_loaded = .true.
      ! load subzone icc rcc
      icc%Ni = Ni ; icc%Nj = Nj ; icc%Nk = Nk
      call icc%alloc
      allocate(rcc_buf(1-icc%gc(1):icc%Ni+icc%gc(2), 1-icc%gc(3):icc%Nj+icc%gc(4), 1-icc%gc(5):icc%Nk+icc%gc(6)))
      read(file_unit)(((rcc_buf(i,j,k),i=i1,i2),j=j1,j2),k=k1,k2)
      icc%rcc = real(rcc_buf, kind=R4P)
      icc%is_loaded = .true.
      call sol%compute_aux(grd=grd, icc=icc, patch=patch, RE=RE, rFR2=rFR2, zfs=zfs)
      if (present(is_cell_centered)) then
         if (.not.is_cell_centered) call sol%interpolate_at_nodes
      endif
   endassociate
   self%is_loaded = .true.
   endsubroutine load_solution
endmodule xview_block_esz_object
