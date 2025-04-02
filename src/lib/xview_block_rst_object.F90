!< xview, block (rst) class definition.
module xview_block_rst_object
!< xview, block (rst) class definition.

use penf
use vecfor

use xview_block_grd_object
use xview_block_icc_object

implicit none
private
public :: block_rst_object

type :: block_rst_object
   !< **Block** class for *file sol*.
   integer(I4P)              :: Ni=0                                        !< Number of cells in i direction.
   integer(I4P)              :: Nj=0                                        !< Number of cells in j direction.
   integer(I4P)              :: Nk=0                                        !< Number of cells in k direction.
   integer(I4P)              :: gc(6)=[2, 2, 2, 2, 2, 2]                    !< Number of ghost cells.
   type(vector), allocatable :: momentum(:,:,:)                             !< Momentum field.
   real(R8P),    allocatable :: pressure(:,:,:)                             !< Pressure field.
   real(R8P),    allocatable :: level_set(:,:,:)                            !< Level set field.
   real(R8P),    allocatable :: level_set_zero(:,:,:)                       !< Zero value of level set field.
   real(R8P),    allocatable :: viscosity(:,:,:)                            !< Viscosity field.
   real(R8P),    allocatable :: turbulent_viscosity(:,:,:)                  !< Turbulent viscosity field.
   real(R8P),    allocatable :: turbulent_kinetic_energy(:,:,:)             !< Turbulent kinetic energy field.
   real(R8P),    allocatable :: turbulent_kinetic_energy_dissipation(:,:,:) !< Turbulent kinetic energy dissipation field.
   real(R8P),    allocatable :: lambda2(:,:,:)                              !< Variable to identify vortices (lambda 2).
   real(R8P),    allocatable :: qfactor(:,:,:)                              !< Variable to identify vortices (q factor).
   real(R8P),    allocatable :: liutex(:,:,:)                               !< Variable to identify vortices (liutex).
   real(R8P),    allocatable :: helicity(:,:,:)                             !< Helicity.
   type(vector), allocatable :: vorticity(:,:,:)                            !< Vorticity.
   type(vector), allocatable :: grad_p(:,:,:)                               !< Pressure Gradient.
   real(R8P),    allocatable :: div2LT(:,:,:)                               !< Double divergence of Lighthill tensor.
   real(R8P),    allocatable :: k_ratio(:,:,:)                              !< Kinetic energies ratio.
   real(R8P),    allocatable :: yplus(:,:,:)                                !< Estimation of y+.
   type(vector), allocatable :: tau(:,:,:)                                  !< Tau wall.
   real(R8P),    allocatable :: div_tau(:,:,:)                              !< Divergence of Tau wall.
   type(vector), allocatable :: force_hydrostatic(:,:,:)                    !< Hydrostic part of forces.
   type(vector), allocatable :: force_pressure(:,:,:)                       !< Pressure part of forces.
   type(vector), allocatable :: force_viscous(:,:,:)                        !< Viscous part of forces.
   type(vector), allocatable :: torque_hydrostatic(:,:,:)                   !< Hydrostic part of torques.
   type(vector), allocatable :: torque_pressure(:,:,:)                      !< Pressure part of torques.
   type(vector), allocatable :: torque_viscous(:,:,:)                       !< Viscous part of torques.
   logical                   :: is_level_set=.false.                        !< Use Level Set model.
   logical                   :: is_zeroeq=.false.                           !< Use *zero* equations turbulence model.
   logical                   :: is_oneeq=.false.                            !< Use *one* equations turbulence model.
   logical                   :: is_twoeq=.false.                            !< Use *two* equations turbulence model.
   logical                   :: has_lambda2=.false.                         !< Solution has lamda2 field.
   logical                   :: has_qfactor=.false.                         !< Solution has qfactor field.
   logical                   :: has_helicity=.false.                        !< Solution has helicity field.
   logical                   :: has_vorticity=.false.                       !< Solution has vorticity field.
   logical                   :: has_div2LT=.false.                          !< Solution has double divergence of Lighthill tensor.
   logical                   :: has_grad_p=.false.                          !< Solution has pressure gradient field.
   logical                   :: has_k_ratio=.false.                         !< Solution has kinetic energy ratio.
   logical                   :: has_yplus=.false.                           !< Solution has y+ field.
   logical                   :: has_tau=.false.                             !< Solution has tau field.
   logical                   :: has_div_tau=.false.                         !< Solution has divergence of tau field.
   logical                   :: has_loads=.false.                           !< Solution has loads (forces and torques).
   logical                   :: is_loaded=.false.                           !< Flag for checking if the block is loaded.
   contains
      ! public methods
      procedure, pass(self) :: destroy                      !< Destroy dynamic memory.
      procedure, pass(self) :: alloc                        !< Allocate dynamic memory.
      procedure, pass(self) :: compute_aux                  !< Compute auxiliary varibles.
      procedure, pass(self) :: compute_div2LT               !< Compute double divergence of Lighthill tensor.
      procedure, pass(self) :: compute_kinetic_energy_ratio !< Compute ratio of kinetic energies.
      procedure, pass(self) :: compute_loads                !< Compute loads (forces and torques) on patches.
      procedure, pass(self) :: compute_tau                  !< Compute Tau wall and its divergence on patches.
      procedure, pass(self) :: compute_vorticity            !< Compute vorticity related varibles.
      procedure, pass(self) :: compute_grad_p               !< Compute pressure gradient.
      procedure, pass(self) :: compute_yplus                !< Compute yplus on patches.
      procedure, pass(self) :: has_level_set                !< Return true if block has level set functions.
      procedure, pass(self) :: has_turb_kinetic_energy      !< Return true if block has turbulent kinetic energy field.
      procedure, pass(self) :: has_turb_kinetic_energy_diss !< Return true if block has turb. k. energy dissipation field.
      procedure, pass(self) :: has_turb_viscosity           !< Return true if block has turbulent viscosity field.
      procedure, pass(self) :: has_viscosity                !< Return true if block has viscosity field.
      procedure, pass(self) :: init                         !< Initialize block.
      procedure, pass(self) :: interpolate_at_nodes         !< Interpolate solution at nodes.
      procedure, pass(self) :: load_dimensions              !< Load block dimensions from file.
      procedure, pass(self) :: load_solution                !< Load block solution from file.
      procedure, pass(self) :: save_dimensions              !< Save block dimensions into file.
      procedure, pass(self) :: save_solution                !< Save block solution into file.
endtype block_rst_object

contains
   ! public methods
   elemental subroutine destroy(self)
   !< Destroy dynamic memory.
   class(block_rst_object), intent(inout) :: self !< Block data.

   self%Ni = 0
   self%Nj = 0
   self%Nk = 0
   self%gc = 2
   if (allocated(self%momentum))                             deallocate(self%momentum)
   if (allocated(self%pressure))                             deallocate(self%pressure)
   if (allocated(self%level_set))                            deallocate(self%level_set)
   if (allocated(self%level_set_zero))                       deallocate(self%level_set_zero)
   if (allocated(self%viscosity))                            deallocate(self%viscosity)
   if (allocated(self%turbulent_viscosity))                  deallocate(self%turbulent_viscosity)
   if (allocated(self%turbulent_kinetic_energy))             deallocate(self%turbulent_kinetic_energy)
   if (allocated(self%turbulent_kinetic_energy_dissipation)) deallocate(self%turbulent_kinetic_energy_dissipation)
   if (allocated(self%lambda2))                              deallocate(self%lambda2)
   if (allocated(self%qfactor))                              deallocate(self%qfactor)
   if (allocated(self%liutex))                               deallocate(self%liutex)
   if (allocated(self%helicity))                             deallocate(self%helicity)
   if (allocated(self%vorticity))                            deallocate(self%vorticity)
   if (allocated(self%grad_p))                               deallocate(self%grad_p)
   if (allocated(self%div2LT))                               deallocate(self%div2LT)
   if (allocated(self%k_ratio))                              deallocate(self%k_ratio)
   if (allocated(self%yplus))                                deallocate(self%yplus)
   if (allocated(self%tau))                                  deallocate(self%tau)
   if (allocated(self%div_tau))                              deallocate(self%div_tau)
   if (allocated(self%force_hydrostatic))                    deallocate(self%force_hydrostatic)
   if (allocated(self%force_pressure))                       deallocate(self%force_pressure)
   if (allocated(self%force_viscous))                        deallocate(self%force_viscous)
   if (allocated(self%torque_hydrostatic))                   deallocate(self%torque_hydrostatic)
   if (allocated(self%torque_pressure))                      deallocate(self%torque_pressure)
   if (allocated(self%torque_viscous))                       deallocate(self%torque_viscous)
   self%is_level_set  =.false.
   self%is_zeroeq     =.false.
   self%is_oneeq      =.false.
   self%is_twoeq      =.false.
   self%has_lambda2   =.false.
   self%has_qfactor   =.false.
   self%has_liutex    =.false.
   self%has_helicity  =.false.
   self%has_vorticity =.false.
   self%has_yplus     =.false.
   self%has_tau       =.false.
   self%has_div_tau   =.false.
   self%has_div2LT    =.false.
   self%has_grad_p    =.false.
   self%has_k_ratio   =.false.
   self%has_loads     =.false.
   self%is_loaded     =.false.
   endsubroutine destroy

   elemental subroutine alloc(self)
   !< Allocate dynamic memory.
   class(block_rst_object), intent(inout) :: self !< Block data.

   associate(Ni => self%Ni, Nj => self%Nj, Nk => self%Nk, gc => self%gc)
   allocate(self%momentum(1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6)))
   allocate(self%pressure(1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6)))
   if (self%is_level_set) then
      allocate(self%level_set(1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6)))
      allocate(self%level_set_zero(1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6)))
   endif
   if (self%is_zeroeq.or.self%is_oneeq.or.self%is_twoeq) then
      allocate(self%viscosity(1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6)))
   endif
   if (self%is_oneeq) allocate(self%turbulent_viscosity(1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6)))
   if (self%is_twoeq) then
      allocate(self%turbulent_kinetic_energy(1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6)))
      allocate(self%turbulent_kinetic_energy_dissipation(1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6)))
   endif
   ! auxiliary variables
   if (self%has_lambda2  ) allocate(self%lambda2(  1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6)))
   if (self%has_qfactor  ) allocate(self%qfactor(  1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6)))
   if (self%has_liutex   ) allocate(self%liutex(   1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6)))
   if (self%has_helicity ) allocate(self%helicity( 1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6)))
   if (self%has_vorticity) allocate(self%vorticity(1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6)))
   if (self%has_div2LT   ) allocate(self%div2LT(   1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6)))
   if (self%has_grad_p   ) allocate(self%grad_p(   1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6)))
   if (self%has_k_ratio  ) allocate(self%k_ratio(  1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6)))
   if (self%has_yplus    ) allocate(self%yplus(    0-gc(1):Ni+gc(2), 0-gc(3):Nj+gc(4), 0-gc(5):Nk+gc(6)))
   if (self%has_tau      ) allocate(self%tau(      0-gc(1):Ni+gc(2), 0-gc(3):Nj+gc(4), 0-gc(5):Nk+gc(6)))
   if (self%has_div_tau  ) allocate(self%div_tau(  0-gc(1):Ni+gc(2), 0-gc(3):Nj+gc(4), 0-gc(5):Nk+gc(6)))
   if (self%has_loads) then
      allocate(self%force_hydrostatic( 0-gc(1):Ni+gc(2), 0-gc(3):Nj+gc(4), 0-gc(5):Nk+gc(6)))
      allocate(self%force_pressure(    0-gc(1):Ni+gc(2), 0-gc(3):Nj+gc(4), 0-gc(5):Nk+gc(6)))
      allocate(self%force_viscous(     0-gc(1):Ni+gc(2), 0-gc(3):Nj+gc(4), 0-gc(5):Nk+gc(6)))
      allocate(self%torque_hydrostatic(0-gc(1):Ni+gc(2), 0-gc(3):Nj+gc(4), 0-gc(5):Nk+gc(6)))
      allocate(self%torque_pressure(   0-gc(1):Ni+gc(2), 0-gc(3):Nj+gc(4), 0-gc(5):Nk+gc(6)))
      allocate(self%torque_viscous(    0-gc(1):Ni+gc(2), 0-gc(3):Nj+gc(4), 0-gc(5):Nk+gc(6)))
   endif
   endassociate
   endsubroutine alloc

   subroutine compute_aux(self, grd, icc, patch, RE, rFR2, zfs)
   !< Compute auxiliary varibles.
   class(block_rst_object), intent(inout)        :: self  !< Block data.
   type(block_grd_object),  intent(in)           :: grd   !< Block grd data.
   type(block_icc_object),  intent(in)           :: icc   !< Block icc data.
   integer(I4P),            intent(in), optional :: patch !< Patch boundary conditions.
   real(R8P),               intent(in), optional :: RE    !< Reynolds number.
   real(R8P),               intent(in), optional :: rFR2  !< 1/(Froude number)^2.
   real(R8P),               intent(in), optional :: zfs   !< Z quote of free surface.

   if ((self%has_lambda2).or.(self%has_qfactor).or.(self%has_liutex).or.(self%has_helicity).or.(self%has_vorticity))  then
      call self%compute_vorticity(grd)
   endif
   if (self%has_div2LT)  call self%compute_div2LT(grd)
   if (self%has_grad_p)  call self%compute_grad_p(grd)
   if (self%has_k_ratio) call self%compute_kinetic_energy_ratio(grd)
   if (present(patch)) then
      if (self%has_yplus)                   call self%compute_yplus(grd=grd, icc=icc, patch=patch, RE=RE)
      if (self%has_tau.or.self%has_div_tau) call self%compute_tau(  grd=grd, icc=icc, patch=patch, RE=RE)
      if (self%has_loads)                   call self%compute_loads(grd=grd, icc=icc, patch=patch, RE=RE, rFR2=rFR2, zfs=zfs)
   endif
   endsubroutine compute_aux

   subroutine compute_div2LT(self, grd)
   !< Compute double divergence of Lighthill tensor.
   class(block_rst_object), intent(inout) :: self          !< Block data.
   type(block_grd_object),  intent(in)    :: grd           !< Block grd data.
   real(R8P), allocatable                 :: G(:,:,:,:,:)  !< Gradient.
   type(vector), allocatable              :: divLT(:,:,:)  !< Divergence of LT.
   integer(I4P)                           :: i,j,k         !< Counter.

   associate(Ni=>self%Ni, Nj=>self%Nj, Nk=>self%Nk, gc=>self%gc,                 &
             NFiS=>grd%NFiS, NFjS=>grd%NFjS, NFkS=>grd%NFkS, volume=>grd%volume, &
             momentum=>self%momentum,div2LT=>self%div2LT)
   allocate(G(1:3,1:3,0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
   allocate(divLT(    1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))

   call grd%compute_gradient(var=momentum, gv=G, tensor=.true.)
   do k=0,Nk+1
      do j=0,Nj+1
         do i=0,Ni+1
            divLT(i,j,k)%x = G(1,1,i,j,k)+G(1,2,i,j,k)+G(1,3,i,j,k) ! dvar%x.var%x/dx + dvar%x.var%y/dy + dvar%x.var%z/dz
            divLT(i,j,k)%y = G(2,1,i,j,k)+G(2,2,i,j,k)+G(2,3,i,j,k) ! dvar%y.var%x/dx + dvar%y.var%y/dy + dvar%y.var%z/dz
            divLT(i,j,k)%z = G(3,1,i,j,k)+G(3,2,i,j,k)+G(3,3,i,j,k) ! dvar%z.var%x/dx + dvar%z.var%y/dy + dvar%z.var%z/dz
         enddo
      enddo
   enddo
   call grd%compute_gradient(var=divLT, gv=G)
   do k=0,Nk+1
      do j=0,Nj+1
         do i=0,Ni+1
            div2LT(i,j,k) = G(1,1,i,j,k)+G(2,2,i,j,k)+G(3,3,i,j,k)  ! dvar%x/dx + dvar%y/dy + dvar%z/dz
         enddo
      enddo
   enddo
   endassociate
   endsubroutine compute_div2LT

   subroutine compute_grad_p(self, grd)
   !< Compute pressure gradient.
   class(block_rst_object), intent(inout) :: self !< Block data.
   type(block_grd_object),  intent(in)    :: grd  !< Block grd data.

   call grd%compute_gradient(var=self%pressure, gv=self%grad_p)
   endsubroutine compute_grad_p

   subroutine compute_kinetic_energy_ratio(self, grd, vel0)
   !< Compute ratio of modelled kinetic energy over total kinetic energy.
   class(block_rst_object), intent(inout)        :: self         !< Block data.
   type(block_grd_object),  intent(in)           :: grd          !< Block grd data.
   type(vector),            intent(in), optional :: vel0         !< Undisturbed velocity.
   type(vector)                                  :: vel0_        !< Undisturbed velocity, local var.
   real(R8P), allocatable                        :: l(:,:,:,:,:) !< Laplacian tensor of velocity.
   type(vector), allocatable                     :: vf(:,:,:)    !< Taylor expansion of filtered velocity.
   real(R8P), allocatable                        :: div(:,:,:)   !< Divergence of filtered velocity.
   real(R8P)                                     :: kmod         !< Modelled kinetic energy.
   real(R8P)                                     :: kres         !< Resolved kinetic energy.
   integer(I4P)                                  :: i,j,k        !< Counter.
   real(R8P)                                     :: delta        !< Coefficient.
   real(R8P), parameter                          :: EPS12=1d-12  !< Tolerances.

   vel0_ = 0._R8P ; if (present(vel0)) vel0_ = vel0
   associate(Ni=>grd%Ni, Nj=>grd%Nj, Nk=>grd%Nk, gc=>grd%gc, volume=>grd%volume, &
             momentum=>self%momentum,k_ratio=>self%k_ratio)
   allocate(vf(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
   call grd%compute_laplacian_tensor(var=momentum, lv=l)
   ! compute Taylor expansion of filtered velocity
   do k=0,Nk+1
      do j=0,Nj+1
         do i=0,Ni+1
            delta = max(EPS12, volume(i,j,k)**4._R8P/3._R8P)/12._R8P
            vf(i,j,k)%x = momentum(i,j,k)%x + delta * (l(1,1,i,j,k) + l(1,2,i,j,k) + l(1,3,i,j,k)) ! d2u/dx2 + d2u/dy2 + d2u/dz2
            vf(i,j,k)%y = momentum(i,j,k)%y + delta * (l(2,1,i,j,k) + l(2,2,i,j,k) + l(2,3,i,j,k)) ! d2v/dx2 + d2v/dy2 + d2v/dz2
            vf(i,j,k)%z = momentum(i,j,k)%z + delta * (l(3,1,i,j,k) + l(3,2,i,j,k) + l(3,3,i,j,k)) ! d2w/dx2 + d2w/dy2 + d2w/dz2
         enddo
      enddo
   enddo
   call grd%compute_divergence(var=vf, div=div)
   ! compute kinetic energy ratio
   do k=0,Nk+1
      do j=0,Nj+1
         do i=0,Ni+1
            delta = max(EPS12, volume(i,j,k)**4._R8P/3._R8P)/24._R8P
            kmod = delta * div(i,j,k) * div(i,j,k)
            kmod = div(i,j,k) * div(i,j,k)
            kres = 0.5_R8P * ((vf(i,j,k)%x - vel0_%x) * (vf(i,j,k)%x - vel0_%x) + &
                              (vf(i,j,k)%y - vel0_%y) * (vf(i,j,k)%y - vel0_%y) + &
                              (vf(i,j,k)%z - vel0_%z) * (vf(i,j,k)%z - vel0_%z))
            k_ratio(i,j,k) = kmod / (kmod + kres)
         enddo
      enddo
   enddo
   endassociate
   endsubroutine compute_kinetic_energy_ratio

   subroutine compute_loads(self, grd, icc, patch, RE, rFR2, zfs)
   !< Compute loads (forces and torques) on patches.
   class(block_rst_object), intent(inout) :: self    !< Block data.
   type(block_grd_object),  intent(in)    :: grd     !< Block grd data.
   type(block_icc_object),  intent(in)    :: icc     !< Block icc data.
   integer(I4P),            intent(in)    :: patch   !< Patch bc to be found.
   real(R8P),               intent(in)    :: RE      !< Reynolds number.
   real(R8P),               intent(in)    :: rFR2    !< 1/(Froude number)^2.
   real(R8P),               intent(in)    :: zfs     !< Z quote of free surface.
   type(vector)                           :: NdS     !< Normal "tilde" for viscous part of forces (or distance for y+).
   type(vector)                           :: fco     !< Pacth center coordinates.
   integer(I4P)                           :: i,j,k,p !< Counter.

   if (.not.allocated(grd%patches_extents)) return
   do p=1, size(grd%patches_extents, dim=1)
      associate(face=>grd%patches_extents(p,0),                                   &
                ci1=>grd%patches_extents(p,1 ), ci2=>grd%patches_extents(p,2 ),   &
                cj1=>grd%patches_extents(p,3 ), cj2=>grd%patches_extents(p,4 ),   &
                ck1=>grd%patches_extents(p,5 ), ck2=>grd%patches_extents(p,6 ),   &
                ni1=>grd%patches_extents(p,7 ), ni2=>grd%patches_extents(p,8 ),   &
                nj1=>grd%patches_extents(p,9 ), nj2=>grd%patches_extents(p,10),   &
                nk1=>grd%patches_extents(p,11), nk2=>grd%patches_extents(p,12),   &
                nodes=>grd%nodes, NFiS=>grd%NFiS, NFjS=>grd%NFjS, NFkS=>grd%NFkS, &
                NFi=>grd%NFi, NFj=>grd%NFj, NFk=>grd%NFk,                         &
                volume=>grd%volume,                                               &
                iicc=>icc%icc,ticc=>icc%tcc,                                      &
                is_level_set=>self%is_level_set, f0=>self%level_set_zero,         &
                momentum=>self%momentum, pressure=>self%pressure,                 &
                force_hydrostatic=>self%force_hydrostatic,                        &
                force_pressure=>self%force_pressure,                              &
                force_viscous=>self%force_viscous,                                &
                torque_hydrostatic=>self%torque_hydrostatic,                      &
                torque_pressure=>self%torque_pressure,                            &
                torque_viscous=>self%torque_viscous)
      select case(face)
      case(1,2)
         do k=ck1,ck2
            do j=cj1,cj2
               ! patch center coordinates
               fco = 0.25_R8P*(nodes(ni1,j,k)+nodes(ni1,j-1,k)+nodes(ni1,j,k-1)+nodes(ni1,j-1,k-1))
               ! check if this is an active cell
               if (ticc(ni1-1+face,j,k)/=patch.or.iicc(ci1,j,k)/=0) cycle
               ! hydrostatic part
               if (is_level_set) then
                  if (f0(ni1-1+face,j,k)>0._R8P) cycle
                  force_hydrostatic(ni1,j,k) = (fco%z - zfs)*rFR2*NFiS(ni1,j,k)
               endif
               ! pressure part of forces
               force_pressure(ni1,j,k) = -(1.5_R8P*pressure(ci1,j,k) - 0.5_R8P*pressure(ci1+3-2*face,j,k))*NFiS(ni1,j,k)
               ! viscous part of forces
               NdS = (2._R8P*NFiS(ni1,j,k) + NFiS(ni1+1,j,k) + NFiS(ni1-1,j,k))/(2._R8P*(volume(ni1,j,k)+volume(ni1+1,j,k)))
               force_viscous(ni1,j,k) = (momentum(ci1-1+face,j,k) - momentum(ci1-2+face,j,k))*(NFi(ni1,j,k).dot.NdS)/RE
               if (face==2) then
                  force_hydrostatic(ni1,j,k) = -force_hydrostatic(ni1,j,k)
                  force_pressure(   ni1,j,k) = -force_pressure(   ni1,j,k)
                  force_viscous(    ni1,j,k) = -force_viscous(    ni1,j,k)
               endif
               ! torques
               torque_hydrostatic(ni1,j,k) = fco.cross.force_hydrostatic(ni1,j,k)
               torque_pressure(   ni1,j,k) = fco.cross.force_pressure(   ni1,j,k)
               torque_viscous(    ni1,j,k) = fco.cross.force_viscous(    ni1,j,k)
            enddo
         enddo
      case(3,4)
         do k=ck1,ck2
            do i=ci1,ci2
               ! patch center coordinates
               fco = 0.25_R8P*(nodes(i,nj1,k)+nodes(i-1,nj1,k)+nodes(i,nj1,k-1)+nodes(i-1,nj1,k-1))
               ! check if this is an active cell
               if (ticc(i,nj1-1+(face-2),k)/=patch.or.iicc(i,cj1,k)/=0) cycle
               ! hydrostatic part
               if (is_level_set) then
                 if (f0(i,nj1-1+(face-2),k)>0._R8P) cycle
                 force_hydrostatic(i,nj1,k) = (fco%z - zfs)*rFR2*NFjS(i,nj1,k)
               endif
               ! pressure part of forces
               force_pressure(i,nj1,k) = -(1.5_R_P*pressure(i,cj1,k) - 0.5_R_P*pressure(i,cj1+3-2*(face-2),k))*NFjS(i,nj1,k)
               ! viscous part of forces
               NdS = (2._R8P*NFjS(i,nj1,k) + NFjS(i,nj1+1,k) + NFjS(i,nj1-1,k))/(2._R8P*(volume(i,nj1,k)+volume(i,nj1+1,k)))
               force_viscous(i,nj1,k) = (momentum(i,cj1-1+(face-2),k) - momentum(i,cj1-2+(face-2),k))*(NFj(i,nj1,k).dot.NdS)/RE
               if (face==4) then
                  force_hydrostatic(i,nj1,k) = -force_hydrostatic(i,nj1,k)
                  force_pressure(   i,nj1,k) = -force_pressure(   i,nj1,k)
                  force_viscous(    i,nj1,k) = -force_viscous(    i,nj1,k)
               endif
               ! torques
               torque_hydrostatic(i,nj1,k) = fco.cross.force_hydrostatic(i,nj1,k)
               torque_pressure(   i,nj1,k) = fco.cross.force_pressure(   i,nj1,k)
               torque_viscous(    i,nj1,k) = fco.cross.force_viscous(    i,nj1,k)
            enddo
         enddo
      case(5,6)
         do j=cj1,cj2
            do i=ci1,ci2
               ! patch center coordinates
               fco = 0.25_R8P*(nodes(i,j,nk1)+nodes(i-1,j,nk1)+nodes(i,j-1,nk1)+nodes(i-1,j-1,nk1))
               ! check if this is an active cell
               if (ticc(i,j,nk1-1+(face-4))/=patch.or.iicc(i,j,ck1)/=0) cycle
               ! hydrostatic part
               if (is_level_set) then
                  if (f0(i,j,nk1-1+(face-4))>0._R8P) cycle
                  force_hydrostatic(i,j,nk1) = (fco%z - zfs)*rFR2*NFkS(i,j,nk1)
               endif
               ! pressure part of forces
               force_pressure(i,j,nk1) = -(1.5_R_P*pressure(i,j,ck1) - 0.5_R_P*pressure(i,j,ck1+3-2*(face-4)))*NFkS(i,j,nk1)
               ! viscous part of forces
               NdS = (2._R8P*NFkS(i,j,nk1) + NFkS(i,j,nk1+1) + NFkS(i,j,nk1-1))/(2._R8P*(volume(i,j,nk1)+volume(i,j,nk1+1)))
               force_viscous(i,j,nk1) = (momentum(i,j,ck1-1+(face-4)) - momentum(i,j,ck1-2+(face-4)))*(NFk(i,j,nk1).dot.NdS)/RE
               if (face==6) then
                  force_hydrostatic(i,j,nk1) = -force_hydrostatic(i,j,nk1)
                  force_pressure(   i,j,nk1) = -force_pressure(   i,j,nk1)
                  force_viscous(    i,j,nk1) = -force_viscous(    i,j,nk1)
               endif
               ! torques
               torque_hydrostatic(i,j,nk1) = fco.cross.force_hydrostatic(i,j,nk1)
               torque_pressure(   i,j,nk1) = fco.cross.force_pressure(   i,j,nk1)
               torque_viscous(    i,j,nk1) = fco.cross.force_viscous(    i,j,nk1)
            enddo
         enddo
      endselect
      endassociate
   enddo
   endsubroutine compute_loads

   subroutine compute_tau(self, grd, icc, patch, RE)
   !< Compute tau wall and its divergence.
   class(block_rst_object), intent(inout) :: self                                !< Block data.
   type(block_grd_object),  intent(in)    :: grd                                 !< Block grd data.
   type(block_icc_object),  intent(in)    :: icc                                 !< Block icc data.
   integer(I4P),            intent(in)    :: patch                               !< Patch bc to be found.
   real(R8P),               intent(in)    :: RE                                  !< Reynolds number.
   type(vector)                           :: NdS                                 !< Normal "tilde" for viscous part of forces.
   integer(I4P)                           :: i,j,k,p                             !< Counter.
   type(vector)                           :: cn,cs,cc,cw,ce                      !< Cell centers.
   type(vector)                           :: e_zet_n,e_zet_s,e_eta_e,e_eta_w     !< Curvilinear corrdinate versors.
   type(vector)                           :: tw_n,tw_s,tw_e,tw_w                 !< Cell faces value of the ambiguos vector.
   real(R8P)                              :: tw_zet_n,tw_zet_s,tw_eta_e,tw_eta_w !< Curvilinear components of the ambiguos vector.
   real(R8P)                              :: ds                                  !< Face area.

   if (.not.allocated(grd%patches_extents)) return
   do p=1, size(grd%patches_extents, dim=1)
      associate(face=>grd%patches_extents(p,0),                                   &
                ci1=>grd%patches_extents(p,1 ), ci2=>grd%patches_extents(p,2 ),   &
                cj1=>grd%patches_extents(p,3 ), cj2=>grd%patches_extents(p,4 ),   &
                ck1=>grd%patches_extents(p,5 ), ck2=>grd%patches_extents(p,6 ),   &
                ni1=>grd%patches_extents(p,7 ), ni2=>grd%patches_extents(p,8 ),   &
                nj1=>grd%patches_extents(p,9 ), nj2=>grd%patches_extents(p,10),   &
                nk1=>grd%patches_extents(p,11), nk2=>grd%patches_extents(p,12),   &
                nodes=>grd%nodes, NFiS=>grd%NFiS, NFjS=>grd%NFjS, NFkS=>grd%NFkS, &
                NFi=>grd%NFi, NFj=>grd%NFj, NFk=>grd%NFk,                         &
                volume=>grd%volume,                                               &
                iicc=>icc%icc,ticc=>icc%tcc,                                      &
                is_level_set=>self%is_level_set,                                  &
                momentum=>self%momentum, f0=>self%level_set_zero, &
                tau=>self%tau, div_tau=>self%div_tau)
      select case(face)
      case(1,2)
         do k=ck1,ck2
            do j=cj1,cj2
               ! check if this is an active cell
               if (ticc(ni1-1+face,j,k)/=patch) cycle
               if (iicc(ci1,j,k)/=0) cycle
               if (is_level_set) then
                  if (f0(ni1-1+face,j,k)>0._R8P) cycle
               endif
               ! viscous part of forces
               NdS = (2._R8P*NFiS(ni1,j,k) + NFiS(ni1+1,j,k) + NFiS(ni1-1,j,k))/(2._R8P*(volume(ni1,j,k)+volume(ni1+1,j,k)))
               tau(ci1,j,k) = (momentum(ci1-1+face,j,k) - momentum(ci1-2+face,j,k))*(NFi(ni1,j,k).dot.NdS)/RE
               tau(ci1,j,k) = (tau(ci1,j,k).ortho.(NFi(ni1,j,k)))
               if (face==2) tau(ci1,j,k) = -tau(ci1,j,k)
            enddo
         enddo
         ! compute divergence
         do k=ck1,ck2
            do j=cj1,cj2
               ! check if this is an active cell
               if (ticc(ni1-1+face,j,k)/=patch) cycle
               if (iicc(ci1,j,k)/=0) cycle
               if (is_level_set) then
                  if (f0(ni1-1+face,j,k)>0._R8P) cycle
               endif
               ! Compute face center coordinates
               cn = 0.25_R8P*(nodes(ni1,j  ,k  ) + nodes(ni1,j  ,k+1) + nodes(ni1,j-1,k+1) + nodes(ni1,j-1,k  ))
               cs = 0.25_R8P*(nodes(ni1,j  ,k-2) + nodes(ni1,j  ,k-1) + nodes(ni1,j-1,k-1) + nodes(ni1,j-1,k-2))
               cc = 0.25_R8P*(nodes(ni1,j  ,k-1) + nodes(ni1,j  ,k  ) + nodes(ni1,j-1,k  ) + nodes(ni1,j-1,k-1))
               cw = 0.25_R8P*(nodes(ni1,j-1,k-1) + nodes(ni1,j-1,k  ) + nodes(ni1,j-2,k  ) + nodes(ni1,j-2,k-1))
               ce = 0.25_R8P*(nodes(ni1,j+1,k-1) + nodes(ni1,j+1,k  ) + nodes(ni1,j  ,k  ) + nodes(ni1,j  ,k-1))

               e_zet_n = (cn-cc)/normL2(cn-cc)
               e_zet_s = (cc-cs)/normL2(cc-cs)
               e_eta_e = (ce-cc)/normL2(ce-cc)
               e_eta_w = (cc-cw)/normL2(cc-cw)

               tw_n = 0.5_R8P*(tau(ci1,j  ,k+1) - tau(ci1,j  ,k  ))
               tw_s = 0.5_R8P*(tau(ci1,j  ,k  ) - tau(ci1,j  ,k-1))
               tw_e = 0.5_R8P*(tau(ci1,j+1,k  ) - tau(ci1,j  ,k  ))
               tw_w = 0.5_R8P*(tau(ci1,j  ,k  ) - tau(ci1,j-1,k  ))

               tw_zet_n = tw_n.dot.e_zet_n
               tw_zet_s = tw_s.dot.e_zet_s
               tw_eta_e = tw_e.dot.e_eta_e
               tw_eta_w = tw_w.dot.e_eta_w

               div_tau(ci1,j,k) = (tw_zet_n - tw_zet_s) + (tw_eta_e - tw_eta_w)
            enddo
         enddo
      case(3,4)
         do k=ck1,ck2
            do i=ci1,ci2
               ! check if this is an active cell
               if (ticc(i,nj1-1+(face-2),k)/=patch) cycle
               if (iicc(i,cj1,k)/=0) cycle
               if (is_level_set) then
                  if (f0(i,nj1-1+(face-2),k)>0._R8P) cycle
               endif
               ! viscous part of forces
               NdS = (2._R8P*NFjS(i,nj1,k) + NFjS(i,nj1+1,k) + NFjS(i,nj1-1,k))/(2._R8P*(volume(i,nj1,k)+volume(i,nj1+1,k)))
               tau(i,cj1,k) = (momentum(i,cj1-1+(face-2),k) - momentum(i,cj1-2+(face-2),k))*(NFj(i,nj1,k).dot.NdS)/RE
               tau(i,cj1,k) = (tau(i,cj1,k).ortho.(NFj(i,nj1,k)))
               if (face==4) tau(i,cj1,k) = -tau(i,cj1,k)
            enddo
         enddo
         ! compute divergence
         do k=ck1,ck2
            do i=ci1,ci2
               ! check if this is an active cell
               if (ticc(i,nj1-1+(face-2),k)/=patch) cycle
               if (iicc(i,cj1,k)/=0) cycle
               if (is_level_set) then
                  if (f0(i,nj1-1+(face-2),k)>0._R8P) cycle
               endif
               ! compute face center coordinates
               cn = 0.25_R8P*(nodes(i  ,nj1,k  ) + nodes(i  ,nj1,k+1) + nodes(i-1,nj1,k+1) + nodes(i-1,nj1,k  ))
               cs = 0.25_R8P*(nodes(i  ,nj1,k-2) + nodes(i  ,nj1,k-1) + nodes(i-1,nj1,k-1) + nodes(i-1,nj1,k-2))
               cc = 0.25_R8P*(nodes(i  ,nj1,k-1) + nodes(i  ,nj1,k  ) + nodes(i-1,nj1,k  ) + nodes(i-1,nj1,k-1))
               cw = 0.25_R8P*(nodes(i-1,nj1,k-1) + nodes(i-1,nj1,k  ) + nodes(i-2,nj1,k  ) + nodes(i-2,nj1,k-1))
               ce = 0.25_R8P*(nodes(i+1,nj1,k-1) + nodes(i+1,nj1,k  ) + nodes(i  ,nj1,k  ) + nodes(i  ,nj1,k-1))

               e_zet_n = (cn-cc)/normL2(cn-cc)
               e_zet_s = (cc-cs)/normL2(cc-cs)
               e_eta_e = (ce-cc)/normL2(ce-cc)
               e_eta_w = (cc-cw)/normL2(cc-cw)

               tw_n = 0.5_R8P*(tau(i  ,cj1,k+1) - tau(i  ,cj1,k  ))
               tw_s = 0.5_R8P*(tau(i  ,cj1,k  ) - tau(i  ,cj1,k-1))
               tw_e = 0.5_R8P*(tau(i+1,cj1,k  ) - tau(i  ,cj1,k  ))
               tw_w = 0.5_R8P*(tau(i  ,cj1,k  ) - tau(i-1,cj1,k  ))

               tw_zet_n = tw_n.dot.e_zet_n
               tw_zet_s = tw_s.dot.e_zet_s
               tw_eta_e = tw_e.dot.e_eta_e
               tw_eta_w = tw_w.dot.e_eta_w

               div_tau(i,cj1,k) = (tw_zet_n - tw_zet_s) + (tw_eta_e - tw_eta_w)
            enddo
         enddo
      case(5,6)
         do j=cj1,cj2
            do i=ci1,ci2
               ! check if this is an active cell
               if (ticc(i,j,nk1-1+(face-4))/=patch) cycle
               if (iicc(i,j,ck1)/=0) cycle
               if (is_level_set) then
                  if (f0(i,j,nk1-1+(face-4))>0._R8P) cycle
               endif
               ! viscous part of forces
               NdS = (2._R8P*NFkS(i,j,nk1) + NFkS(i,j,nk1+1) + NFkS(i,j,nk1-1))/(2._R8P*(volume(i,j,nk1)+volume(i,j,nk1+1)))
               tau(i,j,ck1) = (momentum(i,j,ck1-1+(face-4)) - momentum(i,j,ck1-2+(face-4)))*(NFk(i,j,nk1).dot.NdS)/RE
               tau(i,j,ck1) = (tau(i,j,ck1).ortho.(NFk(i,j,nk1)))
               if (face==6) tau(i,j,ck1) = -tau(i,j,ck1)
            enddo
         enddo
         ! compute of divergence
         do k=ck1,ck2
            do i=ci1,ci2
               ! check if this is an active cell
               if (ticc(i,j,nk1-1+(face-4))/=patch) cycle
               if (iicc(i,j,ck1)/=0) cycle
               if (is_level_set) then
                  if (f0(i,j,nk1-1+(face-4))>0._R8P) cycle
               endif
               ! compute face center coordinates
               cn = 0.25_R8P*(nodes(i  ,j  ,nk1) + nodes(i  ,j+1,nk1) + nodes(i-1,j+1,nk1) + nodes(i-1,j  ,nk1))
               cs = 0.25_R8P*(nodes(i  ,j-2,nk1) + nodes(i  ,j-1,nk1) + nodes(i-1,j-1,nk1) + nodes(i-1,j-2,nk1))
               cc = 0.25_R8P*(nodes(i  ,j-1,nk1) + nodes(i  ,j  ,nk1) + nodes(i-1,j  ,nk1) + nodes(i-1,j-1,nk1))
               cw = 0.25_R8P*(nodes(i-1,j-1,nk1) + nodes(i-1,j  ,nk1) + nodes(i-2,j  ,nk1) + nodes(i-2,j-1,nk1))
               ce = 0.25_R8P*(nodes(i+1,j-1,nk1) + nodes(i+1,j  ,nk1) + nodes(i  ,j  ,nk1) + nodes(i  ,j-1,nk1))

               e_zet_n = (cn-cc)/normL2(cn-cc)
               e_zet_s = (cc-cs)/normL2(cc-cs)
               e_eta_e = (ce-cc)/normL2(ce-cc)
               e_eta_w = (cc-cw)/normL2(cc-cw)

               tw_n = 0.5_R8P*(tau(i  ,j+1,ck1) - tau(i  ,j  ,ck1))
               tw_s = 0.5_R8P*(tau(i  ,j  ,ck1) - tau(i  ,j-1,ck1))
               tw_e = 0.5_R8P*(tau(i+1,j  ,ck1) - tau(i  ,j  ,ck1))
               tw_w = 0.5_R8P*(tau(i  ,j  ,ck1) - tau(i-1,j  ,ck1))

               tw_zet_n = tw_n.dot.e_zet_n
               tw_zet_s = tw_s.dot.e_zet_s
               tw_eta_e = tw_e.dot.e_eta_e
               tw_eta_w = tw_w.dot.e_eta_w

               div_tau(i,j,ck1) = (tw_zet_n - tw_zet_s) + (tw_eta_e - tw_eta_w)
            enddo
         enddo
      endselect
      endassociate
   enddo
   endsubroutine compute_tau

   subroutine compute_vorticity(self, grd)
   !< Compute vorticity related varibles.
   class(block_rst_object), intent(inout) :: self                              !< Block data.
   type(block_grd_object),  intent(in)    :: grd                               !< Block grd data.
   real(R8P), allocatable                 :: gv(:,:,:,:,:)                     !< Gradient of velocity.
   real(R8P)                              :: c(0:2),emin,emax,eval,fval,mu     !< Dummy reals.
   type(vector)                           :: um                                !< Dummy vector.
   real(R8P), dimension(3,3)              :: SO,S,O                            !< Matrices.
   integer(I4P)                           :: i,j,k,ii,jj,kk,iter               !< Counter.
   real(R8P), parameter                   :: EPS6=1d-6, EPS9=1d-9, EPS12=1d-12 !< Tolerances.
   ! Variables for liutex method
   real(R8P)                              :: norm_r_star                       !< Dummy variable for liutex method.
   real(R8P)                              :: w_dot_r                           !< Dummy variable for liutex method.
   real(R8P), dimension(3)                :: eig_vec_r, r_star                 !< Dummy variable for liutex method.
   real(R8P), dimension(3,3)              :: r_liutex                          !< Dummy variable for liutex method.
   real(R8P)                              :: aa, b                             !< Dummy variable for liutex method.
   real(R8P)                              :: delta, delta1, delta2, delta3     !< Dummy variable for liutex method.
   real(R8P)                              :: eig3r                             !< Dummy variable for liutex method.
   complex(R8P)                           :: eig1c, eig2c                      !< Dummy variable for liutex method.
   associate(Ni=>grd%Ni,Nj=>grd%Nj,Nk=>grd%Nk,gc=>grd%gc, &
             momentum=>self%momentum,lambda2=>self%lambda2,qfactor=>self%qfactor,helicity=>self%helicity,vorticity=>self%vorticity)
   ! print '(A)', 'compute vorticity'
   call grd%compute_gradient(var=momentum, gv=gv)
   if (self%has_helicity.or.self%has_vorticity) then
      ! compute vorticity vector
      do k=0,Nk+1
         do j=0,Nj+1
            do i=0,Ni+1
               vorticity%x =  (gv(3,2,i,j,k) - gv(2,3,i,j,k))
               vorticity%y = -(gv(1,3,i,j,k) - gv(3,1,i,j,k))
               vorticity%z =  (gv(2,1,i,j,k) - gv(1,2,i,j,k))
            enddo
         enddo
      enddo
   endif
   if (self%has_helicity) then
      ! print '(A)', 'compute helicity'
      ! compute helicity
      do k=0,Nk+1
         do j=0,Nj+1
            do i=0,Ni+1
               helicity(i,j,k) = momentum(i,j,k).dot.vorticity(i,j,k)/(max(EPS12,normL2(momentum(i,j,k))*normL2(vorticity(i,j,k))))
            enddo
         enddo
      enddo
   endif
   if (self%has_lambda2.or.self%has_qfactor.or.self%has_liutex) then
      ! if (self%has_lambda2) print '(A)', 'compute lambda2'
      ! if (self%has_qfactor) print '(A)', 'compute qfactor'
      do k=0,Nk+1
         do j=0,Nj+1
            do i=0,Ni+1
               ! tensor S^2 + O^2
               S = 0._R8P
               O = 0._R8P
               do kk=1,3
                  do jj=1,3
                     S(jj,kk) = 0.5_R8P*(gv(jj,kk,i,j,k)+gv(kk,jj,i,j,k))
                     O(jj,kk) = 0.5_R8P*(gv(jj,kk,i,j,k)-gv(kk,jj,i,j,k))
                  enddo
               enddo
               SO = matmul(S,S) + matmul(O,O)

               if (self%has_qfactor) then
                  qfactor(i,j,k) = 0.5_R8P*(dot_product(O(1,:),O(1,:)) + &
                                            dot_product(O(2,:),O(2,:)) + &
                                            dot_product(O(3,:),O(3,:)) - &
                                           (dot_product(S(1,:),S(1,:)) + &
                                            dot_product(S(2,:),S(2,:)) + &
                                            dot_product(S(3,:),S(3,:))))
               endif

               if (self%has_lambda2) then
                  ! coefficients of characterist polynomial: lamda^3 + c(2)*lamda^2 + c(1)*lamda + c(0) = 0
                  c(2)          = -(SO(1,1) + SO(2,2) + SO(3,3))
                  c(1)          =   SO(1,1)*SO(2,2) + SO(1,1)*SO(3,3) + SO(2,2)*SO(3,3) - SO(2,3)**2 - SO(1,3)**2 - SO(1,2)**2
                  c(0)          =   SO(1,1)*SO(2,3)**2 + SO(2,2)*SO(1,3)**2 + SO(3,3)*SO(1,2)**2 - 2._R8P*SO(2,3)*SO(1,3)*SO(1,2) &
                                   -SO(1,1)*SO(2,2)*SO(3,3)
                  ! compute second eigenvalue of characteristic polynomial
                  mu   = sqrt(c(2)**2 - 3._R8P*c(1))
                  emin = (-c(2)-mu)/3._R8P
                  emax = (-c(2)+mu)/3._R8P
                  do iter=1,100
                     eval = 0.5_R8P*(emin+emax)
                     fval = eval**3 + c(2)*eval**2 + c(1)*eval + c(0)
                     if (fval<0._R8P) then
                        emax = eval
                     else
                        emin = eval
                     end if
                     if (abs(fval)<eps9 .and.((emax-emin)/eval)<eps6) exit
                  end do
                  lambda2(i,j,k) = eval
               endif

               if (self%has_liutex) then
                  tt            = matmul(gv(:,:,i,j,k),gv(:,:,i,j,k))
                  c(0)          = -(gv(1,1,i,j,k)*(gv(2,2,i,j,k)*gv(3,3,i,j,k)-gv(2,3,i,j,k)*gv(3,2,i,j,k)) &
                                  - gv(1,2,i,j,k)*(gv(2,1,i,j,k)*gv(3,3,i,j,k)-gv(2,3,i,j,k)*gv(3,1,i,j,k)) &
                                  + gv(1,3,i,j,k)*(gv(2,1,i,j,k)*gv(3,2,i,j,k)-gv(2,2,i,j,k)*gv(3,1,i,j,k)))
                  c(2)          = -(gv(1,1,i,j,k) + gv(2,2,i,j,k) + gv(3,3,i,j,k))
                  c(1)          = -0.5_R8P*(tt(1,1) + tt(2,2) + tt(3,3) - (c(2))**2)
                  delta         = 18._R8P*c(1)*c(2)*c(0) - 4._R8P*c(2)**3 * c(0) + c(1)**2 * c(2)**2 - 4._R8P*c(1)**3 - 27._R8P*c(0)**2
                  delta         = -delta/108._R8P
                  liutex(i,j,k) = 0._R8P  
                  if (delta>0._R8P) then
                     mu = c(2)**2 - 3._R8P*c(1)/9._R8P
                     t  = (2._R8P*c(2)**3 - 9._R8P*c(1)*c(2) + 27._R8P*c(0)) / 54._R8P
                     aa = -sign(1._R8P, t) * ( abs(t) + sqrt(delta) )**(1._R8P/3._R8P)
                     if (aa==0._R8P) then
                         b = 0._R8P
                     else
                         b = mu / aa
                     endif
                     eig1c  = cmplx(-0.5_R8P*(aa+b) - c(2)/3._R8P, 0.5_R8P*sqrt(3._R8P)*(aa-b))
                     eig2c  = cmplx(real(eig1c), -aimag(eig1c)) 
                     eig3r  = aa + b - c(2)/3._R8P
                     delta1 = (gv(1,1,i,j,k) - eig3r) * (gv(2,2,i,j,k) - eig3r) - gv(2,1,i,j,k)*gv(1,2,i,j,k)
                     delta2 = (gv(2,2,i,j,k) - eig3r) * (gv(3,3,i,j,k) - eig3r) - gv(2,3,i,j,k)*gv(3,2,i,j,k)
                     delta3 = (gv(1,1,i,j,k) - eig3r) * (gv(3,3,i,j,k) - eig3r) - gv(1,3,i,j,k)*gv(3,1,i,j,k)

                     if (delta1 == 0._R8P .and. delta2 == 0._R8P .and. delta3 == 0._R8P) then
                        stop
                     end if
                  
                     if (abs(delta1) >= abs(delta2) .and. abs(delta1) >= abs(delta3)) then
                     
                        r_star(1) = (-(gv(2,2,i,j,k)-eig3r)*a(1,3,i,j,k) + gv(1,2,i,j,k)*gv(2,3,i,j,k)) / delta1
                        r_star(2) = (  gv(2,1,i,j,k)*gv(1,3,i,j,k) - (gv(1,1,i,j,k)-eig3r)*gv(2,3,i,j,k)) / delta1
                        r_star(3) = 1._R8P
                     
                     else if (abs(delta2) >= abs(delta1) .and. abs(delta2) >= abs(delta3)) then
                     
                        r_star(1) = 1._R8P
                        r_star(2) = (-(gv(3,3,i,j,k)-eig3r)*gv(2,1,i,j,k) + gv(2,3,i,j,k)*gv(3,1,i,j,k))/ delta2
                        r_star(3) = (  gv(3,2,i,j,k)*gv(2,1,i,j,k)-(gv(2,2,i,j,k)-eig3r)*gv(3,1,i,j,k)) / delta2
                     
                     else if (abs(delta3) >= abs(delta1) .and. abs(delta3) >= abs(delta2)) then
                     
                        r_star(1) = (-(gv(3,3)-eig3r)*gv(1,2)+gv(1,3,i,j,k)*gv(3,2,i,j,k))/delta3
                        r_star(2) = 1._R8P
                        r_star(3) = ( gv(3,1,i,j,k)*gv(1,2,i,j,k) - (gv(1,1,i,j,k)-eig3r)*gv(3,2,i,j,k))/delta3
                     
                     else
                         stop
                     end if

                     !eig_vec_r(1) = (gv(1,2,i,j,k)* gv(2,3,i,j,k)-gv(1,1,i,j,k)*(gv(2,2,i,j,k)-eig3r))
                     !eig_vec_r(2) = (gv(1,3,i,j,k)* gv(2,1,i,j,k)-gv(1,1,i,j,k)* gv(2,3,i,j,k))
                     !eig_vec_r(3) = (gv(1,1,i,j,k)*(gv(2,2,i,j,k)-eig3r)-gv(1,2,i,j,k)*gv(2,1,i,j,k))
                     !eig_vec_r    = eig_vec_r / norm2(eig_vec_r)
                     norm_r_star = sqrt(r_star(1)**2 + r_star(2)**2 + r_star(3)**2)
                     eig_vec_r(1) = r_star(1) / norm_r_star
                     eig_vec_r(2) = r_star(2) / norm_r_star
                     eig_vec_r(3) = r_star(3) / norm_r_star
                     aa           = dot_product(vorticity(i,j,k), eig_vec_r)
                     if (aa<=0._R8P) then
                        r_liutex = - eig_vec_r
                     else
                        r_liutex = + eig_vec_r
                     endif
                     w_dot_r       = dot_product(vorticity(i,j,k), r_liutex)
                     liutex(i,j,k) = w_dot_r - sqrt(w_dot_r**2 - 4._R8P*aimag(eig1c)**2)
                  endif
               endif
            enddo
         enddo
      enddo
   endif
   endassociate
   endsubroutine compute_vorticity

   subroutine compute_yplus(self, grd, icc, patch, RE)
   !< Compute yplus on patches.
   class(block_rst_object), intent(inout) :: self    !< Block data.
   type(block_grd_object),  intent(in)    :: grd     !< Block grd data.
   type(block_icc_object),  intent(in)    :: icc     !< Block icc data.
   integer(I4P),            intent(in)    :: patch   !< Patch bc to be found.
   real(R8P),               intent(in)    :: RE      !< Reynolds number.
   type(vector)                           :: NdS     !< Normal "tilde" for viscous part of forces (or distance for y+).
   integer(I4P)                           :: i,j,k,p !< Counter.
   real(R8P)                              :: twall   !< Shear stress.
   real(R8P)                              :: utau    !< Friction velocity.
   real(R8P)                              :: ds      !< Face area.
   integer(I4P)                           :: i1,i2   !< Counter.
   integer(I4P)                           :: j1,j2   !< Counter.
   integer(I4P)                           :: k1,k2   !< Counter.

   if (.not.allocated(grd%patches_extents)) return
   do p=1, size(grd%patches_extents, dim=1)
      associate(face=>grd%patches_extents(p,0),                                   &
                ci1=>grd%patches_extents(p,1 ), ci2=>grd%patches_extents(p,2 ),   &
                cj1=>grd%patches_extents(p,3 ), cj2=>grd%patches_extents(p,4 ),   &
                ck1=>grd%patches_extents(p,5 ), ck2=>grd%patches_extents(p,6 ),   &
                ni1=>grd%patches_extents(p,7 ), ni2=>grd%patches_extents(p,8 ),   &
                nj1=>grd%patches_extents(p,9 ), nj2=>grd%patches_extents(p,10),   &
                nk1=>grd%patches_extents(p,11), nk2=>grd%patches_extents(p,12),   &
                nodes=>grd%nodes, NFiS=>grd%NFiS, NFjS=>grd%NFjS, NFkS=>grd%NFkS, &
                iicc=>icc%icc,ticc=>icc%tcc,                                      &
                is_level_set=>self%is_level_set,                                  &
                momentum=>self%momentum, f0=>self%level_set_zero, yplus=>self%yplus)
      select case(face)
      case(1,2)
         !$omp parallel default(none) &
         !$omp private(i,j,k,NdS)     &
         !$omp shared(ni1,ni2,nj1,nj2,nk1,nk2,ci1,ci2,cj1,cj2,ck1,ck2,face,patch,nodes,NFiS,iicc,ticc,momentum,f0,RE,yplus)
         !$omp do
         i1=ni1; i2=i1+1
         if (face==2) i2=i1-1
         do k=ck1,ck2
            do j=cj1,cj2
               ! checking if this is an active cell
               if (ticc(ni1-1+face,j,k)/=patch) cycle
               if (iicc(ci1,j,k)/=0) cycle
               if (is_level_set) then
                  if (f0(ni1-1+face,j,k)>0._R8P) cycle
               endif

               ! distance from the patch
               NdS = 0.25_R8P*(nodes(i1,j,k)+nodes(i1,j-1,k)+nodes(i1,j,k-1)+nodes(i1,j-1,k-1)) - &
                     0.25_R8P*(nodes(i2,j,k)+nodes(i2,j-1,k)+nodes(i2,j,k-1)+nodes(i2,j-1,k-1))
               NdS = 0.5_R8P*NdS
               ds = sqrt( normL2(NFiS(ni1,j,k)) )
               twall = normL2((momentum(ci1-1+face,j,k)-momentum(ci1-2+face,j,k)).ortho.NFiS(ni1,j,k))/abs(Nds.dot.NFiS(ni1,j,k))/RE
               utau  = sqrt(twall)
               yplus(ci1,j,k) = RE*utau*( abs(Nds.dot.NFiS(ni1,j,k)/ds) )
            enddo
         enddo
      case(3,4)
         !$omp parallel default(none) &
         !$omp private(i,j,k,NdS)     &
         !$omp shared(ni1,ni2,nj1,nj2,nk1,nk2,ci1,ci2,cj1,cj2,ck1,ck2,face,patch,nodes,NFjS,iicc,ticc,momentum,f0,RE,yplus)
         !$omp do
         j1=nj1; j2=j1+1
         if (face==4) j2=j1-1
         do k=ck1,ck2
            do i=ci1,ci2
               ! checking if this is an active cell
               if (ticc(i,nj1-1+(face-2),k)/=patch) cycle
               if (iicc(i,cj1,k)/=0) cycle
               if (is_level_set) then
                  if (f0(i,nj1-1+(face-2),k)>0._R8P) cycle
               endif

               ! distance from the patch
               NdS = 0.25_R8P*(nodes(i,j1,k)+nodes(i-1,j1,k)+nodes(i,j1,k-1)+nodes(i-1,j1,k-1)) - &
                     0.25_R8P*(nodes(i,j2,k)+nodes(i-1,j2,k)+nodes(i,j2,k-1)+nodes(i-1,j2,k-1))
               NdS = 0.5_R8P*NdS
               ds = sqrt( normL2(NFjS(i,nj1,k)) )
               twall = normL2((momentum(i,cj1-1+(face-2),k) - momentum(i,cj1-2+(face-2),k)).ortho.NFjS(i,nj1,k))/ &
                                   abs(Nds.dot.NFjS(i,nj1,k))/RE
               utau  = sqrt(twall)
               yplus(i,cj1,k) = RE*utau*( abs(Nds.dot.NFjS(i,nj1,k)/ds) )
            enddo
         enddo
      case(5,6)
         !$omp parallel default(none) &
         !$omp private(i,j,k,NdS)     &
         !$omp shared(ni1,ni2,nj1,nj2,nk1,nk2,ci1,ci2,cj1,cj2,ck1,ck2,face,patch,nodes,NFkS,iicc,ticc,momentum,f0,RE,yplus)
         !$omp do
         k1=nk1; k2=k1+1
         if (face==6) k2=k1-1
         do j=cj1,cj2
            do i=ci1,ci2
               ! checking if this is an active cell
               if (ticc(i,j,nk1-1+(face-4))/=patch) cycle
               if (iicc(i,j,ck1)/=0) cycle
               if (is_level_set) then
                   if (f0(i,j,nk1-1+(face-4))>0._R8P) cycle
               endif

               ! distance from the patch
               NdS = 0.25_R8P*(nodes(i,j,k1)+nodes(i-1,j,k1)+nodes(i,j-1,k1)+nodes(i-1,j-1,k1)) - &
                     0.25_R8P*(nodes(i,j,k2)+nodes(i-1,j,k2)+nodes(i,j-1,k2)+nodes(i-1,j-1,k2))
               NdS = 0.5_R8P*NdS
               ds = sqrt( normL2(NFkS(i,j,nk1)) )
               twall = normL2((momentum(i,j,ck1-1+(face-4)) - momentum(i,j,ck1-2+(face-4))).ortho.NFkS(i,j,nk1))/ &
                               abs(Nds.dot.NFkS(i,j,nk1))/RE
               utau  = sqrt(twall)
               yplus(i,j,ck1) = RE*utau*( abs(Nds.dot.NFkS(i,j,nk1)/ds) )
            enddo
         enddo
      endselect
      endassociate
   enddo
   endsubroutine compute_yplus

   elemental function has_level_set(self)
   !< Return true if solution has level set functions.
   class(block_rst_object), intent(in) :: self          !< Block data.
   logical                             :: has_level_set !< Result of check.

   has_level_set = allocated(self%level_set).and.allocated(self%level_set_zero)
   endfunction has_level_set

   elemental function has_turb_kinetic_energy(self)
   !< Return true if solution has turbulent kinetic energy field.
   class(block_rst_object), intent(in) :: self                    !< Block data.
   logical                             :: has_turb_kinetic_energy !< Result of check.

   has_turb_kinetic_energy = allocated(self%turbulent_kinetic_energy)
   endfunction has_turb_kinetic_energy

   elemental function has_turb_kinetic_energy_diss(self)
   !< Return true if solution has turbulent kinetic energy dissipation field.
   class(block_rst_object), intent(in) :: self                         !< Block data.
   logical                             :: has_turb_kinetic_energy_diss !< Result of check.

   has_turb_kinetic_energy_diss = allocated(self%turbulent_kinetic_energy_dissipation)
   endfunction has_turb_kinetic_energy_diss

   elemental function has_turb_viscosity(self)
   !< Return true if solution has turbulent viscosity field.
   class(block_rst_object), intent(in) :: self               !< Block data.
   logical                             :: has_turb_viscosity !< Result of check.

   has_turb_viscosity = allocated(self%turbulent_viscosity)
   endfunction has_turb_viscosity

   elemental function has_viscosity(self)
   !< Return true if solution has viscosity field.
   class(block_rst_object), intent(in) :: self          !< Block data.
   logical                             :: has_viscosity !< Result of check.

   has_viscosity = allocated(self%viscosity)
   endfunction has_viscosity

   subroutine init(self, Ni, Nj, Nk, gc, is_level_set, is_zeroeq, is_oneeq, is_twoeq,             &
                   has_lambda2, has_qfactor, has_helicity, has_vorticity, has_div2LT, has_grad_p, &
                   has_k_ratio,                                                                   &
                   has_yplus, has_tau, has_div_tau, has_loads)
   !< Initialize block.
   class(block_rst_object), intent(inout)        :: self          !< Block data.
   integer(I4P),            intent(in)           :: Ni            !< Number of cells in i direction.
   integer(I4P),            intent(in)           :: Nj            !< Number of cells in j direction.
   integer(I4P),            intent(in)           :: Nk            !< Number of cells in k direction.
   integer(I4P),            intent(in), optional :: gc(6)         !< Number of ghost cells.
   logical,                 intent(in), optional :: is_level_set  !< Flag for level set function presence.
   logical,                 intent(in), optional :: is_zeroeq     !< Use *zero* equations turbulence model.
   logical,                 intent(in), optional :: is_oneeq      !< Use *one* equations turbulence model.
   logical,                 intent(in), optional :: is_twoeq      !< Use *two* equations turbulence model.
   logical,                 intent(in), optional :: has_lambda2   !< Solution has lamda2 field.
   logical,                 intent(in), optional :: has_qfactor   !< Solution has qfactor field.
   logical,                 intent(in), optional :: has_helicity  !< Solution has helicity field.
   logical,                 intent(in), optional :: has_vorticity !< Solution has vorticity field.
   logical,                 intent(in), optional :: has_div2LT    !< Solution has double divergence of Lighthill tensor.
   logical,                 intent(in), optional :: has_grad_p    !< Solution has pressure gradient.
   logical,                 intent(in), optional :: has_liutex    !< Solution has liutex field.
   logical,                 intent(in), optional :: has_k_ratio   !< Solution has kinetic energy ratio.
   logical,                 intent(in), optional :: has_yplus     !< Solution has y+ field.
   logical,                 intent(in), optional :: has_tau       !< Solution has tau field.
   logical,                 intent(in), optional :: has_div_tau   !< Solution has divergence of tau field.
   logical,                 intent(in), optional :: has_loads     !< Solution has loads (forces and torques).

   call self%destroy
   self%Ni = Ni
   self%Nj = Nj
   self%Nk = Nk
   if (present(gc            )) self%gc            = gc
   if (present(is_level_set  )) self%is_level_set  = is_level_set
   if (present(is_zeroeq     )) self%is_zeroeq     = is_zeroeq
   if (present(is_oneeq      )) self%is_oneeq      = is_oneeq
   if (present(is_twoeq      )) self%is_twoeq      = is_twoeq
   if (present(has_lambda2   )) self%has_lambda2   = has_lambda2
   if (present(has_qfactor   )) self%has_qfactor   = has_qfactor
   if (present(has_liutex    )) self%has_liutex    = has_liutex
   if (present(has_helicity  )) self%has_helicity  = has_helicity
   if (present(has_vorticity )) self%has_vorticity = has_vorticity
   if (present(has_div2LT    )) self%has_div2LT    = has_div2LT
   if (present(has_grad_p    )) self%has_grad_p    = has_grad_p
   if (present(has_k_ratio   )) self%has_k_ratio   = has_k_ratio
   if (present(has_yplus     )) self%has_yplus     = has_yplus
   if (present(has_tau       )) self%has_tau       = has_tau
   if (present(has_div_tau   )) self%has_div_tau   = has_div_tau
   if (present(has_loads     )) self%has_loads     = has_loads
   call self%alloc
   endsubroutine init

   subroutine interpolate_at_nodes(self)
   !< Interpolate solution at nodes.
   class(block_rst_object), intent(inout) :: self       !< Block data.
   real(R8P), allocatable                 :: tmp(:,:,:) !< Temporary buffer.

   tmp = self%momentum%x ; call interpolate_R8P(var=tmp, vari=self%momentum%x)
   tmp = self%momentum%y ; call interpolate_R8P(var=tmp, vari=self%momentum%y)
   tmp = self%momentum%z ; call interpolate_R8P(var=tmp, vari=self%momentum%z)
   tmp = self%pressure   ; call interpolate_R8P(var=tmp, vari=self%pressure  )
   if (allocated(self%level_set)) then
      tmp = self%level_set      ; call interpolate_R8P(var=tmp, vari=self%level_set     )
      tmp = self%level_set_zero ; call interpolate_R8P(var=tmp, vari=self%level_set_zero)
   endif
   if (allocated(self%viscosity)) then
      tmp = self%viscosity ; call interpolate_R8P(var=tmp, vari=self%viscosity)
   endif
   if (allocated(self%turbulent_viscosity)) then
      tmp = self%turbulent_viscosity ; call interpolate_R8P(var=tmp, vari=self%turbulent_viscosity)
   endif
   if (allocated(self%turbulent_kinetic_energy)) then
      tmp = self%turbulent_kinetic_energy ; call interpolate_R8P(var=tmp, vari=self%turbulent_kinetic_energy)
   endif
   if (allocated(self%turbulent_kinetic_energy_dissipation)) then
      tmp = self%turbulent_kinetic_energy_dissipation ; call interpolate_R8P(var=tmp,vari=self%turbulent_kinetic_energy_dissipation)
   endif
   if (allocated(self%lambda2)) then
      tmp = self%lambda2 ; call interpolate_R8P(var=tmp, vari=self%lambda2)
   endif
   if (allocated(self%qfactor)) then
      tmp = self%qfactor ; call interpolate_R8P(var=tmp, vari=self%qfactor)
   endif
   if (allocated(self%liutex)) then
      tmp = self%liutex ; call interpolate_R8P(var=tmp, vari=self%liutex)
   endif
   if (allocated(self%helicity)) then
      tmp = self%helicity ; call interpolate_R8P(var=tmp, vari=self%helicity)
   endif
   if (allocated(self%vorticity)) then
      tmp = self%vorticity%x ; call interpolate_R8P(var=tmp, vari=self%vorticity%x)
      tmp = self%vorticity%y ; call interpolate_R8P(var=tmp, vari=self%vorticity%y)
      tmp = self%vorticity%z ; call interpolate_R8P(var=tmp, vari=self%vorticity%z)
   endif
   if (allocated(self%div2LT)) then
      tmp = self%div2LT ; call interpolate_R8P(var=tmp, vari=self%div2LT)
   endif
   if (allocated(self%k_ratio)) then
      tmp = self%k_ratio ; call interpolate_R8P(var=tmp, vari=self%k_ratio)
   endif
   contains
      subroutine interpolate_R8P(var, vari)
      !< Interpolate a real(R8P) scalar.
      real(R8P), intent(in)  :: var (0:,0:,0:) !< Cell centered variable.
      real(R8P), intent(out) :: vari(0:,0:,0:) !< Node centered interpolated variable.
      integer(I4P)           :: i, j, k        !< Counters.

      do k=0, self%Nk
         do j=0, self%Nj
            do i=0, self%Ni
               vari(i,j,k) = var(i+1,j+1,k+1)  + var(i,j+1,k+1) &
                           + var(i+1,j  ,k+1)  + var(i,j,  k+1) &
                           + var(i+1,j+1,k  )  + var(i,j+1,k  ) &
                           + var(i+1,j  ,k  )  + var(i,j  ,k  )
               vari(i,j,k) = 0.125_R8P * vari(i,j,k)
            enddo
         enddo
      enddo
      !$OMP END PARALLEL
      endsubroutine interpolate_R8P
   endsubroutine interpolate_at_nodes

   subroutine load_dimensions(self, file_unit, is_level_set, is_zeroeq, is_oneeq, is_twoeq,                     &
                              compute_lambda2,compute_qfactor,compute_helicity,compute_vorticity,compute_div2LT,&
                              compute_grad_p,computex_liutex, compute_k_ratio,                                                   &
                              compute_yplus,compute_tau,compute_div_tau,compute_loads)
   !< Load block dimensions from file.
   !<
   !< @note The solution file must be already open and the current record-index must be at the proper block dimensions record.
   class(block_rst_object), intent(inout)        :: self              !< Block data.
   integer(I4P),            intent(in)           :: file_unit         !< Logical file unit.
   logical,                 intent(in), optional :: is_level_set      !< Flag for level set function presence.
   logical,                 intent(in), optional :: is_zeroeq         !< Use *zero* equations turbulence model.
   logical,                 intent(in), optional :: is_oneeq          !< Use *one* equations turbulence model.
   logical,                 intent(in), optional :: is_twoeq          !< Use *two* equations turbulence model.
   logical,                 intent(in), optional :: compute_lambda2   !< Compute lamda2 field.
   logical,                 intent(in), optional :: compute_qfactor   !< Compute qfactor field.
   logical,                 intent(in), optional :: compute_helicity  !< Compute helicity field.
   logical,                 intent(in), optional :: compute_vorticity !< Compute vorticity field.
   logical,                 intent(in), optional :: compute_div2LT    !< Compute double divergence of Lighthill tensor.
   logical,                 intent(in), optional :: compute_grad_p    !< Compute gradient pressure.
   logical,                 intent(in), optional :: compute_liutex    !< Compute liutex.
   logical,                 intent(in), optional :: compute_k_ratio   !< Compute kinetic energy ratio.
   logical,                 intent(in), optional :: compute_yplus     !< Compute y+ field.
   logical,                 intent(in), optional :: compute_tau       !< Compute tau field.
   logical,                 intent(in), optional :: compute_div_tau   !< Compute divergence of tau field.
   logical,                 intent(in), optional :: compute_loads     !< Compute loads (forces and torques).
   integer(I4P)                                  :: Ni                !< Number of cells in i direction.
   integer(I4P)                                  :: Nj                !< Number of cells in j direction.
   integer(I4P)                                  :: Nk                !< Number of cells in k direction.
   integer(I4P)                                  :: gc(6)             !< Number of ghost cells.

   gc = -1
   read(file_unit, end=10, err=10) Ni, Nj, Nk, gc
   10 continue
   if (all(gc >= 0)) then
      call self%init(Ni=Ni,Nj=Nj,Nk=Nk,gc=gc,                                                               &
                     is_level_set=is_level_set,is_zeroeq=is_zeroeq,is_oneeq=is_oneeq,is_twoeq=is_twoeq,     &
                     has_lambda2=compute_lambda2,has_qfactor=compute_qfactor,has_helicity=compute_helicity, &
                     has_vorticity=compute_vorticity,has_div2LT=compute_div2LT,has_grad_p=compute_grad_p,   &
                     has_liutex=compute_liutex, has_k_ratio=compute_k_ratio,                                &
                     has_yplus=compute_yplus,has_tau=compute_tau,has_div_tau=compute_div_tau,has_loads=compute_loads)
   else
      call self%init(Ni=Ni,Nj=Nj,Nk=Nk,                                                                     &
                     is_level_set=is_level_set,is_zeroeq=is_zeroeq,is_oneeq=is_oneeq,is_twoeq=is_twoeq,     &
                     has_lambda2=compute_lambda2,has_qfactor=compute_qfactor,has_helicity=compute_helicity, &
                     has_vorticity=compute_vorticity,has_div2LT=compute_div2LT,has_grad_p=compute_grad_p,   &
                     has_liutex=compute_liutex, has_k_ratio=compute_k_ratio,                                &
                     has_yplus=compute_yplus,has_tau=compute_tau,has_div_tau=compute_div_tau,has_loads=compute_loads)
   endif
   endsubroutine load_dimensions

   subroutine load_solution(self, file_unit, grd, icc, is_cell_centered, patch, RE, rFR2, zfs)
   !< Load block solution from file.
   !<
   !< @note The solution file must be already open and the current record-index must be at the proper block nodes record.
   class(block_rst_object), intent(inout)        :: self              !< Block data.
   integer(I4P),            intent(in)           :: file_unit         !< Logical file unit.
   type(block_grd_object),  intent(in)           :: grd               !< Block grd data.
   type(block_icc_object),  intent(in)           :: icc               !< Block icc data.
   logical,                 intent(in), optional :: is_cell_centered  !< Define variables at cell centers or nodes.
   integer(I4P),            intent(in), optional :: patch             !< Patch boundary conditions.
   real(R8P),               intent(in), optional :: RE                !< Reynolds number.
   real(R8P),               intent(in), optional :: rFR2              !< 1/(Froude number)^2.
   real(R8P),               intent(in), optional :: zfs               !< Z quote of free surface.
   integer(I4P)                                  :: i                 !< Counter.
   integer(I4P)                                  :: j                 !< Counter.
   integer(I4P)                                  :: k                 !< Counter.

   associate(Ni => self%Ni, Nj => self%Nj, Nk => self%Nk, gc => self%gc)
      read(file_unit)(((self%momentum(i, j, k)%x, i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      read(file_unit)(((self%momentum(i, j, k)%y, i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      read(file_unit)(((self%momentum(i, j, k)%z, i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      read(file_unit)(((self%pressure(i, j, k)  , i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      if (self%is_level_set) then
         read(file_unit)(((self%level_set(i, j, k), i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
         read(file_unit)(((self%level_set_zero(i, j, k), i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      endif
      if (self%is_zeroeq.or.self%is_oneeq.or.self%is_twoeq) then
         read(file_unit)(((self%viscosity(i, j, k), i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      endif
      if (self%is_oneeq) read(file_unit)(((self%turbulent_viscosity(i, j, k), i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      if (self%is_twoeq) then
         read(file_unit)(((self%turbulent_kinetic_energy(i, j, k), i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
         read(file_unit)(((self%turbulent_kinetic_energy_dissipation(i, j, k), i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      endif
   endassociate
   self%is_loaded = .true.
   call self%compute_aux(grd=grd, icc=icc, patch=patch, RE=RE, rFR2=rFR2, zfs=zfs)
   if (present(is_cell_centered)) then
      if (.not.is_cell_centered) call self%interpolate_at_nodes
   endif
   endsubroutine load_solution

   subroutine save_dimensions(self, file_unit, save_gc)
   !< Save block dimensions into file.
   !<
   !< @note The solution file must be already open and the current record-index must be at the proper block dimensions record.
   class(block_rst_object), intent(in)           :: self      !< Block data.
   integer(I4P),            intent(in)           :: file_unit !< Logical file unit.
   logical,                 intent(in), optional :: save_gc   !< Activate ghost cells saving.

   if (present(save_gc)) then
      if (save_gc) then
         write(file_unit) self%Ni, self%Nj, self%Nk, self%gc
      else
         write(file_unit) self%Ni, self%Nj, self%Nk
      endif
   else
      write(file_unit) self%Ni, self%Nj, self%Nk
   endif
   endsubroutine save_dimensions

   subroutine save_solution(self, file_unit)
   !< Save block solution into file.
   !<
   !< @note The solution file must be already open and the current record-index must be at the proper block nodes record.
   class(block_rst_object), intent(in) :: self      !< Block data.
   integer(I4P),            intent(in) :: file_unit !< Logical file unit.
   integer(I4P)                        :: i         !< Counter.
   integer(I4P)                        :: j         !< Counter.
   integer(I4P)                        :: k         !< Counter.

   associate(Ni => self%Ni, Nj => self%Nj, Nk => self%Nk, gc => self%gc)
      write(file_unit)(((self%momentum(i, j, k)%x, i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      write(file_unit)(((self%momentum(i, j, k)%y, i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      write(file_unit)(((self%momentum(i, j, k)%z, i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      write(file_unit)(((self%pressure(i, j, k)  , i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      if (self%is_level_set) then
         write(file_unit)(((self%level_set(i, j, k), i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
         write(file_unit)(((self%level_set_zero(i, j, k), i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      endif
      if (self%is_zeroeq.or.self%is_oneeq.or.self%is_twoeq) then
         write(file_unit)(((self%viscosity(i, j, k), i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      endif
      if (self%is_oneeq) write(file_unit)(((self%turbulent_viscosity(i, j, k), i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      if (self%is_twoeq) then
         write(file_unit)(((self%turbulent_kinetic_energy(i, j, k), i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
         write(file_unit)(((self%turbulent_kinetic_energy_dissipation(i, j, k), i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      endif
   endassociate
   endsubroutine save_solution
endmodule xview_block_rst_object
