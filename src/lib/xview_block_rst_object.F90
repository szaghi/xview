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
   logical                   :: is_level_set=.false.                        !< Use Level Set model.
   logical                   :: is_zeroeq=.false.                           !< Use *zero* equations turbulence model.
   logical                   :: is_oneeq=.false.                            !< Use *one* equations turbulence model.
   logical                   :: is_twoeq=.false.                            !< Use *two* equations turbulence model.
   type(vector), allocatable :: momentum(:,:,:)                             !< Momentum field.
   real(R8P),    allocatable :: pressure(:,:,:)                             !< Pressure field.
   real(R8P),    allocatable :: level_set(:,:,:)                            !< Level set field.
   real(R8P),    allocatable :: level_set_zero(:,:,:)                       !< Zero value of level set field.
   real(R8P),    allocatable :: viscosity(:,:,:)                            !< Viscosity field.
   real(R8P),    allocatable :: turbulent_viscosity(:,:,:)                  !< Turbulent viscosity field.
   real(R8P),    allocatable :: turbulent_kinetic_energy(:,:,:)             !< Turbulent kinetic energy field.
   real(R8P),    allocatable :: turbulent_kinetic_energy_dissipation(:,:,:) !< Turbulent kinetic energy dissipation field.
   logical                   :: has_aux=.false.                             !< Compute auxiliary variables.
   real(R8P),    allocatable :: lambda2(:,:,:)                              !< Variable to identify vortices (lambda 2).
   real(R8P),    allocatable :: qfactor(:,:,:)                              !< Variable to identify vortices (q factor).
   real(R8P),    allocatable :: helicity(:,:,:)                             !< Helicity.
   type(vector), allocatable :: vorticity(:,:,:)                            !< Vorticity.
   real(R8P),    allocatable :: yplus(:,:,:)                                !< Estimation of y+.
   type(vector), allocatable :: tau(:,:,:)                                  !< Tau wall.
   real(R8P),    allocatable :: div_tau(:,:,:)                              !< Divergence of Tau wall.
   logical                   :: is_loaded=.false.                           !< Flag for checking if the block is loaded.
   contains
      ! public methods
      procedure, pass(self) :: destroy                      !< Destroy dynamic memory.
      procedure, pass(self) :: alloc                        !< Allocate dynamic memory.
      procedure, pass(self) :: compute_aux                  !< Compute auxiliary varibles.
      procedure, pass(self) :: compute_yplus                !< Compute yplus on patches.
      procedure, pass(self) :: compute_tau                  !< Compute Tau wall and its divergence on patches.
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
   elemental subroutine destroy(self)
   !< Destroy dynamic memory.
   class(block_rst_object), intent(inout) :: self !< Block data.

   self%Ni = 0
   self%Nj = 0
   self%Nk = 0
   self%gc = 2
   self%is_level_set = .false.
   self%is_zeroeq = .false.
   self%is_oneeq = .false.
   self%is_twoeq = .false.
   if (allocated(self%momentum))                             deallocate(self%momentum)
   if (allocated(self%pressure))                             deallocate(self%pressure)
   if (allocated(self%level_set))                            deallocate(self%level_set)
   if (allocated(self%level_set_zero))                       deallocate(self%level_set_zero)
   if (allocated(self%viscosity))                            deallocate(self%viscosity)
   if (allocated(self%turbulent_viscosity))                  deallocate(self%turbulent_viscosity)
   if (allocated(self%turbulent_kinetic_energy))             deallocate(self%turbulent_kinetic_energy)
   if (allocated(self%turbulent_kinetic_energy_dissipation)) deallocate(self%turbulent_kinetic_energy_dissipation)
   self%has_aux = .false.
   if (allocated(self%lambda2))   deallocate(self%lambda2)
   if (allocated(self%qfactor))   deallocate(self%qfactor)
   if (allocated(self%helicity))  deallocate(self%helicity)
   if (allocated(self%vorticity)) deallocate(self%vorticity)
   if (allocated(self%yplus))     deallocate(self%yplus)
   if (allocated(self%tau))       deallocate(self%tau)
   if (allocated(self%div_tau))   deallocate(self%div_tau)
   self%is_loaded = .false.
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
   if (self%has_aux) then
      allocate(self%lambda2(  1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6))) ; self%lambda2   = 0._R8P
      allocate(self%qfactor(  1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6))) ; self%qfactor   = 0._R8P
      allocate(self%helicity( 1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6))) ; self%helicity  = 0._R8P
      allocate(self%vorticity(1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6))) ; self%vorticity = 0._R8P
      allocate(self%yplus(    0-gc(1):Ni+gc(2), 0-gc(3):Nj+gc(4), 0-gc(5):Nk+gc(6))) ; self%yplus     = 0._R8P
      allocate(self%tau(      0-gc(1):Ni+gc(2), 0-gc(3):Nj+gc(4), 0-gc(5):Nk+gc(6))) ; self%tau       = 0._R8P
      allocate(self%div_tau(  0-gc(1):Ni+gc(2), 0-gc(3):Nj+gc(4), 0-gc(5):Nk+gc(6))) ; self%div_tau   = 0._R8P
   endif
   endassociate
   endsubroutine alloc

   subroutine compute_aux(self, grd)
   !< Compute auxiliary varibles.
   class(block_rst_object), intent(inout) :: self                          !< Block data.
   type(block_grd_object),  intent(in)    :: grd                           !< Block grd data.
   real(R8P), allocatable                 :: Fi(:,:,:,:,:)                 !< Fluxes i direction.
   real(R8P), allocatable                 :: Fj(:,:,:,:,:)                 !< Fluxes j direction.
   real(R8P), allocatable                 :: Fk(:,:,:,:,:)                 !< Fluxes k direction.
   type(vector)                           :: um                            !< Dummy vector variables.
   real(R8P)                              :: c(0:2),emin,emax,eval,fval,mu !< Dummy reals.
   real(R8P), dimension(3,3)              :: IDEN,G,S,O                    !< Matrices.
   integer(I4P)                           :: i,j,k,ii,jj,kk,iter           !< Counter.
   real(R8P), parameter                   :: EPS6=1d-6, EPS9=1d-9          !< Tolerances.

   associate(Ni=>self%Ni, Nj=>self%Nj, Nk=>self%Nk, gc=>self%gc,                 &
             NFiS=>grd%NFiS, NFjS=>grd%NFjS, NFkS=>grd%NFkS, volume=>grd%volume, &
             momentum=>self%momentum,lambda2=>self%lambda2,qfactor=>self%qfactor,helicity=>self%helicity,vorticity=>self%vorticity)
   allocate(Fi(1:3,1:3,0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
   allocate(Fj(1:3,1:3,0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
   allocate(Fk(1:3,1:3,0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))

   Fi   = 0._R8P
   Fj   = 0._R8P
   Fk   = 0._R8P
   IDEN = 0._R8P
   do i=1, 3
      IDEN(i,i) = 1._R8P
   enddo
   ! extrapolate momentum from inner cells to ghost cells
   do k=1,Nk
      do j=1,Nj
         momentum(1- gc(1),j,k) = 2._R8P*momentum(0,   j,k)-momentum(1, j,k)
         momentum(Ni+gc(2),j,k) = 2._R8P*momentum(Ni+1,j,k)-momentum(Ni,j,k)
      enddo
   enddo
   do k=1,Nk
      do i=0,Ni+1
         momentum(i,1- gc(3),k) = 2._R8P*momentum(i,0,   k)-momentum(i,1, k)
         momentum(i,Nj+gc(4),k) = 2._R8P*momentum(i,Nj+1,k)-momentum(i,Nj,k)
      enddo
   enddo
   do j=0,Nj+1
      do i=0,Ni+1
         momentum(i,j,1- gc(5)) = 2._R8P*momentum(i,j,0   )-momentum(i,j,1 )
         momentum(i,j,Nk+gc(6)) = 2._R8P*momentum(i,j,Nk+1)-momentum(i,j,Nk)
      enddo
   enddo
   ! compute fluxes
   !$omp parallel default(none) &
   !$omp private(i,j,k,um)      &
   !$omp shared(Ni,Nj,Nk,Fi,Fj,Fk,NFiS,NFjS,NFkS,momentum)
   !$omp do
   do k=0,Nk+1
      do j=0,Nj+1
         do i=-1,Ni+1
            um = 0.5_R8P*(momentum(i,j,k)+momentum(i+1,j,k))
            Fi(1,1,i,j,k) = um%x*NFiS(i,j,k)%x
            Fi(1,2,i,j,k) = um%x*NFiS(i,j,k)%y
            Fi(1,3,i,j,k) = um%x*NFiS(i,j,k)%z
            Fi(2,1,i,j,k) = um%y*NFiS(i,j,k)%x
            Fi(2,2,i,j,k) = um%y*NFiS(i,j,k)%y
            Fi(2,3,i,j,k) = um%y*NFiS(i,j,k)%z
            Fi(3,1,i,j,k) = um%z*NFiS(i,j,k)%x
            Fi(3,2,i,j,k) = um%z*NFiS(i,j,k)%y
            Fi(3,3,i,j,k) = um%z*NFiS(i,j,k)%z
         enddo
      enddo
   enddo
   !$omp do
   do k=0,Nk+1
      do j=-1,Nj+1
         do i=0,Ni+1
            um = 0.5_R8P*(momentum(i,j,k)+momentum(i,j+1,k))
            Fj(1,1,i,j,k) = um%x*NFjS(i,j,k)%x
            Fj(1,2,i,j,k) = um%x*NFjS(i,j,k)%y
            Fj(1,3,i,j,k) = um%x*NFjS(i,j,k)%z
            Fj(2,1,i,j,k) = um%y*NFjS(i,j,k)%x
            Fj(2,2,i,j,k) = um%y*NFjS(i,j,k)%y
            Fj(2,3,i,j,k) = um%y*NFjS(i,j,k)%z
            Fj(3,1,i,j,k) = um%z*NFjS(i,j,k)%x
            Fj(3,2,i,j,k) = um%z*NFjS(i,j,k)%y
            Fj(3,3,i,j,k) = um%z*NFjS(i,j,k)%z
         enddo
      enddo
   enddo
   !$omp do
   do k=-1,Nk+1
      do j=0,Nj+1
         do i=0,Ni+1
            um = 0.5_R8P*(momentum(i,j,k)+momentum(i,j,k+1))
            Fk(1,1,i,j,k) = um%x*NFkS(i,j,k)%x
            Fk(1,2,i,j,k) = um%x*NFkS(i,j,k)%y
            Fk(1,3,i,j,k) = um%x*NFkS(i,j,k)%z
            Fk(2,1,i,j,k) = um%y*NFkS(i,j,k)%x
            Fk(2,2,i,j,k) = um%y*NFkS(i,j,k)%y
            Fk(2,3,i,j,k) = um%y*NFkS(i,j,k)%z
            Fk(3,1,i,j,k) = um%z*NFkS(i,j,k)%x
            Fk(3,2,i,j,k) = um%z*NFkS(i,j,k)%y
            Fk(3,3,i,j,k) = um%z*NFkS(i,j,k)%z
         enddo
      enddo
   enddo
   !$omp end parallel
   do k=0,Nk+1
      do j=0,Nj+1
         do i=0,Ni+1
            ! compute velocity gradient
            do jj=1,3
               do ii=1,3
                  G(ii,jj) = FI(ii,jj,i,j,k)-FI(ii,jj,i-1,j,k)&
                           + FJ(ii,jj,i,j,k)-FJ(ii,jj,i,j-1,k)&
                           + FK(ii,jj,i,j,k)-FK(ii,jj,i,j,k-1)
                  G(ii,jj) = G(ii,jj)/max(eps6*eps6,volume(i,j,k))
               enddo
            enddo

            ! compute vorticity vector
            um%x =   (G(3,2) - G(2,3))
            um%y = - (G(1,3) - G(3,1)) ! control signs!!!!
            um%z =   (G(2,1) - G(1,2))

            ! tensor S^2 + O^2 (saved in G)
            S = 0._R8P
            O = 0._R8P
            do kk=1,3
               do jj=1,3
                  S(jj,kk) = 0.5_R8P*(G(jj,kk)+G(kk,jj))
                  O(jj,kk) = 0.5_R8P*(G(jj,kk)-G(kk,jj))
               enddo
            enddo
            G = matmul(S,S) + matmul(O,O)

            ! coefficients of characterist polynomial: lamda^3 + c(2)*lamda^2 + c(1)*lamda + c(0) = 0
            c(2) = -(G(1,1) + G(2,2) + G(3,3))
            c(1) = G(1,1)*G(2,2) + G(1,1)*G(3,3) + G(2,2)*G(3,3) - G(2,3)**2 - G(1,3)**2 - G(1,2)**2
            c(0) = G(1,1)*G(2,3)**2 + G(2,2)*G(1,3)**2 + G(3,3)*G(1,2)**2 - 2._R8P*G(2,3)*G(1,3)*G(1,2) - G(1,1)*G(2,2)*G(3,3)

            ! compute second eigenvalue of characteristic polynomial
            mu = sqrt(c(2)**2 - 3._R8P*c(1))
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

            lambda2(  i,j,k) = eval
            qfactor(  i,j,k) = 0.5_R8P*(dot_product(O(1,:),O(1,:)) + dot_product(O(2,:),O(2,:)) + dot_product(O(3,:),O(3,:)) - &
                                       (dot_product(S(1,:),S(1,:)) + dot_product(S(2,:),S(2,:)) + dot_product(S(3,:),S(3,:))))
            helicity( i,j,k) = momentum(i,j,k).dot.um / (max(eps6*eps6,normL2(momentum(i,j,k))*normL2(um)))
            vorticity(i,j,k) = um
         enddo
      enddo
   enddo
   endassociate
   endsubroutine compute_aux

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

   subroutine init(self, Ni, Nj, Nk, gc, is_level_set, is_zeroeq, is_oneeq, is_twoeq, is_aux_to_compute)
   !< Initialize block.
   class(block_rst_object), intent(inout)        :: self              !< Block data.
   integer(I4P),            intent(in)           :: Ni                !< Number of cells in i direction.
   integer(I4P),            intent(in)           :: Nj                !< Number of cells in j direction.
   integer(I4P),            intent(in)           :: Nk                !< Number of cells in k direction.
   integer(I4P),            intent(in), optional :: gc(6)             !< Number of ghost cells.
   logical,                 intent(in), optional :: is_level_set      !< Flag for level set function presence.
   logical,                 intent(in), optional :: is_zeroeq         !< Use *zero* equations turbulence model.
   logical,                 intent(in), optional :: is_oneeq          !< Use *one* equations turbulence model.
   logical,                 intent(in), optional :: is_twoeq          !< Use *two* equations turbulence model.
   logical,                 intent(in), optional :: is_aux_to_compute !< Flag to allocate also metrics arrays.

   call self%destroy
   self%Ni = Ni
   self%Nj = Nj
   self%Nk = Nk
   if (present(gc               )) self%gc           = gc
   if (present(is_level_set     )) self%is_level_set = is_level_set
   if (present(is_zeroeq        )) self%is_zeroeq    = is_zeroeq
   if (present(is_oneeq         )) self%is_oneeq     = is_oneeq
   if (present(is_twoeq         )) self%is_twoeq     = is_twoeq
   if (present(is_aux_to_compute)) self%has_aux      = is_aux_to_compute
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
   if (allocated(self%helicity)) then
      tmp = self%helicity ; call interpolate_R8P(var=tmp, vari=self%helicity)
   endif
   if (allocated(self%vorticity)) then
      tmp = self%vorticity%x ; call interpolate_R8P(var=tmp, vari=self%vorticity%x)
      tmp = self%vorticity%y ; call interpolate_R8P(var=tmp, vari=self%vorticity%y)
      tmp = self%vorticity%z ; call interpolate_R8P(var=tmp, vari=self%vorticity%z)
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

   subroutine load_dimensions(self, file_unit, is_level_set, is_zeroeq, is_oneeq, is_twoeq, is_aux_to_compute)
   !< Load block dimensions from file.
   !<
   !< @note The solution file must be already open and the current record-index must be at the proper block dimensions record.
   class(block_rst_object), intent(inout)        :: self              !< Block data.
   integer(I4P),            intent(in)           :: file_unit         !< Logical file unit.
   logical,                 intent(in), optional :: is_level_set      !< Flag for level set function presence.
   logical,                 intent(in), optional :: is_zeroeq         !< Use *zero* equations turbulence model.
   logical,                 intent(in), optional :: is_oneeq          !< Use *one* equations turbulence model.
   logical,                 intent(in), optional :: is_twoeq          !< Use *two* equations turbulence model.
   logical,                 intent(in), optional :: is_aux_to_compute !< Flag to allocate also metrics arrays.
   integer(I4P)                                  :: Ni                !< Number of cells in i direction.
   integer(I4P)                                  :: Nj                !< Number of cells in j direction.
   integer(I4P)                                  :: Nk                !< Number of cells in k direction.
   integer(I4P)                                  :: gc(6)             !< Number of ghost cells.

   gc = -1
   read(file_unit, end=10, err=10) Ni, Nj, Nk, gc
   10 continue
   if (all(gc >= 0)) then
      call self%init(Ni=Ni,Nj=Nj,Nk=Nk,gc=gc,                                                          &
                     is_level_set=is_level_set,is_zeroeq=is_zeroeq,is_oneeq=is_oneeq,is_twoeq=is_twoeq,&
                     is_aux_to_compute=is_aux_to_compute)
   else
      call self%init(Ni=Ni,Nj=Nj,Nk=Nk,                                                                &
                     is_level_set=is_level_set,is_zeroeq=is_zeroeq,is_oneeq=is_oneeq,is_twoeq=is_twoeq,&
                     is_aux_to_compute=is_aux_to_compute)
   endif
   endsubroutine load_dimensions

   subroutine load_solution(self, file_unit, is_cell_centered, is_aux_to_compute, grd)
   !< Load block solution from file.
   !<
   !< @note The solution file must be already open and the current record-index must be at the proper block nodes record.
   class(block_rst_object), intent(inout)        :: self              !< Block data.
   integer(I4P),            intent(in)           :: file_unit         !< Logical file unit.
   logical,                 intent(in), optional :: is_cell_centered  !< Define variables at cell centers or nodes.
   logical,                 intent(in), optional :: is_aux_to_compute !< Flag to allocate also metrics arrays.
   type(block_grd_object),  intent(in), optional :: grd               !< Block grd data.
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
   if (present(is_aux_to_compute)) then
      if (is_aux_to_compute.and.present(grd)) call self%compute_aux(grd)
   endif
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
