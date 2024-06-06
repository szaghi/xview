!< xview, block (grd) class definition.
module xview_block_grd_object
!< xview, block (grd) class definition.

use penf
use vecfor

implicit none
private
public :: block_grd_object

type :: block_grd_object
   !< Block (grd) class definition.
   integer(I4P)              :: Ni=0                     !< Number of cells in i direction.
   integer(I4P)              :: Nj=0                     !< Number of cells in j direction.
   integer(I4P)              :: Nk=0                     !< Number of cells in k direction.
   integer(I4P)              :: gc(6)=[2, 2, 2, 2, 2, 2] !< Number of ghost cells.
   type(vector), allocatable :: nodes(:,:,:)             !< Nodes coordinates.
   type(vector), allocatable :: centers(:,:,:)           !< Centers coordinates.
   type(vector), allocatable :: extents(:)               !< Box extents, [min, max].
   type(vector), allocatable :: sub_extents(:,:)         !< sub block extents, [min, max]: block eight.
   integer(I4P), allocatable :: sub_ijk_extents(:,:,:)   !< sub block ijk extents, [min, max]: block eight.
   type(vector), allocatable :: NFi(:,:,:)               !< Face i normal versor.
   type(vector), allocatable :: NFj(:,:,:)               !< Face j normal versor.
   type(vector), allocatable :: NFk(:,:,:)               !< Face k normal versor.
   type(vector), allocatable :: NFiS(:,:,:)              !< Face i normal versor with surface area module.
   type(vector), allocatable :: NFjS(:,:,:)              !< Face j normal versor with surface area module.
   type(vector), allocatable :: NFkS(:,:,:)              !< Face k normal versor with surface area module.
   real(R8P),    allocatable :: Si(:,:,:)                !< Face i area.
   real(R8P),    allocatable :: Sj(:,:,:)                !< Face j area.
   real(R8P),    allocatable :: Sk(:,:,:)                !< Face k area.
   real(R8P),    allocatable :: volume(:,:,:)            !< Volumes of cells.
   integer(I4P), allocatable :: patches_extents(:,:)     !< Patches extents, [np, 13]. The second index means
                                                         !<+ 0  => patch face (1,2,3,4,5,6);
                                                         !<+ 1  => cell i-min;
                                                         !<+ 2  => cell i-max;
                                                         !<+ 3  => cell j-min;
                                                         !<+ 4  => cell j-max;
                                                         !<+ 5  => cell k-min;
                                                         !<+ 6  => cell k-max;
                                                         !<+ 7  => node i-min;
                                                         !<+ 8  => node i-max;
                                                         !<+ 9  => node j-min;
                                                         !<+ 10 => node j-max;
                                                         !<+ 11 => node k-min;
                                                         !<+ 12 => node k-max;
   logical                   :: is_loaded=.false.        !< Flag for checking if the block is loaded.
   contains
      ! public methods
      procedure, pass(self) :: alloc                   !< Allocate dynamic memory.
      procedure, pass(self) :: compute_metrics         !< Compute metrics.
      procedure, pass(self) :: compute_patches_extents !< Compute patches extents for given bc or whole block extents.
      procedure, pass(self) :: correct_metrics_bc      !< Correct metrics of Boundary Conditions (ghost) cells.
      procedure, pass(self) :: destroy                 !< Destroy dynamic memory.
      procedure, pass(self) :: load_dimensions         !< Load block dimensions from file.
      procedure, pass(self) :: load_nodes              !< Load block nodes from file.
      procedure, pass(self) :: save_dimensions         !< Save block dimensions into file.
      procedure, pass(self) :: save_nodes              !< Save block nodes into file.
      procedure, pass(self) :: traslate                !< Traslate block nodes by a given traslation vector.
endtype block_grd_object

contains
   elemental subroutine alloc(self, is_centers_to_allocate, is_extents_to_allocate, is_metrics_to_allocate)
   !< Allocate dynamic memory.
   class(block_grd_object), intent(inout)        :: self                   !< Block data.
   logical,                 intent(in), optional :: is_centers_to_allocate !< Flag to allocate also centers array.
   logical,                 intent(in), optional :: is_extents_to_allocate !< Flag to allocate also extents arrays.
   logical,                 intent(in), optional :: is_metrics_to_allocate !< Flag to allocate also metrics arrays.

   associate(Ni=>self%Ni, Nj=>self%Nj, Nk=>self%Nk, gc=>self%gc)
   allocate(self%nodes(0-gc(1):Ni+gc(2), 0-gc(3):Nj+gc(4), 0-gc(5):Nk+gc(6)))
   if (present(is_centers_to_allocate)) then
      if (is_centers_to_allocate) allocate(self%centers(1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6)))
   endif
   if (present(is_extents_to_allocate)) then
      if (is_extents_to_allocate) then
         allocate(self%extents(1:2)) ! min-max
         allocate(self%sub_extents(1:8, 1:2)) ! eight parts, min-max
         allocate(self%sub_ijk_extents(1:8, 1:3, 1:2))  ! eight parts, ijk, min-max
      endif
   endif
   if (present(is_metrics_to_allocate)) then
      if (is_metrics_to_allocate) then
         allocate(self%NFi   (0-gc(1):Ni+gc(2), 0-gc(3):Nj+gc(4), 0-gc(5):Nk+gc(6))) ; self%NFi    = 0._R8P
         allocate(self%NFj   (0-gc(1):Ni+gc(2), 0-gc(3):Nj+gc(4), 0-gc(5):Nk+gc(6))) ; self%NFj    = 0._R8P
         allocate(self%NFk   (0-gc(1):Ni+gc(2), 0-gc(3):Nj+gc(4), 0-gc(5):Nk+gc(6))) ; self%NFk    = 0._R8P
         allocate(self%NFiS  (0-gc(1):Ni+gc(2), 0-gc(3):Nj+gc(4), 0-gc(5):Nk+gc(6))) ; self%NFiS   = 0._R8P
         allocate(self%NFjS  (0-gc(1):Ni+gc(2), 0-gc(3):Nj+gc(4), 0-gc(5):Nk+gc(6))) ; self%NFjS   = 0._R8P
         allocate(self%NFkS  (0-gc(1):Ni+gc(2), 0-gc(3):Nj+gc(4), 0-gc(5):Nk+gc(6))) ; self%NFkS   = 0._R8P
         allocate(self%Si    (0-gc(1):Ni+gc(2), 0-gc(3):Nj+gc(4), 0-gc(5):Nk+gc(6))) ; self%Si     = 0._R8P
         allocate(self%Sj    (0-gc(1):Ni+gc(2), 0-gc(3):Nj+gc(4), 0-gc(5):Nk+gc(6))) ; self%Sj     = 0._R8P
         allocate(self%Sk    (0-gc(1):Ni+gc(2), 0-gc(3):Nj+gc(4), 0-gc(5):Nk+gc(6))) ; self%Sk     = 0._R8P
         allocate(self%volume(1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6))) ; self%volume = 0._R8P
      endif
   endif
   endassociate
   endsubroutine alloc

   subroutine compute_metrics(self)
   !< Compute block metrics.
   class(block_grd_object), intent(inout) :: self                !< Block data.
   type(vector)                           :: NFS, s1, s2, nd, db !< Dummy vectors.
   real(R8P)                              :: signi, signj, signk !< Directions of face normals
   real(R8P)                              :: Vx, Vy, Vz          !< Volumes.
   real(R8P)                              :: xp, yp, zp          !< Face coordinates, plus.
   real(R8P)                              :: xm, ym, zm          !< Face coordinates, minus.
   integer(I4P)                           :: i,j,k               !< Counter.

   associate(Ni=>self%Ni, Nj=>self%Nj, Nk=>self%Nk, gc=>self%gc, nodes=>self%nodes, &
             NFi=>self%NFi, NFj=>self%NFj, NFk=>self%NFk,                           &
             NFiS=>self%NFiS, NFjS=>self%NFjS, NFkS=>self%NFkS,                     &
             Si=>self%Si, Sj=>self%Sj, Sk=>self%Sk, volume=>self%volume)
   ! computing faces normals
   ! positioning at the middle of the block
   i = max(1,Ni/2)
   j = max(1,Nj/2)
   k = max(1,Nk/2)
   ! compute direction of i normals
   s1 = nodes(i,j  ,k) - nodes(i,j-1,k-1)
   s2 = nodes(i,j-1,k) - nodes(i,j,  k-1)
   nd = s1.cross.s2
   s1 = 0.25_R8P*(nodes(i,  j,k)+nodes(i,  j-1,k)+nodes(i,  j,k-1)+nodes(i,  j-1,k-1))
   s2 = 0.25_R8P*(nodes(i-1,j,k)+nodes(i-1,j-1,k)+nodes(i-1,j,k-1)+nodes(i-1,j-1,k-1))
   db = s1 - s2
   signi = sign(1._R8P,(nd.dot.db))
   ! compute direction of j normals
   s1 = nodes(i,j,k  ) - nodes(i-1,j,k-1)
   s2 = nodes(i,j,k-1) - nodes(i-1,j,k  )
   nd = s1.cross.s2
   s1 = 0.25_R8P*(nodes(i,j,  k)+nodes(i-1,j,  k)+nodes(i,j,  k-1)+nodes(i-1,j,  k-1))
   s2 = 0.25_R8P*(nodes(i,j-1,k)+nodes(i-1,j-1,k)+nodes(i,j-1,k-1)+nodes(i-1,j-1,k-1))
   db = s1 - s2
   signj = sign(1._R8P,(nd.dot.db))
   ! compute direction of k normals
   s1 = nodes(i,  j,k) - nodes(i-1,j-1,k)
   s2 = nodes(i-1,j,k) - nodes(i,  j-1,k)
   nd = s1.cross.s2
   s1 = 0.25_R8P*(nodes(i,j,k  )+nodes(i-1,j,k  )+nodes(i,j-1,k  )+nodes(i-1,j-1,k  ))
   s2 = 0.25_R8P*(nodes(i,j,k-1)+nodes(i-1,j,k-1)+nodes(i,j-1,k-1)+nodes(i-1,j-1,k-1))
   db = s1 - s2
   signk = sign(1._R8P,(nd.dot.db))
   !$omp parallel default(none)                        &
   !$omp private(i,j,k,NFS,Vx,Vy,Vz,xp,yp,zp,xm,ym,zm) &
   !$omp shared(Ni,Nj,Nk,gc,signi,signj,signk,nodes,Si,Sj,Sk,NFi,NFj,NFk,NFiS,NFjS,NFkS,volume)
   ! compute faces metrics
   !$omp do
   do k=1-gc(5),Nk+gc(6)
      do j=1-gc(3),Nj+gc(4)
         do i=0-gc(1),Ni+gc(2)
            NFS = face_normal4(pt1=nodes(i,j-1,k-1), &
                               pt2=nodes(i,j  ,k-1), &
                               pt3=nodes(i,j  ,k  ), &
                               pt4=nodes(i,j-1,k  ))
            NFS = NFS*signi
            NFiS(i,j,k) = NFS
            NFi (i,j,k) = normalized(NFS)
            Si  (i,j,k) =     normL2(NFS)
         enddo
      enddo
   enddo
   !$omp do
   do k=1-gc(5),Nk+gc(6)
      do j=0-gc(3),Nj+gc(4)
         do i=1-gc(1),Ni+gc(2)
            NFS = face_normal4(pt1=nodes(i-1,j,k-1), &
                               pt2=nodes(i-1,j,k  ), &
                               pt3=nodes(i  ,j,k  ), &
                               pt4=nodes(i  ,j,k-1))
            NFS = NFS*signj
            NFjS(i,j,k) = NFS
            NFj (i,j,k) = normalized(NFS)
            Sj  (i,j,k) =     normL2(NFS)
         enddo
      enddo
   enddo
   !$omp do
   do k=0-gc(5),Nk+gc(6)
     do j=1-gc(3),Nj+gc(4)
       do i=1-gc(1),Ni+gc(2)
         NFS = face_normal4(pt1 = nodes(i-1,j-1,k), &
                            pt2 = nodes(i  ,j-1,k), &
                            pt3 = nodes(i  ,j  ,k), &
                            pt4 = nodes(i-1,j  ,k))
         NFS = NFS*signk
         NFkS(i,j,k) = NFS
         NFk (i,j,k) = normalized(NFS)
         Sk  (i,j,k) =     normL2(NFS)
       enddo
     enddo
   enddo
   ! compute cells volumes
   !$omp do
   do k=1-gc(5),Nk+gc(6)
      do j=1-gc(3),Nj+gc(4)
         do i=1-gc(1),Ni+gc(2)
            Vx = 0._R8P
            Vy = 0._R8P
            Vz = 0._R8P

            xp = 0.25_R8P*(nodes(i  ,j  ,k  )%x + nodes(i  ,j  ,k-1)%x + &
                           nodes(i  ,j-1,k  )%x + nodes(i  ,j-1,k-1)%x)
            yp = 0.25_R8P*(nodes(i  ,j  ,k  )%y + nodes(i  ,j  ,k-1)%y + &
                           nodes(i  ,j-1,k  )%y + nodes(i  ,j-1,k-1)%y)
            zp = 0.25_R8P*(nodes(i  ,j  ,k  )%z + nodes(i  ,j  ,k-1)%z + &
                           nodes(i  ,j-1,k  )%z + nodes(i  ,j-1,k-1)%z)
            xm = 0.25_R8P*(nodes(i-1,j  ,k  )%x + nodes(i-1,j  ,k-1)%x + &
                           nodes(i-1,j-1,k  )%x + nodes(i-1,j-1,k-1)%x)
            ym = 0.25_R8P*(nodes(i-1,j  ,k  )%y + nodes(i-1,j  ,k-1)%y + &
                           nodes(i-1,j-1,k  )%y + nodes(i-1,j-1,k-1)%y)
            zm = 0.25_R8P*(nodes(i-1,j  ,k  )%z + nodes(i-1,j  ,k-1)%z + &
                           nodes(i-1,j-1,k  )%z + nodes(i-1,j-1,k-1)%z)

            Vx = Vx + xp*NFiS(i,j,k)%x - xm*NFiS(i-1,j,k)%x
            Vy = Vy + yp*NFiS(i,j,k)%y - ym*NFiS(i-1,j,k)%y
            Vz = Vz + zp*NFiS(i,j,k)%z - zm*NFiS(i-1,j,k)%z

            xp = 0.25_R8P*(nodes(i  ,j  ,k  )%x + nodes(i  ,j  ,k-1)%x + &
                           nodes(i-1,j  ,k  )%x + nodes(i-1,j  ,k-1)%x)
            yp = 0.25_R8P*(nodes(i  ,j  ,k  )%y + nodes(i  ,j  ,k-1)%y + &
                           nodes(i-1,j  ,k  )%y + nodes(i-1,j  ,k-1)%y)
            zp = 0.25_R8P*(nodes(i  ,j  ,k  )%z + nodes(i  ,j  ,k-1)%z + &
                           nodes(i-1,j  ,k  )%z + nodes(i-1,j  ,k-1)%z)
            xm = 0.25_R8P*(nodes(i  ,j-1,k  )%x + nodes(i  ,j-1,k-1)%x + &
                           nodes(i-1,j-1,k  )%x + nodes(i-1,j-1,k-1)%x)
            ym = 0.25_R8P*(nodes(i  ,j-1,k  )%y + nodes(i  ,j-1,k-1)%y + &
                           nodes(i-1,j-1,k  )%y + nodes(i-1,j-1,k-1)%y)
            zm = 0.25_R8P*(nodes(i  ,j-1,k  )%z + nodes(i  ,j-1,k-1)%z + &
                           nodes(i-1,j-1,k  )%z + nodes(i-1,j-1,k-1)%z)

            Vx = Vx + xp*NFjS(i,j,k)%x - xm*NFjS(i,j-1,k)%x
            Vy = Vy + yp*NFjS(i,j,k)%y - ym*NFjS(i,j-1,k)%y
            Vz = Vz + zp*NFjS(i,j,k)%z - zm*NFjS(i,j-1,k)%z

            xp = 0.25_R8P*(nodes(i  ,j  ,k  )%x + nodes(i  ,j-1,k  )%x + &
                           nodes(i-1,j  ,k  )%x + nodes(i-1,j-1,k  )%x)
            yp = 0.25_R8P*(nodes(i  ,j  ,k  )%y + nodes(i  ,j-1,k  )%y + &
                           nodes(i-1,j  ,k  )%y + nodes(i-1,j-1,k  )%y)
            zp = 0.25_R8P*(nodes(i  ,j  ,k  )%z + nodes(i  ,j-1,k  )%z + &
                           nodes(i-1,j  ,k  )%z + nodes(i-1,j-1,k  )%z)
            xm = 0.25_R8P*(nodes(i  ,j  ,k-1)%x + nodes(i  ,j-1,k-1)%x + &
                           nodes(i-1,j  ,k-1)%x + nodes(i-1,j-1,k-1)%x)
            ym = 0.25_R8P*(nodes(i  ,j  ,k-1)%y + nodes(i  ,j-1,k-1)%y + &
                           nodes(i-1,j  ,k-1)%y + nodes(i-1,j-1,k-1)%y)
            zm = 0.25_R8P*(nodes(i  ,j  ,k-1)%z + nodes(i  ,j-1,k-1)%z + &
                           nodes(i-1,j  ,k-1)%z + nodes(i-1,j-1,k-1)%z)

            Vx = Vx + xp*NFkS(i,j,k)%x - xm*NFkS(i,j,k-1)%x
            Vy = Vy + yp*NFkS(i,j,k)%y - ym*NFkS(i,j,k-1)%y
            Vz = Vz + zp*NFkS(i,j,k)%z - zm*NFkS(i,j,k-1)%z

            volume(i,j,k) = max(Vx,Vy,Vz)
         enddo
      enddo
   enddo
   !$omp end parallel
   endassociate
   endsubroutine compute_metrics

   pure subroutine compute_patches_extents(self, tcc, patch, offset)
   !< Compute the patches extents of given patch boundary conditions or the use the whole block.
   class(block_grd_object),   intent(inout)        :: self               !< Block data.
   integer(I4P),              intent(in), optional :: tcc(1-self%gc(1):,&
                                                          1-self%gc(3):,&
                                                          1-self%gc(5):) !< Cells type.
   integer(I4P),              intent(in), optional :: patch              !< Patch bc to be found.
   integer(I4P),              intent(in), optional :: offset             !< Offset from patch.
   integer(I4P)                                    :: np                 !< Number of patches found.
   integer(I4P)                                    :: offset_            !< Offset from patch, local variable.
   integer(I4P), parameter                         :: ci1=1,  ci2=2,  &
                                                      cj1=3,  cj2=4,  &
                                                      ck1=5,  ck2=6,  &
                                                      ni1=7,  ni2=8,  &
                                                      nj1=9,  nj2=10, &
                                                      nk1=11, nk2=12   !< Named indexes.

   associate(Ni=>self%Ni, Nj=>self%Nj, Nk=>self%Nk)
   if (allocated(self%patches_extents)) deallocate(self%patches_extents)
   if ((.not.present(tcc)).or.(.not.present(patch))) then
      ! whole block
      allocate(self%patches_extents(1:1, 0:12))
      self%patches_extents(1, 0 ) = -999 ! whole block, no patches
      self%patches_extents(1, 1 ) = 1 ; self%patches_extents(1, 2 ) = Ni
      self%patches_extents(1, 3 ) = 1 ; self%patches_extents(1, 4 ) = Nj
      self%patches_extents(1, 5 ) = 1 ; self%patches_extents(1, 6 ) = Nk
      self%patches_extents(1, 7 ) = 0 ; self%patches_extents(1, 8 ) = Ni
      self%patches_extents(1, 9 ) = 0 ; self%patches_extents(1, 10) = Nj
      self%patches_extents(1, 11) = 0 ; self%patches_extents(1, 12) = Nk
      return
   elseif ((present(tcc)).and.(present(patch))) then
      ! find the queried patches into the block
      offset_ = 0 ; if (present(offset)) offset_ = offset
      np = 0
      if (any(tcc(0   , :   , :   ) == patch).or.any(tcc(0   , :   , :   ) == patch+10_I4P)) np = np + 1
      if (any(tcc(Ni+1, :   , :   ) == patch).or.any(tcc(Ni+1, :   , :   ) == patch+10_I4P)) np = np + 1
      if (any(tcc(:   , 0   , :   ) == patch).or.any(tcc(:   , 0   , :   ) == patch+10_I4P)) np = np + 1
      if (any(tcc(:   , Nj+1, :   ) == patch).or.any(tcc(:   , Nj+1, :   ) == patch+10_I4P)) np = np + 1
      if (any(tcc(:   , :   , 0   ) == patch).or.any(tcc(:   , :   , 0   ) == patch+10_I4P)) np = np + 1
      if (any(tcc(:   , :   , Nk+1) == patch).or.any(tcc(:   , :   , Nk+1) == patch+10_I4P)) np = np + 1
      if (np > 0) then
         allocate(self%patches_extents(1:np, 0:12))
         np = 0
         if (any(tcc(0,:,:) == patch).or.any(tcc(0,:,:) == patch+10_I4P)) then
            np = np + 1
            self%patches_extents(np, 0) = 1
            self%patches_extents(np, ci1) = 1 + offset_ ; self%patches_extents(np, ci2) = 1 + offset_
            self%patches_extents(np, cj1) = 1           ; self%patches_extents(np, cj2) = Nj
            self%patches_extents(np, ck1) = 1           ; self%patches_extents(np, ck2) = Nk
            self%patches_extents(np, ni1) = 0 + offset_ ; self%patches_extents(np, ni2) = 0 + offset_
            self%patches_extents(np, nj1) = 0           ; self%patches_extents(np, nj2) = Nj
            self%patches_extents(np, nk1) = 0           ; self%patches_extents(np, nk2) = Nk
         endif
         if (any(tcc(Ni+1,:,:) == patch).or.any(tcc(Ni+1,:,:) == patch+10_I4P)) then
            np = np + 1
            self%patches_extents(np, 0) = 2
            self%patches_extents(np, ci1) = Ni - offset_ ; self%patches_extents(np, ci2) = Ni - offset_
            self%patches_extents(np, cj1) = 1            ; self%patches_extents(np, cj2) = Nj
            self%patches_extents(np, ck1) = 1            ; self%patches_extents(np, ck2) = Nk
            self%patches_extents(np, ni1) = Ni - offset_ ; self%patches_extents(np, ni2) = Ni - offset_
            self%patches_extents(np, nj1) = 0            ; self%patches_extents(np, nj2) = Nj
            self%patches_extents(np, nk1) = 0            ; self%patches_extents(np, nk2) = Nk
         endif
         if (any(tcc(:,0,:) == patch).or.any(tcc(:,0,:) == patch+10_I4P)) then
            np = np + 1
            self%patches_extents(np, 0) = 3
            self%patches_extents(np, ci1) = 1           ; self%patches_extents(np, ci2) = Ni
            self%patches_extents(np, cj1) = 1 + offset_ ; self%patches_extents(np, cj2) = 1 + offset_
            self%patches_extents(np, ck1) = 1           ; self%patches_extents(np, ck2) = Nk
            self%patches_extents(np, ni1) = 0           ; self%patches_extents(np, ni2) = Ni
            self%patches_extents(np, nj1) = 0 + offset_ ; self%patches_extents(np, nj2) = 0 + offset_
            self%patches_extents(np, nk1) = 0           ; self%patches_extents(np, nk2) = Nk
         endif
         if (any(tcc(:,Nj+1,:) == patch).or.any(tcc(:,Nj+1,:) == patch+10_I4P)) then
            np = np + 1
            self%patches_extents(np, 0) = 4
            self%patches_extents(np, ci1) = 1            ; self%patches_extents(np, ci2) = Ni
            self%patches_extents(np, cj1) = Nj - offset_ ; self%patches_extents(np, cj2) = Nj - offset_
            self%patches_extents(np, ck1) = 1            ; self%patches_extents(np, ck2) = Nk
            self%patches_extents(np, ni1) = 0            ; self%patches_extents(np, ni2) = Ni
            self%patches_extents(np, nj1) = Nj - offset_ ; self%patches_extents(np, nj2) = Nj - offset_
            self%patches_extents(np, nk1) = 0            ; self%patches_extents(np, nk2) = Nk
         endif
         if (any(tcc(:,:,0) == patch).or.any(tcc(:,:,0) == patch+10_I4P)) then
            np = np + 1
            self%patches_extents(np, 0) = 5
            self%patches_extents(np, ci1) = 1           ; self%patches_extents(np, ci2) = Ni
            self%patches_extents(np, cj1) = 1           ; self%patches_extents(np, cj2) = Nj
            self%patches_extents(np, ck1) = 1 + offset_ ; self%patches_extents(np, ck2) = 1 + offset_
            self%patches_extents(np, ni1) = 0           ; self%patches_extents(np, ni2) = Ni
            self%patches_extents(np, nj1) = 0           ; self%patches_extents(np, nj2) = Nj
            self%patches_extents(np, nk1) = 0 + offset_ ; self%patches_extents(np, nk2) = 0 + offset_
         endif
         if (any(tcc(:,:,Nk+1) == patch).or.any(tcc(:,:,Nk+1) == patch+10_I4P)) then
            np = np + 1
            self%patches_extents(np, 0) = 6
            self%patches_extents(np, ci1) = 1            ; self%patches_extents(np, ci2) = Ni
            self%patches_extents(np, cj1) = 1            ; self%patches_extents(np, cj2) = Nj
            self%patches_extents(np, ck1) = Nk - offset_ ; self%patches_extents(np, ck2) = Nk - offset_
            self%patches_extents(np, ni1) = 0            ; self%patches_extents(np, ni2) = Ni
            self%patches_extents(np, nj1) = 0            ; self%patches_extents(np, nj2) = Nj
            self%patches_extents(np, nk1) = Nk - offset_ ; self%patches_extents(np, nk2) = Nk - offset_
         endif
      endif
   endif
   endassociate
   endsubroutine compute_patches_extents

   subroutine correct_metrics_bc(self, tcc)
   !< Correct metrics of Boundary Conditions (ghost) cells.
   class(block_grd_object), intent(inout) :: self               !< Block data.
   integer(I4P),            intent(in)    :: tcc(1-self%gc(1):,&
                                                 1-self%gc(3):,&
                                                 1-self%gc(5):) !< Cells type.
   logical                                :: bc_correct         !< Flag for checking bc.
   logical                                :: bc_wall            !< Flag for checking wall bc.
   real(R8P)                              :: tm                 !< Tangential metrics parameter (-1 for wall-type bc).
   real(R8P)                              :: sn                 !< Normal     metrics coefficient correction.
   integer(I4P)                           :: i, j, k            !< Counter.
   integer(I4P), parameter                :: wall=1             !< Wall boundary condition.
   integer(I4P), parameter                :: simmetry=2         !< Simmetry boundary condition.
   integer(I4P), parameter                :: movingwall=10      !< Moving wall boundary condition.
   integer(I4P), parameter                :: passivewall=11     !< Passive wall boundary condition.

   associate(Ni=>self%Ni, Nj=>self%Nj, Nk=>self%Nk, gc=>self%gc, &
             NFi=>self%NFi, NFj=>self%NFj, NFk=>self%NFk,        &
             NFiS=>self%NFiS, NFjS=>self%NFjS, NFkS=>self%NFkS,  &
             volume=>self%volume)
   !$omp parallel default(none)                  &
   !$omp private(i,j,k,bc_correct,bc_wall,tm,sn) &
   !$omp shared(Ni,Nj,Nk,tcc,NFi,NFj,NFk,NFiS,NFjS,NFkS,volume)
   ! left i
   !$omp do
   do k=1,Nk
      do j=1,Nj
         i = tcc(0,j,k)
         bc_correct = ((i<0).OR.(volume(0,j,k)<(0.2_R_P*volume(1,j,k))))
         bc_wall    = ((i==wall).OR.(i==simmetry).OR.(i==movingwall).OR.(i==passivewall))
         tm = 1._R_P
         if (bc_wall) tm = -1._R_P
         if (bc_correct) then
            ! normal metrics
            sn = 2._R_P*(NFiS(1,j,k).dot.NFi(0,j,k))
            NFiS( -1,j,  k  ) = -NFiS(1,j,k)+sn*NFi(0,j,k)
            ! tangential metrics
            NFjS(  0,j  ,k  ) = tm*NFjS(1,j  ,k  )
            NFjS(  0,j-1,k  ) = tm*NFjS(1,j-1,k  )
            NFjS(  0,j  ,k-1) = tm*NFjS(1,j  ,k-1)
            NFjS(  0,j-1,k-1) = tm*NFjS(1,j-1,k-1)

            NFkS(  0,j  ,k  ) = tm*NFkS(1,j  ,k  )
            NFkS(  0,j-1,k  ) = tm*NFkS(1,j-1,k  )
            NFkS(  0,j  ,k-1) = tm*NFkS(1,j  ,k-1)
            NFkS(  0,j-1,k-1) = tm*NFkS(1,j-1,k-1)
            ! volume
            volume(0,j,  k  ) = volume(1,j,k)
         endif
      enddo
   enddo
   ! right i
   !$omp do
   do k=1,Nk
      do j=1,Nj
         i = tcc(Ni+1,j,k)
         bc_correct = ((i<0).OR.(volume(Ni+1,j,k)<(0.2_R_P*volume(Ni,j,k))))
         bc_wall    = ((i==wall).OR.(i==simmetry).OR.(i==movingwall).OR.(i==passivewall))
         tm = 1._R_P
         if (bc_wall) tm = -1._R_P
         if (bc_correct) then
            ! normal metrics
            sn = 2._R_P*(NFiS(Ni-1,j,k).dot.NFi(Ni,j,k))
            NFiS(  Ni+1,j,  k  ) = -NFiS(Ni-1,j,k)+sn*NFi(Ni,j,k)
            ! tangential metrics
            NFjS(  Ni+1,j  ,k  ) = tm*NFjS(Ni,j  ,k  )
            NFjS(  Ni+1,j-1,k  ) = tm*NFjS(Ni,j-1,k  )
            NFjS(  Ni+1,j  ,k-1) = tm*NFjS(Ni,j  ,k-1)
            NFjS(  Ni+1,j-1,k-1) = tm*NFjS(Ni,j-1,k-1)

            NFkS(  Ni+1,j  ,k  ) = tm*NFkS(Ni,j  ,k  )
            NFkS(  Ni+1,j-1,k  ) = tm*NFkS(Ni,j-1,k  )
            NFkS(  Ni+1,j  ,k-1) = tm*NFkS(Ni,j  ,k-1)
            NFkS(  Ni+1,j-1,k-1) = tm*NFkS(Ni,j-1,k-1)
            ! volume
            volume(Ni+1,j,  k  ) = volume(Ni,j,k)
         endif
      enddo
   enddo
   ! left j
   !$omp do
   do k=1,Nk
      do i=1,Ni
         j = tcc(i,0,k)
         bc_correct = ((j<0).OR.(volume(i,0,k)<(0.2_R_P*volume(i,1,k))))
         bc_wall    = ((j==wall).OR.(j==simmetry).OR.(j==movingwall).OR.(j==passivewall))
         tm = 1._R_P
         if (bc_wall) tm = -1._R_P
         if (bc_correct) then
            ! normal metrics
            sn = 2._R_P*(NFjS(i,1,k).dot.NFj(i,0,k))
            NFjS(  i, -1,k  ) = -NFjS(i,1,k)+sn*NFj(i,0,k)
            ! tangential metrics
            NFiS(  i  ,0,k  ) = tm*NFiS(i  ,1,k  )
            NFiS(  i-1,0,k  ) = tm*NFiS(i-1,1,k  )
            NFiS(  i  ,0,k-1) = tm*NFiS(i  ,1,k-1)
            NFiS(  i-1,0,k-1) = tm*NFiS(i-1,1,k-1)

            NFkS(  i  ,0,k  ) = tm*NFkS(i  ,1,k  )
            NFkS(  i-1,0,k  ) = tm*NFkS(i-1,1,k  )
            NFkS(  i  ,0,k-1) = tm*NFkS(i  ,1,k-1)
            NFkS(  i-1,0,k-1) = tm*NFkS(i-1,1,k-1)
            ! volume
            volume(i,  0,k  ) = volume(i,1,k)
         endif
      enddo
   enddo
   ! right j
   !$omp do
   do k=1,Nk
      do i=1,Ni
         j = tcc(i,Nj+1,k)
         bc_correct = ((j<0).OR.(volume(i,Nj+1,k)<(0.2_R_P*volume(i,Nj,k))))
         bc_wall    = ((j==wall).OR.(j==simmetry).OR.(j==movingwall).OR.(j==passivewall))
         tm = 1._R_P
         if (bc_wall) tm = -1._R_P
         if (bc_correct) then
            ! normal metrics
            sn = 2._R_P*(NFjS(i,Nj-1,k).dot.NFj(i,Nj,k))
            NFjS(  i,Nj+1,  k  ) = -NFjS(i,Nj-1,k)+sn*NFj(i,Nj,k)
            ! tangential metrics
            NFiS(  i  ,Nj+1,k  ) = tm*NFiS(i  ,Nj,k  )
            NFiS(  i-1,Nj+1,k  ) = tm*NFiS(i-1,Nj,k  )
            NFiS(  i  ,Nj+1,k-1) = tm*NFiS(i  ,Nj,k-1)
            NFiS(  i-1,Nj+1,k-1) = tm*NFiS(i-1,Nj,k-1)

            NFkS(  i  ,Nj+1,k  ) = tm*NFkS(i  ,Nj,k  )
            NFkS(  i-1,Nj+1,k  ) = tm*NFkS(i-1,Nj,k  )
            NFkS(  i  ,Nj+1,k-1) = tm*NFkS(i  ,Nj,k-1)
            NFkS(  i-1,Nj+1,k-1) = tm*NFkS(i-1,Nj,k-1)
            ! volume
            volume(i,Nj+1,  k  ) = volume(i,Nj,k)
         endif
      enddo
   enddo
   ! left k
   !$omp do
   do j=1,Nj
      do i=1,Ni
         k = tcc(i,j,0)
         bc_correct = ((k<0).OR.(volume(i,j,0)<(0.2_R_P*volume(i,j,1))))
         bc_wall    = ((k==wall).OR.(k==simmetry).OR.(k==movingwall).OR.(k==passivewall))
         tm = 1._R_P
         if (bc_wall) tm = -1._R_P
         if (bc_correct) then
            ! normal metrics
            sn = 2._R_P*(NFkS(i,j,1).dot.NFk(i,j,0))
            NFkS(  i,  j, -1) = -NFkS(i,j,1)+sn*NFk(i,j,0)
            ! tangential metrics
            NFiS(  i  ,j  ,0) = tm*NFiS(i  ,j  ,1)
            NFiS(  i-1,j  ,0) = tm*NFiS(i-1,j  ,1)
            NFiS(  i  ,j-1,0) = tm*NFiS(i  ,j-1,1)
            NFiS(  i-1,j-1,0) = tm*NFiS(i-1,j-1,1)

            NFjS(  i  ,j  ,0) = tm*NFjS(i  ,j  ,1)
            NFjS(  i-1,j  ,0) = tm*NFjS(i-1,j  ,1)
            NFjS(  i  ,j-1,0) = tm*NFjS(i  ,j-1,1)
            NFjS(  i-1,j-1,0) = tm*NFjS(i-1,j-1,1)
            ! volume
            volume(i,  j,  0) = volume(i,j,1)
         endif
      enddo
   enddo
   ! right k
   !$omp do
   do j=1,Nj
      do i=1,Ni
         k = tcc(i,j,Nk+1)
         bc_correct = ((k<0).OR.(volume(i,j,Nk+1)<(0.2_R_P*volume(i,j,Nk))))
         bc_wall    = ((k==wall).OR.(k==simmetry).OR.(k==movingwall).OR.(k==passivewall))
         tm = 1._R_P
         if (bc_wall) tm = -1._R_P
         if (bc_correct) then
            ! normal metrics
            sn = 2._R_P*(NFkS(i,j,Nk-1).dot.NFk(i,j,Nk))
            NFkS(  i,  j,  Nk+1) = -NFkS(i,j,Nk-1)+sn*NFk(i,j,Nk)
            ! tangential metrics
            NFiS(  i  ,j  ,Nk+1) = tm*NFiS(i  ,j  ,Nk)
            NFiS(  i-1,j  ,Nk+1) = tm*NFiS(i-1,j  ,Nk)
            NFiS(  i  ,j-1,Nk+1) = tm*NFiS(i  ,j-1,Nk)
            NFiS(  i-1,j-1,Nk+1) = tm*NFiS(i-1,j-1,Nk)

            NFjS(  i  ,j  ,Nk+1) = tm*NFjS(i  ,j  ,Nk)
            NFjS(  i-1,j  ,Nk+1) = tm*NFjS(i-1,j  ,Nk)
            NFjS(  i  ,j-1,Nk+1) = tm*NFjS(i  ,j-1,Nk)
            NFjS(  i-1,j-1,Nk+1) = tm*NFjS(i-1,j-1,Nk)
            ! volume
            volume(i,  j,  Nk+1) = volume(i,j,Nk)
         endif
      enddo
   enddo
   !$omp end parallel
   endassociate
   endsubroutine correct_metrics_bc

   elemental subroutine destroy(self)
   !< Destroy dynamic memory.
   class(block_grd_object), intent(inout) :: self !< Block data.

   self%Ni = 0
   self%Nj = 0
   self%Nk = 0
   self%gc = 2
   if (allocated(self%nodes)) deallocate(self%nodes)
   if (allocated(self%centers)) deallocate(self%centers)
   if (allocated(self%extents)) deallocate(self%extents)
   if (allocated(self%sub_extents)) deallocate(self%sub_extents)
   if (allocated(self%sub_ijk_extents)) deallocate(self%sub_ijk_extents)
   if (allocated(self%NFi)) deallocate(self%NFi )
   if (allocated(self%NFj)) deallocate(self%NFj )
   if (allocated(self%NFk)) deallocate(self%NFk )
   if (allocated(self%NFiS)) deallocate(self%NFiS)
   if (allocated(self%NFjS)) deallocate(self%NFjS)
   if (allocated(self%NFkS)) deallocate(self%NFkS)
   if (allocated(self%Si)) deallocate(self%Si)
   if (allocated(self%Sj)) deallocate(self%Sj)
   if (allocated(self%Sk)) deallocate(self%Sk)
   if (allocated(self%volume)) deallocate(self%volume)
   if (allocated(self%patches_extents)) deallocate(self%patches_extents)
   self%is_loaded = .false.
   endsubroutine destroy

   subroutine load_dimensions(self, file_unit, is_centers_to_allocate, is_extents_to_allocate, is_metrics_to_allocate)
   !< Load block dimensions from file.
   !<
   !< @note The grd file must be already open and the current record-index must be at the proper block dimensions record.
   class(block_grd_object), intent(inout)        :: self                   !< Block data.
   integer(I4P),            intent(in)           :: file_unit              !< Logical unit of grd file.
   logical,                 intent(in), optional :: is_centers_to_allocate !< Flag to allocate also centers array.
   logical,                 intent(in), optional :: is_extents_to_allocate !< Flag to allocate also extents arrays.
   logical,                 intent(in), optional :: is_metrics_to_allocate !< Flag to allocate also metrics arrays.

   call self%destroy
   read(file_unit, end=10, err=10) self%Ni, self%Nj, self%Nk, self%gc
   10 call self%alloc(is_centers_to_allocate=is_centers_to_allocate, &
                      is_extents_to_allocate=is_extents_to_allocate, &
                      is_metrics_to_allocate=is_metrics_to_allocate)
   endsubroutine load_dimensions

   subroutine load_nodes(self, file_unit)
   !< Load block nodes from file.
   !<
   !< @note The grd file must be already open and the current record-index must be at the proper block nodes record.
   !<
   !< @note If *centers* and *extents* are allocated they are computed from nodes values.
   class(block_grd_object), intent(inout) :: self      !< Block data.
   integer(I4P),            intent(in)    :: file_unit !< Logical unit of grd file.
   integer(I4P)                           :: i         !< Counter.
   integer(I4P)                           :: j         !< Counter.
   integer(I4P)                           :: k         !< Counter.

   associate(Ni => self%Ni, Nj => self%Nj, Nk => self%Nk, gc => self%gc, nodes => self%nodes)
      read(file_unit)(((nodes(i, j, k)%x, i=0-gc(1), Ni + gc(2)), j=0-gc(3), Nj + gc(4)), k=0-gc(5), Nk + gc(6))
      read(file_unit)(((nodes(i, j, k)%y, i=0-gc(1), Ni + gc(2)), j=0-gc(3), Nj + gc(4)), k=0-gc(5), Nk + gc(6))
      read(file_unit)(((nodes(i, j, k)%z, i=0-gc(1), Ni + gc(2)), j=0-gc(3), Nj + gc(4)), k=0-gc(5), Nk + gc(6))
   endassociate
   self%is_loaded = .true.

   if (allocated(self%centers)) then
      do k=1 - self%gc(5), self%Nk + self%gc(6)
         do j=1 - self%gc(3), self%Nj + self%gc(4)
            do i=1 - self%gc(1), self%Ni + self%gc(2)
               self%centers(i, j, k) = (self%nodes(i,   j,   k  ) + &
                                        self%nodes(i-1, j,   k  ) + &
                                        self%nodes(i  , j-1, k  ) + &
                                        self%nodes(i  , j  , k-1) + &
                                        self%nodes(i-1, j-1, k-1) + &
                                        self%nodes(i  , j-1, k-1) + &
                                        self%nodes(i-1, j  , k-1) + &
                                        self%nodes(i-1, j-1, k  )) * 0.125_R8P
            enddo
         enddo
      enddo
   endif

   if (allocated(self%extents)) self%extents = compute_extents(i_extents=[0, self%Ni], &
                                                               j_extents=[0, self%Nj], &
                                                               k_extents=[0, self%Nk])

   if (allocated(self%sub_extents).and.allocated(self%sub_ijk_extents)) then
      self%sub_ijk_extents(1, 1, :) = [0, self%Ni/2]
      self%sub_ijk_extents(1, 2, :) = [0, self%Nj/2]
      self%sub_ijk_extents(1, 3, :) = [0, self%Nk/2]
      self%sub_ijk_extents(2, 1, :) = [self%Ni/2, self%Ni  ]
      self%sub_ijk_extents(2, 2, :) = [0        , self%Nj/2]
      self%sub_ijk_extents(2, 3, :) = [0        , self%Nk/2]
      self%sub_ijk_extents(3, 1, :) = [0        , self%Ni/2]
      self%sub_ijk_extents(3, 2, :) = [self%Nj/2, self%Nj  ]
      self%sub_ijk_extents(3, 3, :) = [0        , self%Nk/2]
      self%sub_ijk_extents(4, 1, :) = [self%Ni/2, self%Ni  ]
      self%sub_ijk_extents(4, 2, :) = [self%Nj/2, self%Nj  ]
      self%sub_ijk_extents(4, 3, :) = [0        , self%Nk/2]
      self%sub_ijk_extents(5, 1, :) = [0        , self%Ni/2]
      self%sub_ijk_extents(5, 2, :) = [0        , self%Nj/2]
      self%sub_ijk_extents(5, 3, :) = [self%Nk/2, self%Nk  ]
      self%sub_ijk_extents(6, 1, :) = [self%Ni/2, self%Ni  ]
      self%sub_ijk_extents(6, 2, :) = [0        , self%Nj/2]
      self%sub_ijk_extents(6, 3, :) = [self%Nk/2, self%Nk  ]
      self%sub_ijk_extents(7, 1, :) = [0        , self%Ni/2]
      self%sub_ijk_extents(7, 2, :) = [self%Nj/2, self%Nj  ]
      self%sub_ijk_extents(7, 3, :) = [self%Nk/2, self%Nk  ]
      self%sub_ijk_extents(8, 1, :) = [self%Ni/2, self%Ni]
      self%sub_ijk_extents(8, 2, :) = [self%Nj/2, self%Nj]
      self%sub_ijk_extents(8, 3, :) = [self%Nk/2, self%Nk]

      self%sub_extents(1, :) = compute_extents(i_extents=[0, self%Ni/2], &
                                               j_extents=[0, self%Nj/2], &
                                               k_extents=[0, self%Nk/2])
      self%sub_extents(2, :) = compute_extents(i_extents=[self%Ni/2, self%Ni  ], &
                                               j_extents=[0        , self%Nj/2], &
                                               k_extents=[0        , self%Nk/2])
      self%sub_extents(3, :) = compute_extents(i_extents=[0        , self%Ni/2], &
                                               j_extents=[self%Nj/2, self%Nj  ], &
                                               k_extents=[0        , self%Nk/2])
      self%sub_extents(4, :) = compute_extents(i_extents=[self%Ni/2, self%Ni  ], &
                                               j_extents=[self%Nj/2, self%Nj  ], &
                                               k_extents=[0        , self%Nk/2])
      self%sub_extents(5, :) = compute_extents(i_extents=[0        , self%Ni/2], &
                                               j_extents=[0        , self%Nj/2], &
                                               k_extents=[self%Nk/2, self%Nk  ])
      self%sub_extents(6, :) = compute_extents(i_extents=[self%Ni/2, self%Ni  ], &
                                               j_extents=[0        , self%Nj/2], &
                                               k_extents=[self%Nk/2, self%Nk  ])
      self%sub_extents(7, :) = compute_extents(i_extents=[0        , self%Ni/2], &
                                               j_extents=[self%Nj/2, self%Nj  ], &
                                               k_extents=[self%Nk/2, self%Nk  ])
      self%sub_extents(8, :) = compute_extents(i_extents=[self%Ni/2, self%Ni], &
                                               j_extents=[self%Nj/2, self%Nj], &
                                               k_extents=[self%Nk/2, self%Nk])
   endif

   contains
      function compute_extents(i_extents, j_extents, k_extents) result(extents)
      !< Compute (sub)block extents provided indexes extents.
      integer(I4P), intent(in) :: i_extents(2) !< Index "i" extents (min, max) of the (sub)block.
      integer(I4P), intent(in) :: j_extents(2) !< Index "j" extents (min, max) of the (sub)block.
      integer(I4P), intent(in) :: k_extents(2) !< Index "k" extents (min, max) of the (sub)block.
      type(vector)             :: extents(2)   !< (Sub)block extents.

      extents(1)%x = minval(self%nodes(i_extents(1):i_extents(2), j_extents(1):j_extents(2), k_extents(1):k_extents(2))%x)
      extents(2)%x = maxval(self%nodes(i_extents(1):i_extents(2), j_extents(1):j_extents(2), k_extents(1):k_extents(2))%x)

      extents(1)%y = minval(self%nodes(i_extents(1):i_extents(2), j_extents(1):j_extents(2), k_extents(1):k_extents(2))%y)
      extents(2)%y = maxval(self%nodes(i_extents(1):i_extents(2), j_extents(1):j_extents(2), k_extents(1):k_extents(2))%y)

      extents(1)%z = minval(self%nodes(i_extents(1):i_extents(2), j_extents(1):j_extents(2), k_extents(1):k_extents(2))%z)
      extents(2)%z = maxval(self%nodes(i_extents(1):i_extents(2), j_extents(1):j_extents(2), k_extents(1):k_extents(2))%z)
      endfunction compute_extents
   endsubroutine load_nodes

   subroutine save_dimensions(self, file_unit, save_gc)
   !< Save block dimensions into file.
   !<
   !< @note The grd file must be already open and the current record-index must be at the proper block dimensions record.
   class(block_grd_object), intent(in)           :: self      !< Block data.
   integer(I4P),            intent(in)           :: file_unit !< Logical unit of grd file.
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

   subroutine save_nodes(self, file_unit)
   !< Save block nodes into file.
   !<
   !< @note The grd file must be already open and the current record-index must be at the proper block nodes record.
   class(block_grd_object), intent(in) :: self      !< Block data.
   integer(I4P),            intent(in) :: file_unit !< Logical unit of grd file.
   integer(I4P)                        :: i         !< Counter.
   integer(I4P)                        :: j         !< Counter.
   integer(I4P)                        :: k         !< Counter.

   associate(Ni => self%Ni, Nj => self%Nj, Nk => self%Nk, gc => self%gc, nodes => self%nodes)
      write(file_unit)(((nodes(i, j, k)%x, i=0-gc(1), Ni + gc(2)), j=0-gc(3), Nj + gc(4)), k=0-gc(5), Nk + gc(6))
      write(file_unit)(((nodes(i, j, k)%y, i=0-gc(1), Ni + gc(2)), j=0-gc(3), Nj + gc(4)), k=0-gc(5), Nk + gc(6))
      write(file_unit)(((nodes(i, j, k)%z, i=0-gc(1), Ni + gc(2)), j=0-gc(3), Nj + gc(4)), k=0-gc(5), Nk + gc(6))
   endassociate
   endsubroutine save_nodes

   elemental subroutine traslate(self, traslation)
   !< Traslate block nodes by a given traslation vector.
   class(block_grd_object), intent(inout) :: self       !< Block data.
   type(vector),            intent(in)    :: traslation !< Traslation vector.

   if (allocated(self%nodes)) self%nodes = self%nodes + traslation
   endsubroutine traslate
endmodule xview_block_grd_object
