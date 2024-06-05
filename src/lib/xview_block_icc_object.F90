!< xview, block (icc) class definition.
module xview_block_icc_object
!< xview, block (icc) class definition.

use penf

implicit none
private
public :: block_icc_object

type :: block_icc_object
   !< Block (icc) class definition.
   integer(I4P)              :: Ni=0                     !< Number of cells in i direction.
   integer(I4P)              :: Nj=0                     !< Number of cells in j direction.
   integer(I4P)              :: Nk=0                     !< Number of cells in k direction.
   integer(I4P)              :: gc(6)=[2, 2, 2, 2, 2, 2] !< Number of ghost cells.
   integer(I4P), allocatable :: icc(:,:,:)               !< Cell centered icc values.
   real(R4P),    allocatable :: rcc(:,:,:)               !< Cell/nodes centered rcc values.
   integer(I4P), allocatable :: tcc(:,:,:)               !< Cell/nodes centered tcc values.
   logical                   :: is_loaded=.false.        !< Flag for checking if the block is loaded.
   contains
      ! public methods
      procedure, pass(self) :: alloc                !< Allocate dynamic memory.
      procedure, pass(self) :: compute_cc           !< Compute cells/nodes centered rcc and tcc from unstructured rcc.
      procedure, pass(self) :: destroy              !< Destroy dynamic memory.
      procedure, pass(self) :: interpolate_at_nodes !< Interpolate rcc at nodes.
      procedure, pass(self) :: load_dimensions      !< Load block dimensions from file.
      procedure, pass(self) :: load_icc             !< Load block icc from file.
      procedure, pass(self) :: get_patches_extents  !< Return the patches extents of given patch boundary conditions.
      procedure, pass(self) :: save_dimensions      !< Save block dimensions into file.
      procedure, pass(self) :: save_icc             !< Save block icc into file.
endtype block_icc_object

contains
   elemental subroutine alloc(self)
   !< Allocate dynamic memory.
   class(block_icc_object), intent(inout) :: self !< Block data.

   allocate(self%icc(1-self%gc(1):self%Ni+self%gc(2), 1-self%gc(3):self%Nj+self%gc(4), 1-self%gc(5):self%Nk+self%gc(6)))
   allocate(self%rcc(1-self%gc(1):self%Ni+self%gc(2), 1-self%gc(3):self%Nj+self%gc(4), 1-self%gc(5):self%Nk+self%gc(6)))
   allocate(self%tcc(1-self%gc(1):self%Ni+self%gc(2), 1-self%gc(3):self%Nj+self%gc(4), 1-self%gc(5):self%Nk+self%gc(6)))
   endsubroutine alloc

   subroutine compute_cc(self, rcc, is_cell_centered)
   !< Compute cells/nodes centered rcc and tcc from unstructured rcc.
   class(block_icc_object), intent(inout)        :: self             !< Block data.
   real(R4P),               intent(in)           :: rcc(1:)          !< Unstructured rcc.
   logical,                 intent(in), optional :: is_cell_centered !< Define variables at cell centers or nodes.
   integer(I4P)                                  :: i, j, k          !< Counter.

   self%rcc = 0._R8P
   do k=1, self%Nk
      do j=1, self%Nj
         do i=1, self%Ni
            if (self%icc(i, j, k) > 0) self%rcc(i,j,k) = nint(rcc(self%icc(i,j,k)))*1._R4P
            if (self%rcc(i, j, k) > 28._R4P) self%rcc(i, j, k) = 0._R4P
         enddo
      enddo
   enddo
   self%tcc = 0
   do k=0, self%Nk+1
      do j=0, self%Nj+1
         do i=0, self%Ni+1
            if (self%icc(i,j,k)>0) self%tcc(i,j,k) = abs(nint(rcc(self%icc(i,j,k))))
         enddo
      enddo
   enddo
   if (present(is_cell_centered)) then
      if (.not.is_cell_centered) call self%interpolate_at_nodes
   endif
   endsubroutine compute_cc

   elemental subroutine destroy(self)
   !< Destroy dynamic memory.
   class(block_icc_object), intent(inout) :: self !< Block data.

   self%Ni = 0
   self%Nj = 0
   self%Nk = 0
   self%gc = 2
   if (allocated(self%icc)) deallocate(self%icc)
   if (allocated(self%rcc)) deallocate(self%rcc)
   if (allocated(self%tcc)) deallocate(self%tcc)
   self%is_loaded = .false.
   endsubroutine destroy

   pure subroutine interpolate_at_nodes(self)
   !< Compute rcc at nodes.
   class(block_icc_object), intent(inout) :: self       !< Block data.
   real(R4P), allocatable                 :: tmp(:,:,:) !< Temporary storage.
   integer(I4P)                           :: i, j, k    !< Counter.

   allocate(tmp(0:self%Ni,0:self%Nj,0:self%Nk))
   self%rcc(     0     ,      :     ,      :     ) = self%rcc(     1 ,      : ,      : )
   self%rcc(self%Ni + 1,      :     ,      :     ) = self%rcc(self%Ni,      : ,      : )
   self%rcc(     :     ,      0     ,      :     ) = self%rcc(     : ,      1 ,      : )
   self%rcc(     :     , self%Nj + 1,      :     ) = self%rcc(     : , self%Nj,      : )
   self%rcc(     :     ,      :     ,      0     ) = self%rcc(     : ,      : ,      1 )
   self%rcc(     :     ,      :     , self%Nk + 1) = self%rcc(     : ,      : , self%Nk)
   do k=0, self%Nk
      do j=0, self%Nj
         do i=0, self%Ni
            tmp(i, j, k) = max(0._R4P, maxval(self%rcc(i:i + 1, j:j + 1, k:k + 1)))
         enddo
      enddo
   enddo
   self%rcc(0:self%Ni,0:self%Nj,0:self%Nk) = tmp
   endsubroutine interpolate_at_nodes

   subroutine load_dimensions(self, file_unit)
   !< Load block dimensions from file.
   !<
   !< @note The icc file must be already open and the current record-index must be at the proper block dimensions record.
   class(block_icc_object), intent(inout) :: self      !< Block data.
   integer(I4P),            intent(in)    :: file_unit !< Logical unit of icc file.

   call self%destroy
   read(file_unit, end=10, err=10) self%Ni, self%Nj, self%Nk, self%gc
   10 call self%alloc
   endsubroutine load_dimensions

   subroutine load_icc(self, file_unit)
   !< Load block icc from file.
   !<
   !< @note The icc file must be already open and the current record-index must be at the proper block icc record.
   class(block_icc_object), intent(inout) :: self      !< Block data.
   integer(I4P),            intent(in)    :: file_unit !< Logical unit of icc file.
   integer(I4P)                           :: i         !< Counter.
   integer(I4P)                           :: j         !< Counter.
   integer(I4P)                           :: k         !< Counter.

   associate(Ni => self%Ni, Nj => self%Nj, Nk => self%Nk, gc => self%gc, icc => self%icc)
      read(file_unit)(((icc(i, j, k), i=1-gc(1), Ni + gc(2)), j=1-gc(3), Nj + gc(4)), k=1-gc(5), Nk + gc(6))
   endassociate
   self%is_loaded = .true.
   endsubroutine load_icc

   pure subroutine get_patches_extents(self, patch, patches_extents, offset)
   !< Return the patches extents of given patch boundary conditions.
   class(block_icc_object),   intent(in)           :: self                 !< Block data.
   integer(I4P),              intent(in)           :: patch                !< Patch bc to be found.
   integer(I4P), allocatable, intent(inout)        :: patches_extents(:,:) !< Patches extents, [np, 13]. The second index means
                                                                           !+ 0  => patch face (1,2,3,4,5,6);
                                                                           !+ 1  => cell i-min;
                                                                           !+ 2  => cell i-max;
                                                                           !+ 3  => cell j-min;
                                                                           !+ 4  => cell j-max;
                                                                           !+ 5  => cell k-min;
                                                                           !+ 6  => cell k-max;
                                                                           !+ 7  => node i-min;
                                                                           !+ 8  => node i-max;
                                                                           !+ 9  => node j-min;
                                                                           !+ 10 => node j-max;
                                                                           !+ 11 => node k-min;
                                                                           !+ 12 => node k-max;
   integer(I4P),              intent(in), optional :: offset               !< Offset from patch.
   integer(I4P)                                    :: np                   !< Number of patches found.
   integer(I4P)                                    :: offset_              !< Offset from patch, local variable.
   integer(I4P), parameter                         :: ci1=1,  ci2=2,  &
                                                      cj1=3,  cj2=4,  &
                                                      ck1=5,  ck2=6,  &
                                                      ni1=7,  ni2=8,  &
                                                      nj1=9,  nj2=10, &
                                                      nk1=11, nk2=12       !< Named indexes.

   if (allocated(patches_extents)) deallocate(patches_extents)
   if (.not.allocated(self%tcc)) return
   offset_ = 0 ; if (present(offset)) offset_ = offset
   np = 0
   associate(tcc=>self%tcc, Ni=>self%Ni, Nj=>self%Nj, Nk=>self%Nk)
      if (any(tcc(0   , :   , :   ) == patch).or.any(tcc(0   , :   , :   ) == patch+10_I4P)) np = np + 1
      if (any(tcc(Ni+1, :   , :   ) == patch).or.any(tcc(Ni+1, :   , :   ) == patch+10_I4P)) np = np + 1
      if (any(tcc(:   , 0   , :   ) == patch).or.any(tcc(:   , 0   , :   ) == patch+10_I4P)) np = np + 1
      if (any(tcc(:   , Nj+1, :   ) == patch).or.any(tcc(:   , Nj+1, :   ) == patch+10_I4P)) np = np + 1
      if (any(tcc(:   , :   , 0   ) == patch).or.any(tcc(:   , :   , 0   ) == patch+10_I4P)) np = np + 1
      if (any(tcc(:   , :   , Nk+1) == patch).or.any(tcc(:   , :   , Nk+1) == patch+10_I4P)) np = np + 1
      if (np > 0) then
         allocate(patches_extents(1:np, 0:12))
         np = 0
         if (any(tcc(0,:,:) == patch).or.any(tcc(0,:,:) == patch+10_I4P)) then
            np = np + 1
            patches_extents(np, 0) = 1
            patches_extents(np, ci1) = 1 + offset_ ; patches_extents(np, ci2) = 1 + offset_
            patches_extents(np, cj1) = 1           ; patches_extents(np, cj2) = Nj
            patches_extents(np, ck1) = 1           ; patches_extents(np, ck2) = Nk
            patches_extents(np, ni1) = 0 + offset_ ; patches_extents(np, ni2) = 0 + offset_
            patches_extents(np, nj1) = 0           ; patches_extents(np, nj2) = Nj
            patches_extents(np, nk1) = 0           ; patches_extents(np, nk2) = Nk
         endif
         if (any(tcc(Ni+1,:,:) == patch).or.any(tcc(Ni+1,:,:) == patch+10_I4P)) then
            np = np + 1
            patches_extents(np, 0) = 2
            patches_extents(np, ci1) = Ni - offset_ ; patches_extents(np, ci2) = Ni - offset_
            patches_extents(np, cj1) = 1            ; patches_extents(np, cj2) = Nj
            patches_extents(np, ck1) = 1            ; patches_extents(np, ck2) = Nk
            patches_extents(np, ni1) = Ni - offset_ ; patches_extents(np, ni2) = Ni - offset_
            patches_extents(np, nj1) = 0            ; patches_extents(np, nj2) = Nj
            patches_extents(np, nk1) = 0            ; patches_extents(np, nk2) = Nk
         endif
         if (any(tcc(:,0,:) == patch).or.any(tcc(:,0,:) == patch+10_I4P)) then
            np = np + 1
            patches_extents(np, 0) = 3
            patches_extents(np, ci1) = 1           ; patches_extents(np, ci2) = Ni
            patches_extents(np, cj1) = 1 + offset_ ; patches_extents(np, cj2) = 1 + offset_
            patches_extents(np, ck1) = 1           ; patches_extents(np, ck2) = Nk
            patches_extents(np, ni1) = 0           ; patches_extents(np, ni2) = Ni
            patches_extents(np, nj1) = 0 + offset_ ; patches_extents(np, nj2) = 0 + offset_
            patches_extents(np, nk1) = 0           ; patches_extents(np, nk2) = Nk
         endif
         if (any(tcc(:,Nj+1,:) == patch).or.any(tcc(:,Nj+1,:) == patch+10_I4P)) then
            np = np + 1
            patches_extents(np, 0) = 4
            patches_extents(np, ci1) = 1            ; patches_extents(np, ci2) = Ni
            patches_extents(np, cj1) = Nj - offset_ ; patches_extents(np, cj2) = Nj - offset_
            patches_extents(np, ck1) = 1            ; patches_extents(np, ck2) = Nk
            patches_extents(np, ni1) = 0            ; patches_extents(np, ni2) = Ni
            patches_extents(np, nj1) = Nj - offset_ ; patches_extents(np, nj2) = Nj - offset_
            patches_extents(np, nk1) = 0            ; patches_extents(np, nk2) = Nk
         endif
         if (any(tcc(:,:,0) == patch).or.any(tcc(:,:,0) == patch+10_I4P)) then
            np = np + 1
            patches_extents(np, 0) = 5
            patches_extents(np, ci1) = 1           ; patches_extents(np, ci2) = Ni
            patches_extents(np, cj1) = 1           ; patches_extents(np, cj2) = Nj
            patches_extents(np, ck1) = 1 + offset_ ; patches_extents(np, ck2) = 1 + offset_
            patches_extents(np, ni1) = 0           ; patches_extents(np, ni2) = Ni
            patches_extents(np, nj1) = 0           ; patches_extents(np, nj2) = Nj
            patches_extents(np, nk1) = 0 + offset_ ; patches_extents(np, nk2) = 0 + offset_
         endif
         if (any(tcc(:,:,Nk+1) == patch).or.any(tcc(:,:,Nk+1) == patch+10_I4P)) then
            np = np + 1
            patches_extents(np, 0) = 6
            patches_extents(np, ci1) = 1            ; patches_extents(np, ci2) = Ni
            patches_extents(np, cj1) = 1            ; patches_extents(np, cj2) = Nj
            patches_extents(np, ck1) = Nk - offset_ ; patches_extents(np, ck2) = Nk - offset_
            patches_extents(np, ni1) = 0            ; patches_extents(np, ni2) = Ni
            patches_extents(np, nj1) = 0            ; patches_extents(np, nj2) = Nj
            patches_extents(np, nk1) = Nk - offset_ ; patches_extents(np, nk2) = Nk - offset_
         endif
      endif
   endassociate
   endsubroutine get_patches_extents

   subroutine save_dimensions(self, file_unit, save_gc)
   !< Save block dimensions into file.
   !<
   !< @note The icc file must be already open and the current record-index must be at the proper block dimensions record.
   class(block_icc_object), intent(in)           :: self      !< Block data.
   integer(I4P),            intent(in)           :: file_unit !< Logical unit of icc file.
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

   subroutine save_icc(self, file_unit)
   !< Save block icc into file.
   !<
   !< @note The icc file must be already open and the current record-index must be at the proper block icc record.
   class(block_icc_object), intent(in) :: self      !< Block data.
   integer(I4P),            intent(in) :: file_unit !< Logical unit of icc file.
   integer(I4P)                        :: i         !< Counter.
   integer(I4P)                        :: j         !< Counter.
   integer(I4P)                        :: k         !< Counter.

   associate(Ni => self%Ni, Nj => self%Nj, Nk => self%Nk, gc => self%gc, icc => self%icc)
      write(file_unit)(((icc(i, j, k), i=1-gc(1), Ni + gc(2)), j=1-gc(3), Nj + gc(4)), k=1-gc(5), Nk + gc(6))
   endassociate
   endsubroutine save_icc
endmodule xview_block_icc_object
