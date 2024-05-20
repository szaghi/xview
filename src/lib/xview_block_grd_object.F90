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
   logical                   :: is_loaded=.false.        !< Flag for checking if the block is loaded.
   contains
      ! public methods
      procedure, pass(self) :: destroy         !< Destroy dynamic memory.
      procedure, pass(self) :: alloc           !< Allocate dynamic memory.
      procedure, pass(self) :: load_dimensions !< Load block dimensions from file.
      procedure, pass(self) :: load_nodes      !< Load block nodes from file.
      procedure, pass(self) :: save_dimensions !< Save block dimensions into file.
      procedure, pass(self) :: save_nodes      !< Save block nodes into file.
      procedure, pass(self) :: traslate        !< Traslate block nodes by a given traslation vector.
endtype block_grd_object

contains
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
   self%is_loaded = .false.
   endsubroutine destroy

   elemental subroutine alloc(self, is_centers_to_allocate, is_extents_to_allocate)
   !< Allocate dynamic memory.
   class(block_grd_object), intent(inout) :: self                   !< Block data.
   logical,          intent(in), optional :: is_centers_to_allocate !< Flag to allocate also centers array.
   logical,          intent(in), optional :: is_extents_to_allocate !< Flag to allocate also extents arrays.

   associate(Ni => self%Ni, Nj => self%Nj, Nk => self%Nk, gc => self%gc)
      allocate(self%nodes(0-gc(1):Ni+gc(2), 0-gc(3):Nj+gc(4), 0-gc(5):Nk+gc(6)))
      if (present(is_centers_to_allocate)) then
         if (is_centers_to_allocate) allocate(self%centers(1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6)))
      endif
   endassociate
   if (present(is_extents_to_allocate)) then
      if (is_extents_to_allocate) then
         allocate(self%extents(1:2)) ! min-max
         allocate(self%sub_extents(1:8, 1:2)) ! eight parts, min-max
         allocate(self%sub_ijk_extents(1:8, 1:3, 1:2))  ! eight parts, ijk, min-max
      endif
   endif
   endsubroutine alloc

   subroutine load_dimensions(self, file_unit, is_centers_to_allocate, is_extents_to_allocate)
   !< Load block dimensions from file.
   !<
   !< @note The grd file must be already open and the current record-index must be at the proper block dimensions record.
   class(block_grd_object), intent(inout)        :: self                   !< Block data.
   integer(I4P),            intent(in)           :: file_unit              !< Logical unit of grd file.
   logical,                 intent(in), optional :: is_centers_to_allocate !< Flag to allocate also centers array.
   logical,                 intent(in), optional :: is_extents_to_allocate !< Flag to allocate also extents arrays.

   call self%destroy
   read(file_unit, end=10, err=10) self%Ni, self%Nj, self%Nk, self%gc
   10 call self%alloc(is_centers_to_allocate=is_centers_to_allocate, is_extents_to_allocate=is_extents_to_allocate)
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
