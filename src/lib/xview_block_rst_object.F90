!< xview, block (rst) class definition.
module xview_block_rst_object
!< xview, block (rst) class definition.

use penf
use vecfor

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
   logical                   :: is_zero_eq=.false.                          !< Use *zero* equations turbulence model.
   logical                   :: is_one_eq=.false.                           !< Use *one* equations turbulence model.
   logical                   :: is_two_eq=.false.                           !< Use *two* equations turbulence model.
   type(vector), allocatable :: momentum(:,:,:)                             !< Momentum field.
   real(R8P),    allocatable :: pressure(:,:,:)                             !< Pressure field.
   real(R8P),    allocatable :: level_set(:,:,:)                            !< Level set field.
   real(R8P),    allocatable :: level_set_zero(:,:,:)                       !< Zero value of level set field.
   real(R8P),    allocatable :: viscosity(:,:,:)                            !< Viscosity field.
   real(R8P),    allocatable :: turbulent_viscosity(:,:,:)                  !< Turbulent viscosity field.
   real(R8P),    allocatable :: turbulent_kinetic_energy(:,:,:)             !< Turbulent kinetic energy field.
   real(R8P),    allocatable :: turbulent_kinetic_energy_dissipation(:,:,:) !< Turbulent kinetic energy dissipation field.
   logical                   :: is_loaded=.false.                           !< Flag for checking if the block is loaded.
   contains
      ! public methods
      procedure, pass(self) :: destroy                      !< Destroy dynamic memory.
      procedure, pass(self) :: alloc                        !< Allocate dynamic memory.
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
   self%is_zero_eq = .false.
   self%is_one_eq = .false.
   self%is_two_eq = .false.
   if (allocated(self%momentum))                             deallocate(self%momentum)
   if (allocated(self%pressure))                             deallocate(self%pressure)
   if (allocated(self%level_set))                            deallocate(self%level_set)
   if (allocated(self%level_set_zero))                       deallocate(self%level_set_zero)
   if (allocated(self%viscosity))                            deallocate(self%viscosity)
   if (allocated(self%turbulent_viscosity))                  deallocate(self%turbulent_viscosity)
   if (allocated(self%turbulent_kinetic_energy))             deallocate(self%turbulent_kinetic_energy)
   if (allocated(self%turbulent_kinetic_energy_dissipation)) deallocate(self%turbulent_kinetic_energy_dissipation)
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
      if (self%is_zero_eq.or.self%is_one_eq.or.self%is_two_eq) then
         allocate(self%viscosity(1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6)))
      endif
      if (self%is_one_eq) allocate(self%turbulent_viscosity(1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6)))
      if (self%is_two_eq) then
         allocate(self%turbulent_kinetic_energy(1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6)))
         allocate(self%turbulent_kinetic_energy_dissipation(1-gc(1):Ni+gc(2), 1-gc(3):Nj+gc(4), 1-gc(5):Nk+gc(6)))
      endif
   endassociate
   endsubroutine alloc

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

   subroutine init(self, Ni, Nj, Nk, gc, is_level_set, turbulence_model)
   !< Initialize block.
   class(block_rst_object), intent(inout)        :: self             !< Block data.
   integer(I4P),            intent(in)           :: Ni               !< Number of cells in i direction.
   integer(I4P),            intent(in)           :: Nj               !< Number of cells in j direction.
   integer(I4P),            intent(in)           :: Nk               !< Number of cells in k direction.
   integer(I4P),            intent(in), optional :: gc(6)            !< Number of ghost cells.
   logical,                 intent(in), optional :: is_level_set     !< Flag for level set function presence.
   character(*),            intent(in), optional :: turbulence_model !< Turbulence model: 'zero_eq', 'one_eq', 'two_eq'.

   call self%destroy
   if (present(is_level_set)) self%is_level_set = is_level_set
   if (present(turbulence_model)) then
      select case(trim(adjustl(turbulence_model)))
      case('zero_eq', 'ZERO_EQ')
         self%is_zero_eq = .true.
      case('one_eq', 'ONE_EQ')
         self%is_one_eq = .true.
      case('two_eq', 'TWO_EQ')
         self%is_two_eq = .true.
      case default
         ! do nothing: no turbulence model is assumed, no other turbulence members are allocated
      endselect
   endif
   self%Ni = Ni
   self%Nj = Nj
   self%Nk = Nk
   if (present(gc)) self%gc = gc
   10 call self%alloc
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

   subroutine load_dimensions(self, file_unit, is_level_set, turbulence_model)
   !< Load block dimensions from file.
   !<
   !< @note The solution file must be already open and the current record-index must be at the proper block dimensions record.
   class(block_rst_object), intent(inout)        :: self             !< Block data.
   integer(I4P),            intent(in)           :: file_unit        !< Logical file unit.
   logical,                 intent(in), optional :: is_level_set     !< Flag for level set function presence.
   character(*),            intent(in), optional :: turbulence_model !< Turbulence model: 'zero_eq', 'one_eq', 'two_eq'.
   integer(I4P)                                  :: Ni               !< Number of cells in i direction.
   integer(I4P)                                  :: Nj               !< Number of cells in j direction.
   integer(I4P)                                  :: Nk               !< Number of cells in k direction.
   integer(I4P)                                  :: gc(6)            !< Number of ghost cells.

   gc = -1
   read(file_unit, end=10, err=10) Ni, Nj, Nk, gc
   10 continue
   if (all(gc >= 0)) then
      call self%init(Ni=Ni, Nj=Nj, Nk=Nk, gc=gc, is_level_set=is_level_set, turbulence_model=turbulence_model)
   else
      call self%init(Ni=Ni, Nj=Nj, Nk=Nk, is_level_set=is_level_set, turbulence_model=turbulence_model)
   endif
   endsubroutine load_dimensions

   subroutine load_solution(self, file_unit, is_cell_centered)
   !< Load block solution from file.
   !<
   !< @note The solution file must be already open and the current record-index must be at the proper block nodes record.
   class(block_rst_object), intent(inout)        :: self             !< Block data.
   integer(I4P),            intent(in)           :: file_unit        !< Logical file unit.
   logical,                 intent(in), optional :: is_cell_centered !< Define variables at cell centers or nodes.
   integer(I4P)                                  :: i                !< Counter.
   integer(I4P)                                  :: j                !< Counter.
   integer(I4P)                                  :: k                !< Counter.

   associate(Ni => self%Ni, Nj => self%Nj, Nk => self%Nk, gc => self%gc)
      read(file_unit)(((self%momentum(i, j, k)%x, i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      read(file_unit)(((self%momentum(i, j, k)%y, i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      read(file_unit)(((self%momentum(i, j, k)%z, i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      read(file_unit)(((self%pressure(i, j, k)  , i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      if (self%is_level_set) then
         read(file_unit)(((self%level_set(i, j, k), i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
         read(file_unit)(((self%level_set_zero(i, j, k), i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      endif
      if (self%is_zero_eq.or.self%is_one_eq.or.self%is_two_eq) then
         read(file_unit)(((self%viscosity(i, j, k), i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      endif
      if (self%is_one_eq) read(file_unit)(((self%turbulent_viscosity(i, j, k), i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      if (self%is_two_eq) then
         read(file_unit)(((self%turbulent_kinetic_energy(i, j, k), i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
         read(file_unit)(((self%turbulent_kinetic_energy_dissipation(i, j, k), i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      endif
   endassociate
   self%is_loaded = .true.
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
      if (self%is_zero_eq.or.self%is_one_eq.or.self%is_two_eq) then
         write(file_unit)(((self%viscosity(i, j, k), i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      endif
      if (self%is_one_eq) write(file_unit)(((self%turbulent_viscosity(i, j, k), i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      if (self%is_two_eq) then
         write(file_unit)(((self%turbulent_kinetic_energy(i, j, k), i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
         write(file_unit)(((self%turbulent_kinetic_energy_dissipation(i, j, k), i=0, Ni + 1), j=0, Nj + 1), k=0, Nk + 1)
      endif
   endassociate
   endsubroutine save_solution
endmodule xview_block_rst_object
