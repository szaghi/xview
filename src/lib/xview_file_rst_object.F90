!< xview, file rst class definition.
module xview_file_rst_object
!< xview, file rst class definition.
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf

use xview_block_rst_object
use xview_file_object
use xview_file_grd_object
use xview_file_icc_object

implicit none
private
public :: file_rst_object

type, extends(file_object) :: file_rst_object
   !< File rst class definition.
   integer(I4P)                        :: blocks_number=0 !< Number of blocks contained into the file.
   type(block_rst_object), allocatable :: blocks(:)       !< Blocks contained into the sol file.
   real(R8P)                           :: time=0._R8P     !< Time of solution.
   contains
      ! plublic methods
      procedure, pass(self) :: destroy   !< Destroy dynamic memory.
      procedure, pass(self) :: alloc     !< Allocate dynamic memory.
      procedure, pass(self) :: load_file !< Load file.
      procedure, pass(self) :: save_file !< Save file.
endtype file_rst_object

contains
   ! plublic methods
   elemental subroutine destroy(self)
   !< Destroy dynamic memory.
   class(file_rst_object), intent(inout) :: self !< File data.
   integer(I4P)                          :: b    !< Counter.

   call self%file_object%destroy
   self%blocks_number = 0
   if (allocated(self%blocks)) then
      do b=lbound(self%blocks, dim=1), ubound(self%blocks, dim=1)
         call self%blocks(b)%destroy
      enddo
      deallocate(self%blocks)
   endif
   self%time = 0._R8P
   endsubroutine destroy

   elemental subroutine alloc(self)
   !< Allocate dynamic memory.
   class(file_rst_object), intent(inout) :: self !< File data.

   allocate(self%blocks(1:self%blocks_number))
   endsubroutine alloc

   subroutine load_file(self,file_grd,file_icc,filename,patch,RE,rFR2,zfs,                                                 &
                        is_level_set,is_zeroeq,is_oneeq,is_twoeq,is_cell_centered,                                         &
                        compute_lambda2,compute_qfactor,compute_helicity,compute_vorticity,compute_div2LT,compute_k_ratio, &
                        compute_yplus,compute_tau,compute_div_tau,compute_loads,verbose)
   !< Load file.
   class(file_rst_object), intent(inout)        :: self              !< File data.
   type(file_grd_object),  intent(in)           :: file_grd          !< File grd data.
   type(file_icc_object),  intent(in)           :: file_icc          !< File icc data.
   character(*),           intent(in)           :: filename          !< File name of geo file.
   integer(I4P),           intent(in), optional :: patch             !< Patch boundary conditions.
   real(R8P),              intent(in), optional :: RE                !< Reynolds number.
   real(R8P),              intent(in), optional :: rFR2              !< 1/(Froude number)^2.
   real(R8P),              intent(in), optional :: zfs               !< Z quote of free surface.
   logical,                intent(in), optional :: is_level_set      !< Flag for level set function presence.
   logical,                intent(in), optional :: is_zeroeq         !< Use *zero* equations turbulence model.
   logical,                intent(in), optional :: is_oneeq          !< Use *one* equations turbulence model.
   logical,                intent(in), optional :: is_twoeq          !< Use *two* equations turbulence model.
   logical,                intent(in), optional :: is_cell_centered  !< Define variables at cell centers or nodes.
   logical,                intent(in), optional :: compute_lambda2   !< Compute lamda2 field.
   logical,                intent(in), optional :: compute_qfactor   !< Compute qfactor field.
   logical,                intent(in), optional :: compute_helicity  !< Compute helicity field.
   logical,                intent(in), optional :: compute_vorticity !< Compute vorticity field.
   logical,                intent(in), optional :: compute_div2LT    !< Compute double divergence of Lighthill tensor.
   logical,                intent(in), optional :: compute_k_ratio   !< Compute kinetic energy ratio.
   logical,                intent(in), optional :: compute_yplus     !< Compute y+ field.
   logical,                intent(in), optional :: compute_tau       !< Compute tau field.
   logical,                intent(in), optional :: compute_div_tau   !< Compute divergence of tau field.
   logical,                intent(in), optional :: compute_loads     !< Compute loads (forces and torques).
   logical,                intent(in), optional :: verbose           !< Activate verbose mode.
   logical                                      :: verbose_          !< Activate verbose mode, local variable.
   integer(I4P)                                 :: file_unit         !< Logical file unit of geo file.
   integer(I4P)                                 :: b                 !< Counter.

   call self%destroy
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   associate(blocks_number=>file_grd%blocks_number)
   self%filename = trim(adjustl(filename))
   if (self%is_file_present()) then
      self%blocks_number = blocks_number
      call self%alloc
      open(newunit=file_unit, file=trim(adjustl(filename)), form='unformatted', action='read')
      read(file_unit) self%time
      do b=1, self%blocks_number
         call self%blocks(b)%load_dimensions(file_unit=file_unit,                 &
                                             is_level_set=is_level_set,           &
                                             is_zeroeq=is_zeroeq,                 &
                                             is_oneeq=is_oneeq,                   &
                                             is_twoeq=is_twoeq,                   &
                                             compute_lambda2=compute_lambda2,     &
                                             compute_qfactor=compute_qfactor,     &
                                             compute_helicity=compute_helicity,   &
                                             compute_vorticity=compute_vorticity, &
                                             compute_div2LT=compute_div2LT,       &
                                             compute_k_ratio=compute_k_ratio,     &
                                             compute_yplus=compute_yplus,         &
                                             compute_tau=compute_tau,             &
                                             compute_div_tau=compute_div_tau,     &
                                             compute_loads=compute_loads)
      enddo
      do b=1, self%blocks_number
         call self%blocks(b)%load_solution(file_unit=file_unit,               &
                                           grd=file_grd%blocks(b),            &
                                           icc=file_icc%blocks(b),            &
                                           is_cell_centered=is_cell_centered, &
                                           patch=patch,                       &
                                           RE=RE, rFR2=rFR2, zfs=zfs)
      enddo
      close(file_unit)
      self%is_loaded = .true.
   else
      if (verbose_) write(stderr, "(A)") 'warning: file "'//trim(adjustl(filename))//'" not found!'
      self%is_loaded = .false.
   endif
   endassociate
   endsubroutine load_file

   subroutine save_file(self, filename)
   !< Save file.
   class(file_rst_object), intent(in)           :: self      !< File data.
   character(*),           intent(in), optional :: filename  !< File name.
   integer(I4P)                                 :: file_unit !< Logical file unit.
   integer(I4P)                                 :: b         !< Counter.

   if (present(filename)) then
      open(newunit=file_unit, file=trim(adjustl(filename)), form='unformatted', action='write')
   elseif (self%filename%is_allocated()) then
      open(newunit=file_unit, file=trim(adjustl(self%filename%chars())), form='unformatted', action='write')
   else
      error stop 'error: nor "filename" neither "self%filename" have been specified for "file_rst_object%save_file" method!'
   endif
   write(file_unit) self%time
   do b=1, self%blocks_number
      call self%blocks(b)%save_dimensions(file_unit=file_unit)
   enddo
   do b=1, self%blocks_number
      call self%blocks(b)%save_solution(file_unit=file_unit)
   enddo
   close(file_unit)
   endsubroutine save_file
endmodule xview_file_rst_object
