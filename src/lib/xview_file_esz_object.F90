!< xview, file extracted subzone class definition.
module xview_file_esz_object
!< xview, file extracted subzone class definition.
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf

use xview_block_esz_object
use xview_file_object
use xview_file_ptc_object

implicit none
private
public :: file_esz_object

type, extends(file_object) :: file_esz_object
   !< File rst class definition.
   type(block_esz_object) :: block_esz !< Extracted block subzone.
   contains
      ! plublic methods
      procedure, pass(self) :: destroy   !< Destroy dynamic memory.
      procedure, pass(self) :: load_file !< Load file.
endtype file_esz_object

contains
   ! plublic methods
   elemental subroutine destroy(self)
   !< Destroy dynamic memory.
   class(file_esz_object), intent(inout) :: self !< File data.

   call self%file_object%destroy
   call self%block_esz%destroy
   endsubroutine destroy

   subroutine load_file(self,filename,file_ptc,is_cell_centered,patch,RE,rFR2,zfs,is_level_set,is_zeroeq,is_oneeq,is_twoeq, &
                        compute_metrics,compute_lambda2,compute_qfactor,compute_helicity,compute_vorticity,                 &
                        compute_div2LT,compute_grad_p,compute_k_ratio,compute_yplus,compute_tau,compute_div_tau,            &
                        compute_loads,verbose)
   !< Load file.
   class(file_esz_object), intent(inout)        :: self              !< File data.
   character(*),           intent(in)           :: filename          !< File name of geo file.
   type(file_ptc_object),  intent(in)           :: file_ptc          !< File patches.
   logical,                intent(in), optional :: is_cell_centered  !< Define variables at cell centers or nodes.
   integer(I4P),           intent(in), optional :: patch             !< Patch boundary conditions.
   real(R8P),              intent(in), optional :: RE                !< Reynolds number.
   real(R8P),              intent(in), optional :: rFR2              !< 1/(Froude number)^2.
   real(R8P),              intent(in), optional :: zfs               !< Z quote of free surface.
   logical,                intent(in), optional :: is_level_set      !< Flag for level set function presence.
   logical,                intent(in), optional :: is_zeroeq         !< Use *zero* equations turbulence model.
   logical,                intent(in), optional :: is_oneeq          !< Use *one* equations turbulence model.
   logical,                intent(in), optional :: is_twoeq          !< Use *two* equations turbulence model.
   logical,                intent(in), optional :: compute_metrics   !< Compute metrics.
   logical,                intent(in), optional :: compute_lambda2   !< Compute lamda2 field.
   logical,                intent(in), optional :: compute_qfactor   !< Compute qfactor field.
   logical,                intent(in), optional :: compute_helicity  !< Compute helicity field.
   logical,                intent(in), optional :: compute_vorticity !< Compute vorticity field.
   logical,                intent(in), optional :: compute_grad_p    !< Compute gradient pressure.
   logical,                intent(in), optional :: compute_div2LT    !< Compute double divergence of Lighthill tensor.
   logical,                intent(in), optional :: compute_k_ratio   !< Compute kinetic energy ratio.
   logical,                intent(in), optional :: compute_yplus     !< Compute y+ field.
   logical,                intent(in), optional :: compute_tau       !< Compute tau field.
   logical,                intent(in), optional :: compute_div_tau   !< Compute divergence of tau field.
   logical,                intent(in), optional :: compute_loads     !< Compute loads (forces and torques).
   logical,                intent(in), optional :: verbose           !< Activate verbose mode.
   logical                                      :: verbose_          !< Activate verbose mode, local variable.
   integer(I4P)                                 :: file_unit         !< Logical file unit.

   call self%destroy
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   self%filename = trim(adjustl(filename))
   if (self%is_file_present()) then
      open(newunit=file_unit, file=trim(adjustl(filename)), form='unformatted', action='read')
      call self%block_esz%load_solution(file_unit=file_unit,                 &
                                        file_ptc=file_ptc,                   &
                                        patch=patch,                         &
                                        is_cell_centered=is_cell_centered,   &
                                        RE=RE, rFR2=rFR2, zfs=zfs,           &
                                        is_level_set=is_level_set,           &
                                        is_zeroeq=is_zeroeq,                 &
                                        is_oneeq=is_oneeq,                   &
                                        is_twoeq=is_twoeq,                   &
                                        compute_lambda2=compute_lambda2,     &
                                        compute_qfactor=compute_qfactor,     &
                                        compute_helicity=compute_helicity,   &
                                        compute_vorticity=compute_vorticity, &
                                        compute_div2LT=compute_div2LT,       &
                                        compute_grad_p=compute_grad_p,       &
                                        compute_k_ratio=compute_k_ratio,     &
                                        compute_yplus=compute_yplus,         &
                                        compute_tau=compute_tau,             &
                                        compute_div_tau=compute_div_tau,     &
                                        compute_loads=compute_loads)
      close(file_unit)
      self%is_loaded = .true.
   else
      if (verbose_) write(stderr, "(A)") 'warning: file "'//trim(adjustl(filename))//'" not found!'
      self%is_loaded = .false.
   endif
   endsubroutine load_file
endmodule xview_file_esz_object
