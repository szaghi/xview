!< xview test: Xall restart and running-average files are read back field by field.
program xview_test_xall_rst
!< xview test: Xall restart and running-average files are read back field by field.
!<
!< On the blocks of the xnavis-level-set grid, synthetic files are written in the layout of Xall's memsol.F90 and read
!< with the Xall reader at cell centers. Each record of variable v of block b holds a distinct field, so a record read
!< in the wrong place fails the check. Four cases:
!<
!< 1. restart, one-equation model: p, u, v, w, nu~, viscosity, then the previous time level after all blocks;
!< 2. restart, two-equation model with level set: p, u, v, w, k, eps, level set, viscosity, previous time level;
!< 3. running averages, one-equation model: means of p, u, v, w, nu~, then <uu>, <vv>, <ww>, <uv>, <uw>, <vw>; the
!<    viscosity is rebuilt from the mean nu~ (Spalart-Allmaras) and k_res = (1 + 2 + 3)/2 by construction;
!< 4. running averages, two-equation model: the k, eps means precede the second moments, the viscosity is dropped.
!<
!< Run from the repository root; the scratch file is written in the working directory and deleted.
use penf
use xview_file_grd_object
use xview_file_icc_object
use xview_file_rst_object

implicit none
character(*), parameter :: ipath='src/tests/xnavis-level-set/input/' !< Path to input files.
character(*), parameter :: tmp_file='xview_test_xall_rst.tmp'        !< Scratch restart file.
integer(I4P), parameter :: gc=2                                      !< Ghost cells assumed for Xall files.
integer(I4P), parameter :: iv_p=1, iv_u=2, iv_v=3, iv_w=4            !< Variable ids of the mean flow.
integer(I4P), parameter :: iv_t1=5, iv_t2=6                          !< Variable ids of the turbulence variables.
integer(I4P), parameter :: iv_ls=7, iv_visc=8                        !< Variable ids of level set and viscosity.
integer(I4P), parameter :: iv_cross=9                                !< Variable id of the unused cross moments.
real(R8P),    parameter :: RE=1000._R8P                              !< Reynolds number.
real(R8P),    parameter :: tramin=10._R8P                            !< Start time of the running averages.
real(R8P),    parameter :: time=12.5_R8P                             !< Time of the solution.
real(R8P),    parameter :: dt=0.01_R8P                               !< Time step.
type(file_grd_object)   :: file_grd                                  !< File grd.
type(file_icc_object)   :: file_icc                                  !< File icc.
type(file_rst_object)   :: file_rst                                  !< File rst.
integer(I4P)            :: file_unit                                 !< Logical unit of the scratch file.

call file_grd%load_file(filename=ipath//'CC/cc.02.grd.p002')
call file_icc%load_file(filename=ipath//'CC/cc.02.p002', file_grd=file_grd, correct_metrics=.false., &
                        is_cell_centered=.true.)
if (file_grd%blocks_number < 2) error stop 'error: the test needs at least two blocks'

call write_file(is_ra=.false., nturb=1, is_level_set=.false.)
call file_rst%load_file(file_grd=file_grd, file_icc=file_icc, filename=tmp_file, is_level_set=.false., &
                        is_zeroeq=.false., is_oneeq=.true., is_twoeq=.false., is_xall=.true., is_cell_centered=.true.)
call check(case_name='Xall restart, one-equation', is_ra=.false., nturb=1, is_level_set=.false.)

call write_file(is_ra=.false., nturb=2, is_level_set=.true.)
call file_rst%load_file(file_grd=file_grd, file_icc=file_icc, filename=tmp_file, is_level_set=.true., &
                        is_zeroeq=.false., is_oneeq=.false., is_twoeq=.true., is_xall=.true., is_cell_centered=.true.)
call check(case_name='Xall restart, two-equation, level set', is_ra=.false., nturb=2, is_level_set=.true.)

call write_file(is_ra=.true., nturb=1, is_level_set=.false.)
call file_rst%load_file(file_grd=file_grd, file_icc=file_icc, filename=tmp_file, RE=RE, is_level_set=.false., &
                        is_zeroeq=.false., is_oneeq=.true., is_twoeq=.false., is_xall=.true., is_xall_ra=.true., &
                        is_cell_centered=.true.)
call check(case_name='Xall running averages, one-equation', is_ra=.true., nturb=1, is_level_set=.false.)

call write_file(is_ra=.true., nturb=2, is_level_set=.false.)
call file_rst%load_file(file_grd=file_grd, file_icc=file_icc, filename=tmp_file, RE=RE, is_level_set=.false., &
                        is_zeroeq=.false., is_oneeq=.false., is_twoeq=.true., is_xall=.true., is_xall_ra=.true., &
                        is_cell_centered=.true.)
call check(case_name='Xall running averages, two-equation', is_ra=.true., nturb=2, is_level_set=.false.)

open(newunit=file_unit, file=tmp_file, form='unformatted')
close(file_unit, status='delete')

contains
   pure function field(b, v, Ni, Nj, Nk) result(f)
   !< Field of variable v on block b, cells 1-gc..N+gc: the integer part is v, the hundredths are b, the cells vary below.
   integer(I4P), intent(in) :: b, v, Ni, Nj, Nk !< Block, variable id, block dimensions.
   real(R8P), allocatable   :: f(:,:,:)         !< Field.
   integer(I4P)             :: i, j, k          !< Counters.

   allocate(f(1-gc:Ni+gc, 1-gc:Nj+gc, 1-gc:Nk+gc))
   do k=1-gc, Nk+gc
      do j=1-gc, Nj+gc
         do i=1-gc, Ni+gc
            f(i,j,k) = v + 0.01_R8P * b + 1.e-4_R8P * mod(i + 3 + 7 * (j + 3) + 13 * (k + 3), 10)
         enddo
      enddo
   enddo
   endfunction field

   subroutine write_file(is_ra, nturb, is_level_set)
   !< Write the scratch file in the layout of Xall's memsol.F90.
   logical,      intent(in) :: is_ra        !< Running averages, otherwise restart.
   integer(I4P), intent(in) :: nturb        !< Number of turbulence variables (1 or 2).
   logical,      intent(in) :: is_level_set !< Level set present (restart only).
   integer(I4P)             :: b, v         !< Counters.

   open(newunit=file_unit, file=tmp_file, form='unformatted', action='write', status='replace')
   if (is_ra) then
      write(file_unit) tramin, time, dt, 42_I4P
   else
      write(file_unit) time, dt
   endif
   do b=1, file_grd%blocks_number
      write(file_unit) file_grd%blocks(b)%Ni, file_grd%blocks(b)%Nj, file_grd%blocks(b)%Nk
   enddo
   do b=1, file_grd%blocks_number
      associate(Ni=>file_grd%blocks(b)%Ni, Nj=>file_grd%blocks(b)%Nj, Nk=>file_grd%blocks(b)%Nk)
      do v=iv_p, iv_w + nturb
         write(file_unit) field(b, v, Ni, Nj, Nk)
      enddo
      if (is_ra) then
         write(file_unit) field(b, iv_u, Ni, Nj, Nk)**2 + 1._R8P
         write(file_unit) field(b, iv_v, Ni, Nj, Nk)**2 + 2._R8P
         write(file_unit) field(b, iv_w, Ni, Nj, Nk)**2 + 3._R8P
         do v=1, 3
            write(file_unit) field(b, iv_cross, Ni, Nj, Nk)
         enddo
      else
         if (is_level_set) write(file_unit) field(b, iv_ls, Ni, Nj, Nk)
         write(file_unit) field(b, iv_visc, Ni, Nj, Nk)
      endif
      endassociate
   enddo
   if (.not.is_ra) then ! previous time level, after all blocks: it must not be read
      do b=1, file_grd%blocks_number
         associate(Ni=>file_grd%blocks(b)%Ni, Nj=>file_grd%blocks(b)%Nj, Nk=>file_grd%blocks(b)%Nk)
         do v=iv_p, iv_w + nturb
            write(file_unit) -field(b, v, Ni, Nj, Nk)
         enddo
         endassociate
      enddo
   endif
   close(file_unit)
   endsubroutine write_file

   subroutine check(case_name, is_ra, nturb, is_level_set)
   !< Check every field read back against the one written.
   character(*), intent(in) :: case_name    !< Case name.
   logical,      intent(in) :: is_ra        !< Running averages, otherwise restart.
   integer(I4P), intent(in) :: nturb        !< Number of turbulence variables (1 or 2).
   logical,      intent(in) :: is_level_set !< Level set present.
   real(R8P), allocatable   :: chi(:,:,:)   !< nu~/nu.
   real(R8P)                :: err_max      !< Maximum error of the derived fields.
   integer(I4P)             :: b            !< Counter.

   if (.not.file_rst%is_loaded) call fail(case_name, 'file not loaded')
   if (file_rst%time /= time) call fail(case_name, 'time')
   err_max = 0._R8P
   do b=1, file_rst%blocks_number
      associate(blk=>file_rst%blocks(b), Ni=>file_grd%blocks(b)%Ni, Nj=>file_grd%blocks(b)%Nj, Nk=>file_grd%blocks(b)%Nk)
      if (any(blk%gc /= gc)) call fail(case_name, 'ghost cells')
      if (any(blk%pressure   /= field(b, iv_p, Ni, Nj, Nk))) call fail(case_name, 'p')
      if (any(blk%momentum%x /= field(b, iv_u, Ni, Nj, Nk))) call fail(case_name, 'u')
      if (any(blk%momentum%y /= field(b, iv_v, Ni, Nj, Nk))) call fail(case_name, 'v')
      if (any(blk%momentum%z /= field(b, iv_w, Ni, Nj, Nk))) call fail(case_name, 'w')
      if (nturb == 1) then
         if (any(blk%turbulent_viscosity /= field(b, iv_t1, Ni, Nj, Nk))) call fail(case_name, 'nu~')
      else
         if (any(blk%turbulent_kinetic_energy /= field(b, iv_t1, Ni, Nj, Nk))) call fail(case_name, 'k')
         if (any(blk%turbulent_kinetic_energy_dissipation /= field(b, iv_t2, Ni, Nj, Nk))) call fail(case_name, 'eps')
      endif
      if (is_ra) then
         if (.not.allocated(blk%k_res)) call fail(case_name, 'k_res not allocated')
         err_max = max(err_max, maxval(abs(blk%k_res - 3._R8P)))
         if (nturb == 1) then
            if (.not.allocated(blk%viscosity)) call fail(case_name, 'viscosity not allocated')
            chi = field(b, iv_t1, Ni, Nj, Nk) * RE
            err_max = max(err_max, maxval(abs(blk%viscosity * RE - (1._R8P + chi * chi**3 / (chi**3 + 7.1_R8P**3)))))
         else
            if (allocated(blk%viscosity)) call fail(case_name, 'viscosity not dropped')
         endif
      else
         if (allocated(blk%k_res)) call fail(case_name, 'k_res allocated')
         if (is_level_set) then
            if (any(blk%level_set      /= field(b, iv_ls, Ni, Nj, Nk))) call fail(case_name, 'level set')
            if (any(blk%level_set_zero /= field(b, iv_ls, Ni, Nj, Nk))) call fail(case_name, 'zero level set')
         endif
         if (any(blk%viscosity /= field(b, iv_visc, Ni, Nj, Nk))) call fail(case_name, 'viscosity')
      endif
      endassociate
   enddo
   if (err_max > 1.e-9_R8P) call fail(case_name, 'k_res or rebuilt viscosity, error '//trim(str(err_max)))
   print '(A)', case_name//': '//trim(str(file_rst%blocks_number, no_sign=.true.))//' blocks read back'
   endsubroutine check

   subroutine fail(case_name, what)
   !< Stop with the case and the field that failed.
   character(*), intent(in) :: case_name !< Case name.
   character(*), intent(in) :: what      !< Failed field.

   error stop 'error: '//case_name//': '//what
   endsubroutine fail
endprogram xview_test_xall_rst
