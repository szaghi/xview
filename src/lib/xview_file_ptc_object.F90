!< xview, file patches class definition.
module xview_file_ptc_object
!< xview, file grd class definition.
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf
use stringifor

use xview_block_grd_object
use xview_block_rst_object
use xview_file_object

implicit none
private
public :: patches_extents_object
public :: file_ptc_object

type :: patches_extents_object
   !< Patches extents class.
   integer(I4P)              :: b=0              !< Block id (global numeration).
   integer(I4P)              :: patches_number=0 !< Number of patches, np.
   integer(I4P), allocatable :: extents(:,:)     !< Patches extents [1:np,0:6].
   character(6), allocatable :: dir(:)           !< Patches X/R/THETA direction respect IJK [1:np].
   integer(I4P), allocatable :: id(:)            !< Patches ID [1:np].
endtype patches_extents_object

type, extends(file_object) :: file_ptc_object
   !< **File patches** class, parse, manipulate and emit patches file.
   integer(I4P)                              :: blocks_number=0 !< Number of blocks, nb.
   type(patches_extents_object), allocatable :: patches(:)      !< Patches extents [1:nb].
   contains
      ! plublic methods
      procedure, pass(self) :: destroy                      !< Destroy dynamic memory.
      procedure, pass(self) :: alloc                        !< Allocate dynamic memory.
      procedure, pass(self) :: load_file                    !< Load file.
      procedure, pass(self) :: save_file                    !< Save file.
      procedure, pass(self) :: save_block_file_aero_patches !< Save block patches in aeroacustic format.
endtype file_ptc_object

contains
   ! plublic methods
   elemental subroutine destroy(self)
   !< Destroy dynamic memory.
   class(file_ptc_object), intent(inout) :: self !< File data.
   integer(I4P)                          :: b    !< Counter.

   call self%file_object%destroy
   do b=1, self%blocks_number
      if (allocated(self%patches(b)%extents)) deallocate(self%patches(b)%extents)
      if (allocated(self%patches(b)%dir    )) deallocate(self%patches(b)%dir    )
   enddo
   self%blocks_number = 0
   if (allocated(self%patches)) deallocate(self%patches)
   endsubroutine destroy

   elemental subroutine alloc(self)
   !< Allocate dynamic memory.
   class(file_ptc_object), intent(inout) :: self !< File data.

   allocate(self%patches(1:self%blocks_number))
   endsubroutine alloc

   subroutine load_file(self, filename, verbose)
   !< Load file.
   class(file_ptc_object), intent(inout)        :: self      !< File data.
   character(*),           intent(in)           :: filename  !< File name.
   logical,                intent(in), optional :: verbose   !< Activate verbose mode.
   logical                                      :: verbose_  !< Activate verbose mode, local variable.
   integer(I4P)                                 :: file_unit !< Logical file unit.
   integer(I4P)                                 :: b, p, np  !< Counter.

   call self%destroy
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   self%filename = trim(adjustl(filename))
   if (self%is_file_present()) then
      open(newunit=file_unit, file=trim(adjustl(filename)), form='formatted', action='read')
      read(file_unit, *) self%blocks_number
      call self%alloc
      do b=1, self%blocks_number
         read(file_unit, *) self%patches(b)%b, np
         self%patches(b)%patches_number = np
         allocate(self%patches(b)%extents(1:np,0:6))
         allocate(self%patches(b)%dir(1:np))
         allocate(self%patches(b)%id(1:np))
         do p=1, np
            read(file_unit, *) self%patches(b)%extents(p,0:6), self%patches(b)%dir(p), self%patches(b)%id(p)
         enddo
      enddo
      self%is_loaded = .true.
   else
      if (verbose_) write(stderr, "(A)") 'error: file "'//trim(adjustl(filename))//'" not found!'
      self%is_loaded = .false.
   endif
   endsubroutine load_file

   subroutine save_file(self, filename)
   !< Save file.
   class(file_ptc_object), intent(in)           :: self      !< File data.
   character(*),           intent(in), optional :: filename  !< File name.
   integer(I4P)                                 :: file_unit !< Logical file unit.
   integer(I4P)                                 :: b         !< Counter.

   ! if (present(filename)) then
   !    open(newunit=file_unit, file=trim(adjustl(filename)), form='unformatted', action='write')
   ! elseif (self%filename%is_allocated()) then
   !    open(newunit=file_unit, file=trim(adjustl(self%filename%chars())), form='unformatted', action='write')
   ! else
   !    error stop 'error: nor "filename" neither "self%filename" have been specified for "file_grd_object%save_file" method!'
   ! endif
   ! write(file_unit) self%blocks_number
   ! do b=1, self%blocks_number
   !    call self%blocks(b)%save_dimensions(file_unit=file_unit)
   ! enddo
   ! do b=1, self%blocks_number
   !    call self%blocks(b)%save_nodes(file_unit=file_unit)
   ! enddo
   ! close(file_unit)
   endsubroutine save_file

   subroutine save_block_file_aero_patches(self, grd, sol, files_blocks)
   !< Save block patches in aeroacustic format.
   class(file_ptc_object), intent(in) :: self             !< File data.
   type(block_grd_object), intent(in) :: grd              !< Grid of block.
   type(block_rst_object), intent(in) :: sol              !< Solution of block.
   type(string),           intent(in) :: files_blocks(1:) !< Files blocks group.
   type(string), allocatable          :: tokens(:)        !< String tokens buffer.
   integer(I4P)                       :: file_g_unit      !< Geometry unit file.
   integer(I4P)                       :: file_vp_unit     !< Velocity-pressure unit file.
   integer(I4P)                       :: error            !< Error status.
   integer(I4P)                       :: nx1,nx2          !< First and last node x indexes.
   integer(I4P)                       :: nr1,nr2          !< First and last node r indexes.
   integer(I4P)                       :: nt1,nt2          !< First and last node t indexes.
   integer(I4P)                       :: p, i, j, k       !< Counter.

   if (.not.grd%is_loaded) error stop 'error: cannot save aeropatches block file because passed grd block has not been loaded'
   if (.not.sol%is_loaded) error stop 'error: cannot save aeropatches block file because passed sol block has not been loaded'
   if (.not.allocated(grd%patches_extents)) return ! queried a patch bc that is not present into this block
   do p=1, size(grd%patches_extents, dim=1) ! loop over pathces or whole block
      associate(nci1=>grd%patches_extents(p, 1 ), nci2=>grd%patches_extents(p, 2 ), &
                ncj1=>grd%patches_extents(p, 3 ), ncj2=>grd%patches_extents(p, 4 ), &
                nck1=>grd%patches_extents(p, 5 ), nck2=>grd%patches_extents(p, 6 ), &
                nni1=>grd%patches_extents(p, 7 ), nni2=>grd%patches_extents(p, 8 ), &
                nnj1=>grd%patches_extents(p, 9 ), nnj2=>grd%patches_extents(p, 10), &
                nnk1=>grd%patches_extents(p, 11), nnk2=>grd%patches_extents(p, 12))
         call files_blocks(3)%split(tokens=tokens,sep='_')
         open(newunit=file_g_unit,  file='patches/'//tokens(2)//'-patch-g-' //trim(strz(grd%patches_id(p), nz_pad=3))//'.dat')
         open(newunit=file_vp_unit, file='patches/'//tokens(2)//'-patch-vp-'//trim(strz(grd%patches_id(p), nz_pad=3))//'.dat')
         write(file_g_unit,'(A)') 'VARIABLES="X" "Y" "Z"'
         write(file_vp_unit,'(A)') 'VARIABLES="U" "V" "W" "P"'
         select case(trim(adjustl(grd%patches_dir(p))))
         case('tx-r')
            ! geometry
            nt1 = nni1 ; nt2 = nni2
            nx1 = nnj1 ; nx2 = nnj2
            nr1 = nnk1 ; nr2 = nnk2
            write(file_g_unit,'(A)') 'ZONE I='//trim(str(n=nt2-nt1+1))//&
                                         ' J='//trim(str(n=nx2-nx1+1))//&
                                         ' K='//trim(str(n=nr2-nr1+1))//' F=POINT'
            do k=nr2, nr1, -1
               do j=nx1, nx2
                  do i=nt1, nt2
                     write(file_g_unit,'(A)') trim(str(n=grd%nodes(i,j,k)%x))//' '//&
                                              trim(str(n=grd%nodes(i,j,k)%y))//' '//&
                                              trim(str(n=grd%nodes(i,j,k)%z))
                  enddo
               enddo
            enddo
            ! pressure/velocity
            nt1 = nci1 ; nt2 = nci2
            nx1 = ncj1 ; nx2 = ncj2
            nr1 = nck1 ; nr2 = nck2
            write(file_vp_unit,'(A)') 'ZONE I='//trim(str(n=nt2-nt1+1))//&
                                          ' J='//trim(str(n=nx2-nx1+1))//&
                                          ' K='//trim(str(n=nr2-nr1+1))//' F=POINT'
            do k=nr2, nr1, -1
               do j=nx1, nx2
                  do i=nt1, nt2
                     write(file_vp_unit,'(A)') trim(str(n=sol%momentum(i,j,k)%x))//' '//&
                                               trim(str(n=sol%momentum(i,j,k)%y))//' '//&
                                               trim(str(n=sol%momentum(i,j,k)%z))//' '//&
                                               trim(str(n=sol%pressure(i,j,k)))
                  enddo
               enddo
            enddo
         case('xt-r')
            ! geometry
            nt1 = nnj1 ; nt2 = nnj2
            nx1 = nni1 ; nx2 = nni2
            nr1 = nnk1 ; nr2 = nnk2
            write(file_g_unit,'(A)') 'ZONE I='//trim(str(n=nt2-nt1+1))//&
                                         ' J='//trim(str(n=nx2-nx1+1))//&
                                         ' K='//trim(str(n=nr2-nr1+1))//' F=POINT'
            do k=nr2, nr1, -1
               do j=nx1, nx2
                  do i=nt1, nt2
                     write(file_g_unit,'(A)') trim(str(n=grd%nodes(j,i,k)%x))//' '//&
                                              trim(str(n=grd%nodes(j,i,k)%y))//' '//&
                                              trim(str(n=grd%nodes(j,i,k)%z))
                  enddo
               enddo
            enddo
            ! pressure/velocity
            nt1 = ncj1 ; nt2 = ncj2
            nx1 = nci1 ; nx2 = nci2
            nr1 = nck1 ; nr2 = nck2
            write(file_vp_unit,'(A)') 'ZONE I='//trim(str(n=nt2-nt1+1))//&
                                          ' J='//trim(str(n=nx2-nx1+1))//&
                                          ' K='//trim(str(n=nr2-nr1+1))//' F=POINT'
            do k=nr2, nr1, -1
               do j=nx1, nx2
                  do i=nt1, nt2
                     write(file_vp_unit,'(A)') trim(str(n=sol%momentum(j,i,k)%x))//' '//&
                                               trim(str(n=sol%momentum(j,i,k)%y))//' '//&
                                               trim(str(n=sol%momentum(j,i,k)%z))//' '//&
                                               trim(str(n=sol%pressure(j,i,k)))
                  enddo
               enddo
            enddo
         case('txr')
            ! geometry
            nt1 = nni1 ; nt2 = nni2
            nx1 = nnj1 ; nx2 = nnj2
            nr1 = nnk1 ; nr2 = nnk2
            write(file_g_unit,'(A)') 'ZONE I='//trim(str(n=nt2-nt1+1))//&
                                         ' J='//trim(str(n=nx2-nx1+1))//&
                                         ' K='//trim(str(n=nr2-nr1+1))//' F=POINT'
            do k=nr1, nr2
               do j=nx1, nx2
                  do i=nt1, nt2
                     write(file_g_unit,'(A)') trim(str(n=grd%nodes(i,j,k)%x))//' '//&
                                              trim(str(n=grd%nodes(i,j,k)%y))//' '//&
                                              trim(str(n=grd%nodes(i,j,k)%z))
                  enddo
               enddo
            enddo
            ! pressure/velocity
            nt1 = nci1 ; nt2 = nci2
            nx1 = ncj1 ; nx2 = ncj2
            nr1 = nck1 ; nr2 = nck2
            write(file_vp_unit,'(A)') 'ZONE I='//trim(str(n=nt2-nt1+1))//&
                                          ' J='//trim(str(n=nx2-nx1+1))//&
                                          ' K='//trim(str(n=nr2-nr1+1))//' F=POINT'
            do k=nr1, nr2
               do j=nx1, nx2
                  do i=nt1, nt2
                     write(file_vp_unit,'(A)') trim(str(n=sol%momentum(i,j,k)%x))//' '//&
                                               trim(str(n=sol%momentum(i,j,k)%y))//' '//&
                                               trim(str(n=sol%momentum(i,j,k)%z))//' '//&
                                               trim(str(n=sol%pressure(i,j,k)))
                  enddo
               enddo
            enddo
         case('aaaa')
            ! geometry
            nt1 = nni1 ; nt2 = nni2
            nx1 = nnj1 ; nx2 = nnj2
            nr1 = nnk1 ; nr2 = nnk2
            write(file_g_unit,'(A)') 'ZONE I='//trim(str(n=nt2-nt1+1))//&
                                         ' J='//trim(str(n=nx2-nx1+1))//&
                                         ' K='//trim(str(n=nr2-nr1+1))//' F=POINT'
            do k=nr1, nr2
               do j=nx1, nx2
                  do i=nt1, nt2
                     write(file_g_unit,'(A)') trim(str(n=grd%nodes(i,j,k)%x))//' '//&
                                              trim(str(n=grd%nodes(i,j,k)%y))//' '//&
                                              trim(str(n=grd%nodes(i,j,k)%z))
                  enddo
               enddo
            enddo
            ! y x z no
            ! pressure/velocity
            nt1 = nci1 ; nt2 = nci2
            nx1 = ncj1 ; nx2 = ncj2
            nr1 = nck1 ; nr2 = nck2
            write(file_vp_unit,'(A)') 'ZONE I='//trim(str(n=nt2-nt1+1))//&
                                          ' J='//trim(str(n=nx2-nx1+1))//&
                                          ' K='//trim(str(n=nr2-nr1+1))//' F=POINT'
            do k=nr1, nr2
               do j=nx1, nx2
                  do i=nt1, nt2
                     write(file_vp_unit,'(A)') trim(str(n=sol%momentum(i,j,k)%x))//' '//&
                                               trim(str(n=sol%momentum(i,j,k)%y))//' '//&
                                               trim(str(n=sol%momentum(i,j,k)%z))//' '//&
                                               trim(str(n=sol%pressure(i,j,k)))
                  enddo
               enddo
            enddo
         case('bbbb')
            ! geometry
            nt1 = nni1 ; nt2 = nni2
            nx1 = nnj1 ; nx2 = nnj2
            nr1 = nnk1 ; nr2 = nnk2
            write(file_g_unit,'(A)') 'ZONE I='//trim(str(n=nt2-nt1+1))//&
                                         ' J='//trim(str(n=nx2-nx1+1))//&
                                         ' K='//trim(str(n=nr2-nr1+1))//' F=POINT'
            do k=nr1, nr2
               do j=nx1, nx2
                  do i=nt1, nt2
                     write(file_g_unit,'(A)') trim(str(n=grd%nodes(i,j,k)%x))//' '//&
                                              trim(str(n=grd%nodes(i,j,k)%y))//' '//&
                                              trim(str(n=grd%nodes(i,j,k)%z))
                  enddo
               enddo
            enddo
            ! y x z no
            ! pressure/velocity
            nt1 = nci1 ; nt2 = nci2
            nx1 = ncj1 ; nx2 = ncj2
            nr1 = nck1 ; nr2 = nck2
            write(file_vp_unit,'(A)') 'ZONE I='//trim(str(n=nt2-nt1+1))//&
                                          ' J='//trim(str(n=nx2-nx1+1))//&
                                          ' K='//trim(str(n=nr2-nr1+1))//' F=POINT'
            do k=nr1, nr2
               do j=nx1, nx2
                  do i=nt1, nt2
                     write(file_vp_unit,'(A)') trim(str(n=sol%momentum(i,j,k)%x))//' '//&
                                               trim(str(n=sol%momentum(i,j,k)%y))//' '//&
                                               trim(str(n=sol%momentum(i,j,k)%z))//' '//&
                                               trim(str(n=sol%pressure(i,j,k)))
                  enddo
               enddo
            enddo
         case('cccc')
            ! geometry
            nt1 = nni1 ; nt2 = nni2
            nx1 = nnj1 ; nx2 = nnj2
            nr1 = nnk1 ; nr2 = nnk2
            write(file_g_unit,'(A)') 'ZONE I='//trim(str(n=nt2-nt1+1))//&
                                         ' J='//trim(str(n=nx2-nx1+1))//&
                                         ' K='//trim(str(n=nr2-nr1+1))//' F=POINT'
            do k=nr1, nr2
               do j=nx1, nx2
                  do i=nt1, nt2
                     write(file_g_unit,'(A)') trim(str(n=grd%nodes(i,j,k)%x))//' '//&
                                              trim(str(n=grd%nodes(i,j,k)%y))//' '//&
                                              trim(str(n=grd%nodes(i,j,k)%z))
                  enddo
               enddo
            enddo
            ! y x z no
            ! pressure/velocity
            nt1 = nci1 ; nt2 = nci2
            nx1 = ncj1 ; nx2 = ncj2
            nr1 = nck1 ; nr2 = nck2
            write(file_vp_unit,'(A)') 'ZONE I='//trim(str(n=nt2-nt1+1))//&
                                          ' J='//trim(str(n=nx2-nx1+1))//&
                                          ' K='//trim(str(n=nr2-nr1+1))//' F=POINT'
            do k=nr1, nr2
               do j=nx1, nx2
                  do i=nt1, nt2
                     write(file_vp_unit,'(A)') trim(str(n=sol%momentum(i,j,k)%x))//' '//&
                                               trim(str(n=sol%momentum(i,j,k)%y))//' '//&
                                               trim(str(n=sol%momentum(i,j,k)%z))//' '//&
                                               trim(str(n=sol%pressure(i,j,k)))
                  enddo
               enddo
            enddo
         case('xtr')
            ! geometry
            nt1 = nnj1 ; nt2 = nnj2
            nx1 = nni1 ; nx2 = nni2
            nr1 = nnk1 ; nr2 = nnk2
            write(file_g_unit,'(A)') 'ZONE I='//trim(str(n=nt2-nt1+1))//&
                                         ' J='//trim(str(n=nx2-nx1+1))//&
                                         ' K='//trim(str(n=nr2-nr1+1))//' F=POINT'
            do k=nr1, nr2
               do j=nx1, nx2
                  do i=nt1, nt2
                     write(file_g_unit,'(A)') trim(str(n=grd%nodes(j,i,k)%x))//' '//&
                                              trim(str(n=grd%nodes(j,i,k)%y))//' '//&
                                              trim(str(n=grd%nodes(j,i,k)%z))
                  enddo
               enddo
            enddo
            ! pressure/velocity
            nt1 = ncj1 ; nt2 = ncj2
            nx1 = nci1 ; nx2 = nci2
            nr1 = nck1 ; nr2 = nck2
            write(file_vp_unit,'(A)') 'ZONE I='//trim(str(n=nt2-nt1+1))//&
                                          ' J='//trim(str(n=nx2-nx1+1))//&
                                          ' K='//trim(str(n=nr2-nr1+1))//' F=POINT'
            do k=nr1, nr2
               do j=nx1, nx2
                  do i=nt1, nt2
                     write(file_vp_unit,'(A)') trim(str(n=sol%momentum(j,i,k)%x))//' '//&
                                               trim(str(n=sol%momentum(j,i,k)%y))//' '//&
                                               trim(str(n=sol%momentum(j,i,k)%z))//' '//&
                                               trim(str(n=sol%pressure(j,i,k)))
                  enddo
               enddo
            enddo
         case('ijk')
            ! geometry
            nt1 = nni1 ; nt2 = nni2
            nx1 = nnj1 ; nx2 = nnj2
            nr1 = nnk1 ; nr2 = nnk2
            write(file_g_unit,'(A)') 'ZONE I='//trim(str(n=nt2-nt1+1))//&
                                         ' J='//trim(str(n=nx2-nx1+1))//&
                                         ' K='//trim(str(n=nr2-nr1+1))//' F=POINT'
            do k=nr1, nr2
               do j=nx1, nx2
                  do i=nt1, nt2
                     write(file_g_unit,'(A)') trim(str(n=grd%nodes(i,j,k)%x))//' '//&
                                              trim(str(n=grd%nodes(i,j,k)%y))//' '//&
                                              trim(str(n=grd%nodes(i,j,k)%z))
                  enddo
               enddo
            enddo
            ! pressure/velocity
            nt1 = nci1 ; nt2 = nci2
            nx1 = ncj1 ; nx2 = ncj2
            nr1 = nck1 ; nr2 = nck2
            write(file_vp_unit,'(A)') 'ZONE I='//trim(str(n=nt2-nt1+1))//&
                                          ' J='//trim(str(n=nx2-nx1+1))//&
                                          ' K='//trim(str(n=nr2-nr1+1))//' F=POINT'
            do k=nr1, nr2
               do j=nx1, nx2
                  do i=nt1, nt2
                     write(file_vp_unit,'(A)') trim(str(n=sol%momentum(i,j,k)%x))//' '//&
                                               trim(str(n=sol%momentum(i,j,k)%y))//' '//&
                                               trim(str(n=sol%momentum(i,j,k)%z))//' '//&
                                               trim(str(n=sol%pressure(i,j,k)))
                  enddo
               enddo
            enddo
         endselect
         close(file_g_unit)
         close(file_vp_unit)
      endassociate
   enddo
   endsubroutine
endmodule xview_file_ptc_object
