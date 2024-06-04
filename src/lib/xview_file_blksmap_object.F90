!< xview, blocks map ini file class definition.
module xview_file_blksmap_object
!< xview, blocks map ini file class definition.
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf
use stringifor

use xview_file_object

implicit none
private
public :: file_blksmap_object

type, extends(file_object) :: file_blksmap_object
   !< Blocks map ini file class definition.
   integer(I4P)              :: groups_number=0_I4P !< Number of blocks groups.
   type(string), allocatable :: blocks_map(:,:)     !< Blocks groups map.
   contains
      ! public methods
      procedure, pass(self) :: description      !< Return pretty-printed object description.
      procedure, pass(self) :: destroy          !< Destroy dynamic memory.
      procedure, pass(self) :: get_files_blocks !< Return list of blocks per files accordingly blocks groups map.
      procedure, pass(self) :: load_file        !< Load mb.par file data.
endtype file_blksmap_object

contains
   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(file_blksmap_object), intent(in) :: self             !< File data.
   character(len=:), allocatable          :: desc             !< Description.
   character(len=1), parameter            :: NL=new_line('a') !< New line character.
   integer(I4P)                           :: g                !< Counter.

   if (self%is_loaded) then
      desc =       self%filename//' parsed data:'//NL
      desc = desc//'  Groups number: '//trim(str(self%groups_number))
      do g=1, self%groups_number
         desc = desc//NL//'    group  :     '//self%blocks_map(1,g)
         desc = desc//NL//'    blocks :     '//self%blocks_map(2,g)
      enddo
   else
      desc = 'warning: file "'//self%filename//'" not loaded!'
   endif
   endfunction description

   elemental subroutine destroy(self)
   !< Destroy dynamic memory.
   class(file_blksmap_object), intent(inout) :: self !< File data.

   call self%file_object%destroy
   self%groups_number = 0_I4P
   deallocate(self%blocks_map)
   endsubroutine destroy

   subroutine get_files_blocks(self, files_blocks)
   !< Return list of blocks per files accordingly blocks groups map.
   class(file_blksmap_object), intent(inout) :: self                !< File data.
   type(string), allocatable,  intent(inout) :: files_blocks(:,:)   !< Files blocks grouping.
   type(string), allocatable                 :: blocks(:)           !< Blocks list of each group.
   type(string)                              :: buf                 !< Buffer string.
   type(string), allocatable                 :: tok1(:), tok2(:)    !< Tokens string.
   integer(I4P)                              :: b, bb, bbb, bbbb, g !< Counter.

   if (self%is_loaded) then
      ! linearize blocks list given grouped per files
      if (allocated(files_blocks)) then
         buf = ''
         do g=1, size(files_blocks, dim=2)
            buf = buf//' '//files_blocks(1,g)
         enddo
         call buf%split(tokens=blocks, sep=' ')
         deallocate(files_blocks)
         allocate(files_blocks(3,self%groups_number))
         do g=1, self%groups_number
            files_blocks(1,g) = ''
            files_blocks(2,g) = ''
            files_blocks(3,g) = self%blocks_map(1,g)
         enddo
      endif
      ! re-group blocks accordingly map ini
      do b=1, size(blocks, dim=1)
         ! find current block number (global)
         call blocks(b)%split(tokens=tok1, sep='-blk_')
         call tok1(2)%split(tokens=tok2, sep='-grp_')
         bb = tok2(1)%to_number(1_I4P)
         ! search into blocks map for current block
         do g=1, self%groups_number
            call self%blocks_map(2,g)%split(tokens=tok1, sep=' ')
            do bbb=1, size(tok1, dim=1)
               bbbb = tok1(bbb)%to_number(1_I4P)
               if (bb==bbbb) then
                  files_blocks(1,g) = files_blocks(1,g)//' '//blocks(b)
                  files_blocks(2,g) = files_blocks(2,g)//' '//blocks(b)%basename(strip_last_extension=.true.)
               endif
            enddo
         enddo
      enddo

   endif
   endsubroutine get_files_blocks

   subroutine load_file(self, path, filename, verbose)
   !< Load mb.par file data.
   class(file_blksmap_object), intent(inout)        :: self          !< File data.
   character(*),               intent(in)           :: path          !< Path to mb.par.
   character(*),               intent(in)           :: filename      !< File name of mb.par, optionally customizable.
   logical,                    intent(in), optional :: verbose       !< Activate verbose mode.
   logical                                          :: verbose_      !< Activate verbose mode, local variable.
   type(string)                                     :: ini_str       !< Blocks map ini file parsing string.
   type(string), allocatable                        :: ini_groups(:) !< Ini file groups.
   integer(I4P)                                     :: g             !< Counter.

   verbose_  = .false. ; if (present(verbose))  verbose_  = verbose
   self%filename = trim(adjustl(path))//trim(adjustl(filename))

   self%is_loaded = .false.
   if (self%is_file_present()) then
      call ini_str%read_file(file=self%filename%chars()) ! read ini file as a single stream
      call ini_str%split(tokens=ini_groups, sep='[')     ! get ini file groups
      ! parse ini rows (aka groups)
      self%groups_number = size(ini_groups, dim=1)
      allocate(self%blocks_map(2,self%groups_number))
      do g=1, self%groups_number
         call parse_group(row=ini_groups(g), group=self%blocks_map(:,g))
      enddo
      self%is_loaded = .true.
   else
      if (verbose_) write(stderr, "(A)") 'warning: file "'//self%filename//'" not found!'
   endif
   endsubroutine load_file

   ! non TBP methods
   subroutine parse_group(row, group)
   !< Get token in given row by column (n) and columns separator.
   type(string), intent(inout) :: row        !< Row to be parsed.
   type(string), intent(inout) :: group(3)   !< Group data.
   type(string), allocatable   :: columns(:) !< Columns.

   row = row%replace(old=new_line('a'), new='')
   call row%split(tokens=columns, sep='=')
   group(1) = trim(adjustl(columns(1)%chars()))
   group(1) = group(1)%replace(old=']blocks', new='')
   group(2) = trim(adjustl(columns(2)%chars()))
   endsubroutine parse_group
endmodule xview_file_blksmap_object
