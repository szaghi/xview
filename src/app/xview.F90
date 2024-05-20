!< xview.
program xview

use xview_mbpar_object
use xview_procinput_object
use xview_ui_object

implicit none

type(ui_object)        :: ui        !< User Interface.
type(procinput_object) :: procinput !< proc.input, Xnavis input.
type(mbpar_object)     :: mbpar     !< mb.par, Xnavis input.

call ui%parse_cli
if (ui%cli%is_passed(switch='--procinput')) then
   call procinput%load_file(path=ui%ipath, filename=ui%procinput_filename)
   print '(A)', procinput%description()
endif
if (ui%cli%is_passed(switch='--mbpar')) then
   call mbpar%load_file(path=ui%ipath, filename=ui%mbpar_filename)
   print '(A)', mbpar%description()
endif
endprogram xview
