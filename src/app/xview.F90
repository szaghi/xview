!< xview.
program xview

use xview_mbpar_object
use xview_ui_object

implicit none

type(mbpar_object) :: mbpar !< mb.par, Xnavis input.
type(ui_object)    :: ui    !< User Interface.

call ui%parse_cli
if (ui%cli%is_passed(switch='--mbpar')) then
   call mbpar%load(path=ui%mbpar_path, filename=ui%mbpar_filename)
   print '(A)', mbpar%description()
endif
endprogram xview
