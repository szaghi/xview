!< xview.
program xview

use xview_file_grd_object
use xview_file_icc_object
use xview_file_rst_object
use xview_file_mbpar_object
use xview_file_procinput_object
use xview_ui_object

implicit none

type(ui_object)             :: ui             !< User Interface.
type(file_grd_object)       :: file_grd       !< File grd.
type(file_icc_object)       :: file_icc       !< File icc.
type(file_rst_object)       :: file_rst       !< File rst.
type(file_procinput_object) :: file_procinput !< proc.input.
type(file_mbpar_object)     :: file_mbpar     !< mb.par.

call ui%parse_cli
if (ui%cli%is_passed(switch='--procinput')) then
   call file_procinput%load_file(path=ui%ipath, filename=ui%procinput_filename)
   print '(A)', file_procinput%description()
endif
if (ui%cli%is_passed(switch='--mbpar')) then
   call file_mbpar%load_file(path=ui%ipath, filename=ui%mbpar_filename)
   print '(A)', file_mbpar%description()
endif
endprogram xview
