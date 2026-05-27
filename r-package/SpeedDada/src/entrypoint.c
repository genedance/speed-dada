/* extendr entrypoint — do not edit manually.
   extendr_module! generates R_init_SpeedDada_extendr; R calls R_init_SpeedDada on load. */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_init_SpeedDada_extendr(DllInfo *info);

void R_init_SpeedDada(DllInfo *info) {
    R_init_SpeedDada_extendr(info);
}
