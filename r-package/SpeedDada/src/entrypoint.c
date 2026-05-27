/* extendr entrypoint — do not edit manually.
   extendr_module! generates R_init_dada2rs_extendr; R calls R_init_dada2rs on load. */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_init_dada2rs_extendr(DllInfo *info);

void R_init_dada2rs(DllInfo *info) {
    R_init_dada2rs_extendr(info);
}
