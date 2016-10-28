#ifndef pf_io_h
#define pf_io_h

#include "fda.h"

void pf_write_frame(FDA *fda, rvec const *x, gmx_mtop_t *top_global);
void pf_open(FDA *fda);
void pf_close(FDA *fda);

void pf_save_and_write_scalar_time_averages(FDA *fda, rvec const *x, gmx_mtop_t *top_global);
void pf_write_scalar_time_averages(FDA *fda);

#endif  /* pf_io_h */
