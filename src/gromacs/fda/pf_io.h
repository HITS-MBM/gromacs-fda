#ifndef pf_io_h
#define pf_io_h

#include "types/pf_array.h"

void pf_write_frame(t_pf_global *pf_global, rvec *x, gmx_mtop_t *top_global);
void pf_open(t_pf_global *pf_global);
void pf_close(t_pf_global *pf_global);

void pf_save_and_write_scalar_time_averages(t_pf_global *pf_global, rvec *x, gmx_mtop_t *top_global);
void pf_write_scalar_time_averages(t_pf_global *pf_global);

#endif  /* pf_io_h */
