#ifndef SR_H
#define SR_H

#include <stdbool.h>
#include "lndsr.h"
#include "lut.h"
#include "error.h"

typedef struct {
  bool first[NBAND_SR_MAX];
  int sr_min[NBAND_SR_MAX];
  int sr_max[NBAND_SR_MAX];
  long nfill[NBAND_SR_MAX];
  long nsatu[NBAND_SR_MAX];
  long nout_range[NBAND_SR_MAX];
} Sr_stats_t;


bool Sr(Lut_t *lut, int nsamp, int il, atmos_t *atmos_coef,
        atmos_t *interpol_atmos_coef, uint16_t **line_in, uint16_t *qa_line,
        uint16_t **line_out, Sr_stats_t *sr_stats);
#endif
