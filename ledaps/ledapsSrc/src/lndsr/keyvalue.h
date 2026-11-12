#ifndef KEYVALUE_H
#define KEYVALUE_H

#include "lndsr.h"

/* Satellite type values */

const Key_string_t Sat_string[SAT_MAX] = {
  {(int)SAT_LANDSAT_1, "LANDSAT_1"},
  {(int)SAT_LANDSAT_2, "LANDSAT_2"},
  {(int)SAT_LANDSAT_3, "LANDSAT_3"},
  {(int)SAT_LANDSAT_4, "LANDSAT_4"},
  {(int)SAT_LANDSAT_5, "LANDSAT_5"},
  {(int)SAT_LANDSAT_7, "LANDSAT_7"}
};

/* Instrument type values */

const Key_string_t Inst_string[INST_MAX] = {
  {(int)INST_MSS, "MSS"},
  {(int)INST_TM,  "TM"},
  {(int)INST_ETM, "ETM"},
};

/* World Reference System (WRS) type values */

const Key_string_t Wrs_string[WRS_MAX] = {
  {(int)WRS_1, "1"},
  {(int)WRS_2, "2"},
};

/* Satellite type values */

const Key_string_t Ozsrc_string[OZSRC_MAX] = {
  {(int)OZSRC_NIMBUS7,         "NIMBUS7"},
  {(int)OZSRC_METEOR3,         "METEOR3"},
  {(int)OZSRC_EARTHPROBE,      "EARTHPROBE"},
  {(int)OZSRC_NIMBUS7_FILL,    "NIMBUS7_FILL"},
  {(int)OZSRC_METEOR3_FILL,    "METEOR3_FILL"},
  {(int)OZSRC_EARTHPROBE_FILL, "EARTHPROBE_FILL"},
  {(int)OZSRC_FILL,       "FILL"},
};
#endif

