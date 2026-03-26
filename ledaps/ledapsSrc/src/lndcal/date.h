#ifndef DATE_H
#define DATE_H

#include <stdbool.h>

/* Date/time type definition */

#define MAX_DATE_LEN (28)

typedef enum {
  DATE_FORMAT_DATEA_TIME,  /* yyyy-mm-ddThh:mm:ss.ssssssZ" */
  DATE_FORMAT_DATEB_TIME,  /* yyyy-dddThh:mm:ss.ssssssZ" */
  DATE_FORMAT_DATEA,       /* yyyy-mm-dd" */
  DATE_FORMAT_DATEB,       /* yyyy-ddd" */
  DATE_FORMAT_TIME         /* hh:mm:ss.ssssss" */
} Date_format_t;

typedef struct {
  bool fill;
  int year;
  int doy;
  int month;
  int day;
  int hour;
  int minute;
  double second;
  long jday2000;
  double sod;
} Date_t;

bool DateInit(Date_t *this, char *s, Date_format_t iformat);
bool DateDiff(Date_t *d1, Date_t *d2, double *diff);
bool DateCopy(Date_t *this, Date_t *copy);
bool FormatDate(Date_t *this, Date_format_t iformat, char *s);

#endif
