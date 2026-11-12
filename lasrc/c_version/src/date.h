#ifndef DATE_H
#define DATE_H

#include <stdbool.h>
#include "espa_metadata.h"
#include "error_handler.h"

/* Date/time type definition */

#define MAX_DATE_LEN 28

typedef enum {
    DATE_FORMAT_DATEA_TIME,  /* yyyy-mm-ddThh:mm:ss.ssssssZ" */
    DATE_FORMAT_DATEB_TIME,  /* yyyy-dddThh:mm:ss.ssssssZ" */
    DATE_FORMAT_DATEA,       /* yyyy-mm-dd" */
    DATE_FORMAT_DATEB,       /* yyyy-ddd" */
    DATE_FORMAT_TIME         /* hh:mm:ss.ssssss" */
} Date_format_t;

typedef struct {
    bool fill;               /* is this fill? */
    int year;                /* year */
    int doy;                 /* day of year */
    int month;               /* month */
    int day;                 /* day of month */
    int hour;                /* hour */
    int minute;              /* minute */
    double second;           /* second */
    long jday2000;           /* Julian day circa 2000 */
    double sod;              /* seconds of day */
} Date_t;

bool date_init
(
    Date_t *this,           /* I/O: date structure to be initialized */
    char *s,                /* I: date string */
    Date_format_t iformat   /* I: format of the date string */
);

bool date_diff
(
    Date_t *d1,     /* I: first date for comparison */
    Date_t *d2,     /* I: second date for comparison */
    double *diff    /* O: difference between the dates (days) */
);

bool date_copy
(
    Date_t *this,   /* I: date to be copied */
    Date_t *copy    /* O: copy of the date */
);

bool format_date
(
    Date_t *this,             /* I: input date structure to be formatted */
    Date_format_t iformat,    /* I: date format type */
    char *s                   /* O: output date as a character string */
);

#endif
