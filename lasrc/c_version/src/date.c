/*****************************************************************************
FILE: date.c
  
PURPOSE: Contains functions for handling dates and date formatting.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
*****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "date.h"
#include "error.h"

/******************************************************************************
MODULE:  date_init

PURPOSE:  Initializes the date structure from the input character string and
specified date format.

RETURN VALUE:
Type = bool
Value      Description
-----      -----------
false      Error occurred during the initialization
true       Successful completion

NOTES:
******************************************************************************/
bool date_init
(
    Date_t *this,           /* I/O: date structure to be initialized */
    char *s,                /* I: date string */
    Date_format_t iformat   /* I: format of the date string */
)
{
    char FUNC_NAME[] = "date_init";   /* function name */
    char errmsg[STR_SIZE];            /* error message */
    char *date = NULL;   /* date string */
    char *time = NULL;   /* time string */
    bool leap;           /* is this a leap year? */
    int year1;           /* temporary year variable */
    int nday[12] = {31, 29, 31, 30,  31,  30,  31,  31,  30,  31,  30,  31};
                         /* number of days in each month (leap year) */
    int idoy[12] = { 1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336};
                         /* starting day of each month (leap year) */
    int len;             /* string length */
    int jleap;           /* Julian days leap year */
    int idoy_nonleap;    /* day of year for non leap year */
  
    /* Validate the date format */
    this->fill = true;
    if (iformat != DATE_FORMAT_DATEA_TIME  && 
        iformat != DATE_FORMAT_DATEB_TIME  &&
        iformat != DATE_FORMAT_DATEA  &&
        iformat != DATE_FORMAT_DATEB)
    {
        strcpy (errmsg, "Invalid format parameter type");
        error_handler (true, FUNC_NAME, errmsg);
        return (false);
    }

    /* Grab the length of the input date string */
    len = strlen (s);

    /*printf("String = %s\n", s); */

    /* Handle each format, verifying the format itself (i.e. string length and
       other characters) then separate the date and the time based on the
       characteristics of the format.  This will obtain a pointer to the
       date and the time, depending on the format type. */
    if (iformat == DATE_FORMAT_DATEA_TIME)
    {
        /* Validate the string for the current format type */
        if (len < 20 || len > 27) 
        {
            strcpy (errmsg, "Invalid date/time string length");
            error_handler (true, FUNC_NAME, errmsg);
            return (false);
        }

        if (s[10] != 'T' || s[len - 1] != 'Z')
        {
            strcpy (errmsg, "Invalid date/time format");
            error_handler (true, FUNC_NAME, errmsg);
            return (false);
        }

        /* Grab the date and time pointer */
        date = &s[0];
        time = &s[11];
    }
    else if (iformat == DATE_FORMAT_DATEB_TIME)
    {
        /* Validate the string for the current format type */
        if (len < 18 || len > 25) 
        {
            strcpy (errmsg, "Invalid date/time string length");
            error_handler (true, FUNC_NAME, errmsg);
            return (false);
        }

        if (s[8] != 'T' || s[len - 1] != 'Z')
        {
            strcpy (errmsg, "Invalid date/time format");
            error_handler (true, FUNC_NAME, errmsg);
            return (false);
        }

        /* Grab the date and time pointer */
        date = &s[0];
        time = &s[9];
    }
    else if (iformat == DATE_FORMAT_DATEA)
    {
        /* Validate the string for the current format type */
        if (len != 10) 
        {
            strcpy (errmsg, "Invalid date/time string length");
            error_handler (true, FUNC_NAME, errmsg);
            return (false);
        }

        /* Grab the date pointer */
        date = s;
    }
    else if (iformat == DATE_FORMAT_DATEB)
    {
        /* Validate the string for the current format type */
        if (len != 8) 
        {
            strcpy (errmsg, "Invalid date/time string length");
            error_handler (true, FUNC_NAME, errmsg);
            return (false);
        }

        /* Grab the date pointer */
        date = s;
    }
  
    /* Using the date pointer, separate the year, day, month */
    if (iformat == DATE_FORMAT_DATEA_TIME || iformat == DATE_FORMAT_DATEA)
    {
        if (sscanf (date, "%4d-%2d-%2d",
            &this->year, &this->month, &this->day) != 3) 
        {
            strcpy (errmsg, "Invalid date format");
            error_handler (true, FUNC_NAME, errmsg);
            return (false);
        }

        if (this->year < 1900 || this->year > 2400) 
        {
            strcpy (errmsg, "Invalid year");
            error_handler (true, FUNC_NAME, errmsg);
            return (false);
        }

        if (this->month < 1 || this->month > 12) 
        {
            strcpy (errmsg, "Invalid month");
            error_handler (true, FUNC_NAME, errmsg);
            return (false);
        }

        if (this->day < 1 || this->day > nday[this->month-1])
        {
            strcpy (errmsg, "Invalid day of month");
            error_handler (true, FUNC_NAME, errmsg);
            return (false);
        }
        this->doy = this->day + idoy[this->month - 1] - 1;
    }
    else
    {
        if (sscanf (date, "%4d-%3d", &this->year, &this->doy) != 2) 
        {
            strcpy (errmsg, "Invalid date format");
            error_handler (true, FUNC_NAME, errmsg);
            return (false);
        }

        if (this->year < 1900 || this->year > 2400) 
        {
            strcpy (errmsg, "Invalid year");
            error_handler (true, FUNC_NAME, errmsg);
            return (false);
        }

        if (this->doy < 1 || this->doy > 366) 
        {
            strcpy (errmsg, "Invalid day of year");
            error_handler (true, FUNC_NAME, errmsg);
            return (false);
        }
    }
  
    /* Is this a leap year? */
    leap = (bool) (this->year % 4 == 0 &&
        (this->year % 100 != 0 || this->year % 400 == 0)); 

    /* Handle the leap years */
    if (iformat == DATE_FORMAT_DATEA_TIME ||
        iformat == DATE_FORMAT_DATEA)
    {
        if ((this->month == 2) && !leap && (this->day > 28))
        {
            strcpy (errmsg, "Invalid day of month");
            error_handler (true, FUNC_NAME, errmsg);
            return (false);
        }
        if (!leap && (this->month > 2))
            this->doy--;
    }
    else
    {
        /* Determine the month and day for the current year */
        if (leap)
        {
            for (this->month = 0; this->month < 12; this->month++)
                if (this->doy < idoy[this->month])
                    break;
        }
        else
        {
            if (this->doy > 365) 
            {
                strcpy (errmsg, "Invalid day of year");
                error_handler (true, FUNC_NAME, errmsg);
                return (false);
            }
            for (this->month = 0; this->month < 12; this->month++)
            {
                idoy_nonleap = idoy[this->month];
                if (this->month > 1)
                    idoy_nonleap--;
                if (this->doy < idoy_nonleap)
                    break;
            }
        }
    }
  
    /* Convert to Julian days ca. 2000 (1 = Jan. 1, 2000) */
    year1 = this->year - 1900;
    if (year1 > 0)
    {
        jleap = (year1 - 1) / 4;
        if (this->year > 2100)
            jleap -= (this->year - 2001) / 100;
    }
    else
        jleap = 0;

    this->jday2000 = (year1 * 365) + jleap + this->doy;
    this->jday2000 -= 36524;
  
    /* Parse and check the time.  If NULL then set to 0. */
    if (time != NULL)
    {
        if (sscanf (time, "%2d:%2d:%lf", 
            &this->hour, &this->minute, &this->second) != 3)
        {
            strcpy (errmsg, "Invalid time format");
            error_handler (true, FUNC_NAME, errmsg);
            return (false);
        }
    }
    else
    {
        this->hour = 0;
        this->minute = 0;
        this->second = 0.0;
    }

    if (this->hour < 0 || this->hour > 23)
    {
        strcpy (errmsg, "Invalid hour");
        error_handler (true, FUNC_NAME, errmsg);
        return (false);
    }
    if (this->minute < 0 || this->minute > 59)
    {
        strcpy (errmsg, "Invalid minute");
        error_handler (true, FUNC_NAME, errmsg);
        return (false);
    }
    if (this->second < 0.0 || this->second > 59.999999)
    {
        strcpy (errmsg, "Invalid second");
        error_handler (true, FUNC_NAME, errmsg);
        return (false);
    }
  
    /* Convert to seconds of day */
    this->sod = (((this->hour * 60) + this->minute) * 60) + this->second;
    this->fill = false;
  
    return (true);
}


/******************************************************************************
MODULE:  date_diff

PURPOSE:  Computes the difference between two dates.

RETURN VALUE:
Type = bool
Value      Description
-----      -----------
false      Error occurred computing the difference
true       Successful completion

NOTES:
******************************************************************************/
bool date_diff
(
    Date_t *d1,     /* I: first date for comparison */
    Date_t *d2,     /* I: second date for comparison */
    double *diff    /* O: difference between the dates (days) */
)
{
    char FUNC_NAME[] = "date_diff";   /* function name */
    char errmsg[STR_SIZE];            /* error message */

    /* Validate the dates */
    if (d1 == NULL || d2 == NULL)
    {
        strcpy (errmsg, "Invalid date structure");
        error_handler (true, FUNC_NAME, errmsg);
        return (false);
    }
  
    if (d1->fill || d2->fill) 
    {
        strcpy (errmsg, "One or both of the date structures are fill");
        error_handler (true, FUNC_NAME, errmsg);
        return (false);
    }
  
    /* Compute the difference in days */
    *diff = d1->jday2000 - d2->jday2000;
    *diff += (d1->sod - d2->sod) / 86400.0;     /* convert seconds to days */
  
    return (true);
}


/******************************************************************************
MODULE:  date_copy

PURPOSE:  Copies one date structure to another.

RETURN VALUE:
Type = bool
Value      Description
-----      -----------
false      Error occurred during the copy
true       Successful completion

NOTES:
******************************************************************************/
bool date_copy
(
    Date_t *this,   /* I: date to be copied */
    Date_t *copy    /* O: copy of the date */
)
{
    char FUNC_NAME[] = "date_copy";   /* function name */
    char errmsg[STR_SIZE];            /* error message */

    /* Validate the date */
    if (this == NULL || copy == NULL) 
    {
        strcpy (errmsg, "Invalid date structure");
        error_handler (true, FUNC_NAME, errmsg);
        return (false);
    }
  
    /* Copy the date fields */
    copy->fill = this->fill;
    copy->year = this->year;
    copy->doy = this->doy;
    copy->month = this->month;
    copy->day = this->day;
    copy->hour = this->hour;
    copy->minute = this->minute;
    copy->second = this->second;
    copy->jday2000 = this->jday2000;
    copy->sod = this->sod;
  
    return (true);
}


/******************************************************************************
MODULE:  format_date

PURPOSE:  
specified date format.

RETURN VALUE:
Type = bool
Value      Description
-----      -----------
false      Error occurred during the initialization
true       Successful completion

NOTES:
******************************************************************************/
bool format_date
(
    Date_t *this,             /* I: input date structure to be formatted */
    Date_format_t iformat,    /* I: date format type */
    char *s                   /* O: output date as a character string */
)
{
    char FUNC_NAME[] = "format_date";   /* function name */
    char errmsg[STR_SIZE];              /* error message */

    /* Validate the date */
    if (this == NULL)
    {
        strcpy (errmsg, "Invalid date structure");
        error_handler (true, FUNC_NAME, errmsg);
        return (false);
    }
  
    /* Format the date based on the formatting type */
    if (iformat == DATE_FORMAT_DATEA_TIME)
    {
        if (sprintf (s, "%4d-%02d-%02dT%02d:%02d:%09.6fZ", this->year,
            this->month, this->day, this->hour, this->minute, this->second)
            < 0) 
        {
            strcpy (errmsg, "Formatting the date and time");
            error_handler (true, FUNC_NAME, errmsg);
            return (false);
        }
    }
    else if (iformat == DATE_FORMAT_DATEB_TIME)
    {
        if (sprintf (s, "%4d-%03dT%02d:%02d:%09.6fZ", this->year, this->doy, 
            this->hour, this->minute, this->second) < 0) 
        {
            strcpy (errmsg, "Formatting the date and time");
            error_handler (true, FUNC_NAME, errmsg);
            return (false);
        }
    }
    else if (iformat == DATE_FORMAT_DATEA)
    {
        if (sprintf (s, "%4d-%02d-%02d", this->year, this->month, this->day)
            < 0) 
        {
            strcpy (errmsg, "Formatting the date");
            error_handler (true, FUNC_NAME, errmsg);
            return (false);
        }
    }
    else if (iformat == DATE_FORMAT_DATEB)
    {
        if (sprintf (s, "%4d-%03d", this->year, this->doy) < 0) 
        {
            strcpy (errmsg, "Formatting the date");
            error_handler (true, FUNC_NAME, errmsg);
            return (false);
        }
    }
    else if (iformat == DATE_FORMAT_TIME)
    {
        if (sprintf(s, "%02d:%02d:%09.6f", this->hour, this->minute,
            this->second) < 0) 
        {
            strcpy (errmsg, "Formatting the time");
            error_handler (true, FUNC_NAME, errmsg);
            return (false);
        }
    }
    else 
    {
        strcpy (errmsg, "Formatting the time");
        error_handler (true, FUNC_NAME, errmsg);
        return (false);
    }
    
    return (true);
}
