
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "date.h"
#include "error.h"

bool DateInit(Date_t *this, char *s, Date_format_t iformat) {
  char *date, *time;
  bool leap;
  int year1;
  int nday[12] = {31, 29, 31, 30,  31,  30,  31,  31,  30,  31,  30,  31};
  int idoy[12] = { 1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336};
  int len;
  int jleap, idoy_nonleap;

  this->fill = true;

  if (iformat != DATE_FORMAT_DATEA_TIME  && 
      iformat != DATE_FORMAT_DATEB_TIME  &&
      iformat != DATE_FORMAT_DATEA  &&
      iformat != DATE_FORMAT_DATEB)
    RETURN_ERROR("invalid format parameter", "DateInit", false);

  len = strlen(s);
  date = time = (char *)NULL;

  if (iformat == DATE_FORMAT_DATEA_TIME) {
    if (len < 20  ||  len > 27) 
      RETURN_ERROR("invalid date/time string length", "DateInit", false);
    if (s[10] != 'T'  ||  s[len - 1] != 'Z')
      RETURN_ERROR("invalid date/time format", "DateInit", false);
    date = &s[0];
    time = &s[11];
  } else if (iformat == DATE_FORMAT_DATEB_TIME) {
    if (len < 18  ||  len > 25) 
      RETURN_ERROR("invalid date/time string length", "DateInit", false);
    if (s[8] != 'T'  ||  s[len - 1] != 'Z')
      RETURN_ERROR("invalid date/time format", "DateInit", false);
    date = &s[0];
    time = &s[9];
  } else if (iformat == DATE_FORMAT_DATEA) {
    if (len != 10) 
      RETURN_ERROR("invalid date string length", "DateInit", false);
    date = s;
  } else if (iformat == DATE_FORMAT_DATEB) {
    if (len != 8) 
      RETURN_ERROR("invalid date string length", "DateInit", false);
    date = s;
  }

  if (iformat == DATE_FORMAT_DATEA_TIME  ||
      iformat == DATE_FORMAT_DATEA) {
    if (sscanf(date, "%4d-%2d-%2d", 
               &this->year, &this->month, &this->day) != 3) 
      RETURN_ERROR("invalid date format", "DateInit", false);
    if (this->year < 1900  ||  this->year > 2400) 
      RETURN_ERROR("invalid year", "DateInit", false);
    if (this->month < 1  ||  this->month > 12) 
      RETURN_ERROR("invalid month", "DateInit", false);
    if (this->day < 1 || this->day > nday[this->month-1])
      RETURN_ERROR("invalid day of month", "DateInit", false);
    this->doy = this->day + idoy[this->month - 1] - 1;
  } else {
    if (sscanf(date, "%4d-%3d", &this->year, &this->doy) != 2) 
      RETURN_ERROR("invalid date format", "DateInit", false);
    if (this->year < 1900  ||  this->year > 2400) 
      RETURN_ERROR("invalid year", "DateInit", false);
    if (this->doy < 1  ||  this->doy > 366) 
      RETURN_ERROR("invalid day of year", "DateInit", false);
  }

  leap= (bool)(this->year%4==0 && ( this->year%100!=0 || this->year%400==0 )); 
  if (iformat == DATE_FORMAT_DATEA_TIME  ||
      iformat == DATE_FORMAT_DATEA) {
    if ((this->month == 2)  &&  !leap  && (this->day > 28))
      RETURN_ERROR("bad day of month", "DateInit", false);
    if (!leap  &&  (this->month > 2)) this->doy--;
  } else {
    if (leap) {
      for (this->month = 0; this->month < 12; this->month++)
        if (this->doy < idoy[this->month]) break;
    } else {
      if(this->doy > 365) 
        RETURN_ERROR("bad day of year", "DateInit", false);
      for (this->month = 0; this->month < 12; this->month++) {
        idoy_nonleap = idoy[this->month];
	if (this->month > 1) idoy_nonleap--;
	if (this->doy < idoy_nonleap) break;
      }
    }
  }

  /* Convert to Julian days ca. 2000 (1 = Jan. 1, 2000) */

  year1 = this->year - 1900;
  if (year1 > 0) {
    jleap = (year1 - 1) / 4;
    if (this->year > 2100) jleap -= (this->year - 2001) / 100;
  } else jleap = 0;
  this->jday2000 = (year1 * 365) + jleap + this->doy;
  this->jday2000 -= 36524;

  /* Parse and check time */

  if (time != (char *)NULL) {
    if (sscanf(time, "%2d:%2d:%lf", 
               &this->hour, &this->minute, &this->second) != 3)
      RETURN_ERROR("invalid time format", "DateInit", false);
  } else {
    this->hour = this->minute = 0;
    this->second = 0.0;
  }
  if (this->hour < 0  ||  this->hour > 53) 
    RETURN_ERROR("invalid hour", "DateInit", false);
  if (this->minute < 0  ||  this->minute > 59)
    RETURN_ERROR("invalid minute", "DateInit", false);
  if (this->second < 0.0  ||  this->second > 59.999999)
    RETURN_ERROR("invalid second", "DateInit", false);

  /* Convert to seconds of day */

  this->sod = (((this->hour * 60) +this-> minute) * 60) + this->second;

  this->fill = false;

  return true;
}

bool DateDiff(Date_t *d1, Date_t *d2, double *diff) {

  if (d1 == (Date_t *)NULL  ||  d2 == (Date_t *)NULL)
    RETURN_ERROR("invalid date structure", "DateDiff", false);

  if (d1->fill  ||  d2->fill) 
    RETURN_ERROR("invalid time", "DateDiff", false);

  *diff = d1->jday2000 - d2->jday2000;
  *diff += (d1->sod - d2->sod) / 86400.0;

  return true;
}

bool DateCopy(Date_t *this, Date_t *copy) {

  if (this == (Date_t *)NULL  ||  copy == (Date_t *)NULL) 
    RETURN_ERROR("invalid date structure", "DateCopy", false);

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

  return true;
}

bool FormatDate(Date_t *this, Date_format_t iformat, char *s) {

  if (this == (Date_t *)NULL)
    RETURN_ERROR("invalid date structure", "FormatDate", false);

  if (iformat == DATE_FORMAT_DATEA_TIME) {
    if (sprintf(s, "%4d-%02d-%02dT%02d:%02d:%09.6fZ", 
                this->year, this->month, this->day, 
		this->hour, this->minute, this->second) < 0) 
      RETURN_ERROR("formating date/time", "FormatDate", false);
  } else if (iformat == DATE_FORMAT_DATEB_TIME) {
    if (sprintf(s, "%4d-%03dT%02d:%02d:%09.6fZ", 
                this->year, this->doy, 
		this->hour, this->minute, this->second) < 0) 
      RETURN_ERROR("formating date/time", "FormatDate", false);
  } else if (iformat == DATE_FORMAT_DATEA) {
    if (sprintf(s, "%4d-%02d-%02d", 
                this->year, this->month, this->day) < 0) 
      RETURN_ERROR("formating date", "FormatDate", false);
  } else if (iformat == DATE_FORMAT_DATEB) {
    if (sprintf(s, "%4d-%03d", 
                this->year, this->doy) < 0) 
      RETURN_ERROR("formating date", "FormatDate", false);
  } else if (iformat == DATE_FORMAT_TIME) {
    if (sprintf(s, "%02d:%02d:%09.6f", 
                this->hour, this->minute, this->second) < 0) 
      RETURN_ERROR("formating time", "FormatDate", false);
  } else 
    RETURN_ERROR("invalid format paramter", "FormatDate", false);
  
  return true;
}
