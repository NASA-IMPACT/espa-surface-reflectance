#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "read_grib_tools.h"
#include "grib.h"
#include "date.h"

int read_grib_anc
(
    t_ncep_ancillary *anc,
    int datatype
)
{
    FILE *fd = NULL;
    char where[50],tag[50],date[50];
    int i,grib_ret,ny,nx;
    short year,doy,month,day,hour,minute;
    float sec;

    switch (datatype) {
        case TYPE_OZONE_DATA:
            strcpy(tag,OZONE_GRIBTAG);
            strcpy(where,OZONE_GRIBLEVEL);
            break;
        case TYPE_WV_DATA:
            strcpy(tag,WV_GRIBTAG);
            strcpy(where,WV_GRIBLEVEL);
            break;
        case TYPE_SP_DATA:
            strcpy(tag,SP_GRIBTAG);
            strcpy(where,SP_GRIBLEVEL);
            break;
        case TYPE_ATEMP_DATA:
            strcpy(tag,TMP_GRIBTAG);
            strcpy(where,TMP_GRIBLEVEL);
            break;
        case TYPE_UWND_DATA:
            strcpy(tag,UWD_GRIBTAG);
            strcpy(where,UWD_GRIBLEVEL);
            break;
        case TYPE_VWND_DATA:
            strcpy(tag,VWD_GRIBTAG);
            strcpy(where,VWD_GRIBLEVEL);
            break;
        default:
            return -1;
    }        
    anc->latmin=-90;
    anc->latmax=90;
    anc->lonmin=-180;
    anc->lonmax=180;
    anc->deltalat=1;
    anc->deltalon=1;
    
    anc->nbrows=-1;
    anc->nbcols=-1;
    anc->year=-1;
    anc->doy=-1;
    for (i=0;i<anc->nblayers;i++) {
        printf("reading file %s\n",anc->filename[i]);
        if ((fd=fopen(anc->filename[i],"rb")) != NULL) {
            read_grib_date(fd, tag, where, date);
            printf("date=%s\n",date);
            sscanf(date,"%4hd-%2hd-%2hdT%2hd:%2hd:%f",&year,&month,&day,&hour,
                &minute,&sec);
            if (anc->year == -1)
                anc->year=year;
            else if (anc->year != year) {
                fprintf(stderr,"ERROR: inconsistent year in %s\n",
                    anc->filename[i]);
                return (-1);
            }
            doy=getdoy(year,month,day);
            if (anc->doy==-1)
                anc->doy=doy;
            else if (anc->doy != doy) {
                fprintf(stderr,"ERROR: inconsistent day in %s\n",
                    anc->filename[i]);
                return (-1);
            }
            anc->time[i]=sec/3600.+ (float)minute/60.+(float)hour;
            printf("date=%04d-%02d-%02dT%02d:%02d:%09.6f   %03d %09.6f\n",
                year,month,day,hour,minute,sec,anc->doy,anc->time[i]);
            
            grib_ret=read_grib_array(fd, tag, where, &ny, &nx, &(anc->data[i]));
            if (grib_ret != 0)
            {
                fprintf(stderr, "Error %d returned from read_grib_array().\n",
                        grib_ret);
                return -1;
            }
            if (anc->nbrows == -1)
                anc->nbrows = ny;
            else if (anc->nbrows != ny) {
                fprintf(stderr,"ERROR: inconsistent nbrows in %s\n",
                    anc->filename[i]);
                return (-1);
            }
            if (anc->nbcols == -1)
                anc->nbcols = nx;
            else if (anc->nbcols != nx) {
                fprintf(stderr,"ERROR: inconsistent ncols in %s\n",
                    anc->filename[i]);
                return (-1);
            }
            fclose(fd);
        }
        else
            return -1;
    }
    
    return 0;
}

void interpol_spatial_anc
(
    t_ncep_ancillary *anc,    /* I: ancillary structure information */
    float lat,                /* I: latitude */
    float lon,                /* I: longitude */
    float *value              /* O: interpolated anciliary data for this
                                    lat/long location (anc->nblayers values
                                    reside in this array) */
)
{
/* 
  Point order:

    0 ---- 1    +--> sample
    |      |    |
    |      |    v
    2 ---- 3   line

 */

    int j;
    int p0, p1, p2, p3;
    int row, col;
    float dl, ds;
    float grid_line, grid_sample;

    /* Ensure the lat/long values are within the grid. */
    if (lat > anc->latmax)
        lat = anc->latmax;
    else if (lat < anc->latmin)
        lat = anc->latmin;
    if (lon > anc->lonmax)
        lon = anc->lonmax;
    else if (lon < anc->lonmin)
        lon= anc->lonmin;

    /* Set the four corner point indices. */
    grid_line = (anc->latmax - lat)/anc->deltalat;
    grid_sample = (lon - anc->lonmin)/anc->deltalon;
    row = (int)grid_line;
    col = (int)grid_sample;
    p0 = row*anc->nbcols + col;
    p1 = p0 + 1;
    p2 = p0 + anc->nbcols;
    p3 = p2 + 1;

    /* Compute the fractional grid cell offset of the sample point from
       the upper-left corner of the cell. */
    dl = grid_line - row;
    ds = grid_sample - col;

    /* Apply the bilinear interpolation for the sample point. */
    for (j=0;j<anc->nblayers;j++)
    {
        float *data = anc->data[j];
        value[j] = data[p0]
                 + dl*(data[p2] - data[p0])
                 + ds*(data[p1] - data[p0])
                 + dl*ds*(data[p0] - data[p1] - data[p2] + data[p3]);
    }
}

void free_anc_data(t_ncep_ancillary *anc) {
    int i;
    for (i=0;i<anc->nblayers;i++)
        free(anc->data[i]);
}

void print_anc_data(t_ncep_ancillary *anc, char* ancftype)
{
    int i;
    printf("\n*********************************************\n");
    printf("**** anc file %s ***\n",ancftype);
    printf("**** source  = %s ***\n",anc->source);
    printf("**** nblayers= %d ***\n",anc->nblayers);
    printf("**** year,doy=(%d,%d) ***\n",anc->year,anc->doy);
    printf("**** timeres = %f ***\n",anc->timeres);
    printf("**** time    = ");
    for (i=0; i<anc->nblayers; i++)
        printf("%8.1f ",anc->time[i]);
    printf("****\n");
    printf("**** latmin = %f ***\n",anc->latmin);
    printf("**** latmax = %f ***\n",anc->latmax);
    printf("**** deltalat = %f ***\n",anc->deltalat);
    printf("**** lonmin = %f ***\n",anc->lonmin);
    printf("**** lonmax = %f ***\n",anc->lonmax);
    printf("**** deltalon = %f ***\n",anc->deltalon);
    printf("**** nbrows = %d ***\n",anc->nbrows);
    printf("**** nbcols = %d ***\n",anc->nbcols);
    printf("*********************************************\n");
}
