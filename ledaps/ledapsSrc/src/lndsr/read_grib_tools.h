#ifndef READ_GRIB_TOOLS_H
#define READ_GRIB_TOOLS_H

#define TYPE_OZONE_DATA	1
#define TYPE_WV_DATA	2
#define TYPE_SP_DATA	3
#define TYPE_ATEMP_DATA	4
#define TYPE_UWND_DATA	5
#define TYPE_VWND_DATA	6

#define OZONE_GRIBTAG    "TOZNE\0"
#define SP_GRIBTAG    "PRES\0"
#define WV_GRIBTAG       "PWAT\0"
#define TMP_GRIBTAG      "TMP\0"
#define UWD_GRIBTAG      "UGRD\0"
#define VWD_GRIBTAG      "VGRD\0"
#define OZONE_GRIBLEVEL  "atmos col\0"
#define SP_GRIBLEVEL  "sfc\0"
#define WV_GRIBLEVEL     "atmos col\0"
#define TMP_GRIBLEVEL    "2 m above gnd\0"
#define UWD_GRIBLEVEL    "10 m above gnd\0"
#define VWD_GRIBLEVEL    "10 m above gnd\0"

#define GRIB_WV_SP_TIMERES 6.0    /* hours */
#define GRIB_OZONE_TIMERES 24.0   /* hours */

#define MAX_NB_LAYERS 10

typedef struct {
	char source[256];
	short nblayers;
	short year,doy;
	float timeres,time[MAX_NB_LAYERS];
	char filename[MAX_NB_LAYERS][256];
	float *data[MAX_NB_LAYERS];
	float latmin,latmax,deltalat;
	float lonmin,lonmax,deltalon;
	int nbrows,nbcols;
} t_ncep_ancillary;

int read_grib_anc(t_ncep_ancillary *anc,int datatype);
void interpol_spatial_anc(t_ncep_ancillary *anc, float lat, float lon,
                          float *value);
void free_anc_data(t_ncep_ancillary *anc);
void print_anc_data(t_ncep_ancillary *anc, char* ancftype);

#endif
