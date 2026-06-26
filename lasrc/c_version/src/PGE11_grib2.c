/* NODAAC may need to be defined when running this code
   in any non-DAAC processing mode.  This feature was added
   by MODIS SDST to ensure compliance with ESDIS standards upon
   delivering this software to the DAAC. */

#ifdef NODAAC

#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "PGE11_grb2.h"
#include "PGE11_grib2.h"

extern unsigned char *mem_buffer[N_mem_buffers];
extern size_t mem_buffer_size[N_mem_buffers];
extern size_t mem_buffer_pos[N_mem_buffers];

#define UINT2(a, b) ((int)((a << 8) + (b)))
#define GB2_Center(sec) UINT2(sec[1][5], sec[1][6])
#define GB2_Subcenter(sec) UINT2(sec[1][7], sec[1][8])

/* v6.4.11, (23-MAR-21)

Note: grib2-format reads can be performed in one of two ways: 1) with the code in this file
(PGE11_grib2.c) only, or 2) with the libwgrib2.a library, which can be obtained from NOAA
and built with the command "make <makefile> lib".

The code here defaults to using only this file, PGE11_grib2.c, self-contained, and not using
the library.  The routines used by read_grib2_array() and read_grib2_date() have the same
names as those in the library, but are rewritten here, with some modifications.

If we have to use the libwgrib2.a library (some parameter we need can't be read by this code),
we must change the makefiles (MOD_PR09.mk, MOD_PR09DB.mk) like this:

From this:

linux | linux64 ) \
$(MAKE) $(TARGET) -f MOD_PR09.mk \
ADD_CFLAGS="-O0 -Df2cFortran -DLINUX -D_POSIX_SOURCE -DOLDSDSNAMES "\
ADD_FFLAGS="-O0 -D_POSIX_SOURCE " \
LIB="-L$(PGSLIB) -lPGSTK \
-L$(HDFLIB) -lmfhdf -ldf -ljpeg -lsz \
-L$(HDFEOS_LIB) -lhdfeos \
-L$(HDFEOS5_LIB) -lhe5_hdfeos -lGctp \
-L$(HDF5LIB) -lhdf5 -lz -ldl -lgfortran -lm"\
To this:

linux | linux64 ) \
$(MAKE) $(TARGET) -f MOD_PR09.mk \
ADD_CFLAGS="-O0 -fopenmp -DUSE_WGRIB2_LIBRARY -Df2cFortran -DLINUX -D_POSIX_SOURCE -DOLDSDSNAMES "\
ADD_FFLAGS="-O0 -D_POSIX_SOURCE " \
LIB="-L$(PGSLIB) -lPGSTK <whatever path>libwgrib2.a \
-L$(HDFLIB) -lmfhdf -ldf -ljpeg -lsz \
-L$(HDFEOS_LIB) -lhdfeos \
-L$(HDFEOS5_LIB) -lhe5_hdfeos -lGctp \
-L$(HDF5LIB) -lhdf5 -lz -ldl -lgfortran -lm"\

The definition 'USE_WGRIB2_LIBRARY' removes the code I took from the library and wrote here
(with modification) from compilation, and forces the routines in read_grib2_array() and
read_grib2_date() to use the library.


*/

int read_grib2_array(char *filename, char *what, char *where, int *ny, int *nx, float **narray);
int read_grib2_date(char *filename, char *what, char *where, char *date, int *frcst_val);

int read_grib2_array(char *filename, char *what, char *where, int *ny, int *nx, float **narray)
/*
   !C*****************************************************************************
   !Description:  Routine reads data from a GRIB2-format file.  The data to be read
   has the character strings 'what' and 'where' in its internal
   GRIB2-format label.

   !Input Parameters:
   filename: GRIB2-format filename
   what:     pointer to character string of what is to be looked for
   where:    pointer to character string of the height

   Stuff typically used:

   variable               what     where                                              units
   ~~~~~~~~               ~~~~     ~~~~~                                              ~~~~~

   precipitable water     "PWAT"    "entire atmosphere (considered as a single layer)" kg/m^2
   pressure               "PRMSL"   "mean sea level"                                   Pa
   ozone                  "TOZNE"   "entire atmosphere (considered as a single layer)" Dobsons

   Note: the exact syntax used for 'what' and 'where' depends upon the source
   of the GRIB2 file.

   !Output Parameters:
   n_rows:   the number of rows and
   n_cols:       columns of data retrieved
   narray:   array that contains the data

   !Input/Output Parameters:
   narray -- should be null when this routine is called, should be non-null and filled
   when routine returns.

   !Returns:
   1  parameter not found
   0  success

   !Revision History:
   Revision 0.0, 23-MAR-21, Jim Ray (SSAI)

   Original version:
   v1.7.0b1 (8-24-98) Wesley Ebisuzaki (NCEP/NCAR Reanalysis Project).

   !Team-unique Header:

   This software developed by the MODIS Land Science Team for
   the National Aeronautics and Space Administration, Goddard
   Space Flight Center, under contract NAS5-96062

   !References and Credits

   Jim Ray
   Science Systems and Applications Inc.
   NASA's Goddard Space Flight Center Code 619
   Greenbelt MD, 20771
   james.p.ray@nasa.gov, (301) 614-6684

   !Design Notes
   !END***************************************************************************
   */
{
    struct seq_file in_file;
    unsigned char *sec[10], *msg;
    long int pos;
    unsigned long int len;
    int num_submsgs;
    unsigned char *rd_grib2_msg_seq_file(unsigned char **sec, struct seq_file *input, long int *pos,
                                         unsigned long int *len, int *num_submsgs);
    int getName(unsigned char **sec, int mode, char *inv_out, char *name, char *desc, char *unit);
    unsigned int nnpnts;
    int nres, nscan, res, scan;
    int parse_next_msg(unsigned char *sec[10]);
    int parse_1st_msg(unsigned char *sec[10]);
    int Y, M, D, H, MN, S, msg_no, submsg;
    char name[200], desc[200], unit[200];
    char *inv_out;
    int table_4_5, scale_factor, scale_value;
    int level_type1, level_type2;
    float val1, val2;
    int center, subcenter;
    int undef_val1, undef_val2;
    void fixed_surfaces(unsigned char **sec, int *level_type1, float *val1, int *undef_val1, int *level_type2,
                        float *val2, int *undef_val2);
    int level2(int mode, int type1, int undef_val1, float value1, int type2, int undef_val2, float value2, int center,
               int subcenter, char *inv_out);
    int basic_ang, sub_ang;
    double units;
    int sub_angle(unsigned const char *p);
    void flip_data(float *tmpfltarray, float *tmp, int ny, int nx);
    void reverse_data(float *tmpfltarray, float *tmp, int ny, int nx);
    char ASCIIdate[28];
    int i, j, flip, reverse;
    float *tmp;
    int retval = 0;
    char location[50] = "PGE11_grib2.c, read_grib2_array\0";

    flip = reverse = 0;
    j = 1;
    inv_out = (char *)malloc(100 * sizeof(char));

    fopen_file(&in_file, filename, "rb");

    rd_grib2_msg_seq_file(sec, &in_file, &pos, &len, &num_submsgs);

    parse_1st_msg(sec);

    getName(sec, 1, NULL, name, NULL, unit); /* gets name and units of parameters... */

    {
        unsigned int nxTmp, nyTmp;
        get_nxny_(sec, &nxTmp, &nyTmp, &nnpnts, &nres, &nscan); /* gets dimensions of parameter */
        *nx = nxTmp;
        *ny = nyTmp;
    }

    center = GB2_Center(sec);
    subcenter = GB2_Subcenter(sec);
    fixed_surfaces(sec, &level_type1, &val1, &undef_val1, &level_type2, &val2, &undef_val2);
    level2(1, level_type1, undef_val1, val1, level_type2, undef_val2, val2, center, subcenter,
           inv_out); /* gets level of parameter */

    // printf("%s -- %s\n", name, inv_out);

    if (((strcmp(name, what) && !strcmp(inv_out, where))) || ((!strcmp(name, what) && strcmp(inv_out, where))) ||
        ((strcmp(name, what) && strcmp(inv_out, where)))) {
        do {
            // printf("%d %s %s %d -- will skip ...\n", j, name, inv_out, pos);
            for (i = 0; i < 100; i++) inv_out[i] = '\0';
            msg = rd_grib2_msg_seq_file(sec, &in_file, &pos, &len, &num_submsgs);
            if (msg == NULL) break;
            j++;

            parse_1st_msg(sec);

            getName(sec, 1, NULL, name, NULL, unit);

            {
                unsigned int nxTmp, nyTmp;
                get_nxny_(sec, &nxTmp, &nyTmp, &nnpnts, &nres, &nscan);
                *nx = nxTmp;
                *ny = nyTmp;
            }

            center = GB2_Center(sec);
            subcenter = GB2_Subcenter(sec);
            fixed_surfaces(sec, &level_type1, &val1, &undef_val1, &level_type2, &val2, &undef_val2);
            level2(1, level_type1, undef_val1, val1, level_type2, undef_val2, val2, center, subcenter, inv_out);

            /*
             * Attempts to read units -- this is what wgrib2 itself does for grids like those
             in our grib2 file -- and like wgrib2, it always comes up 0.000001.  But it
             doesn't seem to matter since the actual read of the data ( unpk_grib() below)
             gives us the data scaled correctly (10-MAR-21)
             */
            /*
               basic_ang = GDS_LatLon_basic_ang(sec[3]);
               sub_ang = GDS_LatLon_sub_ang(sec[3]);
               units = (basic_ang == 0 ?  0.000001 : ((float) basic_ang / (float) sub_ang));
               */

        } while ((strcmp(name, what)) || (strcmp(inv_out, where)));
    }

    if ((!strcmp(name, what)) && (!strcmp(inv_out, where))) {
        // printf("FOUND -- %s -- %s\n", name, inv_out);
        // printf("FOUND %d %d %d %d %d '%s' %s -- %s \n", *nx, *ny, nnpnts, pos, nscan, name, unit, inv_out);

        *narray = (float *)malloc(nnpnts * sizeof(float));

        unpk_grib(sec, *narray); /* gets parameter */

        if (nscan == 1) flip = 1; /* data is typically with south pole on top, north at the bottom */
        reverse = 1; /* data is always "WE", e. g., running from 0 (longitude) to +180, then -180 up to zero --
                        reverse() puts it from -180 to zero to 180. */

        if (flip == 1) {
            tmp = (float *)malloc(nnpnts * sizeof(float));
            flip_data(*narray, tmp, *ny, *nx);
            free(tmp);
        }
        if (reverse == 1) {
            tmp = (float *)malloc(nnpnts * sizeof(float));
            reverse_data(*narray, tmp, *ny, *nx);
            free(tmp);
        }

    } else {
        /*printf("Parameter '%s' at level '%s' not found\n", what, where);*/
        retval = 1;
    }
    fseek(in_file.cfile, 0L, SEEK_SET);
    fclose_file(&in_file);

    return (retval);
}

void flip_data(float *tmpfltarray, float *tmp, int ny, int nx) {
    long k, lj, ljj;
    int i, j, jj;
    /* flip latitude-wise */
    lj = 0L;
    for (i = (ny - 1); i >= 0; i--)
        for (j = 0; j < nx; j++) {
            ljj = i * nx + j;
            tmp[ljj] = tmpfltarray[lj];
            lj++;
        }
    for (lj = 0L; lj < ny * nx; lj++) tmpfltarray[lj] = tmp[lj];
    return;
}

void reverse_data(float *tmpfltarray, float *tmp, int ny, int nx) {
    long k, lj, ljj;
    int i, j, jj;
    /* 15-JUN-99: Reformat to start at longitude = -180 (or so). */
    lj = 0L;
    for (i = 0; i < ny; i++)
        for (j = 0; j < nx; j++) {
            jj = j + (nx / 2) - 1;
            if (jj >= nx) jj -= nx;
            ljj = i * nx + jj;
            tmp[ljj] = tmpfltarray[lj];
            lj++;
        }
    for (lj = 0L; lj < ny * nx; lj++) tmpfltarray[lj] = tmp[lj];
    return;
}

int read_grib2_date(char *filename, char *what, char *where, char *date, int *frcst_val)
/*
   !C*****************************************************************************
   !Description:  Routine reads data from a GRIB2-format file.  The data to be read
   has the character strings 'what' and 'where' in its internal
   GRIB2-format label.

   !Input Parameters:
   filename: GRIB2-format filename
   what:     pointer to character string of what is to be looked for
   where:    pointer to character string of the height; e. g., "atmos col" for whole-column
   data, "850 mb" for data at 850 mb., etc.

   For more information about 'what' amd 'where' syntax, see notes for routine
   read_grib2_array(), above.


   !Output Parameters:
   date      the date of data at 'what' and 'where'.  'date' will be unchanged if data at
   'what' and 'where' is not found; 'date' is in ASCII time code A (28 characters).
   frcst_val how many hours the data is forecast for.


   !Input/Output Parameters:
   none

   !Returns:
   1  parameter not found
   0  success

   !Revision History:
   Revision 0.0, 23-MAR-21, Jim Ray (SSAI)
   1.  Edited "main()" of file wgrib.c into this routine.

   Original version:
   v1.7.0b1 (8-24-98) Wesley Ebisuzaki (NCEP/NCAR Reanalysis Project).

   !Team-unique Header:

   This software developed by the MODIS Land Science Team for
   the National Aeronautics and Space Administration, Goddard
   Space Flight Center, under contract NAS5-96062

   !References and Credits

   Jim Ray
   Science Systems and Applications Inc.
   NASA's Goddard Space Flight Center Code 923
   Greenbelt MD, 20771
   jim@cratmos.gsfc.nasa.gov, (301) 614-6684

   !Design Notes
   !END***************************************************************************
   */
{
    struct seq_file in_file;
    unsigned char *sec[10], *msg;
    long int pos;
    unsigned long int len;
    int num_submsgs;
    unsigned char *rd_grib2_msg_seq_file(unsigned char **sec, struct seq_file *input, long int *pos,
                                         unsigned long int *len, int *num_submsgs);
    int getName(unsigned char **sec, int mode, char *inv_out, char *name, char *desc, char *unit);
    unsigned int nnpnts;
    int nres, nscan, res, scan;
    int parse_next_msg(unsigned char *sec[10]);
    int parse_1st_msg(unsigned char *sec[10]);
    int Y, M, D, H, MN, S, msg_no, submsg;
    char name[200], desc[200], unit[200];
    char *inv_out;
    int table_4_5, scale_factor, scale_value;
    int level_type1, level_type2;
    float val1, val2;
    int center, subcenter;
    int undef_val1, undef_val2;
    void fixed_surfaces(unsigned char **sec, int *level_type1, float *val1, int *undef_val1, int *level_type2,
                        float *val2, int *undef_val2);
    int level2(int mode, int type1, int undef_val1, float value1, int type2, int undef_val2, float value2, int center,
               int subcenter, char *inv_out);
    int basic_ang, sub_ang;
    double units;
    int sub_angle(unsigned const char *p);
    int i, j;
    unsigned int nx, ny;
    int retval = 0;
    char location[50] = "PGE11_grib2.c, read_grib2_date\0";
    int forecast_hours(unsigned char **sec);

    j = 1;
    inv_out = (char *)malloc(100 * sizeof(char));

    fopen_file(&in_file, filename, "rb");

    rd_grib2_msg_seq_file(sec, &in_file, &pos, &len, &num_submsgs);

    parse_1st_msg(sec);

    getName(sec, 1, NULL, name, NULL, unit);          /* gets name and units of parameters... */
    get_nxny_(sec, &nx, &ny, &nnpnts, &nres, &nscan); /* gets dimensions of parameter */

    center = GB2_Center(sec);
    subcenter = GB2_Subcenter(sec);
    fixed_surfaces(sec, &level_type1, &val1, &undef_val1, &level_type2, &val2, &undef_val2);
    level2(1, level_type1, undef_val1, val1, level_type2, undef_val2, val2, center, subcenter,
           inv_out); /* gets level of parameter */

    // printf("%s -- %s\n", name, inv_out);

    if (((strcmp(name, what) && !strcmp(inv_out, where))) || ((!strcmp(name, what) && strcmp(inv_out, where))) ||
        ((strcmp(name, what) && strcmp(inv_out, where)))) {
        do {
            // printf("%d %s %s %d -- will skip ...\n", j, name, inv_out, pos);
            for (i = 0; i < 100; i++) inv_out[i] = '\0';
            msg = rd_grib2_msg_seq_file(sec, &in_file, &pos, &len, &num_submsgs);
            if (msg == NULL) break;
            j++;

            parse_1st_msg(sec);

            getName(sec, 1, NULL, name, NULL, unit);
            get_nxny_(sec, &nx, &ny, &nnpnts, &nres, &nscan);
            center = GB2_Center(sec);
            subcenter = GB2_Subcenter(sec);
            fixed_surfaces(sec, &level_type1, &val1, &undef_val1, &level_type2, &val2, &undef_val2);
            level2(1, level_type1, undef_val1, val1, level_type2, undef_val2, val2, center, subcenter, inv_out);

        } while ((strcmp(name, what)) || (strcmp(inv_out, where)));
    }

    if ((!strcmp(name, what)) && (!strcmp(inv_out, where))) {
        // printf("FOUND -- %s -- %s\n", name, inv_out);
        // printf("FOUND %d %d %d %d %d '%s' %s -- %s \n", nx, ny, nnpnts, pos, nscan, name, unit, inv_out);

        *frcst_val = forecast_hours(sec);

        get_time(sec[1] + 12, &Y, &M, &D, &H, &MN, &S); /* gets date/time of parameter */

        // printf("%d %d %d %d %d %d\n", Y, M, D, H, MN, S);

        sprintf(date, "%4.4d-%2.2d-%2.2dT%2.2d:%2.2d:%2.2d.000000Z", Y, M, D, H, MN, S); /* ASCII time code A */

    } else {
        /*printf("Parameter '%s' at level '%s' not found\n", what, where);*/
        retval = 1;
    }
    fseek(in_file.cfile, 0L, SEEK_SET);
    fclose_file(&in_file);

    return (retval);
}

int forecast_hours(unsigned char **sec) {
    int fcst_time, fcst_unit;
    unsigned char *verf_time;

    /* if not a forecast .. no code 4.4 */
    if (code_table_4_4_location(sec) == NULL) {
        return 0;
    }

    /* forecast time */
    if ((fcst_unit = code_table_4_4(sec)) != 1) {
        printf("Warning: units of forecast time other than 'hour'\n");
    }
    fcst_time = forecast_time_in_units(sec);
    if ((verf_time = stat_proc_verf_time_location(sec)) != NULL) {
        printf("Warning: verf_time is non-null, contact developer\n");
    }

    return (fcst_time);
}

#ifndef USE_WGRIB2_LIBRARY

/* definition of gribtable */

struct gribtable_s gribtable[] = {

    {0, 0, 0, 255, 7, 1, 7, 193, "4LFTX", "Best (4 layer) Lifted Index", "K"},
    {0, 1, 0, 255, 0, 0, 7, 11, "4LFTX", "Best (4 layer) Lifted Index", "K"},
    {0, 0, 0, 255, 7, 1, 3, 197, "5WAVA", "5-Wave Geopotential Height Anomaly", "gpm"},
    {0, 1, 0, 255, 0, 0, 3, 19, "5WAVA", "5-Wave Geopotential Height Anomaly", "gpm"},
    {0, 0, 0, 255, 7, 1, 3, 193, "5WAVH", "5-Wave Geopotential Height", "gpm"},
    {0, 1, 0, 255, 0, 0, 3, 15, "5WAVH", "5-Wave Geopotential Height", "gpm"},
    {0, 1, 0, 255, 0, 0, 20, 106, "AACOEF", "Aerosol Absorption Coefficient", "1/m"},
    {0, 1, 0, 255, 0, 0, 2, 11, "ABSD", "Absolute Divergence", "1/s"},
    {4, 1, 0, 255, 0, 0, 2, 5, "ABSFRQ", "HF Absorption Frequency", "Hz"},
    {0, 1, 0, 255, 0, 0, 1, 18, "ABSH", "Absolute Humidity", "kg/m^3"},
    {4, 1, 0, 255, 0, 0, 2, 6, "ABSRB", "HF Absorption", "dB"},
    {0, 1, 0, 255, 0, 0, 2, 10, "ABSV", "Absolute Vorticity", "1/s"},
    {0, 1, 0, 255, 0, 0, 18, 0, "ACCES", "Air Concentration of Caesium 137", "Bq/m^3"},
    {0, 1, 0, 255, 0, 0, 18, 1, "ACIOD", "Air Concentration of Iodine 131", "Bq/m^3"},
    {2, 0, 0, 255, 7, 1, 0, 228, "ACOND", "Aerodynamic conductance", "m/s"},
    {0, 1, 0, 255, 0, 0, 1, 10, "ACPCP", "Convective Precipitation", "kg/m^2"},
    {0, 0, 0, 255, 7, 1, 1, 224, "ACPCPN", "Convective precipitation (nearest grid point)", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 18, 2, "ACRADP", "Air Concentration of Radioactive Pollutant", "Bq/m^3"},
    {1, 1, 0, 255, 0, 0, 2, 11, "ACWSR", "Attenuation Coefficient of Water with Respect to Solar Radiation", "1/m"},
    {10, 1, 0, 255, 0, 0, 4, 13, "ACWSRD", "Attenuation Coefficient Of Water With Respect to Solar Radiation", "1/m"},
    {0, 1, 0, 255, 0, 0, 20, 105, "AECOEF", "Aerosol Extinction Coefficient", "1/m"},
    {0, 1, 0, 255, 0, 0, 20, 3, "AEMFLX", "Atmosphere Emission Mass Flux", "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 13, 0, "AEROT", "Aerosol Type", "-"},
    {0, 1, 0, 255, 0, 0, 2, 36, "AFRWE", "Amplitude Function for Rossby Wave Envelope for Meridional Wind", "m/s"},
    {0, 1, 0, 255, 0, 0, 20, 50, "AIA", "Amount in Atmosphere", "mol"},
    {0, 1, 0, 255, 0, 0, 18, 10, "AIRCON", "Air Concentration", "Bq/m^3"},
    {2, 0, 0, 255, 7, 1, 0, 208, "AKHS", "Surface exchange coefficients for T and Q divided by delta z", "m/s"},
    {2, 0, 0, 255, 7, 1, 0, 209, "AKMS", "Surface exchange coefficients for U and V divided by delta z", "m/s"},
    {0, 1, 0, 255, 0, 0, 19, 1, "ALBDO", "Albedo", "%"},
    {0, 1, 0, 255, 0, 0, 20, 108, "ALBGRD", "Aerosol Lidar Backscatter from the Ground", "1/m/sr"},
    {0, 1, 0, 255, 0, 0, 20, 107, "ALBSAT", "Aerosol Lidar Backscatter from Satellite", "1/m/sr"},
    {10, 1, 0, 255, 0, 0, 0, 38, "ALCWH", "Altimeter Corrected Wave Height", "m"},
    {0, 1, 0, 255, 0, 0, 20, 110, "ALEGRD", "Aerosol Lidar Extinction from the Ground", "1/m"},
    {0, 1, 0, 255, 0, 0, 20, 109, "ALESAT", "Aerosol Lidar Extinction from Satellite", "1/m"},
    {10, 1, 0, 255, 0, 0, 0, 39, "ALRRC", "Altimeter Range Relative Correction", "-"},
    {0, 1, 0, 255, 0, 0, 3, 11, "ALTS", "Altimeter Setting", "Pa"},
    {10, 1, 0, 255, 0, 0, 0, 37, "ALTWH", "Altimeter Wave Height", "m"},
    {2, 0, 0, 255, 7, 1, 0, 219, "AMIXL", "Asymptotic mixing length scale", "m"},
    {3, 0, 0, 255, 7, 1, 192, 11, "AMSRE10", "Simulated Brightness Temperature for AMSRE on Aqua, Channel 10", "K"},
    {3, 0, 0, 255, 7, 1, 192, 12, "AMSRE11", "Simulated Brightness Temperature for AMSRE on Aqua, Channel 11", "K"},
    {3, 0, 0, 255, 7, 1, 192, 13, "AMSRE12", "Simulated Brightness Temperature for AMSRE on Aqua, Channel 12", "K"},
    {3, 0, 0, 255, 7, 1, 192, 10, "AMSRE9", "Simulated Brightness Temperature for AMSRE on Aqua, Channel 9", "K"},
    {0, 1, 0, 255, 0, 0, 20, 59, "ANCON", "Aerosol Number Concentration", "1/m^3"},
    {3, 1, 0, 255, 0, 0, 1, 23, "ANGCOE", "Angstrom Coefficient", "-"},
    {0, 1, 0, 255, 0, 0, 20, 111, "ANGSTEXP", "Angstrom Exponent", "Numeric"},
    {0, 1, 0, 255, 0, 0, 20, 5, "ANPEMFLX", "Atmosphere Net Production And Emision Mass Flux", "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 20, 4, "ANPMFLX", "Atmosphere Net Production Mass Flux", "kg/m^2/s"},
    {10, 0, 0, 255, 7, 1, 3, 197, "AOHFLX", "Net Air-Ocean Heat Flux", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 3, 21, "AOSGSO", "Angle of Sub-Grid Scale Orography", "rad"},
    {3, 1, 0, 255, 0, 0, 1, 20, "AOT06", "Aerosol Optical Thickness at 0.635 um", "-"},
    {3, 1, 0, 255, 0, 0, 1, 21, "AOT08", "Aerosol Optical Thickness at 0.810 um", "-"},
    {3, 1, 0, 255, 0, 0, 1, 22, "AOT16", "Aerosol Optical Thickness at 1.640 um", "-"},
    {0, 1, 0, 255, 0, 0, 20, 102, "AOTK", "Aerosol Optical Thickness", "Numeric"},
    {0, 1, 0, 255, 0, 0, 1, 8, "APCP", "Total Precipitation", "kg/m^2"},
    {0, 0, 0, 255, 7, 1, 1, 223, "APCPN", "Total precipitation (nearest grid point)", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 0, 21, "APTMP", "Apparent Temperature", "K"},
    {0, 0, 0, 255, 7, 1, 1, 221, "ARAIN", "Liquid precipitation (Rainfall)", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 20, 8, "AREMFLX", "Atmosphere Re-Emission Mass Flux", "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 3, 24, "ASGSO", "Anisotropy of Sub-Grid Scale Orography", "Numeric"},
    {10, 0, 0, 255, 7, 1, 3, 198, "ASHFL", "Assimilative Heat Flux", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 20, 60, "ASNCON", "Aerosol Specific Number Concentration", "1/kg"},
    {0, 1, 0, 255, 0, 0, 1, 29, "ASNOW", "Total Snowfall", "m"},
    {0, 1, 0, 255, 0, 0, 20, 104, "ASYSFK", "Asymmetry Factor", "Numeric"},
    {0, 1, 0, 255, 0, 0, 190, 0, "ATEXT", "Arbitrary Text String", "CCITTIA5"},
    {3, 1, 0, 255, 0, 0, 1, 13, "ATMDIV", "Atmospheric Divergence", "1/s"},
    {0, 1, 0, 255, 0, 0, 20, 101, "ATMTK", "Vertical Visual Range", "m"},
    {2, 0, 0, 255, 7, 1, 3, 201, "AVSFT", "Average Surface Skin Temperature", "K"},
    {2, 0, 0, 255, 7, 1, 3, 200, "BARET", "Bare Soil Surface Skin temperature", "K"},
    {10, 1, 0, 255, 0, 0, 4, 7, "BATHY", "Bathymetry", "m"},
    {10, 1, 0, 255, 0, 0, 0, 44, "BENINX", "Benjamin-Feir Index", "-"},
    {1, 0, 0, 255, 7, 1, 0, 192, "BGRUN", "Baseflow-Groundwater Runoff", "kg/m^2"},
    {1, 1, 0, 255, 0, 0, 0, 5, "BGRUN", "Baseflow-Groundwater Runoff", "kg/m^2"},
    {10, 0, 0, 255, 7, 1, 4, 194, "BKENG", "Barotropic Kinectic Energy", "J/kg"},
    {0, 1, 0, 255, 0, 0, 7, 1, "BLI", "Best Lifted Index (to 500 hPa)", "K"},
    {0, 1, 0, 255, 0, 0, 7, 16, "BLKRN", "Bulk Richardson Number", "Numeric"},
    {0, 1, 0, 255, 0, 0, 2, 20, "BLYDP", "Boundary Layer Dissipation", "W/m^2"},
    {2, 0, 0, 255, 7, 1, 0, 197, "BMIXL", "Blackadars Mixing Length Scale", "m"},
    {2, 1, 0, 255, 0, 0, 0, 14, "BMIXL", "Blackadars Mixing Length Scale", "m"},
    {0, 0, 0, 255, 7, 1, 7, 201, "BNEGELAY", "Bourgoiun Negative Energy Layer (surface to freezing level)", "J/kg"},
    {2, 1, 0, 255, 0, 0, 3, 4, "BOTLST", "Bottom Layer Soil Temperature", "K"},
    {0, 0, 0, 255, 7, 1, 7, 202, "BPOSELAY", "Bourgoiun Positive Energy Layer (2k ft AGL to 400 hPa)", "J/kg"},
    {0, 1, 0, 255, 0, 0, 15, 1, "BREF", "Base Reflectivity", "dB"},
    {3, 1, 0, 255, 0, 0, 1, 27, "BRFLF", "Bidirectional Reflecance Factor", "Numeric"},
    {0, 1, 0, 255, 0, 0, 5, 7, "BRTEMP", "Brightness Temperature", "K"},
    {0, 1, 0, 255, 0, 0, 4, 4, "BRTMP", "Brightness Temperature", "K"},
    {0, 1, 0, 255, 0, 0, 15, 2, "BRVEL", "Base Radial Velocity", "m/s"},
    {0, 1, 0, 255, 0, 0, 15, 0, "BSWID", "Base Spectrum Width", "m/s"},
    {4, 1, 0, 255, 0, 0, 3, 0, "BTOT", "Magnetic Field Magnitude", "T"},
    {4, 1, 0, 255, 0, 0, 3, 1, "BVEC1", "1st Vector Component of Magnetic Field", "T"},
    {4, 1, 0, 255, 0, 0, 3, 2, "BVEC2", "2nd Vector Component of Magnetic Field", "T"},
    {4, 1, 0, 255, 0, 0, 3, 3, "BVEC3", "3rd Vector Component of Magnetic Field", "T"},
    {0, 1, 0, 255, 0, 0, 18, 18, "CAACL", "Column-Averaged Air Concentration in Layer", "Bq/m^3"},
    {4, 1, 0, 255, 0, 0, 8, 4, "CAIIRAD", "CaII-K Radiance", "W/sr/m^2"},
    {0, 0, 0, 255, 7, 1, 7, 206, "CANGLE", "Critical Angle", "degree"},
    {2, 0, 0, 255, 7, 1, 1, 192, "CANL", "Cold Advisory for Newborn Livestock", "-"},
    {0, 1, 0, 255, 0, 0, 7, 6, "CAPE", "Convective Available Potential Energy", "J/kg"},
    {0, 1, 0, 255, 0, 0, 19, 22, "CAT", "Clear Air Turbulence (CAT)", "%"},
    {0, 1, 0, 255, 0, 0, 1, 88, "CATCP", "Categorical Convective Precipitation", "-"},
    {0, 1, 0, 255, 0, 0, 19, 29, "CATEDR", "Clear Air Turbulence (CAT) (Eddy Dissipation Rate)", "m2/3s-1"},
    {0, 1, 0, 255, 0, 0, 20, 63, "CAVEMDL", "Column-Averaged Mass Density in Layer", "kg/m^3"},
    {0, 1, 0, 255, 0, 0, 20, 70, "CBECSLSP",
     "Column-integrated below-cloud scavenging rate by large-scale precipitation", "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 6, 25, "CBHE", "Horizontal Extent of Cumulonimbus (CB)", "%"},
    {0, 1, 0, 255, 0, 0, 20, 67, "CBLCLDSP", "Column-integrated below-cloud scavenging rate by precipitation",
     "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 20, 73, "CBLCSRCP",
     "Column-integrated below-cloud scavenging rate by convective precipitation", "kg/m^2/s"},
    {3, 1, 0, 255, 0, 0, 1, 14, "CBTMP", "Cloudy Brightness Temperature", "K"},
    {3, 1, 0, 255, 0, 0, 1, 98, "CCMPEMRR",
     "Correlation coefficient between MPE rain rates for the co-located IR data and the microwave data rain rates",
     "Numeric"},
    {2, 0, 0, 255, 7, 1, 0, 199, "CCOND", "Canopy Conductance", "m/s"},
    {2, 1, 0, 255, 0, 0, 0, 15, "CCOND", "Canopy Conductance", "m/s"},
    {0, 0, 0, 255, 7, 1, 2, 196, "CD", "Drag Coefficient", "non-dim"},
    {0, 1, 0, 255, 0, 0, 2, 29, "CD", "Drag Coefficient", "Numeric"},
    {0, 1, 0, 255, 0, 0, 6, 7, "CDCA", "Cloud Amount", "%"},
    {0, 1, 0, 255, 0, 0, 6, 11, "CDCB", "Cloud Base", "m"},
    {0, 1, 0, 255, 0, 0, 6, 22, "CDCC", "Cloud Cover", "%"},
    {0, 1, 0, 255, 0, 0, 17, 3, "CDCDLTFD", "Cloud-to-Cloud Lightning Flash Density", "km-2day-1"},
    {0, 1, 0, 255, 0, 0, 6, 23, "CDCIMR", "Cloud Ice Mixing Ratio", "kg/kg"},
    {0, 1, 0, 255, 0, 0, 6, 2, "CDCON", "Convective Cloud Cover", "%"},
    {0, 1, 0, 255, 0, 0, 6, 8, "CDCT", "Cloud Type", "-"},
    {0, 1, 0, 255, 0, 0, 6, 12, "CDCTOP", "Cloud Top", "m"},
    {0, 1, 0, 255, 0, 0, 17, 2, "CDGDLTFD", "Cloud-to-Ground Lightning Flash Density", "km-2day-1"},
    {0, 0, 0, 255, 7, 1, 6, 192, "CDLYR", "Non-Convective Cloud Cover", "%"},
    {0, 1, 0, 255, 0, 0, 6, 14, "CDLYR", "Non-Convective Cloud Cover", "%"},
    {0, 0, 0, 255, 7, 1, 4, 195, "CDUVB", "Clear sky UV-B Downward Solar Flux", "W/m^2"},
    {10, 1, 0, 255, 0, 0, 0, 16, "CDWW", "Coefficient of Drag With Waves", "-"},
    {0, 1, 0, 255, 0, 0, 6, 13, "CEIL", "Ceiling", "m"},
    {0, 0, 0, 255, 7, 1, 5, 197, "CFNLF", "Cloud Forcing Net Long Wave Flux", "W/m^2"},
    {0, 0, 0, 255, 7, 1, 4, 199, "CFNSF", "Cloud Forcing Net Solar Flux", "W/m^2"},
    {0, 0, 0, 255, 7, 1, 1, 193, "CFRZR", "Categorical Freezing Rain", "-"},
    {0, 1, 0, 255, 0, 0, 1, 34, "CFRZR", "Categorical Freezing Rain", "-"},
    {0, 1, 0, 255, 0, 0, 20, 54, "CGDRC", "Chemical Gross Destruction Rate of Concentration", "mol/m^3/s"},
    {0, 1, 0, 255, 0, 0, 20, 53, "CGPRC", "Chemical Gross Production Rate of Concentration", "mol/m^3/s"},
    {10, 1, 0, 255, 0, 0, 3, 2, "CH", "Heat Exchange Coefficient", "-"},
    {0, 1, 0, 255, 0, 0, 18, 17, "CIAIRC", "Column-Integrated Air Concentration", "Bq/m^2"},
    {0, 1, 0, 255, 0, 0, 6, 0, "CICE", "Cloud Ice", "kg/m^2"},
    {0, 0, 0, 255, 7, 1, 19, 206, "CICEL", "Confidence - Ceiling", "-"},
    {0, 0, 0, 255, 7, 1, 1, 194, "CICEP", "Categorical Ice Pellets", "-"},
    {0, 1, 0, 255, 0, 0, 1, 35, "CICEP", "Categorical Ice Pellets", "-"},
    {10, 1, 0, 255, 0, 0, 2, 12, "CICES", "Compressive Ice Strength", "N/m"},
    {0, 0, 0, 255, 7, 1, 19, 208, "CIFLT", "Confidence - Flight Category", "-"},
    {0, 1, 0, 255, 0, 0, 1, 82, "CIMIXR", "Cloud Ice Mixing Ratio", "kg/kg"},
    {0, 1, 0, 255, 0, 0, 7, 7, "CIN", "Convective Inhibition", "J/kg"},
    {0, 1, 0, 255, 0, 0, 20, 66, "CINCLDSP", "Column-integrated in-cloud scavenging rate by precipitation", "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 20, 69, "CINCSLSP", "Column-integrated in-cloud scavenging rate by large-scale precipitation",
     "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 20, 72, "CINCSRCP", "Column-integrated in-cloud scavenging rate by convective precipitation",
     "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 20, 68, "CIRELREP", "Column-integrated release rate from evaporating precipitation",
     "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 20, 74, "CIRERECP", "Column-integrated release rate from evaporating convective precipitation",
     "kg/m^2/s"},
    {2, 1, 0, 255, 0, 0, 3, 22, "CISICE", "Column-Integrated Soil Ice", "kg/m^2"},
    {2, 1, 0, 255, 0, 0, 3, 20, "CISOILM", "Column-Integrated Soil Moisture", "kg/m^2"},
    {2, 1, 0, 255, 0, 0, 0, 23, "CISOILW", "Column-Integrated Soil Water", "kg/m^2"},
    {0, 0, 0, 255, 7, 1, 19, 207, "CIVIS", "Confidence - Visibility", "-"},
    {3, 1, 0, 255, 0, 0, 2, 9, "CLDALB", "Cloud Albedo", "Numeric"},
    {3, 1, 0, 255, 0, 0, 2, 10, "CLDEMISS", "Cloud Emissivity", "Numeric"},
    {3, 1, 0, 255, 0, 0, 2, 8, "CLDIWP", "Cloud Ice Water Path", "kg/m^2"},
    {3, 1, 0, 255, 0, 0, 2, 7, "CLDLWP", "Cloud Liquid Water Path", "kg/m^2"},
    {3, 1, 0, 255, 0, 0, 2, 5, "CLDODEP", "Cloud Optical Depth", "Numeric"},
    {3, 1, 0, 255, 0, 0, 2, 6, "CLDPER", "Cloud Particle Effective Radius", "m"},
    {3, 1, 0, 255, 0, 0, 2, 4, "CLDPHAS", "Cloud Phase", "-"},
    {3, 1, 0, 255, 0, 0, 1, 16, "CLDRAD", "Cloudy Radiance (with respect to wave number)", "W/m/sr"},
    {3, 1, 0, 255, 0, 0, 2, 3, "CLDTYPE", "Cloud Type", "-"},
    {0, 0, 0, 255, 7, 1, 1, 235, "CLLMR", "Cloud Liquid Mixing Ratio", "kg/kg"},
    {0, 1, 0, 255, 0, 0, 1, 22, "CLMR", "Cloud Mixing Ratio", "kg/kg"},
    {3, 1, 0, 255, 0, 0, 0, 7, "CLOUDM", "Cloud Mask", "-"},
    {0, 0, 0, 255, 7, 1, 2, 216, "CNGWDU", "Convective Gravity wave drag zonal acceleration", "m/s^2"},
    {0, 0, 0, 255, 7, 1, 2, 217, "CNGWDV", "Convective Gravity wave drag meridional acceleration", "m/s^2"},
    {0, 0, 0, 255, 7, 1, 3, 209, "CNVDEMF", "Convective detrainment mass flux", "kg/m^2/s"},
    {0, 0, 0, 255, 7, 1, 3, 208, "CNVDMF", "Convective downdraft mass flux", "kg/m^2/s"},
    {0, 0, 0, 255, 7, 1, 0, 196, "CNVHR", "Deep Convective Heating Rate", "K/s"},
    {0, 0, 0, 255, 7, 1, 1, 213, "CNVMR", "Deep Convective Moistening Rate", "kg/kg/s"},
    {0, 0, 0, 255, 7, 1, 2, 212, "CNVU", "Convective zonal momentum mixing acceleration", "m/s^2"},
    {0, 0, 0, 255, 7, 1, 3, 207, "CNVUMF", "Convective updraft mass flux", "kg/m^2/s"},
    {0, 0, 0, 255, 7, 1, 2, 213, "CNVV", "Convective meridional momentum mixing acceleration", "m/s^2"},
    {2, 0, 0, 255, 7, 1, 0, 196, "CNWAT", "Plant Canopy Surface Water", "kg/m^2"},
    {2, 1, 0, 255, 0, 0, 0, 13, "CNWAT", "Plant Canopy Surface Water", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 20, 56, "COAIA", "Changes Of Amount in Atmosphere", "mol/s"},
    {0, 1, 0, 255, 0, 0, 20, 1, "COLMD", "Column-Integrated Mass Density", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 20, 51, "CONAIR", "Concentration In Air", "mol/m^3"},
    {0, 1, 0, 255, 0, 0, 7, 19, "CONAPES", "Convective Available Potential Energy Shear", "m^2/s^2"},
    {0, 0, 0, 255, 7, 1, 1, 216, "CONDP", "Condensation Pressure of Parcali Lifted From Indicate Surface", "Pa"},
    {0, 1, 0, 255, 0, 0, 19, 16, "CONTB", "Contrail Base", "m"},
    {0, 1, 0, 255, 0, 0, 19, 14, "CONTET", "Contrail Engine Type", "-"},
    {0, 1, 0, 255, 0, 0, 19, 13, "CONTI", "Contrail Intensity", "-"},
    {0, 1, 0, 255, 0, 0, 19, 24, "CONTKE", "Convective Turbulent Kinetic Energy", "J/kg"},
    {0, 1, 0, 255, 0, 0, 19, 15, "CONTT", "Contrail Top", "m"},
    {0, 1, 0, 255, 0, 0, 19, 26, "CONVO", "Convective Outlook", "-"},
    {0, 0, 0, 255, 7, 1, 19, 222, "CONVP", "Convection Potential", "-"},
    {0, 0, 0, 255, 7, 1, 192, 6, "COVMM",
     "Covariance between meridional and meridional components of the wind. Defined as [vv]-[v][v], where [] indicates "
     "the mean over the indicated time span.",
     "m^2/s^2"},
    {0, 0, 0, 255, 7, 1, 192, 1, "COVMZ",
     "Covariance between zonal and meridional components of the wind. Defined as [uv]-[u][v], where [] indicates the "
     "mean over the indicated time span.",
     "m^2/s^2"},
    {0, 0, 0, 255, 7, 1, 2, 205, "COVMZ", "Covariance between Meridional and Zonal Components of the wind.", "m^2/s^2"},
    {0, 0, 0, 255, 7, 1, 192, 11, "COVPSPS",
     "Covariance between surface pressure and surface pressure. Defined as [Psfc]-[Psfc][Psfc], where [] indicates the "
     "mean over the indicated time span.",
     "Pa*Pa"},
    {0, 0, 0, 255, 7, 1, 192, 8, "COVQM",
     "Covariance between specific humidity and meridional components of the wind. Defined as [vq]-[v][q], where [] "
     "indicates the mean over the indicated time span.",
     "kg/kg*m/s"},
    {0, 0, 0, 255, 7, 1, 192, 12, "COVQQ",
     "Covariance between specific humidity and specific humidy. Defined as [qq]-[q][q], where [] indicates the mean "
     "over the indicated time span.",
     "kg/kg*kg/kg"},
    {0, 0, 0, 255, 7, 1, 192, 10, "COVQVV",
     "Covariance between specific humidity and vertical components of the wind. Defined as [Omegaq]-[Omega][q], where "
     "[] indicates the mean over the indicated time span.",
     "kg/kg*Pa/s"},
    {0, 0, 0, 255, 7, 1, 192, 7, "COVQZ",
     "Covariance between specific humidity and zonal components of the wind. Defined as [uq]-[u][q], where [] "
     "indicates the mean over the indicated time span.",
     "kg/kg*m/s"},
    {0, 0, 0, 255, 7, 1, 192, 3, "COVTM",
     "Covariance between meridional component of the wind and temperature. Defined as [vT]-[v][T], where [] indicates "
     "the mean over the indicated time span.",
     "K*m/s"},
    {0, 0, 0, 255, 7, 1, 2, 207, "COVTM", "Covariance between Temperature and Meridional Components of the wind.",
     "K*m/s"},
    {0, 0, 0, 255, 7, 1, 192, 14, "COVTT",
     "Covariance between temperature and temperature. Defined as [TT]-[T][T], where [] indicates the mean over the "
     "indicated time span.",
     "K*K"},
    {0, 0, 0, 255, 7, 1, 192, 9, "COVTVV",
     "Covariance between temperature and vertical components of the wind. Defined as [OmegaT]-[Omega][T], where [] "
     "indicates the mean over the indicated time span.",
     "K*Pa/s"},
    {0, 0, 0, 255, 7, 1, 192, 4, "COVTW",
     "Covariance between temperature and vertical component of the wind. Defined as [wT]-[w][T], where [] indicates "
     "the mean over the indicated time span.",
     "K*m/s"},
    {0, 0, 0, 255, 7, 1, 192, 2, "COVTZ",
     "Covariance between zonal component of the wind and temperature. Defined as [uT]-[u][T], where [] indicates the "
     "mean over the indicated time span.",
     "K*m/s"},
    {0, 0, 0, 255, 7, 1, 2, 206, "COVTZ", "Covariance between Temperature and Zonal Components of the wind.", "K*m/s"},
    {0, 0, 0, 255, 7, 1, 192, 13, "COVVVVV",
     "Covariance between vertical and vertical components of the wind. Defined as [OmegaOmega]-[Omega][Omega], where "
     "[] indicates the mean over the indicated time span.",
     "Pa^2/s^2"},
    {0, 0, 0, 255, 7, 1, 192, 5, "COVZZ",
     "Covariance between zonal and zonal components of the wind. Defined as [uu]-[u][u], where [] indicates the mean "
     "over the indicated time span.",
     "m^2/s^2"},
    {0, 1, 0, 255, 0, 0, 1, 39, "CPOFP", "Percent frozen precipitation", "%"},
    {1, 0, 0, 255, 7, 1, 1, 193, "CPOFP", "Percent of Frozen Precipitation", "%"},
    {1, 0, 0, 255, 7, 1, 1, 192, "CPOZP", "Probability of Freezing Precipitation", "%"},
    {1, 1, 0, 255, 0, 0, 1, 0, "CPPOP",
     "Conditional percent precipitation amount fractile for an overall period (encoded as an accumulation)", "kg/m^2"},
    {0, 0, 0, 255, 7, 1, 1, 196, "CPRAT", "Convective Precipitation Rate", "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 1, 37, "CPRAT", "Convective Precipitation Rate", "kg/m^2/s"},
    {0, 0, 0, 255, 7, 1, 1, 192, "CRAIN", "Categorical Rain", "-"},
    {0, 1, 0, 255, 0, 0, 1, 33, "CRAIN", "Categorical Rain", "-"},
    {0, 1, 0, 255, 0, 0, 20, 71, "CRERELSP",
     "Column-integrated release rate from evaporating large-scale precipitation", "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 1, 76, "CRRATE", "Convective Rain Rate", "kg/m^2/s"},
    {4, 1, 0, 255, 0, 0, 2, 9, "CRTFRQ", "Critical Frequency", "Hz"},
    {1, 1, 0, 255, 0, 0, 2, 13, "CSAFC", "Cross Sectional Area of Flow in Channel", "m2"},
    {3, 1, 0, 255, 0, 0, 1, 15, "CSBTMP", "Clear Sky Brightness Temperature", "K"},
    {0, 0, 0, 255, 7, 1, 5, 196, "CSDLF", "Clear Sky Downward Long Wave Flux", "W/m^2"},
    {0, 0, 0, 255, 7, 1, 4, 196, "CSDSF", "Clear Sky Downward Solar Flux", "W/m^2"},
    {3, 1, 0, 255, 0, 0, 2, 0, "CSKPROB", "Clear Sky Probability", "%"},
    {3, 1, 0, 255, 0, 0, 1, 17, "CSKYRAD", "Clear Sky Radiance (with respect to wave number)", "W/m/sr"},
    {0, 0, 0, 255, 7, 1, 1, 195, "CSNOW", "Categorical Snow", "-"},
    {0, 1, 0, 255, 0, 0, 1, 36, "CSNOW", "Categorical Snow", "-"},
    {0, 1, 0, 255, 0, 0, 1, 58, "CSRATE", "Convective Snowfall Rate", "m/s"},
    {0, 1, 0, 255, 0, 0, 1, 55, "CSRWE", "Convective Snowfall Rate Water Equivalent", "kg/m^2/s"},
    {0, 0, 0, 255, 7, 1, 5, 195, "CSULF", "Clear Sky Upward Long Wave Flux", "W/m^2"},
    {0, 0, 0, 255, 7, 1, 4, 198, "CSUSF", "Clear Sky Upward Solar Flux", "W/m^2"},
    {3, 1, 0, 255, 0, 0, 1, 2, "CTOPH", "Cloud Top Height", "m"},
    {3, 1, 0, 255, 0, 0, 1, 3, "CTOPHQI", "Cloud Top Height Quality Indicator", "-"},
    {3, 1, 0, 255, 0, 0, 2, 2, "CTOPRES", "Cloud Top Pressure", "Pa"},
    {3, 1, 0, 255, 0, 0, 2, 1, "CTOPTMP", "Cloud Top Temperature", "K"},
    {0, 1, 0, 255, 0, 0, 19, 21, "CTP", "In-Cloud Turbulence", "%"},
    {0, 0, 0, 255, 7, 1, 6, 194, "CUEFI", "Convective Cloud Efficiency", "non-dim"},
    {0, 1, 0, 255, 0, 0, 6, 16, "CUEFI", "Convective Cloud Efficiency", "Proportion"},
    {0, 1, 0, 255, 0, 0, 6, 6, "CWAT", "Cloud Water", "kg/m^2"},
    {0, 0, 0, 255, 7, 1, 7, 195, "CWDI", "Convective Weather Detection Index", "-"},
    {0, 0, 0, 255, 7, 1, 6, 193, "CWORK", "Cloud Work Function", "J/kg"},
    {0, 1, 0, 255, 0, 0, 6, 15, "CWORK", "Cloud Work Function", "J/kg"},
    {0, 1, 0, 255, 0, 0, 1, 48, "CWP", "Convective Water Precipitation", "kg/m^2"},
    {1, 0, 0, 255, 7, 1, 1, 195, "CWR", "Probability of Wetting Rain, exceeding in 0.10 in a given time period", "%"},
    {10, 0, 0, 255, 7, 1, 4, 195, "DBSS", "Geometric Depth Below Sea Surface", "m"},
    {0, 0, 0, 255, 7, 1, 7, 203, "DCAPE", "Downdraft CAPE", "J/kg"},
    {0, 1, 0, 255, 0, 0, 20, 12, "DDMFLX", "Dry Deposition Mass Flux", "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 3, 30, "DDRATE", "Downdraught Detrainment Rate", "kg/m^3/s"},
    {0, 1, 0, 255, 0, 0, 20, 15, "DDVEL", "Dry deposition velocity", "m/s"},
    {2, 1, 0, 255, 0, 0, 0, 30, "DECFC", "Deciduous Forest Cover", "Proportion"},
    {0, 1, 0, 255, 0, 0, 3, 10, "DEN", "Density", "kg/m^3"},
    {0, 1, 0, 255, 0, 0, 3, 14, "DENALT", "Density Altitude", "m"},
    {0, 1, 0, 255, 0, 0, 0, 7, "DEPR", "Dew Point Depression (or Deficit)", "K"},
    {1, 1, 0, 255, 0, 0, 0, 13, "DEPWSS", "Depth of Water on Soil Surface", "kg/m^2"},
    {10, 1, 0, 255, 0, 0, 2, 2, "DICED", "Direction of Ice Drift", "deg"},
    {4, 1, 0, 255, 0, 0, 4, 2, "DIFEFLUX", "Electron Flux (Differential)", "1/(m^2s*sr*eV)"},
    {4, 1, 0, 255, 0, 0, 4, 4, "DIFIFLUX", "Heavy Ion Flux (Differential)", "1/(m^2s*sr*eV/nuc)"},
    {4, 1, 0, 255, 0, 0, 4, 0, "DIFPFLUX", "Proton Flux (Differential)", "1/(m^2s*sr*eV)"},
    {3, 1, 0, 255, 0, 0, 6, 5, "DIFSOLEX", "Diffuse Solar Exposure", "J/m^2"},
    {3, 1, 0, 255, 0, 0, 6, 4, "DIFSOLIR", "Diffuse Solar Irradiance", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 4, 14, "DIFSWRF", "Diffuse Short Wave Radiation Flux", "W/m^2"},
    {10, 1, 0, 255, 0, 0, 1, 0, "DIRC", "Current Direction", "deg"},
    {2, 1, 0, 255, 0, 0, 3, 14, "DIREC", "Direct Evaporation Cease(Soil Moisture)", "kg/m^3"},
    {10, 1, 0, 255, 0, 0, 0, 10, "DIRPW", "Primary Wave Direction", "deg"},
    {3, 1, 0, 255, 0, 0, 6, 3, "DIRSOLEX", "Direct Solar Exposure", "J/m^2"},
    {3, 1, 0, 255, 0, 0, 6, 2, "DIRSOLIR", "Direct Solar Irradiance", "W/m^2"},
    {10, 1, 0, 255, 0, 0, 0, 12, "DIRSW", "Secondary Wave Direction", "deg"},
    {10, 1, 0, 255, 0, 0, 0, 33, "DIRWTS", "Directional Width of The Total Swell", "-"},
    {10, 1, 0, 255, 0, 0, 0, 32, "DIRWWW", "Directional Width of The Wind Waves", "-"},
    {1, 1, 0, 255, 0, 0, 0, 7, "DISRS", "Discharge from Rivers or Streams", "m^3/s"},
    {0, 1, 0, 255, 0, 0, 3, 6, "DIST", "Geometric Height", "m"},
    {0, 0, 0, 255, 7, 1, 5, 192, "DLWRF", "Downward Long-Wave Rad. Flux", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 5, 3, "DLWRF", "Downward Long-Wave Rad. Flux", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 5, 8, "DLWRFCS", "Downward Long-Wave Radiation Flux, Clear Sky", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 3, 28, "DMFLX", "Downdraught Mass Flux", "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 0, 6, "DPT", "Dew Point Temperature", "K"},
    {2, 1, 0, 255, 0, 0, 4, 8, "DRTCODE", "Drought Code (Canadian Forest Service)", "Numeric"},
    {0, 1, 0, 255, 0, 0, 18, 12, "DRYDEP", "Dry Deposition", "Bq/m^2"},
    {0, 0, 0, 255, 7, 1, 19, 237, "DRYTPROB", "Dry Thunderstorm Probability", "%"},
    {4, 1, 0, 255, 0, 0, 7, 2, "DSKDAY", "Disk Intensity Day", "1/m^2/s"},
    {4, 1, 0, 255, 0, 0, 7, 1, "DSKINT", "Disk Intensity", "1/m^2/s"},
    {4, 1, 0, 255, 0, 0, 7, 3, "DSKNGT", "Disk Intensity Night", "1/m^2/s"},
    {10, 1, 0, 255, 0, 0, 3, 1, "DSLM", "Deviation of Sea Level from Mean", "m"},
    {0, 1, 0, 255, 0, 0, 191, 3, "DSLOBS", "Days Since Last Observation", "day"},
    {10, 1, 0, 255, 0, 0, 191, 3, "DSLOBSO", "Days Since Last Observation", "day"},
    {0, 0, 0, 255, 7, 1, 4, 192, "DSWRF", "Downward Short-Wave Radiation Flux", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 4, 7, "DSWRF", "Downward Short-Wave Radiation Flux", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 4, 52, "DSWRFCS", "Downward Short-Wave Radiation Flux, Clear Sky", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 4, 13, "DSWRFLX", "Direct Short Wave Radiation Flux", "W/m^2"},
    {0, 0, 0, 255, 7, 1, 4, 204, "DTRF", "Downward Total Radiation Flux", "W/m^2"},
    {2, 1, 0, 255, 0, 0, 4, 7, "DUFMCODE", "Duff Moisture Code (Canadian Forest Service)", "Numeric"},
    {0, 0, 0, 255, 7, 1, 4, 194, "DUVB", "UV-B Downward Solar Flux", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 4, 12, "DWUVR", "Downward UV Radiation", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 2, 9, "DZDT", "Vertical Velocity (Geometric)", "m/s"},
    {3, 1, 0, 255, 0, 0, 2, 11, "EAODR", "Effective Absorption Optical Depth Ratio", "Numeric"},
    {3, 1, 0, 255, 0, 0, 5, 5, "EBSDSSTS",
     "Estimated bias Standard Deviation between Sea-Surface Temperature and Standard", "K"},
    {3, 1, 0, 255, 0, 0, 5, 4, "EBSSTSTD", "Estimated bias between Sea-Surface Temperature and Standard", "K"},
    {0, 1, 0, 255, 0, 0, 19, 30, "EDPARM", "Eddy Dissipation Parameter", "m2/3s-1"},
    {0, 1, 0, 255, 0, 0, 1, 138, "EFARCICE", "Effective Aspect Ratio of Cloud Ice", "-"},
    {0, 1, 0, 255, 0, 0, 1, 140, "EFARGRL", "Effective Aspect Ratio of Graupel", "-"},
    {0, 1, 0, 255, 0, 0, 1, 141, "EFARHAIL", "Effective Aspect Ratio of Hail", "-"},
    {0, 1, 0, 255, 0, 0, 1, 137, "EFARRAIN", "Effective Aspect Ratio of Rain", "-"},
    {0, 1, 0, 255, 0, 0, 1, 142, "EFARSIC", "Effective Aspect Ratio of Subgrid Ice Clouds", "-"},
    {0, 1, 0, 255, 0, 0, 1, 139, "EFARSNOW", "Effective Aspect Ratio of Snow", "-"},
    {0, 0, 0, 255, 7, 1, 7, 204, "EFHL", "Effective Storm Relative Helicity", "m^2/s^2"},
    {0, 1, 0, 255, 0, 0, 1, 131, "EFRCICE", "Effective Radius of Cloud Ice", "m"},
    {0, 1, 0, 255, 0, 0, 1, 129, "EFRCWAT", "Effective Radius of Cloud Water", "m"},
    {0, 1, 0, 255, 0, 0, 1, 133, "EFRGRL", "Effective Radius of Graupel", "m"},
    {0, 1, 0, 255, 0, 0, 1, 134, "EFRHAIL", "Effective Radius of Hail", "m"},
    {0, 1, 0, 255, 0, 0, 1, 130, "EFRRAIN", "Effective Radius of Rain", "m"},
    {0, 1, 0, 255, 0, 0, 1, 136, "EFRSICEC", "Effective Radius of Subgrid Ice Clouds", "m"},
    {0, 1, 0, 255, 0, 0, 1, 135, "EFRSLC", "Effective Radius of Subgrid Liquid Clouds", "m"},
    {0, 1, 0, 255, 0, 0, 1, 132, "EFRSNOW", "Effective Radius of Snow", "m"},
    {0, 1, 0, 255, 0, 0, 7, 9, "EHLX", "Energy Helicity Index", "Numeric"},
    {4, 1, 0, 255, 0, 0, 2, 1, "ELCDEN", "Electron Density", "1/m^3"},
    {4, 1, 0, 255, 0, 0, 0, 1, "ELECTMP", "Electron Temperature", "K"},
    {10, 0, 0, 255, 7, 1, 3, 194, "ELEV", "Ocean Surface Elevation Relative to Geoid", "m"},
    {0, 0, 0, 255, 7, 1, 19, 238, "ELLINX", "Ellrod Index", "-"},
    {0, 0, 0, 255, 7, 1, 191, 193, "ELON", "East Longitude (0 to 360)", "deg"},
    {0, 0, 0, 255, 7, 1, 191, 197, "ELONN", "East Longitude (nearest neighbor) (0 to 360)", "deg"},
    {0, 1, 0, 255, 0, 0, 20, 76, "EMISFLX", "Emission Rate", "kg/kg/s"},
    {0, 0, 0, 255, 7, 1, 1, 211, "EMNP", "Evaporation - Precipitation", "cm/day"},
    {0, 1, 0, 255, 0, 0, 0, 3, "EPOT", "Pseudo-Adiabatic Potential Temperature (or Equivalent Potential Temperature)",
     "K"},
    {0, 0, 0, 255, 7, 1, 19, 218, "EPSR", "Radiative emissivity", "-"},
    {10, 0, 0, 255, 7, 1, 3, 252, "EROSNP", "Erosion Occurrence Probability", "%"},
    {1, 1, 0, 255, 0, 0, 0, 3, "ESCT", "Elevation of Snow Covered Terrain", "-"},
    {0, 0, 0, 255, 7, 1, 7, 205, "ESP", "Enhanced Stretching Potential", "Numeric"},
    {3, 1, 0, 255, 0, 0, 1, 0, "ESTP", "Estimated Precipitation", "kg/m^2"},
    {3, 1, 0, 255, 0, 0, 1, 4, "ESTUGRD", "Estimated u-Component of Wind", "m/s"},
    {3, 1, 0, 255, 0, 0, 1, 5, "ESTVGRD", "Estimated v-Component of Wind", "m/s"},
    {0, 1, 0, 255, 0, 0, 2, 32, "ETACVV", "Eta Coordinate Vertical Velocity", "1/s"},
    {10, 0, 0, 255, 7, 1, 3, 250, "ETCWL", "Extra Tropical Storm Surge Combined Surge and Tide", "m"},
    {4, 1, 0, 255, 0, 0, 3, 4, "ETOT", "Electric Field Magnitude", "V/m"},
    {10, 0, 0, 255, 7, 1, 3, 193, "ETSRG", "Extra Tropical Storm Surge", "m"},
    {0, 1, 0, 255, 0, 0, 2, 38, "ETSS", "Eastward Turbulent Surface Stress", "N/m^2*s"},
    {4, 1, 0, 255, 0, 0, 6, 3, "EUVIRR", "Solar EUV Irradiance", "W/m^2"},
    {4, 1, 0, 255, 0, 0, 8, 1, "EUVRAD", "EUV Radiance", "W/sr/m^2"},
    {2, 1, 0, 255, 0, 0, 0, 6, "EVAPT", "Evapotranspiration", "1/kg^2/s"},
    {2, 1, 0, 255, 0, 0, 0, 39, "EVAPTRAT", "Evapotranspiration Rate", "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 1, 79, "EVARATE", "Evaporation Rate", "kg/m^2/s"},
    {2, 0, 0, 255, 7, 1, 3, 198, "EVBS", "Direct Evaporation from Bare Soil", "W/m^2"},
    {2, 0, 0, 255, 7, 1, 0, 229, "EVCW", "Canopy water evaporation", "W/m^2"},
    {4, 1, 0, 255, 0, 0, 3, 5, "EVEC1", "1st Vector Component of Electric Field", "V/m"},
    {4, 1, 0, 255, 0, 0, 3, 6, "EVEC2", "2nd Vector Component of Electric Field", "V/m"},
    {4, 1, 0, 255, 0, 0, 3, 7, "EVEC3", "3rd Vector Component of Electric Field", "V/m"},
    {2, 1, 0, 255, 0, 0, 0, 29, "EVGFC", "Evergreen Forest Cover", "Proportion"},
    {0, 1, 0, 255, 0, 0, 1, 6, "EVP", "Evaporation", "kg/m^2"},
    {2, 0, 0, 255, 7, 1, 0, 213, "EWATR", "Open water evaporation (standing water)", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 2, 39, "EWTPARM", "Eastward Wind Tendency Due to Parameterizations", "m/s^2"},
    {0, 1, 0, 255, 0, 0, 3, 26, "EXPRES", "Exner Pressure", "Numeric"},
    {4, 1, 0, 255, 0, 0, 6, 5, "F107", "F10.7", "W/m^2/Hz"},
    {2, 1, 0, 255, 0, 0, 4, 3, "FBAREA", "Fire Burned Area", "%"},
    {2, 1, 0, 255, 0, 0, 4, 10, "FBUPINX", "Fire Build Up Index (Canadian Forest Service)", "Numeric"},
    {0, 1, 0, 255, 0, 0, 6, 37, "FCONPC", "Fraction of Convective Precipitation Cover", "Proportion"},
    {3, 1, 0, 255, 0, 0, 5, 3, "FDNSSTMP", "Foundation Sea-Surface Temperature", "K"},
    {2, 1, 0, 255, 0, 0, 4, 11, "FDSRTE", "Fire Daily Severity Rating (Canadian Forest Service)", "Numeric"},
    {1, 1, 0, 255, 0, 0, 0, 0, "FFLDG",
     "Flash Flood Guidance (Encoded as an accumulation over a floating subinterval of time between the reference time "
     "and valid time)",
     "kg/m^2"},
    {1, 1, 0, 255, 0, 0, 0, 1, "FFLDRO",
     "Flash Flood Runoff (Encoded as an accumulation over a floating subinterval of time)", "kg/m^2"},
    {2, 1, 0, 255, 0, 0, 4, 6, "FFMCODE", "Fine Fuel Moisture Code (Canadian Forest Service)", "Numeric"},
    {0, 0, 0, 255, 7, 1, 6, 199, "FICE", "Ice fraction of total condensate", "non-dim"},
    {0, 1, 0, 255, 0, 0, 6, 21, "FICE", "Ice fraction of total condensate", "Proportion"},
    {0, 0, 0, 255, 7, 1, 1, 228, "FICEAC", "Flat Ice Accumulation (FRAM)", "kg/m^2"},
    {3, 1, 0, 255, 0, 0, 0, 9, "FIREDI", "Fire Detection Indicator", "-"},
    {2, 1, 0, 255, 0, 0, 4, 1, "FIREODT", "Fire Outlook Due to Dry Thunderstorm", "-"},
    {2, 1, 0, 255, 0, 0, 4, 0, "FIREOLK", "Fire Outlook", "-"},
    {2, 0, 0, 255, 7, 1, 3, 203, "FLDCP", "Field Capacity", "Fraction"},
    {1, 1, 0, 255, 0, 0, 0, 12, "FLDPSW", "Flood Plain Storage of Water", "m3"},
    {0, 0, 0, 255, 7, 1, 19, 205, "FLGHT", "Flight Category", "-"},
    {0, 1, 0, 255, 0, 0, 7, 18, "FLXRN", "Flux Richardson Number", "Numeric"},
    {2, 1, 0, 255, 0, 0, 4, 4, "FOSINDX", "Fosberg Index", "Numeric"},
    {0, 1, 0, 255, 0, 0, 1, 67, "FPRATE", "Freezing Rain Precipitation Rate", "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 6, 32, "FRACCC", "Fraction of Cloud Cover", "Numeric"},
    {0, 0, 0, 255, 7, 1, 1, 202, "FRAIN", "Rain Fraction of Total Liquid Water", "non-dim"},
    {0, 1, 0, 255, 0, 0, 1, 43, "FRAIN", "Rain Fraction of Total Cloud Water", "Proportion"},
    {10, 1, 0, 255, 0, 0, 0, 64, "FREWTSW", "Frequency width of total swell", "-"},
    {10, 1, 0, 255, 0, 0, 0, 63, "FREWWW", "Frequency width of wind waves", "-"},
    {0, 0, 0, 255, 7, 1, 2, 197, "FRICV", "Frictional Velocity", "m/s"},
    {0, 1, 0, 255, 0, 0, 2, 30, "FRICV", "Frictional Velocity", "m/s"},
    {10, 1, 0, 255, 0, 0, 0, 17, "FRICVW", "Friction Velocity", "m/s"},
    {0, 0, 0, 255, 7, 1, 1, 227, "FROZR", "Frozen Rain", "kg/m^2"},
    {2, 1, 0, 255, 0, 0, 3, 24, "FRSTINX", "Frost Index", "kg/day"},
    {0, 0, 0, 255, 7, 1, 1, 225, "FRZR", "Freezing Rain", "kg/m^2"},
    {10, 0, 0, 255, 7, 1, 3, 204, "FRZSPR", "Freezing Spray", "-"},
    {0, 1, 0, 255, 0, 0, 1, 121, "FSNOWC", "Fraction of Snow Cover", "Proportion"},
    {0, 1, 0, 255, 0, 0, 6, 36, "FSTRPC", "Fraction of Stratiform Precipitation Cover", "Proportion"},
    {2, 1, 0, 255, 0, 0, 4, 5, "FWINX", "Fire Weath Index (Canadian Forest Service)", "Numeric"},
    {0, 1, 0, 255, 0, 0, 1, 95, "FZPRATE", "Freezing or Frozen Precipitation Rate", "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 18, 3, "GDCES", "Ground Deposition of Caesium 137", "Bq/m^2"},
    {0, 1, 0, 255, 0, 0, 18, 4, "GDIOD", "Ground Deposition of Iodine 131", "Bq/m^2"},
    {0, 1, 0, 255, 0, 0, 18, 5, "GDRADP", "Ground Deposition of Radioactive Pollutant", "Bq/m^2"},
    {0, 1, 0, 255, 0, 0, 191, 1, "GEOLAT", "Geographical Latitude", "deg"},
    {0, 1, 0, 255, 0, 0, 191, 2, "GEOLON", "Geographical Longitude", "deg"},
    {0, 1, 0, 255, 0, 0, 2, 43, "GEOWD", "Geostrophic Wind Direction", "deg"},
    {0, 1, 0, 255, 0, 0, 2, 44, "GEOWS", "Geostrophic Wind Speed", "m/s"},
    {2, 0, 0, 255, 7, 1, 0, 193, "GFLUX", "Ground Heat Flux", "W/m^2"},
    {2, 1, 0, 255, 0, 0, 0, 10, "GFLUX", "Ground Heat Flux", "W/m^2"},
    {2, 1, 0, 255, 0, 0, 5, 1, "GLACTMP", "Glacier Temperature", "K"},
    {0, 1, 0, 255, 0, 0, 3, 4, "GP", "Geopotential", "m^2/s^2"},
    {0, 1, 0, 255, 0, 0, 3, 9, "GPA", "Geopotential Height Anomaly", "gpm"},
    {0, 1, 0, 255, 0, 0, 1, 75, "GPRATE", "Graupel (Snow Pellets) Prepitation Rate", "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 4, 3, "GRAD", "Global Radiation Flux", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 7, 17, "GRDRN", "Gradient Richardson Number", "Numeric"},
    {0, 1, 0, 255, 0, 0, 1, 32, "GRLE", "Graupel", "kg/kg"},
    {3, 1, 0, 255, 0, 0, 6, 1, "GSOLEXP", "Global Solar Exposure", "J/m^2"},
    {3, 1, 0, 255, 0, 0, 6, 0, "GSOLIRR", "Global Solar Irradiance", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 2, 22, "GUST", "Wind Speed (Gust)", "m/s"},
    {0, 1, 0, 255, 0, 0, 3, 23, "GWD", "Gravity Wave Dissipation", "W/m^2"},
    {0, 0, 0, 255, 7, 1, 2, 210, "GWDU", "Gravity wave drag zonal acceleration", "m/s^2"},
    {0, 0, 0, 255, 7, 1, 2, 211, "GWDV", "Gravity wave drag meridional acceleration", "m/s^2"},
    {1, 1, 0, 255, 0, 0, 0, 9, "GWLOWS", "Group Water Lower Storage", "kg/m^2"},
    {2, 0, 0, 255, 7, 1, 0, 214, "GWREC", "Groundwater recharge", "kg/m^2"},
    {1, 1, 0, 255, 0, 0, 0, 8, "GWUPS", "Group Water Upper Storage", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 1, 31, "HAIL", "Hail", "m"},
    {0, 1, 0, 255, 0, 0, 1, 71, "HAILMXR", "Hail Mixing Ratio", "kg/kg"},
    {0, 1, 0, 255, 0, 0, 1, 73, "HAILPR", "Hail Prepitation Rate", "kg/m^2/s"},
    {0, 0, 0, 255, 7, 1, 19, 198, "HAILPROB", "Hail probability", "%"},
    {4, 1, 0, 255, 0, 0, 8, 2, "HARAD", "H-Alpha Radiance", "W/sr/m^2"},
    {0, 0, 0, 255, 7, 1, 19, 210, "HAVNI", "High-Level aviation interest", "-"},
    {0, 1, 0, 255, 0, 0, 6, 5, "HCDC", "High Cloud Cover", "%"},
    {0, 1, 0, 255, 0, 0, 6, 26, "HCONCB", "Height of Convective Cloud Base", "m"},
    {0, 1, 0, 255, 0, 0, 6, 27, "HCONCT", "Height of Convective Cloud Top", "m"},
    {0, 1, 0, 255, 0, 0, 0, 12, "HEATX", "Heat Index", "K"},
    {4, 1, 0, 255, 0, 0, 8, 6, "HELCOR", "Heliospheric Radiance", "W/sr/m^2"},
    {2, 1, 0, 255, 0, 0, 0, 24, "HFLUX", "Heat Flux", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 3, 5, "HGT", "Geopotential Height", "gpm"},
    {0, 1, 0, 255, 0, 0, 20, 62, "HGTMD", "Height of Mass Density", "m"},
    {0, 0, 0, 255, 7, 1, 3, 211, "HGTN", "Geopotential Height (nearest grid point)", "gpm"},
    {0, 0, 0, 255, 7, 1, 3, 203, "HGTX", "X-gradient of Height", "1/m"},
    {0, 0, 0, 255, 7, 1, 3, 204, "HGTY", "Y-gradient of Height", "1/m"},
    {0, 1, 0, 255, 0, 0, 19, 32, "HIFREL", "Highest Freezing Level", "m"},
    {2, 1, 0, 255, 0, 0, 4, 2, "HINDEX", "Haines Index", "Numeric"},
    {0, 1, 0, 255, 0, 0, 7, 8, "HLCY", "Storm Relative Helicity", "m^2/s^2"},
    {0, 1, 0, 255, 0, 0, 18, 16, "HMXACON", "Height of Maximum of Air Concentration", "m"},
    {0, 0, 0, 255, 7, 1, 3, 196, "HPBL", "Planetary Boundary Layer Height", "m"},
    {0, 1, 0, 255, 0, 0, 3, 18, "HPBL", "Planetary Boundary Layer Height", "m"},
    {4, 1, 0, 255, 0, 0, 2, 8, "HPRIMF", "hF", "m"},
    {0, 0, 0, 255, 7, 1, 19, 196, "HRCONO", "High risk convective outlook", "Categorical"},
    {0, 1, 0, 255, 0, 0, 15, 15, "HSR", "Hybrid Scan Reflectivity", "dB"},
    {0, 1, 0, 255, 0, 0, 15, 16, "HSRHT", "Hybrid Scan Reflectivity Height", "m"},
    {0, 1, 0, 255, 0, 0, 3, 7, "HSTDV", "Standard Deviation of Height", "m"},
    {10, 1, 0, 255, 0, 0, 0, 3, "HTSGW", "Significant Height of Combined Wind Waves and Swell", "m"},
    {0, 1, 0, 255, 0, 0, 3, 3, "ICAHT", "ICAO Standard Atmosphere Reference Height", "m"},
    {10, 1, 0, 255, 0, 0, 2, 0, "ICEC", "Ice Cover", "Proportion"},
    {1, 1, 0, 255, 0, 0, 2, 7, "ICECIL", "Ice Cover", "Proportion"},
    {10, 1, 0, 255, 0, 0, 2, 7, "ICED", "Ice Divergence", "1/s"},
    {10, 1, 0, 255, 0, 0, 2, 6, "ICEG", "Ice Growth Rate", "m/s"},
    {10, 1, 0, 255, 0, 0, 2, 9, "ICEPRS", "Module of Ice Internal Pressure", "Pa*m"},
    {0, 1, 0, 255, 0, 0, 19, 27, "ICESC", "Icing Scenario", "-"},
    {0, 1, 0, 255, 0, 0, 19, 37, "ICESEV", "Icing Severity", "non-dim"},
    {1, 1, 0, 255, 0, 0, 2, 6, "ICETIL", "Ice Temperature", "K"},
    {10, 1, 0, 255, 0, 0, 2, 1, "ICETK", "Ice Thickness", "m"},
    {10, 1, 0, 255, 0, 0, 2, 8, "ICETMP", "Ice Temperature", "K"},
    {0, 1, 0, 255, 0, 0, 19, 7, "ICI", "Icing", "-"},
    {0, 1, 0, 255, 0, 0, 19, 6, "ICIB", "Icing Base", "m"},
    {0, 1, 0, 255, 0, 0, 19, 20, "ICIP", "Icing", "%"},
    {0, 1, 0, 255, 0, 0, 19, 5, "ICIT", "Icing Top", "m"},
    {0, 1, 0, 255, 0, 0, 1, 23, "ICMR", "Ice Water Mixing Ratio", "kg/kg"},
    {0, 0, 0, 255, 7, 1, 19, 233, "ICPRB", "Icing probability", "non-dim"},
    {0, 0, 0, 255, 7, 1, 19, 234, "ICSEV", "Icing Severity", "non-dim"},
    {1, 1, 0, 255, 0, 0, 2, 5, "ICTKIL", "Ice Thickness", "m"},
    {2, 0, 0, 255, 7, 1, 0, 207, "ICWAT", "Ice-free water surface", "%"},
    {0, 1, 0, 255, 0, 0, 1, 20, "ILIQW", "Integrated Liquid Water", "kg/m^2"},
    {10, 1, 0, 255, 0, 0, 0, 27, "IMFTSW", "Inverse Mean Frequency of The Total Swell", "s"},
    {10, 1, 0, 255, 0, 0, 0, 26, "IMFWW", "Inverse Mean Frequency of The Wind Waves", "s"},
    {255, 0, 0, 255, 7, 1, 255, 255, "IMGD", "Image data", "-"},
    {10, 1, 0, 255, 0, 0, 0, 25, "IMWF", "Inverse Mean Wave Frequency", "s"},
    {2, 1, 0, 255, 0, 0, 4, 9, "INFSINX", "Initial Fire Spread Index (Canadian Forest Service)", "Numeric"},
    {4, 1, 0, 255, 0, 0, 4, 3, "INTEFLUX", "Electron Flux (Integral)", "1/(m^2s*sr)"},
    {10, 0, 0, 255, 7, 1, 4, 196, "INTFD", "Interface Depths", "m"},
    {4, 1, 0, 255, 0, 0, 4, 5, "INTIFLUX", "Heavy Ion Flux (iIntegral)", "1/(m^2s*sr)"},
    {4, 1, 0, 255, 0, 0, 4, 1, "INTPFLUX", "Proton Flux (Integral)", "1/(m^2s*sr)"},
    {4, 1, 0, 255, 0, 0, 2, 3, "IONDEN", "Ion Density", "1/m^3"},
    {4, 1, 0, 255, 0, 0, 0, 3, "IONTMP", "Ion Temperature", "K"},
    {0, 1, 0, 255, 0, 0, 1, 68, "IPRATE", "Ice Pellets Precipitation Rate", "kg/m^2/s"},
    {3, 1, 0, 255, 0, 0, 1, 1, "IRRATE", "Instantaneous Rain Rate", "kg/m^2/s"},
    {10, 1, 0, 255, 0, 0, 191, 0, "IRTSEC", "Seconds Prior To Initial Reference Time", "s"},
    {3, 1, 0, 255, 0, 0, 5, 0, "ISSTMP", "Interface Sea-Surface Temperature", "K"},
    {0, 0, 0, 255, 7, 1, 19, 235, "JFWPRB", "Joint Fire Weather Probability", "%"},
    {10, 0, 0, 255, 7, 1, 3, 201, "KENG", "Kinetic Energy", "J/kg"},
    {0, 1, 0, 255, 0, 0, 7, 3, "KOX", "KO Index", "K"},
    {10, 1, 0, 255, 0, 0, 0, 43, "KSSEW", "Kurtosis of The Sea Surface Elevation Due to Waves", "-"},
    {0, 1, 0, 255, 0, 0, 7, 2, "KX", "K Index", "K"},
    {0, 0, 0, 255, 7, 1, 7, 198, "LAI", "Leaf Area Index", "Numeric"},
    {2, 1, 0, 255, 0, 0, 0, 0, "LAND", "Land Cover (0=sea, 1=land)", "Proportion"},
    {1, 1, 0, 255, 0, 0, 2, 8, "LANDIL", "Land Cover (0=water, 1=land)", "Proportion"},
    {2, 0, 0, 255, 7, 1, 0, 218, "LANDN", "Land-sea coverage (nearest neighbor) [land=1,sea=0]", "-"},
    {2, 1, 0, 255, 0, 0, 0, 8, "LANDU", "Land Use", "-"},
    {0, 0, 0, 255, 7, 1, 2, 202, "LAPP", "Latitude of Presure Point", "deg"},
    {0, 1, 0, 255, 0, 0, 0, 8, "LAPR", "Lapse Rate", "K/m"},
    {0, 0, 0, 255, 7, 1, 2, 198, "LAUV", "Latitude of U Wind Component of Velocity", "deg"},
    {0, 0, 0, 255, 7, 1, 19, 209, "LAVNI", "Low-Level aviation interest", "-"},
    {0, 0, 0, 255, 7, 1, 2, 200, "LAVV", "Latitude of V Wind Component of Velocity", "deg"},
    {0, 0, 0, 255, 7, 1, 3, 205, "LAYTH", "Layer Thickness", "m"},
    {0, 1, 0, 255, 0, 0, 6, 3, "LCDC", "Low Cloud Cover", "%"},
    {10, 0, 0, 255, 7, 1, 3, 203, "LCH", "Heat Exchange Coefficient", "-"},
    {2, 1, 0, 255, 0, 0, 0, 28, "LEAINX", "Leaf Area Index", "Numeric"},
    {0, 0, 0, 255, 7, 1, 7, 192, "LFTX", "Surface Lifted Index", "K"},
    {0, 1, 0, 255, 0, 0, 7, 10, "LFTX", "Surface Lifted Index", "K"},
    {0, 1, 0, 255, 0, 0, 0, 10, "LHTFL", "Latent Heat Net Flux", "W/m^2"},
    {0, 0, 0, 255, 7, 1, 1, 229, "LICEAC", "Line Ice Accumulation (FRAM)", "kg/m^2"},
    {0, 0, 0, 255, 7, 1, 13, 195, "LIPMF", "Integrated column particulate matter (fine)", "log10(10^-6g/m^3)"},
    {2, 1, 0, 255, 0, 0, 3, 10, "LIQVSM", "Liquid Volumetric Soil Moisture (Non-Frozen)", "m^3/m^3"},
    {0, 1, 0, 255, 0, 0, 15, 4, "LMAXBR", "Layer Maximum Base Reflectivity", "dB"},
    {4, 1, 0, 255, 0, 0, 7, 0, "LMBINT", "Limb Intensity", "1/m^2/s"},
    {0, 0, 0, 255, 7, 1, 3, 210, "LMH", "Mass Point Model Surface", "-"},
    {0, 0, 0, 255, 7, 1, 2, 218, "LMV", "Velocity Point Model Surface", "-"},
    {0, 0, 0, 255, 7, 1, 2, 203, "LOPP", "Longitude of Presure Point", "deg"},
    {0, 0, 0, 255, 7, 1, 2, 199, "LOUV", "Longitude of U Wind Component of Velocity", "deg"},
    {0, 0, 0, 255, 7, 1, 2, 201, "LOVV", "Longitude of V Wind Component of Velocity", "deg"},
    {2, 1, 0, 255, 0, 0, 3, 3, "LOWLSM", "Lower Layer Soil Moisture", "kg/m^3"},
    {0, 0, 0, 255, 7, 1, 13, 194, "LPMTF", "Particulate matter (fine)", "log10(10^-6g/m^3)"},
    {0, 0, 0, 255, 7, 1, 3, 201, "LPSX", "X-gradient of Log Pressure", "1/m"},
    {0, 0, 0, 255, 7, 1, 3, 202, "LPSY", "Y-gradient of Log Pressure", "1/m"},
    {0, 0, 0, 255, 7, 1, 0, 195, "LRGHR", "Large Scale Condensate Heating Rate", "K/s"},
    {0, 0, 0, 255, 7, 1, 1, 217, "LRGMR", "Large scale moistening rate", "kg/kg/s"},
    {2, 0, 0, 255, 7, 1, 0, 212, "LSOIL", "Liquid soil moisture content (non-frozen)", "kg/m^2"},
    {2, 0, 0, 255, 7, 1, 3, 199, "LSPA", "Land Surface Precipitation Accumulation", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 1, 54, "LSPRATE", "Large Scale Precipitation Rate", "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 1, 77, "LSRRATE", "Large Scale Rain Rate", "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 1, 59, "LSSRATE", "Large Scale Snowfall Rate", "m/s"},
    {0, 1, 0, 255, 0, 0, 1, 56, "LSSRWE", "Large Scale Snowfall Rate Water Equivalent", "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 1, 47, "LSWP", "Large Scale Water Precipitation (Non-Convective)", "kg/m^2"},
    {0, 0, 0, 255, 7, 1, 17, 192, "LTNG", "Lightning", "non-dim"},
    {0, 1, 0, 255, 0, 0, 17, 0, "LTNGSD", "Lightning Strike Density", "1/m^2/s"},
    {0, 1, 0, 255, 0, 0, 17, 1, "LTPINX", "Lightning Potential Index (LPI)", "J/kg"},
    {0, 1, 0, 255, 0, 0, 5, 2, "LWAVR", "Long-Wave Radiation Flux", "W/m^2"},
    {0, 0, 0, 255, 7, 1, 5, 194, "LWHR", "Long-Wave Radiative Heating Rate", "K/s"},
    {0, 1, 0, 255, 0, 0, 4, 5, "LWRAD", "Radiance (with respect to wave number)", "W/m/sr"},
    {2, 1, 0, 255, 0, 0, 3, 23, "LWSNWP", "Liquid Water in Snow Pack", "kg/m^2"},
    {4, 1, 0, 255, 0, 0, 8, 7, "MASK", "Thematic Mask", "Numeric"},
    {0, 1, 0, 255, 0, 0, 6, 38, "MASSDCD", "Mass Density of Cloud Droplets", "kg/m^3"},
    {0, 1, 0, 255, 0, 0, 6, 39, "MASSDCI", "Mass Density of Cloud Ice", "kg/m^3"},
    {0, 1, 0, 255, 0, 0, 20, 0, "MASSDEN", "Mass Density (Concentration)", "kg/m^3"},
    {0, 1, 0, 255, 0, 0, 1, 98, "MASSDG", "Mass Density of Graupel", "kg/m^3"},
    {0, 1, 0, 255, 0, 0, 1, 99, "MASSDH", "Mass Density of Hail", "kg/m^3"},
    {0, 1, 0, 255, 0, 0, 1, 96, "MASSDR", "Mass Density of Rain", "kg/m^3"},
    {0, 1, 0, 255, 0, 0, 1, 97, "MASSDS", "Mass Density of Snow", "kg/m^3"},
    {0, 1, 0, 255, 0, 0, 20, 2, "MASSMR", "Mass Mixing Ratio (Mass Fraction in Air)", "kg/kg"},
    {0, 1, 0, 255, 0, 0, 18, 15, "MAXACON", "Maximum of Air Concentration in Layer", "Bq/m^3"},
    {0, 1, 0, 255, 0, 0, 1, 28, "MAXAH", "Maximum Absolute Humidity", "kg/m^3"},
    {0, 0, 0, 255, 7, 1, 2, 221, "MAXDVV", "Hourly Maximum of Downward Vertical Velocity in the lowest 400hPa", "m/s"},
    {0, 1, 0, 255, 0, 0, 2, 21, "MAXGUST", "Maximum Wind Speed", "m/s"},
    {0, 0, 0, 255, 7, 1, 16, 198, "MAXREF", "Hourly Maximum of Simulated Reflectivity at 1 km AGL", "dB"},
    {0, 1, 0, 255, 0, 0, 1, 27, "MAXRH", "Maximum Relative Humidity", "%"},
    {0, 0, 0, 255, 7, 1, 2, 220, "MAXUVV", "Hourly Maximum of Upward Vertical Velocity in the lowest 400hPa", "m/s"},
    {0, 0, 0, 255, 7, 1, 2, 222, "MAXUW", "U Component of Hourly Maximum 10m Wind Speed", "m/s"},
    {0, 0, 0, 255, 7, 1, 2, 223, "MAXVW", "V Component of Hourly Maximum 10m Wind Speed", "m/s"},
    {10, 1, 0, 255, 0, 0, 0, 24, "MAXWH", "Maximum Individual Wave Height", "m"},
    {0, 1, 0, 255, 0, 0, 6, 4, "MCDC", "Medium Cloud Cover", "%"},
    {0, 1, 0, 255, 0, 0, 1, 26, "MCONV", "Horizontal Moisture Convergence", "kg/kg/s"},
    {0, 1, 0, 255, 0, 0, 6, 40, "MDCCWD", "Mass Density of Convective Cloud Water Droplets", "kg/m^3"},
    {0, 0, 0, 255, 7, 1, 1, 197, "MDIV", "Horizontal Moisture Divergence", "kg/kg/s"},
    {0, 1, 0, 255, 0, 0, 1, 38, "MDIVER", "Horizontal Moisture Divergence", "kg/kg/s"},
    {0, 1, 0, 255, 0, 0, 1, 112, "MDLWGVA",
     "Mass Density of Liquid Water Coating on Graupel Expressed as Mass of Liquid Water per Unit Volume of Air",
     "kg/m^3"},
    {0, 1, 0, 255, 0, 0, 1, 109, "MDLWHVA",
     "Mass Density of Liquid Water Coating on Hail Expressed as Mass of Liquid Water per Unit Volume of Air", "kg/m^3"},
    {0, 1, 0, 255, 0, 0, 1, 115, "MDLWSVA",
     "Mass Density of Liquid Water Coating on Snow Expressed as Mass of Liquid Water per Unit Volume of Air", "kg/m^3"},
    {3, 1, 0, 255, 0, 0, 2, 30, "MEACST", "Measurement cost", "Numeric"},
    {0, 0, 0, 255, 7, 1, 6, 200, "MFLUX", "Convective Cloud Mass Flux", "Pa/s"},
    {0, 0, 0, 255, 7, 1, 2, 193, "MFLX", "Horizontal Momentum Flux", "N/m^2"},
    {0, 1, 0, 255, 0, 0, 2, 26, "MFLX", "Horizontal Momentum Flux", "N/m^2"},
    {0, 1, 0, 255, 0, 0, 0, 14, "MINDPD", "Minimum Dew Point Depression", "K"},
    {0, 0, 0, 255, 7, 1, 1, 198, "MINRH", "Minimum Relative Humidity", "%"},
    {0, 1, 0, 255, 0, 0, 19, 3, "MIXHT", "Mixed Layer Depth", "m"},
    {0, 0, 0, 255, 7, 1, 19, 204, "MIXLY", "Number of mixed layers next to surface", "Integer"},
    {0, 1, 0, 255, 0, 0, 1, 2, "MIXR", "Humidity Mixing Ratio", "kg/kg"},
    {0, 0, 0, 255, 7, 1, 191, 195, "MLYNO", "Model Layer number (From bottom up)", "-"},
    {0, 1, 0, 255, 0, 0, 1, 114, "MMLWGDA",
     "Mass Mixing Ratio of Liquid Water Coating on Graupel Expressed as Mass of Liquid Water per Unit Mass of Dry Air",
     "kg/kg"},
    {0, 1, 0, 255, 0, 0, 1, 111, "MMLWHDA",
     "Mass Mixing Ratio of Liquid Water Coating on Hail Expressed as Mass of Liquid Water per Unit Mass of Dry Air",
     "kg/kg"},
    {0, 1, 0, 255, 0, 0, 1, 117, "MMLWSDA",
     "Mass Mixing Ratio of Liquid Water Coating on Snow Expressed as Mass of Liquid Water per Unit Mass of Dry Air",
     "kg/kg"},
    {0, 1, 0, 255, 0, 0, 2, 6, "MNTSF", "Montgomery Stream Function", "m^2/s^2"},
    {0, 0, 0, 255, 7, 1, 7, 200, "MNUPHL", "Hourly Minimum of Updraft Helicity", "m^2/s^2"},
    {10, 1, 0, 255, 0, 0, 0, 40, "MNWSOW", "10 Metre Neutral Wind Speed Over Waves", "m/s"},
    {0, 1, 0, 255, 0, 0, 20, 64, "MOLRDRYA", "Mole fraction with respect to dry air", "mol/mol"},
    {0, 1, 0, 255, 0, 0, 20, 65, "MOLRWETA", "Mole fraction with respect to wet air", "mol/mol"},
    {10, 1, 0, 255, 0, 0, 191, 1, "MOSF", "Meridional Overturning Stream Function", "m^3/s"},
    {0, 0, 0, 255, 7, 1, 19, 195, "MRCONO", "Moderate risk convective outlook", "Categorical"},
    {0, 0, 0, 255, 7, 1, 3, 192, "MSLET", "MSLP (Eta model reduction)", "Pa"},
    {0, 0, 0, 255, 7, 1, 3, 198, "MSLMA", "MSLP (MAPS System Reduction)", "Pa"},
    {0, 1, 0, 255, 0, 0, 20, 16, "MSSRDRYA", "Mass mixing ratio with respect to dry air", "kg/kg"},
    {0, 1, 0, 255, 0, 0, 20, 17, "MSSRWETA", "Mass mixing ratio with respect to wet air", "kg/kg"},
    {10, 1, 0, 255, 0, 0, 0, 20, "MSSW", "Mean Square Slope of Waves", "-"},
    {2, 0, 0, 255, 7, 1, 0, 194, "MSTAV", "Moisture Availability", "%"},
    {2, 1, 0, 255, 0, 0, 0, 11, "MSTAV", "Moisture Availability", "%"},
    {2, 1, 0, 255, 0, 0, 0, 7, "MTERH", "Model Terrain Height", "m"},
    {10, 1, 0, 255, 0, 0, 4, 1, "MTHA", "Main Thermocline Anomaly", "m"},
    {10, 1, 0, 255, 0, 0, 4, 0, "MTHD", "Main Thermocline Depth", "m"},
    {10, 1, 0, 255, 0, 0, 2, 11, "MVCICEP",
     "Meridional Vector Component of Vertically Integrated Ice Internal Pressure", "Pa*m"},
    {10, 1, 0, 255, 0, 0, 0, 53, "MWDFSWEL", "Mean wave direction of first swell partition", "deg"},
    {10, 1, 0, 255, 0, 0, 0, 41, "MWDIRW", "10 Metre Wind Direction Over Waves", "deg"},
    {10, 1, 0, 255, 0, 0, 0, 54, "MWDSSWEL", "Mean wave direction of second swell partition", "deg"},
    {10, 1, 0, 255, 0, 0, 0, 55, "MWDTSWEL", "Mean wave direction of third swell partition", "deg"},
    {10, 1, 0, 255, 0, 0, 0, 50, "MWPFSWEL", "Mean wave period of first swell partition", "s"},
    {10, 1, 0, 255, 0, 0, 0, 51, "MWPSSWEL", "Mean wave period of second swell partition", "s"},
    {10, 1, 0, 255, 0, 0, 0, 52, "MWPTSWEL", "Mean wave period of third swell partition", "s"},
    {10, 1, 0, 255, 0, 0, 0, 15, "MWSPER", "Mean Period of Combined Wind Waves and Swell", "s"},
    {0, 1, 0, 255, 0, 0, 19, 28, "MWTURB", "Mountain Wave Turbulence (Eddy Dissipation Rate)", "m2/3s-1"},
    {0, 1, 0, 255, 0, 0, 19, 31, "MXEDPRM", "Maximum of Eddy Dissipation Parameter in Layer", "m2/3s-1"},
    {0, 1, 0, 255, 0, 0, 20, 61, "MXMASSD", "Maximum of Mass Density", "kg/m^3"},
    {0, 0, 0, 255, 7, 1, 19, 192, "MXSALB", "Maximum Snow Albedo", "%"},
    {0, 1, 0, 255, 0, 0, 19, 17, "MXSALB", "Maximum Snow Albedosee Note 1", "%"},
    {0, 0, 0, 255, 7, 1, 7, 199, "MXUPHL", "Hourly Maximum of Updraft Helicity over Layer 2km to 5 km AGL", "m^2/s^2"},
    {10, 1, 0, 255, 0, 0, 0, 30, "MZPTSW", "Mean Zero-Crossing Period of The Total Swell", "s"},
    {10, 1, 0, 255, 0, 0, 0, 29, "MZPWW", "Mean Zero-Crossing Period of The Wind Waves", "s"},
    {10, 1, 0, 255, 0, 0, 0, 28, "MZWPER", "Mean Zero-Crossing Wave Period", "s"},
    {0, 0, 0, 255, 7, 1, 4, 202, "NBDSF", "Near IR Beam Downward Solar Flux", "W/m^2"},
    {0, 0, 0, 255, 7, 1, 19, 213, "NBSALB", "Near IR, Black Sky Albedo", "%"},
    {0, 1, 0, 255, 0, 0, 6, 29, "NCCICE", "Number Concentration of Cloud Ice", "1/kg"},
    {0, 0, 0, 255, 7, 1, 1, 207, "NCIP", "Number concentration for ice particles", "non-dim"},
    {0, 1, 0, 255, 0, 0, 6, 28, "NCONCD", "Number Concentration of Cloud Droplets", "1/kg"},
    {0, 1, 0, 255, 0, 0, 1, 9, "NCPCP", "Large-Scale Precipitation (non-convective)", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 6, 31, "NDCICE", "Number Density of Cloud Ice", "1/m^3"},
    {0, 0, 0, 255, 7, 1, 4, 203, "NDDSF", "Near IR Diffuse Downward Solar Flux", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 6, 30, "NDENCD", "Number Density of Cloud Droplets", "1/m^3"},
    {2, 0, 0, 255, 7, 1, 0, 217, "NDVI", "Normalized Difference Vegetation Index", "-"},
    {2, 1, 0, 255, 0, 0, 0, 31, "NDVINX", "Normalized Differential Vegetation Index (NDVI)", "Numeric"},
    {0, 0, 0, 255, 7, 1, 191, 192, "NLAT", "Latitude (-90 to 90)", "deg"},
    {0, 0, 0, 255, 7, 1, 191, 196, "NLATN", "Latitude (nearest neighbor) (-90 to 90)", "deg"},
    {0, 0, 0, 255, 7, 1, 3, 206, "NLGSP", "Natural Log of Surface Pressure", "ln(kPa)"},
    {0, 1, 0, 255, 0, 0, 3, 25, "NLPRES", "Natural Logarithm of Pressure in Pa", "Numeric"},
    {0, 1, 0, 255, 0, 0, 5, 6, "NLWRCS", "Net Long-Wave Radiation Flux, Clear Sky", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 5, 5, "NLWRF", "Net Long-Wave Radiation Flux", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 5, 0, "NLWRS", "Net Long-Wave Radiation Flux (Surface)", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 5, 1, "NLWRT", "Net Long-Wave Radiation Flux (Top of Atmosphere)", "W/m^2"},
    {3, 1, 0, 255, 0, 0, 1, 6, "NPIXU", "Number Of Pixels Used", "Numeric"},
    {0, 1, 0, 255, 0, 0, 4, 9, "NSWRF", "Net Short Wave Radiation Flux", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 4, 11, "NSWRFCS", "Net Short-Wave Radiation Flux, Clear Sky", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 4, 0, "NSWRS", "Net Short-Wave Radiation Flux (Surface)", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 4, 1, "NSWRT", "Net Short-Wave Radiation Flux (Top of Atmosphere)", "W/m^2"},
    {4, 1, 0, 255, 0, 0, 4, 6, "NTRNFLUX", "Cosmic Ray Neutron Flux", "1/h"},
    {0, 1, 0, 255, 0, 0, 2, 37, "NTSS", "Northward Turbulent Surface Stress", "N/m^2*s"},
    {0, 1, 0, 255, 0, 0, 1, 106, "NUMDG", "Number Density of Graupel", "1/m^3"},
    {0, 1, 0, 255, 0, 0, 1, 107, "NUMDH", "Number Density of Hail", "1/m^3"},
    {0, 1, 0, 255, 0, 0, 1, 104, "NUMDR", "Number Density of Rain", "1/m^3"},
    {0, 1, 0, 255, 0, 0, 1, 105, "NUMDS", "Number Density of Snow", "1/m^3"},
    {0, 1, 0, 255, 0, 0, 2, 34, "NWIND", "Normal Wind Component", "m/s"},
    {0, 0, 0, 255, 7, 1, 19, 214, "NWSALB", "Near IR, White Sky Albedo", "%"},
    {10, 1, 0, 255, 0, 0, 0, 19, "NWSTR", "Normalised Waves Stress", "-"},
    {0, 1, 0, 255, 0, 0, 2, 40, "NWTPARM", "Northward Wind Tendency Due to Parameterizations", "m/s^2"},
    {0, 0, 0, 255, 7, 1, 14, 192, "O3MR", "Ozone Mixing Ratio", "kg/kg"},
    {0, 1, 0, 255, 0, 0, 14, 1, "O3MR", "Ozone Mixing Ratio", "kg/kg"},
    {10, 0, 0, 255, 7, 1, 4, 197, "OHC", "Ocean Heat Content", "J/m^2"},
    {0, 0, 0, 255, 7, 1, 2, 215, "OMGALF", "Omega (Dp/Dt) divide by density", "K"},
    {10, 0, 0, 255, 7, 1, 1, 192, "OMLU", "Ocean Mixed Layer U Velocity", "m/s"},
    {10, 0, 0, 255, 7, 1, 1, 193, "OMLV", "Ocean Mixed Layer V Velocity", "m/s"},
    {0, 0, 0, 255, 7, 1, 3, 217, "ORASNW", "Orographic Asymmetry, NW Component", "-"},
    {0, 0, 0, 255, 7, 1, 3, 215, "ORASS", "Orographic Asymmetry, S Component", "-"},
    {0, 0, 0, 255, 7, 1, 3, 216, "ORASSW", "Orographic Asymmetry, SW Component", "-"},
    {0, 0, 0, 255, 7, 1, 3, 214, "ORASW", "Orographic Asymmetry, W Component", "-"},
    {0, 0, 0, 255, 7, 1, 3, 213, "ORCONV", "Orographic Convexity", "-"},
    {0, 0, 0, 255, 7, 1, 3, 221, "ORLSNW", "Orographic Length Scale, NW Component", "-"},
    {0, 0, 0, 255, 7, 1, 3, 219, "ORLSS", "Orographic Length Scale, S Component", "-"},
    {0, 0, 0, 255, 7, 1, 3, 220, "ORLSSW", "Orographic Length Scale, SW Component", "-"},
    {0, 0, 0, 255, 7, 1, 3, 218, "ORLSW", "Orographic Length Scale, W Component", "-"},
    {10, 1, 0, 255, 0, 0, 4, 4, "OVHD", "Ocean Vertical Heat Diffusivity", "m^2/s"},
    {10, 1, 0, 255, 0, 0, 4, 6, "OVMD", "Ocean Vertical Momentum Diffusivity", "m^2/s"},
    {10, 1, 0, 255, 0, 0, 4, 5, "OVSD", "Ocean Vertical Salt Diffusivity", "m^2/s"},
    {10, 0, 0, 255, 7, 1, 3, 253, "OWASHP", "Overwash Occurrence Probability", "%"},
    {0, 0, 0, 255, 7, 1, 14, 194, "OZCAT", "Categorical Ozone Concentration", "non-dim"},
    {0, 0, 0, 255, 7, 1, 14, 193, "OZCON", "Ozone Concentration", "ppb"},
    {0, 0, 0, 255, 7, 1, 14, 200, "OZMAX1", "Ozone Daily Max from 1-hour Average", "ppbV"},
    {0, 0, 0, 255, 7, 1, 14, 201, "OZMAX8", "Ozone Daily Max from 8-hour Average", "ppbV"},
    {10, 0, 0, 255, 7, 1, 3, 196, "P2OMLT", "Ocean Mixed Layer Potential Density (Reference 2000m)", "kg/m^3"},
    {3, 1, 0, 255, 0, 0, 3, 2, "PBINFRC", "Probability of Encountering Instrument Flight Rules Conditions", "%"},
    {3, 1, 0, 255, 0, 0, 3, 1, "PBLIFRC", "Probability of Encountering Low Instrument Flight Rules Conditions", "%"},
    {0, 1, 0, 255, 0, 0, 19, 12, "PBLREG", "Planetary Boundary Layer Regime", "-"},
    {3, 1, 0, 255, 0, 0, 3, 0, "PBMVFRC", "Probability of Encountering Marginal Visual Flight Rules Conditions", "%"},
    {0, 0, 0, 255, 7, 1, 1, 234, "PCPDUR", "Precipitation Duration", "hour"},
    {0, 0, 0, 255, 7, 1, 14, 202, "PDMAX1", "PM 2.5 Daily Max from 1-hour Average", "10^-6g/m^3"},
    {0, 0, 0, 255, 7, 1, 14, 203, "PDMAX24", "PM 2.5 Daily Max from 24-hour Average", "10^-6g/m^3"},
    {10, 1, 0, 255, 0, 0, 0, 11, "PERPW", "Primary Wave Mean Period", "s"},
    {1, 1, 0, 255, 0, 0, 0, 16, "PERRATE", "Percolation Rate", "kg/m^2/s"},
    {10, 1, 0, 255, 0, 0, 0, 13, "PERSW", "Secondary Wave Mean Period", "s"},
    {0, 0, 0, 255, 7, 1, 1, 199, "PEVAP", "Potential Evaporation", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 1, 40, "PEVAP", "Potential Evaporation", "kg/m^2"},
    {0, 0, 0, 255, 7, 1, 1, 200, "PEVPR", "Potential Evaporation Rate", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 1, 41, "PEVPR", "Potential Evaporation Rate", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 4, 10, "PHOTAR", "Photosynthetically Active Radiation", "W/m^2"},
    {3, 1, 0, 255, 0, 0, 0, 8, "PIXST", "Pixel scene type", "-"},
    {0, 1, 0, 255, 0, 0, 7, 0, "PLI", "Parcel Lifted Index (to 500 hPa)", "K"},
    {0, 0, 0, 255, 7, 1, 3, 200, "PLPL", "Pressure of level from which parcel was lifted", "Pa"},
    {4, 1, 0, 255, 0, 0, 2, 0, "PLSMDEN", "Particle Number Density", "1/m^3"},
    {10, 1, 0, 255, 0, 0, 0, 23, "PMAXWH", "Period of Maximum Individual Wave Height", "s"},
    {0, 0, 0, 255, 7, 1, 13, 192, "PMTC", "Particulate matter (coarse)", "10^-6g/m^3"},
    {0, 0, 0, 255, 7, 1, 13, 193, "PMTF", "Particulate matter (fine)", "10^-6g/m^3"},
    {1, 1, 0, 255, 0, 0, 1, 2, "POP", "Probability of 0.01 inch of precipitation (POP)", "%"},
    {2, 0, 0, 255, 7, 1, 3, 197, "POROS", "Soil Porosity", "Proportion"},
    {2, 1, 0, 255, 0, 0, 3, 9, "POROS", "Soil Porosity", "Proportion"},
    {0, 1, 0, 255, 0, 0, 0, 2, "POT", "Potential Temperature", "K"},
    {0, 0, 0, 255, 7, 1, 14, 196, "POZ", "Ozone Production", "kg/kg/s"},
    {0, 0, 0, 255, 7, 1, 14, 199, "POZO", "Ozone Production from Column Ozone Term", "kg/kg/s"},
    {0, 0, 0, 255, 7, 1, 14, 198, "POZT", "Ozone Production from Temperature Term", "kg/kg/s"},
    {10, 1, 0, 255, 0, 0, 0, 36, "PPERTS", "Peak Period of The Total Swell", "s"},
    {10, 1, 0, 255, 0, 0, 0, 35, "PPERWW", "Peak Period of The Wind Waves", "s"},
    {1, 0, 0, 255, 7, 1, 1, 194, "PPFFG", "Probability of precipitation exceeding flash flood guidance values", "%"},
    {0, 0, 0, 255, 7, 1, 1, 231, "PPINDX", "Precipitation Potential Index", "%"},
    {1, 1, 0, 255, 0, 0, 1, 1, "PPOSP",
     "Percent Precipitation in a sub-period of an overall period (encoded as a percent accumulation over the "
     "sub-period)",
     "%"},
    {10, 1, 0, 255, 0, 0, 3, 3, "PRACTSAL", "Practical salinity", "Numeric"},
    {0, 1, 0, 255, 0, 0, 1, 7, "PRATE", "Precipitation Rate", "kg/m^2/s"},
    {4, 1, 0, 255, 0, 0, 0, 4, "PRATMP", "Parallel Temperature", "K"},
    {0, 1, 0, 255, 0, 0, 15, 5, "PREC", "Precipitation", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 3, 0, "PRES", "Pressure", "Pa"},
    {0, 1, 0, 255, 0, 0, 3, 8, "PRESA", "Pressure Anomaly", "Pa"},
    {0, 1, 0, 255, 0, 0, 3, 13, "PRESALT", "Pressure Altitude", "m"},
    {0, 0, 0, 255, 7, 1, 3, 212, "PRESN", "Pressure (nearest grid point)", "Pa"},
    {0, 1, 0, 255, 0, 0, 3, 1, "PRMSL", "Pressure Reduced to MSL", "Pa"},
    {0, 0, 0, 255, 7, 1, 1, 232, "PROBCIP", "Probability Cloud Ice Present", "%"},
    {0, 0, 0, 255, 7, 1, 19, 221, "PROCON", "Probability of Convection", "%"},
    {4, 1, 0, 255, 0, 0, 2, 2, "PROTDEN", "Proton Density", "1/m^3"},
    {4, 1, 0, 255, 0, 0, 0, 2, "PROTTMP", "Proton Temperature", "K"},
    {4, 1, 0, 255, 0, 0, 0, 5, "PRPTMP", "Perpendicular Temperature", "K"},
    {0, 0, 0, 255, 7, 1, 19, 216, "PRSIGSVR", "Total Probability of Extreme Severe Thunderstorms (Days 2,3)", "%"},
    {0, 0, 0, 255, 7, 1, 19, 215, "PRSVR", "Total Probability of Severe Thunderstorms (Days 2,3)", "%"},
    {10, 1, 0, 255, 0, 0, 4, 21, "PRTSAL", "Practical salinity", "psu(numeric)"},
    {0, 1, 0, 255, 0, 0, 19, 36, "PSNOWS", "Presence of Snow Squalls", "-"},
    {0, 1, 0, 255, 0, 0, 3, 2, "PTEND", "Pressure Tendency", "Pa/s"},
    {0, 1, 0, 255, 0, 0, 1, 19, "PTYPE", "Precipitation Type", "-"},
    {0, 0, 0, 255, 7, 1, 2, 219, "PVMWW", "Potential Vorticity (Mass-Weighted)", "1/s/m"},
    {0, 1, 0, 255, 0, 0, 2, 14, "PVORT", "Potential Vorticity", "Km^2/kg/s"},
    {0, 1, 0, 255, 0, 0, 1, 3, "PWAT", "Precipitable Water", "kg/m^2"},
    {10, 1, 0, 255, 0, 0, 0, 46, "PWAVEDIR", "Peak wave direction", "deg"},
    {0, 1, 0, 255, 0, 0, 1, 30, "PWCAT", "Precipitable Water Category", "-"},
    {10, 1, 0, 255, 0, 0, 0, 34, "PWPER", "Peak Wave Period", "s"},
    {0, 0, 0, 255, 7, 1, 1, 226, "PWTHER", "Predominant Weather", "Numeric"},
    {0, 0, 0, 255, 7, 1, 1, 219, "QMAX", "Maximum specific humidity at 2m", "kg/kg"},
    {0, 0, 0, 255, 7, 1, 1, 220, "QMIN", "Minimum specific humidity at 2m", "kg/kg"},
    {2, 0, 0, 255, 7, 1, 0, 215, "QREC", "Flood plain recharge", "kg/m^2"},
    {0, 0, 0, 255, 7, 1, 1, 218, "QZ0", "Specific humidity at top of viscous sublayer", "kg/kg"},
    {2, 0, 0, 255, 7, 1, 3, 202, "RADT", "Effective Radiative Skin Temperature", "K"},
    {3, 1, 0, 255, 0, 0, 1, 8, "RAZA", "Relative Azimuth Angle", "deg"},
    {2, 0, 0, 255, 7, 1, 0, 204, "RCQ", "Humidity parameter in canopy conductance", "Fraction"},
    {2, 1, 0, 255, 0, 0, 0, 20, "RCQ", "Humidity parameter in canopy conductance", "Proportion"},
    {2, 0, 0, 255, 7, 1, 0, 202, "RCS", "Solar parameter in canopy conductance", "Fraction"},
    {2, 1, 0, 255, 0, 0, 0, 18, "RCS", "Solar parameter in canopy conductance", "Proportion"},
    {2, 0, 0, 255, 7, 1, 0, 205, "RCSOL", "Soil moisture parameter in canopy conductance", "Fraction"},
    {2, 1, 0, 255, 0, 0, 0, 21, "RCSOL", "Soil moisture parameter in canopy conductance", "Proportion"},
    {2, 0, 0, 255, 7, 1, 0, 203, "RCT", "Temperature parameter in canopy conductance", "Fraction"},
    {2, 1, 0, 255, 0, 0, 0, 19, "RCT", "Temperature parameter in canopy", "Proportion"},
    {2, 0, 0, 255, 7, 1, 0, 206, "RDRIP", "Rate of water dropping from canopy to ground", "-"},
    {0, 1, 0, 255, 0, 0, 15, 6, "RDSP1", "Radar Spectra (1)", "-"},
    {0, 1, 0, 255, 0, 0, 15, 7, "RDSP2", "Radar Spectra (2)", "-"},
    {0, 1, 0, 255, 0, 0, 15, 8, "RDSP3", "Radar Spectra (3)", "-"},
    {2, 1, 0, 255, 0, 0, 0, 32, "RDVEG", "Root Depth of Vegetation", "m"},
    {0, 0, 0, 255, 7, 1, 16, 196, "REFC", "Composite reflectivity", "dB"},
    {0, 1, 0, 255, 0, 0, 16, 5, "REFC", "Composite reflectivity", "dB"},
    {0, 0, 0, 255, 7, 1, 16, 195, "REFD", "Reflectivity", "dB"},
    {0, 1, 0, 255, 0, 0, 16, 4, "REFD", "Reflectivity", "dB"},
    {0, 0, 0, 255, 7, 1, 16, 194, "REFZC", "Equivalent radar reflectivity factor for parameterized convection",
     "mm^6/m^3"},
    {0, 1, 0, 255, 0, 0, 16, 2, "REFZC", "Equivalent radar reflectivity factor for parameterized convection",
     "mm^6/m^3"},
    {0, 0, 0, 255, 7, 1, 16, 193, "REFZI", "Equivalent radar reflectivity factor for snow", "mm^6/m^3"},
    {0, 1, 0, 255, 0, 0, 16, 1, "REFZI", "Equivalent radar reflectivity factor for snow", "mm^6/m^3"},
    {0, 0, 0, 255, 7, 1, 16, 192, "REFZR", "Equivalent radar reflectivity factor for rain", "mm^6/m^3"},
    {0, 1, 0, 255, 0, 0, 16, 0, "REFZR", "Equivalent radar reflectivity factor for rain", "mm^6/m^3"},
    {0, 1, 0, 255, 0, 0, 2, 13, "RELD", "Relative Divergence", "1/s"},
    {0, 1, 0, 255, 0, 0, 2, 12, "RELV", "Relative Vorticity", "1/s"},
    {0, 0, 0, 255, 7, 1, 16, 197, "RETOP", "Echo Top", "m"},
    {0, 1, 0, 255, 0, 0, 16, 3, "RETOP", "Echo Top", "m"},
    {0, 0, 0, 255, 7, 1, 0, 194, "REV", "Relative Error Variance", "-"},
    {0, 1, 0, 255, 0, 0, 15, 9, "RFCD", "Reflectivity of Cloud Droplets", "dB"},
    {0, 1, 0, 255, 0, 0, 15, 10, "RFCI", "Reflectivity of Cloud Ice", "dB"},
    {0, 1, 0, 255, 0, 0, 15, 13, "RFGRPL", "Reflectivity of Graupel", "dB"},
    {0, 1, 0, 255, 0, 0, 15, 14, "RFHAIL", "Reflectivity of Hail", "dB"},
    {3, 1, 0, 255, 0, 0, 1, 9, "RFL06", "Reflectance in 0.6 Micron Channel", "%"},
    {3, 1, 0, 255, 0, 0, 1, 10, "RFL08", "Reflectance in 0.8 Micron Channel", "%"},
    {3, 1, 0, 255, 0, 0, 1, 11, "RFL16", "Reflectance in 1.6 Micron Channel", "%"},
    {3, 1, 0, 255, 0, 0, 1, 12, "RFL39", "Reflectance in 3.9 Micron Channel", "%"},
    {0, 1, 0, 255, 0, 0, 15, 12, "RFRAIN", "Reflectivity of Rain", "dB"},
    {0, 1, 0, 255, 0, 0, 15, 11, "RFSNOW", "Reflectivity of Snow", "dB"},
    {0, 1, 0, 255, 0, 0, 1, 1, "RH", "Relative Humidity", "%"},
    {0, 1, 0, 255, 0, 0, 1, 94, "RHICE", "Relative Humidity With Respect to Ice", "%"},
    {0, 0, 0, 255, 7, 1, 1, 242, "RHPW", "Relative Humidity with Respect to Precipitable Water", "%"},
    {0, 1, 0, 255, 0, 0, 1, 93, "RHWATER", "Relative Humidity With Respect to Water", "%"},
    {0, 0, 0, 255, 7, 1, 7, 194, "RI", "Richardson Number", "Numeric"},
    {0, 1, 0, 255, 0, 0, 7, 12, "RI", "Richardson Number", "Numeric"},
    {0, 0, 0, 255, 7, 1, 1, 203, "RIME", "Rime Factor", "non-dim"},
    {0, 1, 0, 255, 0, 0, 1, 44, "RIME", "Rime Factor", "Numeric"},
    {10, 1, 0, 255, 0, 0, 1, 4, "RIPCOP", "Rip Current Occurrence Probability", "%"},
    {2, 0, 0, 255, 7, 1, 3, 193, "RLYRS", "Number of Soil Layers in Root Zone", "non-dim"},
    {2, 1, 0, 255, 0, 0, 3, 6, "RLYRS", "Number of Soil Layers in Root Zone", "Numeric"},
    {0, 1, 0, 255, 0, 0, 1, 65, "RPRATE", "Rain Precipitation Rate", "kg/m^2/s"},
    {2, 0, 0, 255, 7, 1, 0, 200, "RSMIN", "Minimal Stomatal Resistance", "s/m"},
    {2, 1, 0, 255, 0, 0, 0, 16, "RSMIN", "Minimal Stomatal Resistance", "s/m"},
    {1, 1, 0, 255, 0, 0, 0, 2, "RSSC", "Remotely Sensed Snow Cover", "-"},
    {0, 0, 0, 255, 7, 1, 191, 194, "RTSEC", "Seconds prior to initial reference time", "s"},
    {10, 0, 0, 255, 7, 1, 3, 206, "RUNUP", "Total Water Level Increase due to Waves", "m"},
    {1, 1, 0, 255, 0, 0, 0, 11, "RVERSW", "River Storage of Water", "m3"},
    {0, 1, 0, 255, 0, 0, 1, 24, "RWMR", "Rain Mixing Ratio", "kg/kg"},
    {0, 1, 0, 255, 0, 0, 18, 14, "SACON", "Specific Activity Concentration", "Bq/kg"},
    {0, 1, 0, 255, 0, 0, 20, 100, "SADEN", "Surface Area Density (Aerosol)", "1/m"},
    {0, 1, 0, 255, 0, 0, 19, 19, "SALBD", "Snow Albedo", "%"},
    {3, 1, 0, 255, 0, 0, 0, 1, "SALBEDO", "Scaled Albedo", "Numeric"},
    {10, 0, 0, 255, 7, 1, 4, 193, "SALIN", "3-D Salinity", "psu"},
    {1, 1, 0, 255, 0, 0, 2, 12, "SALTIL", "Salinity", "kg/kg"},
    {10, 1, 0, 255, 0, 0, 4, 3, "SALTY", "Salinity", "kg/kg"},
    {0, 1, 0, 255, 0, 0, 1, 5, "SATD", "Saturation Deficit", "Pa"},
    {2, 1, 0, 255, 0, 0, 3, 17, "SATOSM", "Saturation Of Soil Moisture", "kg/m^3"},
    {3, 0, 0, 255, 7, 1, 192, 4, "SBC123", "Simulated Brightness Counts for GOES 12, Channel 3", "Byte"},
    {3, 0, 0, 255, 7, 1, 192, 5, "SBC124", "Simulated Brightness Counts for GOES 12, Channel 4", "Byte"},
    {0, 0, 0, 255, 7, 1, 19, 211, "SBSALB", "Visible, Black Sky Albedo", "%"},
    {0, 0, 0, 255, 7, 1, 1, 212, "SBSNO", "Sublimation (evaporation from snow)", "W/m^2"},
    {3, 0, 0, 255, 7, 1, 192, 6, "SBT112", "Simulated Brightness Temperature for GOES 11, Channel 2", "K"},
    {3, 0, 0, 255, 7, 1, 192, 7, "SBT113", "Simulated Brightness Temperature for GOES 11, Channel 3", "K"},
    {3, 0, 0, 255, 7, 1, 192, 8, "SBT114", "Simulated Brightness Temperature for GOES 11, Channel 4", "K"},
    {3, 0, 0, 255, 7, 1, 192, 9, "SBT115", "Simulated Brightness Temperature for GOES 11, Channel 5", "K"},
    {3, 0, 0, 255, 7, 1, 192, 0, "SBT122", "Simulated Brightness Temperature for GOES 12, Channel 2", "K"},
    {3, 0, 0, 255, 7, 1, 192, 1, "SBT123", "Simulated Brightness Temperature for GOES 12, Channel 3", "K"},
    {3, 0, 0, 255, 7, 1, 192, 2, "SBT124", "Simulated Brightness Temperature for GOES 12, Channel 4", "K"},
    {3, 0, 0, 255, 7, 1, 192, 3, "SBT126", "Simulated Brightness Temperature for GOES 12, Channel 6", "K"},
    {3, 0, 0, 255, 7, 1, 192, 23, "SBTA1610", "Simulated Brightness Temperature for ABI GOES-16, Band-10", "K"},
    {3, 0, 0, 255, 7, 1, 192, 24, "SBTA1611", "Simulated Brightness Temperature for ABI GOES-16, Band-11", "K"},
    {3, 0, 0, 255, 7, 1, 192, 25, "SBTA1612", "Simulated Brightness Temperature for ABI GOES-16, Band-12", "K"},
    {3, 0, 0, 255, 7, 1, 192, 26, "SBTA1613", "Simulated Brightness Temperature for ABI GOES-16, Band-13", "K"},
    {3, 0, 0, 255, 7, 1, 192, 27, "SBTA1614", "Simulated Brightness Temperature for ABI GOES-16, Band-14", "K"},
    {3, 0, 0, 255, 7, 1, 192, 28, "SBTA1615", "Simulated Brightness Temperature for ABI GOES-16, Band-15", "K"},
    {3, 0, 0, 255, 7, 1, 192, 29, "SBTA1616", "Simulated Brightness Temperature for ABI GOES-16, Band-16", "K"},
    {3, 0, 0, 255, 7, 1, 192, 20, "SBTA167", "Simulated Brightness Temperature for ABI GOES-16, Band-7", "K"},
    {3, 0, 0, 255, 7, 1, 192, 21, "SBTA168", "Simulated Brightness Temperature for ABI GOES-16, Band-8", "K"},
    {3, 0, 0, 255, 7, 1, 192, 22, "SBTA169", "Simulated Brightness Temperature for ABI GOES-16, Band-9", "K"},
    {3, 0, 0, 255, 7, 1, 192, 39, "SBTA1710", "Simulated Brightness Temperature for ABI GOES-17, Band-10", "K"},
    {3, 0, 0, 255, 7, 1, 192, 40, "SBTA1711", "Simulated Brightness Temperature for ABI GOES-17, Band-11", "K"},
    {3, 0, 0, 255, 7, 1, 192, 41, "SBTA1712", "Simulated Brightness Temperature for ABI GOES-17, Band-12", "K"},
    {3, 0, 0, 255, 7, 1, 192, 42, "SBTA1713", "Simulated Brightness Temperature for ABI GOES-17, Band-13", "K"},
    {3, 0, 0, 255, 7, 1, 192, 43, "SBTA1714", "Simulated Brightness Temperature for ABI GOES-17, Band-14", "K"},
    {3, 0, 0, 255, 7, 1, 192, 44, "SBTA1715", "Simulated Brightness Temperature for ABI GOES-17, Band-15", "K"},
    {3, 0, 0, 255, 7, 1, 192, 45, "SBTA1716", "Simulated Brightness Temperature for ABI GOES-17, Band-16", "K"},
    {3, 0, 0, 255, 7, 1, 192, 36, "SBTA177", "Simulated Brightness Temperature for ABI GOES-17, Band-7", "K"},
    {3, 0, 0, 255, 7, 1, 192, 37, "SBTA178", "Simulated Brightness Temperature for ABI GOES-17, Band-8", "K"},
    {3, 0, 0, 255, 7, 1, 192, 38, "SBTA179", "Simulated Brightness Temperature for ABI GOES-17, Band-9", "K"},
    {3, 0, 0, 255, 7, 1, 192, 55, "SBTAGR10", "Simulated Brightness Temperature for nadir ABI GOES-R, Band-10", "-"},
    {3, 0, 0, 255, 7, 1, 192, 56, "SBTAGR11", "Simulated Brightness Temperature for nadir ABI GOES-R, Band-11", "-"},
    {3, 0, 0, 255, 7, 1, 192, 57, "SBTAGR12", "Simulated Brightness Temperature for nadir ABI GOES-R, Band-12", "-"},
    {3, 0, 0, 255, 7, 1, 192, 58, "SBTAGR13", "Simulated Brightness Temperature for nadir ABI GOES-R, Band-13", "-"},
    {3, 0, 0, 255, 7, 1, 192, 59, "SBTAGR14", "Simulated Brightness Temperature for nadir ABI GOES-R, Band-14", "-"},
    {3, 0, 0, 255, 7, 1, 192, 60, "SBTAGR15", "Simulated Brightness Temperature for nadir ABI GOES-R, Band-15", "-"},
    {3, 0, 0, 255, 7, 1, 192, 61, "SBTAGR16", "Simulated Brightness Temperature for nadir ABI GOES-R, Band-16", "-"},
    {3, 0, 0, 255, 7, 1, 192, 52, "SBTAGR7", "Simulated Brightness Temperature for nadir ABI GOES-R, Band-7", "-"},
    {3, 0, 0, 255, 7, 1, 192, 53, "SBTAGR8", "Simulated Brightness Temperature for nadir ABI GOES-R, Band-8", "-"},
    {3, 0, 0, 255, 7, 1, 192, 54, "SBTAGR9", "Simulated Brightness Temperature for nadir ABI GOES-R, Band-9", "-"},
    {3, 1, 0, 255, 0, 0, 0, 2, "SBTMP", "Scaled Brightness Temperature", "Numeric"},
    {4, 1, 0, 255, 0, 0, 2, 10, "SCINT", "Scintillation", "Numeric"},
    {0, 1, 0, 255, 0, 0, 1, 84, "SCLIWC", "Specific Cloud Ice Water Content", "kg/kg"},
    {0, 1, 0, 255, 0, 0, 1, 83, "SCLLWC", "Specific Cloud Liquid Water Content", "kg/kg"},
    {3, 1, 0, 255, 0, 0, 1, 29, "SCRAD", "Scaled Radiance", "Numeric"},
    {0, 1, 0, 255, 0, 0, 20, 112, "SCTAOTK", "Scattering Aerosol Optical Thickness", "Numeric"},
    {3, 1, 0, 255, 0, 0, 0, 5, "SCTPRES", "Scaled Cloud Top Pressure", "Numeric"},
    {0, 1, 0, 255, 0, 0, 20, 6, "SDDMFLX", "Surface Dry Deposition Mass Flux", "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 1, 61, "SDEN", "Snow Density", "kg/m^3"},
    {3, 1, 0, 255, 0, 0, 1, 99, "SDMPEMRR",
     "Standard deviation between MPE rain rates for the co-located IR data and the microwave data rain rates",
     "Numeric"},
    {0, 1, 0, 255, 0, 0, 3, 20, "SDSGSO", "Standard Deviation of Sub-Grid Scale Orography", "m"},
    {0, 1, 0, 255, 0, 0, 1, 60, "SDWE", "Snow Depth Water Equivalent", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 20, 11, "SEDMFLX", "Sedimentation Mass Flux", "kg/m^2/s"},
    {1, 1, 0, 255, 0, 0, 2, 3, "SEDTK", "Sediment Thickness", "m"},
    {1, 1, 0, 255, 0, 0, 2, 4, "SEDTMP", "Sediment Temperature", "K"},
    {10, 0, 0, 255, 7, 1, 3, 207, "SETUP", "Mean Increase in Water Level due to Waves", "m"},
    {0, 1, 0, 255, 0, 0, 1, 62, "SEVAP", "Snow Evaporation", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 20, 77, "SFCEFLX", "Surface Emission flux", "kg/m^2/s"},
    {2, 1, 0, 255, 0, 0, 0, 1, "SFCR", "Surface Roughness", "m"},
    {2, 0, 0, 255, 7, 1, 0, 216, "SFCRH", "Roughness length for heat", "m"},
    {2, 1, 0, 255, 0, 0, 0, 34, "SFCWRO", "Surface Water Runoff", "kg/m^2"},
    {2, 0, 0, 255, 7, 1, 0, 195, "SFEXC", "Exchange Coefficient", "(kg/m^3)(m/s)"},
    {2, 1, 0, 255, 0, 0, 0, 12, "SFEXC", "Exchange Coefficient", "kg/m^2/s"},
    {1, 1, 0, 255, 0, 0, 0, 10, "SFLORC", "Side Flow into River Channel", "m3s-1m-1"},
    {0, 1, 0, 255, 0, 0, 20, 55, "SFLUX", "Surface Flux", "mol/m^2/s"},
    {1, 1, 0, 255, 0, 0, 2, 9, "SFSAL", "Shape Factor with Respect to Salinity Profile", "-"},
    {10, 1, 0, 255, 0, 0, 4, 11, "SFSALP", "Shape Factor With Respect To Salinity Profile", "-"},
    {1, 1, 0, 255, 0, 0, 2, 10, "SFTMP", "Shape Factor with Respect to Temperature Profile in Thermocline", "-"},
    {10, 1, 0, 255, 0, 0, 4, 12, "SFTMPP", "Shape Factor With Respect To Temperature Profile In Thermocline", "-"},
    {0, 1, 0, 255, 0, 0, 2, 7, "SGCVV", "Sigma Coordinate Vertical Velocity", "1/s"},
    {0, 0, 0, 255, 7, 1, 0, 201, "SHAHR", "Shallow Convective Heating Rate", "K/s"},
    {0, 0, 0, 255, 7, 1, 19, 201, "SHAILPRO", "Significant Hail probability", "%"},
    {0, 0, 0, 255, 7, 1, 1, 214, "SHAMR", "Shallow Convective Moistening Rate", "kg/kg/s"},
    {2, 1, 0, 255, 0, 0, 3, 26, "SHFLX", "Soil Heat Flux", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 0, 11, "SHTFL", "Sensible Heat Net Flux", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 1, 108, "SHTPRM", "Specific Humidity Tendency due to Parameterizations", "kg/kg/s"},
    {0, 1, 0, 255, 0, 0, 7, 13, "SHWINX", "Showalter Index", "K"},
    {10, 1, 0, 255, 0, 0, 2, 3, "SICED", "Speed of Ice Drift", "m/s"},
    {4, 1, 0, 255, 0, 0, 9, 1, "SIGHAL", "Hall Conductivity", "S/m"},
    {4, 1, 0, 255, 0, 0, 9, 2, "SIGPAR", "Parallel Conductivity", "S/m"},
    {4, 1, 0, 255, 0, 0, 9, 0, "SIGPED", "Pedersen Conductivity", "S/m"},
    {0, 0, 0, 255, 7, 1, 19, 217, "SIPD", "Supercooled Large Droplet (SLD) Icingsee Note 2", "-"},
    {0, 1, 0, 255, 0, 0, 0, 17, "SKINT", "Skin Temperature", "K"},
    {3, 1, 0, 255, 0, 0, 5, 1, "SKSSTMP", "Skin Sea-Surface Temperature", "K"},
    {0, 0, 0, 255, 7, 1, 1, 230, "SLACC", "Sleet Accumulation", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 19, 23, "SLDP", "Supercooled Large Droplet Probabilitysee Note 2", "%"},
    {3, 1, 0, 255, 0, 0, 0, 4, "SLFTI", "Scaled Lifted Index", "Numeric"},
    {10, 0, 0, 255, 7, 1, 3, 202, "SLTFL", "Salt Flux", "kg/m^2/s"},
    {2, 0, 0, 255, 7, 1, 3, 194, "SLTYP", "Surface Slope Type", "Index"},
    {0, 1, 0, 255, 0, 0, 6, 34, "SLWTC", "Surface Long Wave Effective Total Cloudiness", "Numeric"},
    {2, 0, 0, 255, 7, 1, 3, 196, "SMDRY", "Direct Evaporation Cease (soil moisture)", "Proportion"},
    {2, 1, 0, 255, 0, 0, 3, 8, "SMDRY", "Direct Evaporation Cease (soil moisture)", "Proportion"},
    {0, 1, 0, 255, 0, 0, 1, 113, "SMLWGMA",
     "Specific Mass of Liquid Water Coating on Graupel Expressed as Mass of Liquid Water per Unit Mass of Moist Air",
     "kg/kg"},
    {0, 1, 0, 255, 0, 0, 1, 110, "SMLWHMA",
     "Specific Mass of Liquid Water Coating on Hail Expressed as Mass of Liquid Water per Unit Mass of Moist Air",
     "kg/kg"},
    {0, 1, 0, 255, 0, 0, 1, 116, "SMLWSMA",
     "Specific Mass of Liquid Water Coating on Snow Expressed as Mass of Liquid Water per Unit Mass of Moist Air",
     "kg/kg"},
    {2, 0, 0, 255, 7, 1, 3, 195, "SMREF", "Transpiration Stress-onset (soil moisture)", "Proportion"},
    {2, 1, 0, 255, 0, 0, 3, 7, "SMREF", "Transpiration Stress-onset (soil moisture)", "Proportion"},
    {0, 0, 0, 255, 7, 1, 19, 193, "SNFALB", "Snow-Free Albedo", "%"},
    {0, 1, 0, 255, 0, 0, 19, 18, "SNFALB", "Snow-Free Albedo", "%"},
    {0, 1, 0, 255, 0, 0, 1, 25, "SNMR", "Snow Mixing Ratio", "kg/kg"},
    {0, 1, 0, 255, 0, 0, 1, 17, "SNOAG", "Snow Age", "day"},
    {0, 1, 0, 255, 0, 0, 1, 14, "SNOC", "Convective Snow", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 1, 11, "SNOD", "Snow Depth", "m"},
    {0, 0, 0, 255, 7, 1, 0, 192, "SNOHF", "Snow Phase Change Heat Flux", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 0, 16, "SNOHF", "Snow Phase Change Heat Flux", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 1, 15, "SNOL", "Large-Scale Snow", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 1, 16, "SNOM", "Snow Melt", "kg/m^2"},
    {0, 0, 0, 255, 7, 1, 1, 208, "SNOT", "Snow temperature", "K"},
    {0, 1, 0, 255, 0, 0, 0, 18, "SNOT", "Snow Temperature (top of snow)", "K"},
    {0, 0, 0, 255, 7, 1, 1, 201, "SNOWC", "Snow Cover", "%"},
    {0, 1, 0, 255, 0, 0, 1, 42, "SNOWC", "Snow Cover", "%"},
    {0, 0, 0, 255, 7, 1, 1, 233, "SNOWLR", "Snow Liquid ratio", "kg/kg"},
    {0, 0, 0, 255, 7, 1, 19, 236, "SNOWLVL", "Snow Level", "m"},
    {0, 0, 0, 255, 7, 1, 1, 222, "SNOWT", "Snow temperature, depth-avg", "K"},
    {2, 1, 0, 255, 0, 0, 3, 25, "SNWDEB", "Snow Depth at Elevation Bands", "kg/m^2"},
    {2, 1, 0, 255, 0, 0, 3, 27, "SOILDEP", "Soil Depth", "m"},
    {2, 1, 0, 255, 0, 0, 3, 21, "SOILICE", "Soil Ice", "kg/m^3"},
    {2, 0, 0, 255, 7, 1, 3, 192, "SOILL", "Liquid Volumetric Soil Moisture (non Frozen)", "Proportion"},
    {2, 1, 0, 255, 0, 0, 3, 5, "SOILL", "Liquid Volumetric Soil Moisture (non-frozen)", "Proportion"},
    {2, 1, 0, 255, 0, 0, 0, 3, "SOILM", "Soil Moisture Content", "kg/m^2"},
    {2, 1, 0, 255, 0, 0, 3, 19, "SOILMOI", "Soil Moisture", "kg/m^3"},
    {2, 1, 0, 255, 0, 0, 3, 15, "SOILP", "Soil Porosity", "m^3/m^3"},
    {2, 1, 0, 255, 0, 0, 3, 18, "SOILTMP", "Soil Temperature", "K"},
    {2, 1, 0, 255, 0, 0, 0, 38, "SOILVIC", "Soil Volumetric Ice Content (Water Equivalent)", "m^3/m^3"},
    {2, 0, 0, 255, 7, 1, 0, 192, "SOILW", "Volumetric Soil Moisture Content", "Fraction"},
    {2, 1, 0, 255, 0, 0, 0, 9, "SOILW", "Volumetric Soil Moisture Content", "Proportion"},
    {2, 1, 0, 255, 0, 0, 0, 22, "SOIL_M", "Soil Moisture", "kg/m^3"},
    {4, 1, 0, 255, 0, 0, 6, 6, "SOLRF", "Solar Radio Emissions", "W/m^2/Hz"},
    {3, 1, 0, 255, 0, 0, 1, 7, "SOLZA", "Solar Zenith Angle", "deg"},
    {2, 1, 0, 255, 0, 0, 3, 0, "SOTYP", "Soil Type", "-"},
    {3, 1, 0, 255, 0, 0, 1, 28, "SPBRT", "Brightness Temperture", "K"},
    {10, 1, 0, 255, 0, 0, 1, 1, "SPC", "Current Speed", "m/s"},
    {4, 1, 0, 255, 0, 0, 6, 4, "SPECIRR", "Solar Spectral Irradiance", "W/m^2/m^-9"},
    {4, 1, 0, 255, 0, 0, 1, 0, "SPEED", "Velocity Magnitude (Speed)", "m/s"},
    {0, 1, 0, 255, 0, 0, 1, 0, "SPFH", "Specific Humidity", "kg/kg"},
    {10, 1, 0, 255, 0, 0, 0, 45, "SPFTR", "Spectral Peakedness Factor", "1/s"},
    {0, 1, 0, 255, 0, 0, 1, 102, "SPNCG", "Specific Number Concentration of Graupel", "1/kg"},
    {0, 1, 0, 255, 0, 0, 1, 103, "SPNCH", "Specific Number Concentration of Hail", "1/kg"},
    {0, 1, 0, 255, 0, 0, 1, 100, "SPNCR", "Specific Number Concentration of Rain", "1/kg"},
    {0, 1, 0, 255, 0, 0, 1, 101, "SPNCS", "Specific Number Concentration of Snow", "1/kg"},
    {0, 1, 0, 255, 0, 0, 1, 66, "SPRATE", "Snow Precipitation Rate", "kg/m^2/s"},
    {4, 1, 0, 255, 0, 0, 2, 7, "SPRDF", "Spread F", "m"},
    {3, 1, 0, 255, 0, 0, 0, 3, "SPWAT", "Scaled Precipitable Water", "Numeric"},
    {3, 1, 0, 255, 0, 0, 0, 0, "SRAD", "Scaled Radiance", "Numeric"},
    {0, 1, 0, 255, 0, 0, 1, 85, "SRAINW", "Specific Rain Water Content", "kg/kg"},
    {0, 0, 0, 255, 7, 1, 19, 194, "SRCONO", "Slight risk convective outlook", "Categorical"},
    {3, 0, 0, 255, 7, 1, 192, 14, "SRFA161", "Simulated Reflectance Factor for ABI GOES-16, Band-1", "-"},
    {3, 0, 0, 255, 7, 1, 192, 15, "SRFA162", "Simulated Reflectance Factor for ABI GOES-16, Band-2", "-"},
    {3, 0, 0, 255, 7, 1, 192, 16, "SRFA163", "Simulated Reflectance Factor for ABI GOES-16, Band-3", "-"},
    {3, 0, 0, 255, 7, 1, 192, 17, "SRFA164", "Simulated Reflectance Factor for ABI GOES-16, Band-4", "-"},
    {3, 0, 0, 255, 7, 1, 192, 18, "SRFA165", "Simulated Reflectance Factor for ABI GOES-16, Band-5", "-"},
    {3, 0, 0, 255, 7, 1, 192, 19, "SRFA166", "Simulated Reflectance Factor for ABI GOES-16, Band-6", "-"},
    {3, 0, 0, 255, 7, 1, 192, 30, "SRFA171", "Simulated Reflectance Factor for ABI GOES-17, Band-1", "-"},
    {3, 0, 0, 255, 7, 1, 192, 31, "SRFA172", "Simulated Reflectance Factor for ABI GOES-17, Band-2", "-"},
    {3, 0, 0, 255, 7, 1, 192, 32, "SRFA173", "Simulated Reflectance Factor for ABI GOES-17, Band-3", "-"},
    {3, 0, 0, 255, 7, 1, 192, 33, "SRFA174", "Simulated Reflectance Factor for ABI GOES-17, Band-4", "-"},
    {3, 0, 0, 255, 7, 1, 192, 34, "SRFA175", "Simulated Reflectance Factor for ABI GOES-17, Band-5", "-"},
    {3, 0, 0, 255, 7, 1, 192, 35, "SRFA176", "Simulated Reflectance Factor for ABI GOES-17, Band-6", "-"},
    {3, 0, 0, 255, 7, 1, 192, 46, "SRFAGR1", "Simulated Reflectance Factor for nadir ABI GOES-R, Band-1", "-"},
    {3, 0, 0, 255, 7, 1, 192, 47, "SRFAGR2", "Simulated Reflectance Factor for nadir ABI GOES-R, Band-2", "-"},
    {3, 0, 0, 255, 7, 1, 192, 48, "SRFAGR3", "Simulated Reflectance Factor for nadir ABI GOES-R, Band-3", "-"},
    {3, 0, 0, 255, 7, 1, 192, 49, "SRFAGR4", "Simulated Reflectance Factor for nadir ABI GOES-R, Band-4", "-"},
    {3, 0, 0, 255, 7, 1, 192, 50, "SRFAGR5", "Simulated Reflectance Factor for nadir ABI GOES-R, Band-5", "-"},
    {3, 0, 0, 255, 7, 1, 192, 51, "SRFAGR6", "Simulated Reflectance Factor for nadir ABI GOES-R, Band-6", "-"},
    {0, 1, 0, 255, 0, 0, 1, 12, "SRWEQ", "Snowfall Rate Water Equivalent", "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 20, 103, "SSALBK", "Single Scattering Albedo", "Numeric"},
    {0, 1, 0, 255, 0, 0, 3, 22, "SSGSO", "Slope of Sub-Grid Scale Orography", "Numeric"},
    {10, 0, 0, 255, 7, 1, 3, 195, "SSHG", "Sea Surface Height Relative to Geoid", "m"},
    {3, 1, 0, 255, 0, 0, 5, 2, "SSKSSTMP", "Sub-Skin Sea-Surface Temperature", "K"},
    {0, 1, 0, 255, 0, 0, 1, 86, "SSNOWW", "Specific Snow Water Content", "kg/kg"},
    {1, 0, 0, 255, 7, 1, 0, 193, "SSRUN", "Storm Surface Runoff", "kg/m^2"},
    {1, 1, 0, 255, 0, 0, 0, 6, "SSRUN", "Storm Surface Runoff", "kg/m^2"},
    {10, 0, 0, 255, 7, 1, 3, 200, "SSST", "Surface Salinity Trend", "psu/day"},
    {3, 1, 0, 255, 0, 0, 0, 6, "SSTMP", "Scaled Skin Temperature", "Numeric"},
    {2, 0, 0, 255, 7, 1, 0, 211, "SSTOR", "Surface water storage", "kg/m^2"},
    {10, 0, 0, 255, 7, 1, 3, 199, "SSTT", "Surface Temperature Trend", "K/day"},
    {0, 1, 0, 255, 0, 0, 6, 35, "SSWTC", "Surface Short Wave Effective Total Cloudiness", "Numeric"},
    {0, 0, 0, 255, 7, 1, 19, 200, "STORPROB", "Significant Tornado probability", "%"},
    {0, 1, 0, 255, 0, 0, 2, 4, "STRM", "Stream Function", "m^2/s"},
    {0, 1, 0, 255, 0, 0, 1, 87, "STRPRATE", "Stratiform Precipitation Rate", "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 6, 24, "SUNS", "Sunshine", "Numeric"},
    {0, 0, 0, 255, 7, 1, 6, 201, "SUNSD", "Sunshine Duration", "s"},
    {0, 1, 0, 255, 0, 0, 6, 33, "SUNSD", "Sunshine Duration", "s"},
    {10, 0, 0, 255, 7, 1, 3, 192, "SURGE", "Hurricane Storm Surge", "m"},
    {0, 0, 0, 255, 7, 1, 19, 220, "SVRTS", "Categorical Severe Thunderstorm", "-"},
    {10, 0, 0, 255, 7, 1, 3, 208, "SWASH", "Time-varying Increase in Water Level due to Waves", "m"},
    {0, 1, 0, 255, 0, 0, 4, 2, "SWAVR", "Short-Wave Radiation Flux", "W/m^2"},
    {10, 1, 0, 255, 0, 0, 0, 7, "SWDIR", "Direction of Swell Waves", "deg"},
    {0, 1, 0, 255, 0, 0, 20, 7, "SWDMFLX", "Surface Wet Deposition Mass Flux", "kg/m^2/s"},
    {10, 1, 0, 255, 0, 0, 0, 8, "SWELL", "Significant Height of Swell Waves", "m"},
    {1, 1, 0, 255, 0, 0, 0, 4, "SWEPON", "Snow Water Equivalent Percent of Normal", "%"},
    {10, 1, 0, 255, 0, 0, 0, 47, "SWHFSWEL", "Significant wave height of first swell partition", "m"},
    {0, 0, 0, 255, 7, 1, 4, 197, "SWHR", "Solar Radiative Heating Rate", "K/s"},
    {10, 1, 0, 255, 0, 0, 0, 48, "SWHSSWEL", "Significant wave height of second swell partition", "m"},
    {10, 1, 0, 255, 0, 0, 0, 49, "SWHTSWEL", "Significant wave height of third swell partition", "m"},
    {0, 0, 0, 255, 7, 1, 19, 202, "SWINDPRO", "Significant Wind probability", "%"},
    {10, 1, 0, 255, 0, 0, 0, 9, "SWPER", "Mean Period of Swell Waves", "s"},
    {3, 0, 0, 255, 7, 1, 1, 194, "SWQI", "Scatterometer Wind Quality", "-"},
    {0, 1, 0, 255, 0, 0, 4, 6, "SWRAD", "Radiance (with respect to wavelength)", "W/m^3/sr"},
    {0, 0, 0, 255, 7, 1, 19, 212, "SWSALB", "Visible, White Sky Albedo", "%"},
    {0, 1, 0, 255, 0, 0, 7, 5, "SX", "Sweat Index", "Numeric"},
    {0, 1, 0, 255, 0, 0, 6, 1, "TCDC", "Total Cloud Cover", "%"},
    {0, 0, 0, 255, 7, 1, 0, 204, "TCHP", "Tropical Cyclone Heat Potential", "J/m^2*K"},
    {0, 1, 0, 255, 0, 0, 1, 81, "TCICON", "Total Column-Integrate Condensate", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 14, 2, "TCIOZ", "Total Column Integrated Ozone", "DU"},
    {0, 1, 0, 255, 0, 0, 1, 64, "TCIWV", "Total Column Integrated Water Vapour", "kg/m^2"},
    {2, 1, 0, 255, 0, 0, 0, 35, "TCLASS", "Tile Class", "-"},
    {0, 0, 0, 255, 7, 1, 1, 209, "TCLSW", "Total column-integrated supercooled liquid water", "kg/m^2"},
    {0, 0, 0, 255, 7, 1, 6, 198, "TCOLC", "Total Column-Integrated Condensate", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 6, 20, "TCOLC", "Total Column-Integrated Condensate", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 1, 74, "TCOLG", "Total Column Integrate Graupel", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 1, 72, "TCOLH", "Total Column Integrate Hail", "kg/m^2"},
    {0, 0, 0, 255, 7, 1, 6, 197, "TCOLI", "Total Column-Integrated Cloud Ice", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 1, 70, "TCOLI", "Total Column Integrate Cloud Ice", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 6, 19, "TCOLIold", "Total Column-Integrated Cloud Ice", "kg/m^2"},
    {0, 0, 0, 255, 7, 1, 1, 210, "TCOLM", "Total column-integrated melting ice", "kg/m^2"},
    {0, 0, 0, 255, 7, 1, 1, 204, "TCOLR", "Total Column Integrated Rain", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 1, 45, "TCOLR", "Total Column Integrated Rain", "kg/m^2"},
    {0, 0, 0, 255, 7, 1, 1, 205, "TCOLS", "Total Column Integrated Snow", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 1, 46, "TCOLS", "Total Column Integrated Snow", "kg/m^2"},
    {0, 0, 0, 255, 7, 1, 6, 196, "TCOLW", "Total Column-Integrated Cloud Water", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 1, 69, "TCOLW", "Total Column Integrate Cloud Water", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 1, 78, "TCOLWA", "Total Column Integrate Water (All components including precipitation)",
     "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 6, 18, "TCOLWold", "Total Column-Integrated Cloud Water", "kg/m^2"},
    {0, 0, 0, 255, 7, 1, 6, 195, "TCOND", "Total Condensate", "kg/kg"},
    {0, 1, 0, 255, 0, 0, 1, 21, "TCOND", "Condensate", "kg/kg"},
    {0, 1, 0, 255, 0, 0, 6, 17, "TCONDold", "Total Condensate", "kg/kg"},
    {10, 0, 0, 255, 7, 1, 3, 242, "TCSRG20", "20% Tropical Cyclone Storm Surge Exceedance", "m"},
    {10, 0, 0, 255, 7, 1, 3, 243, "TCSRG30", "30% Tropical Cyclone Storm Surge Exceedance", "m"},
    {10, 0, 0, 255, 7, 1, 3, 244, "TCSRG40", "40% Tropical Cyclone Storm Surge Exceedance", "m"},
    {10, 0, 0, 255, 7, 1, 3, 245, "TCSRG50", "50% Tropical Cyclone Storm Surge Exceedance", "m"},
    {10, 0, 0, 255, 7, 1, 3, 246, "TCSRG60", "60% Tropical Cyclone Storm Surge Exceedance", "m"},
    {10, 0, 0, 255, 7, 1, 3, 247, "TCSRG70", "70% Tropical Cyclone Storm Surge Exceedance", "m"},
    {10, 0, 0, 255, 7, 1, 3, 248, "TCSRG80", "80% Tropical Cyclone Storm Surge Exceedance", "m"},
    {10, 0, 0, 255, 7, 1, 3, 249, "TCSRG90", "90% Tropical Cyclone Storm Surge Exceedance", "m"},
    {0, 1, 0, 255, 0, 0, 1, 51, "TCWAT",
     "Total Column Water (Vertically integrated total water (vapour+cloud water/ice)", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 0, 20, "TDCHT", "Turbulent Diffusion Coefficient for Heat", "m^2/s"},
    {0, 1, 0, 255, 0, 0, 2, 31, "TDCMOM", "Turbulent Diffusion Coefficient for Momentum", "m^2/s"},
    {2, 1, 0, 255, 0, 0, 0, 36, "TFRCT", "Tile Fraction", "Proportion"},
    {0, 0, 0, 255, 7, 1, 0, 197, "THFLX", "Total Downward Heat Flux at Surface", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 3, 12, "THICK", "Thickness", "m"},
    {0, 1, 0, 255, 0, 0, 6, 10, "THUNC", "Thunderstorm Coverage", "-"},
    {0, 0, 0, 255, 7, 1, 0, 203, "THZ0", "Potential Temperature at Top of Viscous Sublayer", "K"},
    {0, 1, 0, 255, 0, 0, 18, 6, "TIACCP", "Time Integrated Air Concentration of Cesium Pollutant", "Bqs/m^3"},
    {0, 1, 0, 255, 0, 0, 18, 7, "TIACIP", "Time Integrated Air Concentration of Iodine Pollutant", "Bqs/m^3"},
    {0, 1, 0, 255, 0, 0, 18, 8, "TIACRP", "Time Integrated Air Concentration of Radioactive Pollutant", "Bqs/m^3"},
    {10, 0, 0, 255, 7, 1, 3, 251, "TIDE", "Tide", "m"},
    {0, 0, 0, 255, 7, 1, 1, 206, "TIPD", "Total Icing Potential Diagnostic", "non-dim"},
    {0, 1, 0, 255, 0, 0, 19, 11, "TKE", "Turbulent Kinetic Energy", "J/kg"},
    {0, 1, 0, 255, 0, 0, 1, 90, "TKMFLX", "Total Kinematic Moisture Flux", "kg/kg*m/s"},
    {0, 1, 0, 255, 0, 0, 17, 4, "TLGTFD", "Total Lightning Flash Density", "km-2day-1"},
    {0, 1, 0, 255, 0, 0, 0, 4, "TMAX", "Maximum Temperature", "K"},
    {0, 1, 0, 255, 0, 0, 6, 9, "TMAXT", "Thunderstorm Maximum Tops", "m"},
    {0, 1, 0, 255, 0, 0, 0, 5, "TMIN", "Minimum Temperature", "K"},
    {0, 1, 0, 255, 0, 0, 0, 0, "TMP", "Temperature", "K"},
    {0, 1, 0, 255, 0, 0, 0, 9, "TMPA", "Temperature Anomaly", "K"},
    {0, 1, 0, 255, 0, 0, 0, 29, "TMPADV", "Temperature Advection", "K/s"},
    {4, 1, 0, 255, 0, 0, 0, 0, "TMPSWP", "Temperature", "K"},
    {0, 0, 0, 255, 7, 1, 2, 227, "TOA10", "Earliest Reasonable Arrival Time (10% exceedance)", "s"},
    {0, 0, 0, 255, 7, 1, 2, 228, "TOA50", "Most Likely Arrival Time (50% exceedance)", "s"},
    {0, 0, 0, 255, 7, 1, 2, 229, "TOD50", "Most Likely Departure Time (50% exceedance)", "s"},
    {0, 0, 0, 255, 7, 1, 2, 230, "TOD90", "Latest Reasonable Departure Time (90% exceedance)", "s"},
    {0, 0, 0, 255, 7, 1, 19, 197, "TORPROB", "Tornado probability", "%"},
    {0, 1, 0, 255, 0, 0, 7, 4, "TOTALX", "Total Totals Index", "K"},
    {0, 1, 0, 255, 0, 0, 1, 80, "TOTCON", "Total Condensate", "kg/kg"},
    {0, 1, 0, 255, 0, 0, 18, 13, "TOTLWD", "Total Deposition (Wet + Dry)", "Bq/m^2"},
    {0, 0, 0, 255, 7, 1, 14, 197, "TOZ", "Ozone Tendency", "kg/kg/s"},
    {0, 1, 0, 255, 0, 0, 14, 0, "TOZNE", "Total Ozone", "DU"},
    {2, 1, 0, 255, 0, 0, 0, 37, "TPERCT", "Tile Percentage", "%"},
    {0, 0, 0, 255, 7, 1, 19, 219, "TPFI", "Turbulence Potential Forecast Index", "-"},
    {0, 1, 0, 255, 0, 0, 1, 52, "TPRATE", "Total Precipitation Rate", "kg/m^2/s"},
    {0, 0, 0, 255, 7, 1, 2, 231, "TPWDIR", "Tropical Wind Direction", "deg"},
    {0, 0, 0, 255, 7, 1, 2, 232, "TPWSPD", "Tropical Wind Speed", "m/s"},
    {0, 1, 0, 255, 0, 0, 20, 13, "TRANHH", "Transfer From Hydrophobic to Hydrophilic", "kg/kg/s"},
    {2, 0, 0, 255, 7, 1, 0, 230, "TRANS", "Transpiration", "W/m^2"},
    {2, 1, 0, 255, 0, 0, 3, 12, "TRANSO", "Transpiration Stree-Onset(Soil Moisture)", "kg/m^3"},
    {0, 1, 0, 255, 0, 0, 20, 14, "TRSDS", "Transfer From SO2 (Sulphur Dioxide) to SO4 (Sulphate)", "kg/kg/s"},
    {0, 0, 0, 255, 7, 1, 2, 226, "TRWDIR", "Transport Wind Direction", "deg"},
    {0, 0, 0, 255, 7, 1, 2, 225, "TRWSPD", "Transport Wind Speed", "m/s"},
    {0, 0, 0, 255, 7, 1, 0, 200, "TSD1D", "Standard Dev. of IR Temp. over 1x1 deg. area", "K"},
    {0, 1, 0, 255, 0, 0, 191, 0, "TSEC", "Seconds prior to initial reference time (defined in Section 1)", "s"},
    {4, 1, 0, 255, 0, 0, 6, 0, "TSI", "Integrated Solar Irradiance", "W/m^2"},
    {0, 0, 0, 255, 7, 1, 3, 199, "TSLSA", "3-hr pressure tendency (Std. Atmos. Reduction)", "Pa/s"},
    {0, 0, 0, 255, 7, 1, 1, 241, "TSNOW", "Total Snow", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 1, 50, "TSNOWP", "Total Snow Precipitation", "kg/m^2"},
    {2, 1, 0, 255, 0, 0, 0, 2, "TSOIL", "Soil Temperature", "K"},
    {0, 1, 0, 255, 0, 0, 1, 57, "TSRATE", "Total Snowfall Rate", "m/s"},
    {0, 1, 0, 255, 0, 0, 1, 53, "TSRWE", "Total Snowfall Rate Water Equivalent", "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 19, 2, "TSTM", "Thunderstorm Probability", "%"},
    {0, 0, 0, 255, 7, 1, 19, 203, "TSTMC", "Categorical Thunderstorm", "-"},
    {0, 1, 0, 255, 0, 0, 0, 19, "TTCHT", "Turbulent Transfer Coefficient for Heat", "Numeric"},
    {0, 0, 0, 255, 7, 1, 0, 198, "TTDIA", "Temperature Tendency by All Physics", "K/s"},
    {10, 1, 0, 255, 0, 0, 4, 2, "TTHDP", "Transient Thermocline Depth", "m"},
    {0, 1, 0, 255, 0, 0, 0, 23, "TTLWR", "Temperature Tendency due to Long-Wave Radiation", "K/s"},
    {0, 1, 0, 255, 0, 0, 0, 25, "TTLWRCS", "Temperature Tendency due to Long-Wave Radiation, Clear Sky", "K/s"},
    {0, 1, 0, 255, 0, 0, 0, 26, "TTPARM", "Temperature Tendency due to parameterizations", "K/s"},
    {0, 0, 0, 255, 7, 1, 0, 199, "TTPHY", "Temperature Tendency by Non-radiation Physics", "K/s"},
    {0, 0, 0, 255, 7, 1, 0, 193, "TTRAD", "Temperature Tendency by All Radiation", "K/s"},
    {0, 1, 0, 255, 0, 0, 0, 22, "TTSWR", "Temperature Tendency due to Short-Wave Radiation", "K/s"},
    {0, 1, 0, 255, 0, 0, 0, 24, "TTSWRCS", "Temperature Tendency due to Short-Wave Radiation, Clear Sky", "K/s"},
    {0, 1, 0, 255, 0, 0, 19, 10, "TURB", "Turbulence", "-"},
    {0, 1, 0, 255, 0, 0, 19, 9, "TURBB", "Turbulence Base", "m"},
    {0, 1, 0, 255, 0, 0, 19, 8, "TURBT", "Turbulence Top", "m"},
    {0, 1, 0, 255, 0, 0, 1, 49, "TWATP", "Total Water Precipitation", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 2, 35, "TWIND", "Tangential Wind Component", "m/s"},
    {10, 0, 0, 255, 7, 1, 3, 210, "TWLDC", "Total Water Level Above Dune Crest", "m"},
    {10, 0, 0, 255, 7, 1, 3, 209, "TWLDT", "Total Water Level Above Dune Toe", "m"},
    {10, 0, 0, 255, 7, 1, 3, 205, "TWLWAV", "Total Water Level Accounting for Tide, Wind and Waves", "m"},
    {0, 1, 0, 255, 0, 0, 20, 58, "TYAAL", "Total Yearly Average Atmospheric Loss", "mol/s"},
    {0, 1, 0, 255, 0, 0, 20, 57, "TYABA", "-", "mol"},
    {0, 0, 0, 255, 7, 1, 3, 194, "U-GWD", "Zonal Flux of Gravity Wave Stress", "N/m^2"},
    {0, 1, 0, 255, 0, 0, 3, 16, "U-GWD", "Zonal Flux of Gravity Wave Stress", "N/m^2"},
    {10, 0, 0, 255, 7, 1, 1, 194, "UBARO", "Barotropic U velocity", "m/s"},
    {0, 1, 0, 255, 0, 0, 3, 31, "UCLSPRS", "Unbalanced Component of Logarithm of Surface Pressure", "-"},
    {0, 1, 0, 255, 0, 0, 1, 120, "UCSCIW", "Unbalanced Component of Specific Cloud Ice Water content", "kg/kg"},
    {0, 1, 0, 255, 0, 0, 1, 119, "UCSCLW", "Unbalanced Component of Specific Cloud Liquid Water content", "kg/kg"},
    {0, 1, 0, 255, 0, 0, 0, 28, "UCTMP", "Unbalanced Component of Temperature", "K"},
    {0, 1, 0, 255, 0, 0, 3, 29, "UDRATE", "Updraught Detrainment Rate", "kg/m^3/s"},
    {0, 1, 0, 255, 0, 0, 2, 17, "UFLX", "Momentum Flux, U-Component", "N/m^2"},
    {0, 1, 0, 255, 0, 0, 2, 2, "UGRD", "U-Component of Wind", "m/s"},
    {0, 1, 0, 255, 0, 0, 2, 23, "UGUST", "U-Component of Wind (Gust)", "m/s"},
    {0, 1, 0, 255, 0, 0, 2, 41, "UGWIND", "U-Component of Geostrophic Wind", "m/s"},
    {10, 1, 0, 255, 0, 0, 2, 4, "UICE", "U-Component of Ice Drift", "m/s"},
    {0, 1, 0, 255, 0, 0, 1, 91, "UKMFLX", "U-component (zonal) Kinematic Moisture Flux", "kg/kg*m/s"},
    {0, 0, 0, 255, 7, 1, 5, 193, "ULWRF", "Upward Long-Wave Rad. Flux", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 5, 4, "ULWRF", "Upward Long-Wave Rad. Flux", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 3, 27, "UMFLX", "Updraught Mass Flux", "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 1, 118, "UNCSH", "Unbalanced Component of Specific Humidity", "kg/kg"},
    {0, 1, 0, 255, 0, 0, 2, 45, "UNDIV", "Unbalanced Component of Divergence", "1/s"},
    {10, 1, 0, 255, 0, 0, 1, 2, "UOGRD", "U-Component of Current", "m/s"},
    {1, 1, 0, 255, 0, 0, 0, 14, "UPAPCP", "Upstream Accumulated Precipitation", "kg/m^2"},
    {1, 1, 0, 255, 0, 0, 0, 15, "UPASM", "Upstream Accumulated Snow Melt", "kg/m^2"},
    {0, 0, 0, 255, 7, 1, 7, 197, "UPHL", "Updraft Helicity", "m^2/s^2"},
    {0, 1, 0, 255, 0, 0, 7, 15, "UPHL", "Updraft Helicity", "m^2/s^2"},
    {2, 1, 0, 255, 0, 0, 3, 2, "UPLSM", "Upper Layer Soil Moisture", "kg/m^3"},
    {2, 1, 0, 255, 0, 0, 3, 1, "UPLST", "Upper Layer Soil Temperature", "K"},
    {3, 0, 0, 255, 7, 1, 1, 192, "USCT", "Scatterometer Estimated U Wind Component", "m/s"},
    {10, 1, 0, 255, 0, 0, 0, 21, "USSD", "U-component Surface Stokes Drift", "m/s"},
    {0, 0, 0, 255, 7, 1, 2, 194, "USTM", "U-Component Storm Motion", "m/s"},
    {0, 1, 0, 255, 0, 0, 2, 27, "USTM", "U-Component Storm Motion", "m/s"},
    {0, 0, 0, 255, 7, 1, 4, 193, "USWRF", "Upward Short-Wave Radiation Flux", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 4, 8, "USWRF", "Upward Short-Wave Radiation Flux", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 4, 53, "USWRFCS", "Upward Short-Wave Radiation Flux, Clear Sky", "W/m^2"},
    {0, 0, 0, 255, 7, 1, 4, 205, "UTRF", "Upward Total Radiation Flux", "W/m^2"},
    {0, 0, 0, 255, 7, 1, 7, 196, "UVI", "Ultra Violet Index", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 4, 51, "UVI", "UV Index", "W/m^2"},
    {0, 1, 0, 255, 0, 0, 4, 50, "UVIUCS", "UV Index (Under Clear Sky)", "Numeric"},
    {0, 0, 0, 255, 7, 1, 3, 195, "V-GWD", "Meridional Flux of Gravity Wave Stress", "N/m^2"},
    {0, 1, 0, 255, 0, 0, 3, 17, "V-GWD", "Meridional Flux of Gravity Wave Stress", "N/m^2"},
    {0, 0, 0, 255, 7, 1, 19, 232, "VAFTD", "Volcanic Ash Forecast Transport and Dispersion", "log10(kg/m^3)"},
    {0, 1, 0, 255, 0, 0, 1, 4, "VAPP", "Vapour Pressure", "Pa"},
    {10, 0, 0, 255, 7, 1, 1, 195, "VBARO", "Barotropic V velocity", "m/s"},
    {0, 0, 0, 255, 7, 1, 4, 200, "VBDSF", "Visible Beam Downward Solar Flux", "W/m^2"},
    {0, 0, 0, 255, 7, 1, 4, 201, "VDDSF", "Visible Diffuse Downward Solar Flux", "W/m^2"},
    {0, 0, 0, 255, 7, 1, 0, 202, "VDFHR", "Vertical Diffusion Heating rate", "K/s"},
    {0, 0, 0, 255, 7, 1, 1, 215, "VDFMR", "Vertical Diffusion Moistening Rate", "kg/kg/s"},
    {0, 0, 0, 255, 7, 1, 14, 195, "VDFOZ", "Ozone Vertical Diffusion", "kg/kg/s"},
    {0, 0, 0, 255, 7, 1, 2, 208, "VDFUA", "Vertical Diffusion Zonal Acceleration", "m/s^2"},
    {0, 0, 0, 255, 7, 1, 2, 209, "VDFVA", "Vertical Diffusion Meridional Acceleration", "m/s^2"},
    {0, 0, 0, 255, 7, 1, 2, 204, "VEDH", "Vertical Eddy Diffusivity Heat exchange", "m^2/s"},
    {2, 1, 0, 255, 0, 0, 0, 4, "VEG", "Vegetation", "%"},
    {2, 0, 0, 255, 7, 1, 0, 210, "VEGT", "Vegetation canopy temperature", "K"},
    {4, 1, 0, 255, 0, 0, 1, 1, "VEL1", "1st Vector Component of Velocity (Coordinate system dependent)", "m/s"},
    {4, 1, 0, 255, 0, 0, 1, 2, "VEL2", "2nd Vector Component of Velocity (Coordinate system dependent)", "m/s"},
    {4, 1, 0, 255, 0, 0, 1, 3, "VEL3", "3rd Vector Component of Velocity (Coordinate system dependent)", "m/s"},
    {0, 1, 0, 255, 0, 0, 2, 18, "VFLX", "Momentum Flux, V-Component", "N/m^2"},
    {0, 1, 0, 255, 0, 0, 6, 48, "VFRCICE", "Volume Fraction of Cloud Ice Particles", "Numeric"},
    {0, 1, 0, 255, 0, 0, 6, 49, "VFRCIW", "Volume Fraction of Cloud (Ice and/or Water)", "Numeric"},
    {0, 1, 0, 255, 0, 0, 6, 47, "VFRCWD", "Volume Fraction of Cloud Water Droplets", "Numeric"},
    {0, 1, 0, 255, 0, 0, 2, 3, "VGRD", "V-Component of Wind", "m/s"},
    {2, 0, 0, 255, 7, 1, 0, 198, "VGTYP", "Vegetation Type", "Integer(0-13)"},
    {0, 1, 0, 255, 0, 0, 2, 24, "VGUST", "V-Component of Wind (Gust)", "m/s"},
    {0, 1, 0, 255, 0, 0, 2, 42, "VGWIND", "V-Component of Geostrophic Wind", "m/s"},
    {10, 1, 0, 255, 0, 0, 2, 5, "VICE", "V-Component of Ice Drift", "m/s"},
    {0, 1, 0, 255, 0, 0, 15, 3, "VIL", "Vertically-Integrated Liquid Water", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 19, 0, "VIS", "Visibility", "m"},
    {0, 1, 0, 255, 0, 0, 19, 35, "VISBSN", "Visibility Through Blowing Snow", "m"},
    {0, 1, 0, 255, 0, 0, 19, 34, "VISIFOG", "Visibility Through Ice Fog", "m"},
    {0, 1, 0, 255, 0, 0, 19, 33, "VISLFOG", "Visibility Through Liquid Fog", "m"},
    {0, 1, 0, 255, 0, 0, 1, 92, "VKMFLX", "V-component (meridional) Kinematic Moisture Flux", "kg/kg*m/s"},
    {0, 1, 0, 255, 0, 0, 20, 52, "VMXR", "Volume Mixing Ratio (Fraction in Air)", "mol/mol"},
    {10, 1, 0, 255, 0, 0, 1, 3, "VOGRD", "V-Component of Current", "m/s"},
    {3, 1, 0, 255, 0, 0, 4, 4, "VOLACDEM", "Volcanic Ash Cloud Emissity", "Numeric"},
    {3, 1, 0, 255, 0, 0, 4, 7, "VOLACDEN", "Volcanic Ash Column Density", "kg/m^2"},
    {3, 1, 0, 255, 0, 0, 4, 6, "VOLACDOD", "Volcanic Ash Cloud Optical Depth", "Numeric"},
    {3, 1, 0, 255, 0, 0, 4, 3, "VOLACDTH", "Volcanic Ash Cloud Top Height", "m"},
    {3, 1, 0, 255, 0, 0, 4, 2, "VOLACDTP", "Volcanic Ash Cloud Top Pressure", "Pa"},
    {3, 1, 0, 255, 0, 0, 4, 1, "VOLACDTT", "Volcanic Ash Cloud Top Temperature", "K"},
    {3, 1, 0, 255, 0, 0, 4, 5, "VOLAEADR", "Volcanic Ash Effective Absorption Depth Ratio", "Numeric"},
    {3, 1, 0, 255, 0, 0, 4, 8, "VOLAPER", "Volcanic Ash Particle Effective Radius", "m"},
    {3, 1, 0, 255, 0, 0, 4, 0, "VOLAPROB", "Volcanic Ash Probability", "%"},
    {0, 1, 0, 255, 0, 0, 19, 4, "VOLASH", "Volcanic Ash", "-"},
    {2, 1, 0, 255, 0, 0, 3, 13, "VOLDEC", "Volumetric Direct Evaporation Cease(Soil Moisture)", "m^3/m^3"},
    {2, 1, 0, 255, 0, 0, 3, 11, "VOLTSO", "Volumetric Transpiration Stree-Onset(Soil Moisture)", "m^3/m^3"},
    {0, 1, 0, 255, 0, 0, 2, 46, "VORTADV", "Vorticity Advection", "s-2"},
    {0, 1, 0, 255, 0, 0, 2, 5, "VPOT", "Velocity Potential", "m^2/s"},
    {0, 1, 0, 255, 0, 0, 0, 15, "VPTMP", "Virtual Potential Temperature", "K"},
    {0, 0, 0, 255, 7, 1, 2, 224, "VRATE", "Ventilation Rate", "m^2/s"},
    {3, 0, 0, 255, 7, 1, 1, 193, "VSCT", "Scatterometer Estimated V Wind Component", "m/s"},
    {2, 1, 0, 255, 0, 0, 0, 25, "VSOILM", "Volumetric Soil Moisture", "m^3/m^3"},
    {2, 1, 0, 255, 0, 0, 3, 16, "VSOSM", "Volumetric Saturation Of Soil Moisture", "m^3/m^3"},
    {10, 1, 0, 255, 0, 0, 0, 22, "VSSD", "V-component Surface Stokes Drift", "m/s"},
    {0, 0, 0, 255, 7, 1, 2, 195, "VSTM", "V-Component Storm Motion", "m/s"},
    {0, 1, 0, 255, 0, 0, 2, 28, "VSTM", "V-Component Storm Motion", "m/s"},
    {4, 1, 0, 255, 0, 0, 2, 4, "VTEC", "Vertical Electron Content", "1/m^2"},
    {0, 1, 0, 255, 0, 0, 0, 1, "VTMP", "Virtual Temperature", "K"},
    {0, 1, 0, 255, 0, 0, 2, 15, "VUCSH", "Vertical U-Component Shear", "1/s"},
    {0, 1, 0, 255, 0, 0, 2, 16, "VVCSH", "Vertical V-Component Shear", "1/s"},
    {0, 1, 0, 255, 0, 0, 2, 8, "VVEL", "Vertical Velocity (Pressure)", "Pa/s"},
    {2, 1, 0, 255, 0, 0, 0, 27, "VWILTP", "Volumetric Wilting Point", "m^3/m^3"},
    {0, 0, 0, 255, 7, 1, 2, 192, "VWSH", "Vertical Speed Shear", "1/s"},
    {0, 1, 0, 255, 0, 0, 2, 25, "VWSH", "Vertical Speed Shear", "1/s"},
    {10, 1, 0, 255, 0, 0, 4, 17, "WATDENA", "Water Density Anomaly (sigma)", "kg/m^3"},
    {10, 1, 0, 255, 0, 0, 4, 16, "WATERDEN", "Water Density (rho)", "kg/m^3"},
    {10, 1, 0, 255, 0, 0, 4, 19, "WATPDEN", "Water potential density (rho theta)", "kg/m^3"},
    {10, 1, 0, 255, 0, 0, 4, 20, "WATPDENA", "Water potential density anomaly (sigma theta) See Note", "kg/m^3"},
    {10, 1, 0, 255, 0, 0, 4, 18, "WATPTEMP", "Water potential temperature (theta)", "K"},
    {2, 1, 0, 255, 0, 0, 0, 5, "WATR", "Water Runoff", "kg/m^2"},
    {10, 1, 0, 255, 0, 0, 0, 62, "WAVEFREW", "Wave frequency width", "-"},
    {2, 0, 0, 255, 7, 1, 0, 223, "WCCONV", "Water Condensate Flux Convergance (Vertical Int)", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 0, 13, "WCF", "Wind Chill Factor", "K"},
    {2, 0, 0, 255, 7, 1, 0, 221, "WCINC", "Water condensate added by precip assimilation", "kg/m^2"},
    {2, 0, 0, 255, 7, 1, 0, 226, "WCUFLX", "Water Condensate Zonal Flux (Vertical Int)", "kg/m^2"},
    {2, 0, 0, 255, 7, 1, 0, 227, "WCVFLX", "Water Condensate Meridional Flux (Vertical Int)", "kg/m^2"},
    {0, 1, 0, 255, 0, 0, 20, 10, "WDCPMFLX", "Wet Deposition by Convective Precipitation Mass Flux", "kg/m^2/s"},
    {10, 1, 0, 255, 0, 0, 4, 14, "WDEPTH", "Water Depth", "m"},
    {0, 1, 0, 255, 0, 0, 2, 0, "WDIR", "Wind Direction (from which blowing)", "deg"},
    {10, 1, 0, 255, 0, 0, 0, 31, "WDIRW", "Wave Directional Width", "-"},
    {1, 1, 0, 255, 0, 0, 2, 0, "WDPTHIL", "Water Depth", "m"},
    {10, 1, 0, 255, 0, 0, 0, 56, "WDWFSWEL", "Wave directional width of first swell partition", "-"},
    {10, 1, 0, 255, 0, 0, 0, 57, "WDWSSWEL", "Wave directional width of second swell partition", "-"},
    {10, 1, 0, 255, 0, 0, 0, 58, "WDWTSWEL", "Wave directional width of third swell partition", "-"},
    {0, 1, 0, 255, 0, 0, 1, 13, "WEASD", "Water Equivalent of Accumulated Snow Depth", "kg/m^2"},
    {10, 1, 0, 255, 0, 0, 0, 42, "WESP", "Wave Engery Spectrum", "1/m^2*s/rad"},
    {0, 1, 0, 255, 0, 0, 0, 27, "WETBT", "Wet Bulb Temperature", "K"},
    {0, 1, 0, 255, 0, 0, 18, 11, "WETDEP", "Wet Deposition", "Bq/m^2"},
    {0, 1, 0, 255, 0, 0, 20, 75, "WFIREFLX", "Wildfire flux", "kg/m^2/s"},
    {1, 1, 0, 255, 0, 0, 2, 2, "WFRACT", "Water Fraction", "Proportion"},
    {10, 1, 0, 255, 0, 0, 0, 59, "WFWFSWEL", "Wave frequency width of first swell partition", "-"},
    {10, 1, 0, 255, 0, 0, 0, 60, "WFWSSWEL", "Wave frequency width of second swell partition", "-"},
    {10, 1, 0, 255, 0, 0, 0, 61, "WFWTSWEL", "Wave frequency width of third swell partition", "-"},
    {4, 1, 0, 255, 0, 0, 8, 5, "WHTCOR", "White Light Coronagraph Radiance", "W/sr/m^2"},
    {4, 1, 0, 255, 0, 0, 8, 3, "WHTRAD", "White Light Radiance", "W/sr/m^2"},
    {2, 0, 0, 255, 7, 1, 0, 201, "WILT", "Wilting Point", "Fraction"},
    {2, 1, 0, 255, 0, 0, 0, 17, "WILT", "Wilting Point", "Proportion"},
    {2, 1, 0, 255, 0, 0, 0, 26, "WILTPT", "Wilting Point", "kg/m^3"},
    {0, 1, 0, 255, 0, 0, 2, 1, "WIND", "Wind Speed", "m/s"},
    {0, 1, 0, 255, 0, 0, 2, 33, "WINDF", "Wind Fetch", "m"},
    {0, 0, 0, 255, 7, 1, 19, 199, "WINDPROB", "Wind probability", "%"},
    {3, 1, 0, 255, 0, 0, 1, 19, "WINDS", "Wind Speed", "m/s"},
    {0, 1, 0, 255, 0, 0, 19, 25, "WIWW", "Weather", "-"},
    {10, 0, 0, 255, 7, 1, 0, 193, "WLENG", "Wave Length", "-"},
    {0, 1, 0, 255, 0, 0, 20, 9, "WLSMFLX", "Wet Deposition by Large-Scale Precipitation Mass Flux", "kg/m^2/s"},
    {0, 1, 0, 255, 0, 0, 2, 19, "WMIXE", "Wind Mixing Energy", "J"},
    {2, 1, 0, 255, 0, 0, 0, 33, "WROD", "Water Runoff and Drainage", "kg/m^2"},
    {10, 0, 0, 255, 7, 1, 0, 192, "WSTP", "Wave Steepness", "Proportion"},
    {10, 1, 0, 255, 0, 0, 0, 18, "WSTR", "Wave Stress", "N/m^2"},
    {0, 0, 0, 255, 7, 1, 2, 214, "WTEND", "Tendency of vertical velocity", "m/s^2"},
    {10, 1, 0, 255, 0, 0, 3, 0, "WTMP", "Water Temperature", "K"},
    {10, 0, 0, 255, 7, 1, 4, 192, "WTMPC", "3-D Temperature", "degC"},
    {1, 1, 0, 255, 0, 0, 2, 1, "WTMPIL", "Water Temperature", "K"},
    {10, 1, 0, 255, 0, 0, 4, 15, "WTMPSS", "Water Temperature", "K"},
    {2, 0, 0, 255, 7, 1, 0, 222, "WVCONV", "Water Vapor Flux Convergance (Vertical Int)", "kg/m^2"},
    {10, 1, 0, 255, 0, 0, 0, 4, "WVDIR", "Direction of Wind Waves", "deg"},
    {10, 1, 0, 255, 0, 0, 0, 5, "WVHGT", "Significant Height of Wind Waves", "m"},
    {2, 0, 0, 255, 7, 1, 0, 220, "WVINC", "Water vapor added by precip assimilation", "kg/m^2"},
    {10, 1, 0, 255, 0, 0, 0, 6, "WVPER", "Mean Period of Wind Waves", "s"},
    {10, 1, 0, 255, 0, 0, 0, 0, "WVSP1", "Wave Spectra (1)", "-"},
    {10, 1, 0, 255, 0, 0, 0, 1, "WVSP2", "Wave Spectra (2)", "-"},
    {10, 1, 0, 255, 0, 0, 0, 2, "WVSP3", "Wave Spectra (3)", "-"},
    {2, 0, 0, 255, 7, 1, 0, 224, "WVUFLX", "Water Vapor Zonal Flux (Vertical Int)", "kg/m^2"},
    {2, 0, 0, 255, 7, 1, 0, 225, "WVVFLX", "Water Vapor Meridional Flux (Vertical Int)", "kg/m^2"},
    {10, 1, 0, 255, 0, 0, 0, 14, "WWSDIR", "Direction of Combined Wind Waves and Swell", "deg"},
    {4, 1, 0, 255, 0, 0, 6, 1, "XLONG", "Solar X-ray Flux (XRS Long)", "W/m^2"},
    {4, 1, 0, 255, 0, 0, 8, 0, "XRAYRAD", "X-Ray Radiance", "W/sr/m^2"},
    {4, 1, 0, 255, 0, 0, 6, 2, "XSHRT", "Solar X-ray Flux (XRS Short)", "W/m^2"},
    {10, 1, 0, 255, 0, 0, 2, 10, "ZVCICEP", "Zonal Vector Component of Vertically Integrated Ice Internal Pressure",
     "Pa*m"},

    /* END MARKER */
    {-1, -1, -1, -1, -1, -1, -1, -1, NULL, NULL, NULL}

};

/* end of definition of gribtable */

/*
 * rd_grib2_msg_seq_file.c *                              Wesley Ebisuzaki
 *
 * unsigned char rd_grib2_msg(FILE *input, long int pos, *int len)
 *
 * to read grib2 message
 *
 *    msg = rd_grib_msg(input, position, &len)
 *
 *    input is the file
 *    position is the byte location (should be set to zero for first record
 *    msg = location of message
 *
 *    to get next record: *  position = position + len;
 *
 * rd_grib_msg allocates its own buffer which can be seen by the
 * parsing routines
 *
 * 1/2007 cleanup M. Schwarb
 * 6/2009 fix repeated bitmaps W. Ebisuzuaki
 * 1/2011 added seq input W. Ebisuzaki
 */

/* from rd_seq_grib.c */

#define BUFF_ALLOC0 (1024 * 64)

/* ascii values of GRIB */
#define G 71
#define R 82
#define I 73
#define B 66

/* memory buffers */
extern unsigned char *mem_buffer[N_mem_buffers];
extern size_t mem_buffer_size[N_mem_buffers];
extern size_t mem_buffer_pos[N_mem_buffers];

/* lots deleted */
int fopen_file(struct seq_file *file, const char *filename, const char *open_mode) {
    int i;

    file->unget_cnt = 0;
    file->pos = 0;
    file->buffer = NULL;

    file->file_type = NOT_OPEN;
    file->cfile = ffopen(filename, open_mode);
    if (file->cfile == NULL) return 1;
    file->file_type = DISK;
    return 0;
}

void fatal_error2(const char *fmt, const char *string) {
    fprintf(stderr, "\n*** FATAL ERROR: ");
    fprintf(stderr, fmt, string);
    fprintf(stderr, " ***\n\n");
    exit(8);
    return;
}

static int getchar_seq_file(struct seq_file *file) {
    file->pos++;
    if (file->unget_cnt) {
        return file->unget_buf[--file->unget_cnt];
    }
    return fgetc_file(file);
}

static int unget_seq_file(struct seq_file *file, int c) {
    file->pos--;
    file->unget_buf[file->unget_cnt++] = (unsigned char)c;
    return 0;
}

size_t fread_file2(void *ptr, size_t size, size_t nmemb, struct seq_file *file) {
    size_t iret;

    if (file->file_type == DISK || file->file_type == PIPE) {
        iret = fread(ptr, size, nmemb, file->cfile);
        if (ferror(file->cfile)) fprintf(stderr, "\nWARNING: file read error\n");
        return iret;
    }
    return 0;
}

int fgetc_file(struct seq_file *file) {
    int n;
    if (file->file_type == DISK || file->file_type == PIPE) return fgetc(file->cfile);
    if (file->file_type == MEM) {
        n = file->n_mem_buffer;
        if (mem_buffer_size[n] > mem_buffer_pos[n]) {
            return mem_buffer[n][mem_buffer_pos[n]++];
        }
    }
    return EOF;
}

void fclose_file(struct seq_file *file) {
    /* do not close stdin or stdout */
    if (file->cfile == stdin || file->cfile == stdout) return;

    if (file->file_type == PIPE || file->file_type == DISK) ffclose(file->cfile);
    if (file->buffer != NULL) {
        if (file->file_type == PIPE || file->file_type == DISK) {
            free(file->buffer);
            file->buffer = NULL;
            file->buffer_size = 0;
        }
    }
    file->file_type = NOT_OPEN;
    return;
}

unsigned char *rd_grib2_msg_seq_file(unsigned char **sec, struct seq_file *input, long int *pos, unsigned long int *len,
                                     int *num_submsgs) {
    int i, c, c1, c2, c3, c4;
    size_t j, len_grib;
    unsigned char *p, *end_of_msg;

    /* setup grib buffer */
    if (input->buffer == NULL) {
        if ((input->buffer = (unsigned char *)malloc(BUFF_ALLOC0)) == NULL) {
            fatal_error2("not enough memory: rd_grib2_msg", "");
        }
        input->buffer_size = BUFF_ALLOC0;
    }

    /* search for GRIB...2 */

    while (1) {
        c = getchar_seq_file(input);
        if (c == EOF) {
            *len = 0;
            return NULL;
        }
        if (c != G) continue;
        if ((c = getchar_seq_file(input)) != R) {
            unget_seq_file(input, c);
            continue;
        }
        if ((c = getchar_seq_file(input)) != I) {
            unget_seq_file(input, c);
            continue;
        }
        if ((c = getchar_seq_file(input)) != B) {
            unget_seq_file(input, c);
            continue;
        }
        c1 = getchar_seq_file(input);
        c2 = getchar_seq_file(input);
        c3 = getchar_seq_file(input);
        c4 = getchar_seq_file(input);
        if (c4 == 1) {
            fprintf(stderr, "grib1 message ignored (use wgrib)\n");
            continue;
        }
        if (c4 != 2) {
            unget_seq_file(input, c4);
            unget_seq_file(input, c3);
            unget_seq_file(input, c2);
            unget_seq_file(input, c1);
            continue;
        }
        input->buffer[0] = G;
        input->buffer[1] = R;
        input->buffer[2] = I;
        input->buffer[3] = B;
        input->buffer[4] = c1;
        input->buffer[5] = c2;
        input->buffer[6] = c3;
        input->buffer[7] = c4;

        /* fill in the size 8-15, unget buffer is empty */

        for (i = 0; i < 8; i++) {
            input->buffer[8 + i] = c = getchar_seq_file(input);
            if (c == EOF) {
                *len = 0;
                return NULL;
            }
        }
        break;
    }

    *len = len_grib = uint8(input->buffer + 8);

    *pos = input->pos - 16;

    if ((size_t)input->buffer_size < len_grib) {
        input->buffer_size = (len_grib * 5) / 4;
        input->buffer = (unsigned char *)realloc((void *)input->buffer, input->buffer_size);
        if (input->buffer == NULL) fatal_error2("rd_grib2_msg_seq_file: grib message size too large", "");
    }

    j = fread_file2(input->buffer + 16, sizeof(unsigned char), len_grib - 16, input);
    input->pos += j;

    if (j != len_grib - 16) fatal_error2("rd_grib2_msg_seq_file, read outside of file, bad grib file", "");

    sec[0] = input->buffer;
    sec[8] = sec[0] + len_grib - 4;
    if (sec[8][0] != 55 || sec[8][1] != 55 || sec[8][2] != 55 || sec[8][3] != 55) {
        fatal_error2("rd_grib2_msg_seq_file, missing end section ('7777')", "");
    }

    /* scan message for number of submessages and perhaps for errors */
    p = sec[0] + GB2_Sec0_size;
    end_of_msg = sec[0] + len_grib;

    i = 0;
    while (p < sec[8]) {
        if (p[4] == 7) i++;
        if (uint4(p) < 5) fatal_error2("rd_grib2_msg_file: illegal grib: section length, section", "");
        p += uint4(p);
        if (p > end_of_msg) fatal_error2("bad grib format", "");
    }
    if (p != sec[8]) {
        fatal_error2("rd_grib2_msg: illegal format, end section expected", "");
    }
    *num_submsgs = i;

    *len = len_grib;
    return sec[0];
}

/* from parse_msg.c */

/*
 * parse_msg.c public domain 2007                              Wesley Ebisuzaki
 *
 * 1/2007 cleanup M. Schwarb
 * 6/2009 fix repeated bitmaps W. Ebisuzuaki
 * 1/2011 added seq input W. Ebisuzaki
 * 12/2014 W. Ebisuzaki, own file, updated to use sec[]
 */

/*
 * with grib 1, a message = 1 field
 * with grib 2, a message can have more than one field
 *
 * this routine parses a grib2 message that has already been read into buffer
 *
 * parse_1st_msg .. returns 1st message starting at sec[0], fills in sec[]
 *   note sec[9] is used to store pointer to last valid bitmap
 */

int parse_1st_msg(unsigned char **sec) {
    unsigned char *p;

    if (sec[0] == NULL) fatal_error2("parse_1st_msg .. sec[0] == NULL", "");
    sec[2] = sec[6] = NULL;
    sec[9] = NULL; /* last valid bitmap */
    p = sec[0] + 16;

    while (sec[8] - p > 0) {
        if (p[4] > 8) fatal_error2("parse_1st_msg illegal section", "");
        sec[p[4]] = p;

        /* Section 6: bitmap */
        if (p[4] == 6) {
            if (p[5] == 0) {
                sec[9] = p; /* last valid bitmap */
            } else if (p[5] >= 1 && p[5] <= 253) {
                fatal_error2("parse_1st_msg: predefined bitmaps are not handled", "");
            } else if (p[5] == 254) {
                fatal_error2("parse_1st_msg: illegal grib msg, bitmap not defined code, table 6.0=254", "");
            }
        }

        /* last section */
        if (p[4] == 7) {
            return 0;
        }
        p += uint4(p);
    }
    fatal_error2("parse_1st_msg illegally format grib", "");
    return 1;
}

int parse_next_msg(unsigned char **sec) {
    unsigned char *p, *end_of_msg;

    end_of_msg = sec[0] + GB2_MsgLen(sec);
    p = sec[7];
    if (p[4] != 7) {
        fatal_error2("parse_next_msg: parsing error", "");
    }
    p += uint4(p);
    if (p > end_of_msg) fatal_error2("bad grib fill", "");

    while (p < sec[8]) {
        sec[p[4]] = p;

        // code to handle code table 6.0
        if (p[4] == 6) {
            if (p[5] == 0) {
                sec[9] = p; /* last valid bitmap */
            } else if (p[5] >= 1 && p[5] <= 253) {
                fatal_error2("parse_next_msg: predefined bitmaps are not handled", "");
            } else if (p[5] == 254) {
                if ((sec[6] = sec[9]) == NULL) {
                    fatal_error2("parse_1st_msg: illegal grib msg, bitmap not defined code, table 6.0=254", "");
                }
            }
        }
        if (p[4] == 7) {  // end of message .. save on sec[]
            return 0;
        }
        p += uint4(p);
        if (p > end_of_msg) fatal_error2("bad grib fill", "");
    }
    return 1;
}

/* from cname.c */

struct gribtable_s gribtable[];
struct gribtable_s *user_gribtable = NULL;

static struct gribtable_s *search_gribtablex(struct gribtable_s *gribtable, unsigned char **sec);

/*
 * get the name information
 *
 * if inv_out, name, desc, unit == NULL, not used
 *
 * v1.0 Wesley Ebisuzaki 2006
 * v1.1 Wesley Ebisuzaki 4/2007 netcdf support
 * v1.2 Wesley Ebisuzaki 4/2007 multiple table support
 * v1.3 Wesley Ebisuzaki 6/2011 make parameter cat >= 192 local
 * v1.4 Wesley Ebisuzaki 2/2012 fixed search_gribtab for local tables
 * v1.5 Wesley Ebisuzaki 4/2013 gribtab -> gribtable, added user_gribtable
 */

int getName(unsigned char **sec, int mode, char *inv_out, char *name, char *desc, char *unit) {
    int discipline, center, mastertab, localtab, parmcat, parmnum;
    int pdt;
    struct gribtable_s *p;
    const char *p_unit;

    p = NULL;
    if (user_gribtable != NULL) p = search_gribtablex(user_gribtable, sec);

    if (p == NULL) p = search_gribtablex(gribtable, sec);

    p_unit = "unit";
    if (p) {
        p_unit = p->unit;
        pdt = code_table_4_0(sec);
        if (pdt == 5 || pdt == 9) p_unit = "prob";
    }

    if (p) {
        if (name) strcpy(name, p->name);
        if (desc) strcpy(desc, p->desc);
        if (unit) strcpy(unit, p_unit);

        if (inv_out) {
            sprintf(inv_out, "%s", p->name);
            inv_out += strlen(inv_out);
            if (mode) sprintf(inv_out, " %s [%s]", p->desc, p_unit);
        }
    } else {
        discipline = GB2_Discipline(sec);
        center = GB2_Center(sec);
        mastertab = GB2_MasterTable(sec);
        localtab = GB2_LocalTable(sec);
        parmcat = GB2_ParmCat(sec);
        parmnum = GB2_ParmNum(sec);

        if (name) sprintf(name, "var%d_%d_%d", discipline, parmcat, parmnum);
        if (desc) strcpy(desc, "desc");
        if (unit) strcpy(unit, p_unit);

        if (inv_out) {
            if ((parmnum >= 192 && parmnum <= 254) || (parmcat >= 192 && parmcat <= 254) ||
                (discipline >= 192 && discipline <= 254)) {
                sprintf(inv_out, "var discipline=%d center=%d local_table=%d parmcat=%d parm=%d", discipline, center,
                        localtab, parmcat, parmnum);
            } else {
                sprintf(inv_out, "var discipline=%d master_table=%d parmcat=%d parm=%d", discipline, mastertab, parmcat,
                        parmnum);
            }
        }
    }

    return 0;
}

/*
 * search the grib table
 */

static struct gribtable_s *search_gribtablex(struct gribtable_s *p, unsigned char **sec) {
    int discipline, center, mastertab, localtab, parmcat, parmnum;
    int use_local_table;
    static int count = 0;

    if (p == NULL) return NULL;

    discipline = GB2_Discipline(sec);
    center = GB2_Center(sec);
    mastertab = GB2_MasterTable(sec);
    localtab = GB2_LocalTable(sec);
    parmcat = GB2_ParmCat(sec);
    parmnum = GB2_ParmNum(sec);

    use_local_table = (mastertab == 255) ? 1 : 0;
    if ((parmnum >= 192 && parmnum <= 254) || (parmcat >= 192 && parmcat <= 254) ||
        (discipline >= 192 && discipline <= 254))
        use_local_table = 1;

    if (use_local_table == 1 && localtab == 0) {
        if (count++ < 6) fprintf(stderr, "**** ERROR: local table = 0 is not allowed, set to 1 ***\n");
        localtab = 1;
    }
    if (use_local_table == 1 && localtab == 255) {
        fatal_error2("local gribtable is undefined (255)", "");
    }

    if (!use_local_table) {
        for (; p->disc >= 0; p++) {
            if (discipline == p->disc && (mastertab >= p->mtab_low) && (mastertab <= p->mtab_high) &&
                parmcat == p->pcat && parmnum == p->pnum) {
                return p;
            }
        }
    } else {
        //  printf(">> cname local find: disc %d center %d localtab %d pcat %d pnum %d\n", discipline, center, localtab,
        //  parmcat, parmnum);
        for (; p->disc >= 0; p++) {
            if (discipline == p->disc && center == p->cntr && localtab == p->ltab && parmcat == p->pcat &&
                parmnum == p->pnum) {
                return p;
            }
        }
    }
    return NULL;
}

/* from Sec3.c */

extern int nx, ny, res, scan;
// extern unsigned int nx_, ny_, npnts;
extern enum output_order_type output_order;
// extern char *nl;
static void print_stagger(int scan, char *inv_out);

int n_variable_dim2 = 0;
int *variable_dim2 = NULL, *raw_variable_dim2 = NULL;

int get_nxny_(unsigned char **sec, unsigned int *nx, unsigned int *ny, unsigned int *npnts, int *res, int *scan) {
    int grid_template, n_var_dim, i, j, n_octets, center;
    unsigned int npoints, n;
    unsigned char *gds, *p;

    center = GB2_Center(sec);
    grid_template = code_table_3_1(sec);
    *res = flag_table_3_3(sec);
    *scan = flag_table_3_4(sec);
    gds = sec[3];

    switch (grid_template) {
        case 0:
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 10:
        case 12:
        case 20:
        case 30:
        case 31:
        case 40:
        case 41:
        case 42:
        case 43:
        case 44:
        case 90:
        case 110:
        case 140:
        case 204:
            *nx = uint4_missing(gds + 30);
            *ny = uint4_missing(gds + 34);
            break;
        case 51:
        case 52:
        case 53:
        case 101:
        case 130:
        case 50:
            *nx = GB2_Sec3_npts(sec);
            *ny = 1;
            break;  // should calculate for from parameters
        case 120:
            *nx = uint4_missing(gds + 14);  // nx = bin along radials, ny = num radials
            *ny = uint4_missing(gds + 18);
            break;
#ifdef WMO_VALIDATION
        case 60:
            *nx = uint4_missing(gds + 30);
            *ny = uint4_missing(gds + 34);
            break;
#endif
        case 32768:
        case 32769:
            if (center == NCEP) {
                *nx = uint4_missing(gds + 30);
                *ny = uint4_missing(gds + 34);
                break;
            }
            *nx = *ny = 0;
            break;
        case 40110:
            if ((center == JMA1) || (center == JMA2)) {
                *nx = uint4_missing(gds + 30);
                *ny = uint4_missing(gds + 34);
                break;
            }
            *nx = *ny = 0;
            break;
        case 50120:
            if ((center == JMA1) || (center == JMA2)) {
                *nx = uint4_missing(gds + 14);  // Nb number of point from origin
                *ny = uint4_missing(gds + 18);  // Nr number of angles
                break;
            }
            *nx = *ny = 0;
            break;
        default:
            *nx = *ny = 0;
            break;
    }

    n_var_dim = 0;
    if (*nx == 0) n_var_dim = *ny;
    if (*ny == 0) n_var_dim = *nx;
    if (*nx == 0 && *ny == 0) n_var_dim = 0;

    p = NULL;
    if (n_var_dim) {
        switch (grid_template) {
            case 0:
                p = gds + 72;
                break;
            case 1:
                p = gds + 84;
                break;
            case 2:
                p = gds + 84;
                break;
            case 3:
                p = gds + 96;
                break;
            case 10:
                p = gds + 72;
                break;
            case 40:
                p = gds + 72;
                break;
            case 41:
                p = gds + 84;
                break;
            case 42:
                p = gds + 84;
                break;
            case 43:
                p = gds + 96;
                break;
            case 32768:
                if (GB2_Center(sec) == NCEP)
                    p = gds + 72;
                else
                    p = NULL;
                break;
            case 32769:
                if (GB2_Center(sec) == NCEP)
                    p = gds + 80;
                else
                    p = NULL;
                break;
            default:
                p = NULL;
                break;
        }
    }

    /* calculate number of grid points, check with GDS */
    npoints = 0;
    if (n_var_dim) {
        if (n_variable_dim2 != n_var_dim) {
            if (variable_dim2) free(variable_dim2);
            if (raw_variable_dim2) free(raw_variable_dim2);
            variable_dim2 = (int *)malloc(n_var_dim * sizeof(int));
            raw_variable_dim2 = (int *)malloc(n_var_dim * sizeof(int));

            if (variable_dim2 == NULL || raw_variable_dim2 == NULL) fatal_error2("ran out of memory", "");
            n_variable_dim2 = n_var_dim;
        }
        n_octets = (int)gds[10]; /* number of octets per integer */
        for (i = 0; i < n_var_dim; i++) {
            for (n = j = 0; j < n_octets; j++) {
                n = (n << 8) + (int)*p++;
            }
            raw_variable_dim2[i] = variable_dim2[i] = (int)n;
            npoints += n;
        }

        /* convert variable_dim2 to SN order if needed */
        if (*nx == 0 && GDS_Scan_y(*scan) == 0 && output_order == wesn) {
            for (i = 0; i < (int)*ny; i++) {
                variable_dim2[i] = raw_variable_dim2[*ny - 1 - i];
            }
        }
        /* convert variable_dim2 to NS order if needed */
        else if (*nx == 0 && GDS_Scan_y(*scan) != 0 && output_order == wens) {
            for (i = 0; i < (int)*ny; i++) {
                variable_dim2[i] = raw_variable_dim2[*ny - 1 - i];
            }
        }
    } else if (*nx > 0 && *ny > 0)
        npoints = *nx * *ny;
    *npnts = GB2_Sec3_npts(sec);

#ifdef WMO_VALIDATION
    /* if global cubed sphere */
    if (grid_template == 60 && GDS_Gnom_tile(gds) == 0) npoints *= 6;
#endif

    if ((*nx != 0 || *ny != 0) && GB2_Sec3_npts(sec) != npoints && GDS_Scan_staggered_storage(*scan) == 0) {
        fprintf(stderr, "two values for number of points %u (GDS) %u (calculated)\n", GB2_Sec3_npts(sec), npoints);
    }

    /*
       for (i = 0; i < n_var_dim; i++) {
       printf("%d ", variable_dim2[i]);
       }
       */

    return 0;
}

/* from int8.c */
/*
 * various conversion routines
 *
 * 2006: Public Domain Wesley Ebisuzaki
 * 1/2007: uint8 fix Wesley Ebisuzaki
 */

/* routines to return various sized integers from GRIB file */

unsigned int uint2(unsigned char const *p) { return (p[0] << 8) + p[1]; }

unsigned int uint4(unsigned const char *p) { return ((p[0] << 24) + (p[1] << 16) + (p[2] << 8) + p[3]); }

/*
 * if len character is 255 -> return 1 (missing)
 */
int is_missing(unsigned const char *p, int len) {
    int i;
    for (i = 0; i < len; i++) {
        if (p[i] != 255) return 0;
    }
    return 1;
}

/*
 * uint4_missing
 * if missing return 0
 * uint4_missing is only used in Sec3.c where an undefined nx/ny == 0 is a good responce
 */
unsigned int uint4_missing(unsigned const char *p) {
    unsigned int t;

    t = p[0];
    t = t << 8 | p[1];
    t = t << 8 | p[2];
    t = t << 8 | p[3];

    if (t == 0xffffffff) return 0;
    return t;
}

unsigned long int uint8(unsigned const char *p) {
#if (ULONG_MAX == 4294967295UL)
    if (p[0] || p[1] || p[2] || p[3]) {
        fatal_error2(
            "unsigned value (8 byte integer) too large for machine\n"
            "fatal error .. run on 64-bit machine",
            "");
    }
    return ((unsigned long int)p[4] << 24) + ((unsigned long int)p[5] << 16) + ((unsigned long int)p[6] << 8) +
           (unsigned long int)p[7];
#else
    return ((unsigned long int)p[0] << 56) + ((unsigned long int)p[1] << 48) + ((unsigned long int)p[2] << 40) +
           ((unsigned long int)p[3] << 32) + ((unsigned long int)p[4] << 24) + ((unsigned long int)p[5] << 16) +
           ((unsigned long int)p[6] << 8) + (unsigned long int)p[7];
#endif
}

/*  uint_n: converts n bytes to unsigned int */

unsigned int uint_n(unsigned const char *p, int n) {
    unsigned int i;
    i = 0;
    while (n-- > 0) {
        i = (i << 8) + *p++;
    }
    return i;
}

int int1(unsigned const char *p) {
    int i;
    if (*p & 0x80) {
        i = -(*p & 0x7f);
    } else {
        i = (int)*p;
    }
    return i;
}

int int2(unsigned const char *p) {
    int i;
    if (p[0] & 0x80) {
        i = -(((p[0] & 0x7f) << 8) + p[1]);
    } else {
        i = (p[0] << 8) + p[1];
    }
    return i;
}

int int4(unsigned const char *p) {
    int i;
    if (p[0] & 0x80) {
        i = -(((p[0] & 0x7f) << 24) + (p[1] << 16) + (p[2] << 8) + p[3]);
    } else {
        i = (p[0] << 24) + (p[1] << 16) + (p[2] << 8) + p[3];
    }
    return i;
}

/*  int_n: converts n bytes to int */

int int_n(unsigned const char *p, int n) {
    int i, sign;

    if (n == 0) return 0;
    sign = *p;
    i = *p++ & 127;
    while (n-- > 1) {
        i = i * 256 + (int)*p++;
    }
    if (sign & 0x80) i = -i;
    return i;
}

//
// 2's complement integer4 -- normal storage
//
int int4_comp(unsigned const char *p) {
    int i;
    unsigned int j;

    if (p[0] & 0x80) {
        j = (p[0] << 24) + (p[1] << 16) + (p[2] << 8) + p[3];
        j = (j ^ 0xffffffff) + 1;
        i = 0 - j;
    } else {
        i = (p[0] << 24) + (p[1] << 16) + (p[2] << 8) + p[3];
    }
    return i;
}

//
// floating point values are often represented as int * power of 10
//
float scaled2flt(int scale_factor, int scale_value) {
    if (scale_factor == 0) return (float)scale_value;
    if (scale_factor < 0) return scale_value * Int_Power(10.0, -scale_factor);
    return scale_value / Int_Power(10.0, scale_factor);
}

double scaled2dbl(int scale_factor, int scale_value) {
    if (scale_factor == 0) return (float)scale_value;
    if (scale_factor < 0) return scale_value * Int_Power(10.0, -scale_factor);
    return scale_value / Int_Power(10.0, scale_factor);
}

//
// inverse of scaled2flt
//
int flt2scaled(int scale_factor, float value) {
    if (scale_factor == 0) return (int)value;
    if (scale_factor > 0) return (int)(value * Int_Power(10.0, scale_factor));
    return (int)(value / Int_Power(10.0, -scale_factor));
}

//
// best scaled values
//
int best_scaled_value(double val, int *scale_factor, int *scale_value) {
    int n;

    if (isinf(val)) {
        fatal_error2("best_scaled_value: encountered an infinite value", "");
    }

    if (val == 0.0) {
        *scale_factor = *scale_value = 0;
        return 0;
    }

    n = 0;

    // scale for large numbers
    if (fabs(val) > INT_MAX) {
        n = 0;
        while (fabs(val) > INT_MAX) {
            val *= 0.1;
            n--;
        }
        *scale_factor = n;
        *scale_value = floor(val + 0.5);
        return 0;
    }

    while (fabs(val * 10.0) < INT_MAX && (val - floor(val)) != 0.0) {
        /* removed 3/2014 WNE
           if (fabs( floor(val+0.5) - val)  < 0.00001*fabs(val) ) {
         *scale_factor = n;
         *scale_value = floor(val + 0.5);
         return 0;
         }
         */
        n++;
        val *= 10.0;
    }
    *scale_factor = n;
    *scale_value = floor(val + 0.5);
    return 0;
}

void uint8_char(unsigned long int i, unsigned char *p) {
    int j;
    for (j = 0; j < 8; j++) {
        p[7 - j] = i & 255;
        i = i >> 8;
    }
}

void uint_char(unsigned int i, unsigned char *p) {
    p[0] = (i >> 24) & 255;
    p[1] = (i >> 16) & 255;
    p[2] = (i >> 8) & 255;
    p[3] = (i) & 255;
}

void int_char(int i, unsigned char *p) {
    int sign = 0;
    if (i < 0) {
        sign = 128;
        i = -i;
    }
    p[0] = ((i >> 24) & 127) | sign;
    p[1] = (i >> 16) & 255;
    p[2] = (i >> 8) & 255;
    p[3] = (i) & 255;
    return;
}

void uint2_char(unsigned int i, unsigned char *p) {
    p[0] = (i >> 8) & 255;
    p[1] = (i) & 255;
    return;
}

void int2_char(int i, unsigned char *p) {
    int sign = 0;
    if (i < 0) {
        sign = 128;
        i = -i;
    }
    p[0] = ((i >> 8) & 127) | sign;
    p[1] = i & 255;
    return;
}

/*
 * originally nx and ny were int with -1 == variable
 * with the large grib conversions, nx and ny became unsigned int
 * and zero became indicator of a variable size
 *
 * to keep the output the same .. have a function that returns a string variable
 * non-threaded!
 */

char *nx_str(unsigned int nx) {
    static char string[30];
    if (nx == 0) return "-1";
    sprintf(string, "%u", nx);
    return string;
}

char *ny_str(unsigned int ny) {
    static char string[30];
    if (ny == 0) return "-1";
    sprintf(string, "%u", ny);
    return string;
}

/*
 *  return sub_angle
 *    note: 0 -> 1
 *          undefined -> 1e6
 *    documentation does not say that subangle is unsigned, assumed signed
 */

int sub_angle(unsigned const char *p) {
    /* 0 -> 1 */
    if (p[0] == 0 && p[1] == 0 && p[2] == 0 && p[3] == 0) return 1;
    if (p[0] == 255 && p[1] == 255 && p[2] == 255 && p[3] == 255) return 1000000;
    return int4(p);
}

/* from Code_Values.c */

int sub_missing_values(unsigned char **sec, float *missing1, float *missing2) {
    int i, j;
    unsigned char *p;

    i = code_table_5_5(sec);
    if (i < 1 || i > 2) return 0;
    j = code_table_5_1(sec);
    p = sec[5];
    if (j == 0) {  // ieee
        if (p[23] == 255 && p[24] == 255 && p[25] == 255 && p[26] == 255)
            *missing1 = UNDEFINED;
        else
            *missing1 = ieee2flt(p + 23);
        if (i == 2) {
            if (p[27] == 255 && p[28] == 255 && p[29] == 255 && p[30] == 255)
                *missing1 = UNDEFINED;
            else
                *missing2 = ieee2flt(p + 27);
        }
    } else if (j == 1) {  // integer
        if (p[23] == 255 && p[24] == 255 && p[25] == 255 && p[26] == 255)
            *missing1 = UNDEFINED;
        else
            *missing1 = (float)int4(p + 23);
        if (i == 2) {
            if (p[27] == 255 && p[28] == 255 && p[29] == 255 && p[30] == 255)
                *missing1 = UNDEFINED;
            else
                *missing2 = (float)int4(p + 27);
        }
    }
    return i;
}

void fixed_surfaces(unsigned char **sec, int *type1, float *surface1, int *undef_val1, int *type2, float *surface2,
                    int *undef_val2) {
    unsigned char *p1, *p2;
    *undef_val1 = *undef_val2 = 1;
    *surface1 = *surface2 = UNDEFINED;
    *type1 = *type2 = 255;

    p1 = code_table_4_5a_location(sec);
    p2 = code_table_4_5b_location(sec);

    if (p1 != NULL && *p1 != 255) {
        *type1 = *p1;
        if (p1[1] != 255) {
            if (p1[2] != 255 || p1[3] != 255 || p1[4] != 255 || p1[5] != 255) {
                *undef_val1 = 0;
                *surface1 = scaled2flt(INT1(p1[1]), int4(p1 + 2));
            }
        }
    }
    if (p2 != NULL && *p2 != 255) {
        *type2 = *p2;
        if (p2[1] != 255) {
            if (p2[2] != 255 || p2[3] != 255 || p2[4] != 255 || p2[5] != 255) {
                *undef_val2 = 0;
                *surface2 = scaled2flt(INT1(p2[1]), int4(p2 + 2));
            }
        }
    }
    return;
}

/* from Level.c */
/* Levels.c
 *   2006: public domain wesley ebisuzaki
 *   1/2007: cleanup M. Schwarb
 *   1/2007: Caser Tejeda Hernandez found error in meter underground
 *   2/2007: level 11
 *   2/2007: spelling error fixed
 *   9/2008: type 241 added Nick Lott
 */

/*
 * HEADER:200:lev:inv:0:level (code table 4.5)
 */

/* code table 4.5 */

const char *level_table2[192] = {
    /* 0 */ "reserved",
    /* 1 */ "surface",
    /* 2 */ "cloud base",
    /* 3 */ "cloud top",
    /* 4 */ "0C isotherm",
    /* 5 */ "level of adiabatic condensation from sfc",
    /* 6 */ "max wind",
    /* 7 */ "tropopause",
    /* 8 */ "top of atmosphere",
    /* 9 */ "sea bottom",
    /* 10 */ "entire atmosphere",
    /* 11 */ "cumulonimbus base",
    /* 12 */ "cumulonimbus top",
    /* 13 */ "lowest level %g%% integrated cloud cover",
    /* 14 */ "level of free convection",
    /* 15 */ "convection condensation level",
    /* 16 */ "level of neutral buoyancy",
    /* 17 */ "reserved",
    /* 18 */ "reserved",
    /* 19 */ "reserved",
    /* 20 */ "%g K level",
    /* 21 */ "lowest level > %g kg/m^3",
    /* 22 */ "highest level > %g kg/m^3",
    /* 23 */ "lowest level > %g Bq/m^3",
    /* 24 */ "highest level > %g Bg/m^3",
    /* 25 */ "highest level > %g dBZ",
    /* 26 */ "reserved",
    /* 27 */ "reserved",
    /* 28 */ "reserved",
    /* 29 */ "reserved",
    /* 30 */ "reserved",
    /* 31 */ "reserved",
    /* 32 */ "reserved",
    /* 33 */ "reserved",
    /* 34 */ "reserved",
    /* 35 */ "reserved",
    /* 36 */ "reserved",
    /* 37 */ "reserved",
    /* 38 */ "reserved",
    /* 39 */ "reserved",
    /* 40 */ "reserved",
    /* 41 */ "reserved",
    /* 42 */ "reserved",
    /* 43 */ "reserved",
    /* 44 */ "reserved",
    /* 45 */ "reserved",
    /* 46 */ "reserved",
    /* 47 */ "reserved",
    /* 48 */ "reserved",
    /* 49 */ "reserved",
    /* 50 */ "reserved",
    /* 51 */ "reserved",
    /* 52 */ "reserved",
    /* 53 */ "reserved",
    /* 54 */ "reserved",
    /* 55 */ "reserved",
    /* 56 */ "reserved",
    /* 57 */ "reserved",
    /* 58 */ "reserved",
    /* 59 */ "reserved",
    /* 60 */ "reserved",
    /* 61 */ "reserved",
    /* 62 */ "reserved",
    /* 63 */ "reserved",
    /* 64 */ "reserved",
    /* 65 */ "reserved",
    /* 66 */ "reserved",
    /* 67 */ "reserved",
    /* 68 */ "reserved",
    /* 69 */ "reserved",
    /* 70 */ "reserved",
    /* 71 */ "reserved",
    /* 72 */ "reserved",
    /* 73 */ "reserved",
    /* 74 */ "reserved",
    /* 75 */ "reserved",
    /* 76 */ "reserved",
    /* 77 */ "reserved",
    /* 78 */ "reserved",
    /* 79 */ "reserved",
    /* 80 */ "reserved",
    /* 81 */ "reserved",
    /* 82 */ "reserved",
    /* 83 */ "reserved",
    /* 84 */ "reserved",
    /* 85 */ "reserved",
    /* 86 */ "reserved",
    /* 87 */ "reserved",
    /* 88 */ "reserved",
    /* 89 */ "reserved",
    /* 90 */ "reserved",
    /* 91 */ "reserved",
    /* 92 */ "reserved",
    /* 93 */ "reserved",
    /* 94 */ "reserved",
    /* 95 */ "reserved",
    /* 96 */ "reserved",
    /* 97 */ "reserved",
    /* 98 */ "reserved",
    /* 99 */ "reserved",
    /* 100 */ "%g mb",
    /* 101 */ "mean sea level",
    /* 102 */ "%g m above mean sea level",
    /* 103 */ "%g m above ground",
    /* 104 */ "%g sigma level",
    /* 105 */ "%g hybrid level",
    /* 106 */ "%g m underground",
    /* 107 */ "%g K isentropic level",
    /* 108 */ "%g mb above ground",
    /* 109 */ "PV=%g (Km^2/kg/s) surface",
    /* 110 */ "reserved",
    /* 111 */ "%g Eta level",
    /* 112 */ "reserved",
    /* 113 */ "%g logarithmic hybrid level",
    /* 114 */ "snow level",
    /* 115 */ "%g sigma height level",
    /* 116 */ "reserved",
    /* 117 */ "mixed layer depth",
    /* 118 */ "%g hybrid height level",
    /* 119 */ "%g hybrid pressure level",
    /* 120 */ "reserved",
    /* 121 */ "reserved",
    /* 122 */ "reserved",
    /* 123 */ "reserved",
    /* 124 */ "reserved",
    /* 125 */ "reserved",
    /* 126 */ "reserved",
    /* 127 */ "reserved",
    /* 128 */ "reserved",
    /* 129 */ "reserved",
    /* 130 */ "reserved",
    /* 131 */ "reserved",
    /* 132 */ "reserved",
    /* 133 */ "reserved",
    /* 134 */ "reserved",
    /* 135 */ "reserved",
    /* 136 */ "reserved",
    /* 137 */ "reserved",
    /* 138 */ "reserved",
    /* 139 */ "reserved",
    /* 140 */ "reserved",
    /* 141 */ "reserved",
    /* 142 */ "reserved",
    /* 143 */ "reserved",
    /* 144 */ "reserved",
    /* 145 */ "reserved",
    /* 146 */ "reserved",
    /* 147 */ "reserved",
    /* 148 */ "reserved",
    /* 149 */ "reserved",
    /* 150 */ "%g generalized vertical height coordinate",
    /* 151 */ "soil level %g",
    /* 152 */ "reserved",
    /* 153 */ "reserved",
    /* 154 */ "reserved",
    /* 155 */ "reserved",
    /* 156 */ "reserved",
    /* 157 */ "reserved",
    /* 158 */ "reserved",
    /* 159 */ "reserved",
    /* 160 */ "%g m below sea level",
    /* 161 */ "%g m below water surface",
    /* 162 */ "lake or river bottom",
    /* 163 */ "bottom of sediment layer",
    /* 164 */ "bottom of thermally active sediment layer",
    /* 165 */ "bottom of sediment layer penetrated by thermal wave",
    /* 166 */ "maxing layer",
    /* 167 */ "bottom of root zone",
    /* 168 */ "ocean model level %g",
    /* 169 */ "ocean level %g (kg*m-3) density difference to near surface",
    /* 170 */ "ocean level %g (K) pot. temperature difference to near surface",
    /* 171 */ "reserved",
    /* 172 */ "reserved",
    /* 173 */ "reserved",
    /* 174 */ "top surface of ice on sea, lake or river",
    /* 175 */ "top surface of ice, und snow on sea, lake or river",
    /* 176 */ "bottom surface ice on sea, lake or river",
    /* 177 */ "deep soil",
    /* 178 */ "reserved",
    /* 179 */ "top surface of glacier ice and inland ice",
    /* 180 */ "deep inland or glacier ice",
    /* 181 */ "grid tile land fraction as a model surface",
    /* 182 */ "grid tile water fraction as a model surface",
    /* 183 */ "grid tile ice fraction on sea, lake or river as a model surface",
    /* 184 */ "grid tile glacier ice and inland ice fraction as a model surface",
    /* 185 */ "reserved",
    /* 186 */ "reserved",
    /* 187 */ "reserved",
    /* 188 */ "reserved",
    /* 189 */ "reserved",
    /* 190 */ "reserved",
    /* 191 */ "reserved"};

// 192..255
const char *ncep_level_table2[64] = {
    /* 192 */ "reserved",
    /* 193 */ "reserved",
    /* 194 */ "reserved",
    /* 195 */ "reserved",
    /* 196 */ "reserved",
    /* 197 */ "reserved",
    /* 198 */ "reserved",
    /* 199 */ "reserved",
    /* 200 */ "entire atmosphere (considered as a single layer)",
    /* 201 */ "entire ocean (considered as a single layer)",
    /* 202 */ "reserved",
    /* 203 */ "reserved",
    /* 204 */ "highest tropospheric freezing level",
    /* 205 */ "reserved",
    /* 206 */ "grid scale cloud bottom level",
    /* 207 */ "grid scale cloud top level",
    /* 208 */ "reserved",
    /* 209 */ "boundary layer cloud bottom level",
    /* 210 */ "boundary layer cloud top level",
    /* 211 */ "boundary layer cloud layer",
    /* 212 */ "low cloud bottom level",
    /* 213 */ "low cloud top level",
    /* 214 */ "low cloud layer",
    /* 215 */ "cloud ceiling",
    /* 216 */ "reserved",
    /* 217 */ "reserved",
    /* 218 */ "reserved",
    /* 219 */ "reserved",
    /* 220 */ "planetary boundary layer",
    /* 221 */ "layer between two hybrid levels",
    /* 222 */ "middle cloud bottom level",
    /* 223 */ "middle cloud top level",
    /* 224 */ "middle cloud layer",
    /* 225 */ "reserved",
    /* 226 */ "reserved",
    /* 227 */ "reserved",
    /* 228 */ "reserved",
    /* 229 */ "reserved",
    /* 230 */ "reserved",
    /* 231 */ "reserved",
    /* 232 */ "high cloud bottom level",
    /* 233 */ "high cloud top level",
    /* 234 */ "high cloud layer",
    /* 235 */ "%gC ocean isotherm",
    /* 236 */ "layer between two depths below ocean surface",
    /* 237 */ "bottom of ocean mixed layer",
    /* 238 */ "bottom of ocean isothermal layer",
    /* 239 */ "layer ocean surface and 26C ocean isothermal level",
    /* 240 */ "ocean mixed layer",
    /* 241 */ "%g in sequence",
    /* 242 */ "convective cloud bottom level",
    /* 243 */ "convective cloud top level",
    /* 244 */ "convective cloud layer",
    /* 245 */ "lowest level of the wet bulb zero",
    /* 246 */ "maximum equivalent potential temperature level",
    /* 247 */ "equilibrium level",
    /* 248 */ "shallow convective cloud bottom level",
    /* 249 */ "shallow convective cloud top level",
    /* 250 */ "reserved",
    /* 251 */ "deep convective cloud bottom level",
    /* 252 */ "deep convective cloud top level",
    /* 253 */ "lowest bottom level of supercooled liquid water layer",
    /* 254 */ "highest top level of supercooled liquid water layer",
    /* 255 */ "missing"};

// 192..255
const char *kma_level_table2[64] = {
    /* 192 */ "reserved",
    /* 193 */ "reserved",
    /* 194 */ "reserved",
    /* 195 */ "hybrid sigma pressure",
    /* 196 */ "Height from ground KFT",
    /* 197 */ "canopy",
    /* 198 */ "ICAO convective cloud top level",
    /* 199 */ "ICAO convective cloud base level",
    /* 200 */ "entire atmosphere (considered as a single layer)",
    /* 201 */ "entire ocean (considered as a single layer)",
    /* 202 */ "reserved",
    /* 203 */ "reserved",
    /* 204 */ "highest tropospheric freezing level",
    /* 205 */ "reserved",
    /* 206 */ "grid scale cloud bottom level",
    /* 207 */ "grid scale cloud top level",
    /* 208 */ "reserved",
    /* 209 */ "boundary layer cloud bottom level",
    /* 210 */ "boundary layer cloud top level",
    /* 211 */ "boundary layer cloud layer",
    /* 212 */ "low cloud bottom level",
    /* 213 */ "low cloud top level",
    /* 214 */ "low cloud layer",
    /* 215 */ "cloud ceiling",
    /* 216 */ "reserved",
    /* 217 */ "reserved",
    /* 218 */ "reserved",
    /* 219 */ "reserved",
    /* 220 */ "planetary boundary layer",
    /* 221 */ "layer between two hybrid levels",
    /* 222 */ "middle cloud bottom level",
    /* 223 */ "middle cloud top level",
    /* 224 */ "middle cloud layer",
    /* 225 */ "reserved",
    /* 226 */ "reserved",
    /* 227 */ "reserved",
    /* 228 */ "reserved",
    /* 229 */ "reserved",
    /* 230 */ "reserved",
    /* 231 */ "reserved",
    /* 232 */ "high cloud bottom level",
    /* 233 */ "high cloud top level",
    /* 234 */ "high cloud layer",
    /* 235 */ "%gC ocean isotherm",
    /* 236 */ "layer between two depths below ocean surface",
    /* 237 */ "bottom of ocean mixed layer",
    /* 238 */ "bottom of ocean isothermal layer",
    /* 239 */ "layer ocean surface and 26C ocean isothermal level",
    /* 240 */ "ocean mixed layer",
    /* 241 */ "%g in sequence",
    /* 242 */ "convective cloud bottom level",
    /* 243 */ "convective cloud top level",
    /* 244 */ "convective cloud layer",
    /* 245 */ "lowest level of the wet bulb zero",
    /* 246 */ "maximum equivalent potential temperature level",
    /* 247 */ "equilibrium level",
    /* 248 */ "shallow convective cloud bottom level",
    /* 249 */ "shallow convective cloud top level",
    /* 250 */ "reserved",
    /* 251 */ "deep convective cloud bottom level",
    /* 252 */ "deep convective cloud top level",
    /* 253 */ "lowest bottom level of supercooled liquid water layer",
    /* 254 */ "highest top level of supercooled liquid water layer",
    /* 255 */ "missing"};

int level1two(int mode, int type, int undef_val, float value, int center, int subcenter, char *inv_out);
int level2(int mode, int type1, int undef_val1, float value1, int type2, int undef_val2, float value2, int center,
           int subcenter, char *inv_out);

/*
 * level2 is for layers
 */

int level2(int mode, int type1, int undef_val1, float value1, int type2, int undef_val2, float value2, int center,
           int subcenter, char *inv_out) {
    if (type1 == 100 && type2 == 100) {
        sprintf(inv_out, "%g-%g mb", value1 / 100, value2 / 100);
    } else if (type1 == 102 && type2 == 102) {
        sprintf(inv_out, "%g-%g m above mean sea level", value1, value2);
    } else if (type1 == 103 && type2 == 103) {
        sprintf(inv_out, "%g-%g m above ground", value1, value2);
    } else if (type1 == 104 && type2 == 104) {
        sprintf(inv_out, "%g-%g sigma layer", value1, value2);
    } else if (type1 == 105 && type2 == 105) {
        sprintf(inv_out, "%g-%g hybrid layer", value1, value2);
    } else if (type1 == 106 && type2 == 106) {
        /* sprintf(inv_out,"%g-%g m below ground",value1/100,value2/100); removed 1/2007 */
        sprintf(inv_out, "%g-%g m below ground", value1, value2);
    } else if (type1 == 107 && type2 == 107) {
        sprintf(inv_out, "%g-%g K isentropic layer", value1, value2);
    } else if (type1 == 108 && type2 == 108) {
        sprintf(inv_out, "%g-%g mb above ground", value1 / 100, value2 / 100);
    } else if (type1 == 111 && type2 == 111) {
        sprintf(inv_out, "%g-%g Eta layer", value1, value2);
    } else if (type1 == 115 && type2 == 115) {
        sprintf(inv_out, "%g-%g sigma height layer", value1, value2);
    } else if (type1 == 118 && type2 == 118) {
        sprintf(inv_out, "%g-%g hybrid height layer", value1, value2);
    } else if (type1 == 119 && type2 == 119) {
        sprintf(inv_out, "%g-%g hybrid pressure layer", value1, value2);
    } else if (type1 == 160 && type2 == 160) {
        sprintf(inv_out, "%g-%g m below sea level", value1, value2);
    } else if (type1 == 161 && type2 == 161) {
        sprintf(inv_out, "%g-%g m ocean layer", value1, value2);
    } else if (type1 == 1 && type2 == 8) {
        sprintf(inv_out, "atmos col");  // compatible with wgrib
    } else if (type1 == 9 && type2 == 1) {
        sprintf(inv_out, "ocean column");
    } else if (center == NCEP && type1 == 235 && type2 == 235) {
        sprintf(inv_out, "%g-%gC ocean isotherm layer", value1 / 10, value2 / 10);
    } else if (center == NCEP && type1 == 236 && type2 == 236) {  // obsolete
        sprintf(inv_out, "%g-%g m ocean layer", value1 * 10, value2 * 10);
    } else if (type1 == 255 && type2 == 255) {
        sprintf(inv_out, "no_level");
    } else {
        level1two(mode, type1, undef_val1, value1, center, subcenter, inv_out);
        inv_out += strlen(inv_out);
        if (type2 != 255) {
            sprintf(inv_out, " - ");
            inv_out += strlen(inv_out);
            level1two(mode, type2, undef_val2, value2, center, subcenter, inv_out);
        }
    }
    return 0;
}

/*
 * level1 is for a single level (not a layer)
 */

int level1two(int mode, int type, int undef_val, float val, int center, int subcenter, char *inv_out) {
    /* WMO defined levels */

    if (type < 192) {
        if (type == 100 || type == 108) val = val * 0.01;  // Pa -> mb
        sprintf(inv_out, level_table2[type], val);
        return 0;
    }

    // no numeric information
    if (type == 255) return 8;

    /* local table for NCEP */
    if (center == NCEP && type >= 192 && type <= 254) {
        if (type == 235) val *= 0.01;  // C -> 0.1C
        sprintf(inv_out, ncep_level_table2[type - 192], val);
    }

    else if (center == KMA && type >= 192 && type <= 254) {
        if (type == 235) val *= 0.01;  // C -> 0.1C
        sprintf(inv_out, kma_level_table2[type - 192], val);
    }

    else {
        if (undef_val == 0)
            sprintf(inv_out, "local level type %d %g", type, val);
        else
            sprintf(inv_out, "local level type %d", type);
    }

    return 0;
}

/* from unpk.c */

int unpk_grib(unsigned char **sec, float *data) {
    int packing, bitmap_flag, nbits;
    unsigned int ndata, ii;
    unsigned char *mask_pointer, mask;
    unsigned char *ieee, *p;
    float tmp;
    // float reference, tmp;
    double reference;
    double bin_scale, dec_scale, b;

#ifdef USE_PNG
    int width, height;
#endif

#ifdef USE_JASPER
    jas_image_t *image;
    char *opts;
    jas_stream_t *jpcstream;
    jas_image_cmpt_t *pcmpt;
    jas_matrix_t *jas_data;
    int j, k;
#endif

#if (defined USE_PNG || defined USE_AEC)
    unsigned char *c;
#endif

#ifdef USE_AEC
    struct aec_stream strm;
    int status;
    int numBitsNeeded;
    size_t size;
#endif

    packing = code_table_5_0(sec);
    // ndata = (int) GB2_Sec3_npts(sec);
    ndata = GB2_Sec3_npts(sec);
    bitmap_flag = code_table_6_0(sec);

    if (bitmap_flag != 0 && bitmap_flag != 254 && bitmap_flag != 255) fatal_error2("unknown bitmap", "");

    if (packing == 4) {  // ieee
        if (sec[5][11] != 1) fatal_error2("unpk ieee grib file precision not supported", "");

        // ieee depacking -- simple no bitmap
        if (bitmap_flag == 255) {
            for (ii = 0; ii < ndata; ii++) {
                data[ii] = ieee2flt_nan(sec[7] + 5 + ii * 4);
            }
            return 0;
        }
        if (bitmap_flag == 0 || bitmap_flag == 254) {
            mask_pointer = sec[6] + 6;
            ieee = sec[7] + 5;
            mask = 0;
            for (ii = 0; ii < ndata; ii++) {
                if ((ii & 7) == 0) mask = *mask_pointer++;
                if (mask & 128) {
                    data[ii] = ieee2flt_nan(ieee);
                    ieee += 4;
                } else {
                    data[ii] = UNDEFINED;
                }
                mask <<= 1;
            }
            return 0;
        }
        fatal_error2("unknown bitmap", "");
    } else if (packing == 0 || packing == 61) {  // simple grib1 packing  61 -- log preprocessing

        p = sec[5];
        reference = ieee2flt(p + 11);
        bin_scale = Int_Power(2.0, int2(p + 15));
        dec_scale = Int_Power(10.0, -int2(p + 17));
        nbits = p[19];
        b = 0.0;
        if (packing == 61) b = ieee2flt(p + 20);

        if (bitmap_flag != 0 && bitmap_flag != 254 && bitmap_flag != 255) fatal_error2("unknown bitmap", "");

        if (nbits == 0) {
            tmp = reference * dec_scale;
            if (packing == 61) tmp = exp(tmp) - b;  // remove log prescaling
            if (bitmap_flag == 255) {
                for (ii = 0; ii < ndata; ii++) {
                    data[ii] = tmp;
                }
                return 0;
            }
            if (bitmap_flag == 0 || bitmap_flag == 254) {
                mask_pointer = sec[6] + 6;
                mask = 0;
                for (ii = 0; ii < ndata; ii++) {
                    if ((ii & 7) == 0) mask = *mask_pointer++;
                    data[ii] = (mask & 128) ? tmp : UNDEFINED;
                    mask <<= 1;
                }
                return 0;
            }
        }

        mask_pointer = (bitmap_flag == 255) ? NULL : sec[6] + 6;

        unpk_0(data, sec[7] + 5, mask_pointer, nbits, ndata, reference, bin_scale, dec_scale);

        if (packing == 61) {  // remove log prescaling
            // #pragma omp parallel for private(ii) schedule(static)
            for (ii = 0; ii < ndata; ii++) {
                if (DEFINED_VAL(data[ii])) data[ii] = exp(data[ii]) - b;
            }
        }
        return 0;
    }

    else if (packing == 2 || packing == 3) {  // complex
        return unpk_complex(sec, data, ndata);
    } else if (packing == 200) {  // run length
        return unpk_run_length(sec, data, ndata);
    }
#ifdef USE_JASPER
    else if (packing == 40 || packing == 40000) {  // jpeg2000
        p = sec[5];
        reference = ieee2flt(p + 11);
        bin_scale = Int_Power(2.0, int2(p + 15));
        dec_scale = Int_Power(10.0, -int2(p + 17));
        nbits = p[19];

        if (nbits == 0) {
            tmp = reference * dec_scale;
            if (bitmap_flag == 255) {
                for (ii = 0; ii < ndata; ii++) {
                    data[ii] = tmp;
                }
                return 0;
            }
            if (bitmap_flag == 0 || bitmap_flag == 254) {
                mask_pointer = sec[6] + 6;
                ieee = sec[7] + 5;
                mask = 0;
                for (ii = 0; ii < ndata; ii++) {
                    if ((ii & 7) == 0) mask = *mask_pointer++;
                    data[ii] = (mask & 128) ? tmp : UNDEFINED;
                    mask <<= 1;
                }
                return 0;
            }
            fatal_error2("unknown bitmap", "");
        }

        // decode jpeg2000

        image = NULL;
        opts = NULL;
        jpcstream = jas_stream_memopen((char *)sec[7] + 5, (int)GB2_Sec7_size(sec) - 5);
        image = jpc_decode(jpcstream, opts);
        if (image == NULL) fatal_error2("jpeg2000 decoding", "");
        pcmpt = image->cmpts_[0];
        if (image->numcmpts_ != 1) fatal_error2("unpk: Found color image.  Grayscale expected", "");

        jas_data = jas_matrix_create(jas_image_height(image), jas_image_width(image));
        jas_image_readcmpt(image, 0, 0, 0, jas_image_width(image), jas_image_height(image), jas_data);

        // transfer data

        k = ndata - pcmpt->height_ * pcmpt->width_;

        // #pragma omp parallel for private(ii,j)
        for (ii = 0; ii < pcmpt->height_; ii++) {
            for (j = 0; j < pcmpt->width_; j++) {
                //      data[k++] = (((jas_data->rows_[ii][j])*bin_scale)+reference)*dec_scale;
                data[k + j + ii * pcmpt->width_] = (((jas_data->rows_[ii][j]) * bin_scale) + reference) * dec_scale;
            }
        }

        if (bitmap_flag == 0 || bitmap_flag == 254) {
            k = ndata - pcmpt->height_ * pcmpt->width_;
            mask_pointer = sec[6] + 6;
            mask = 0;
            for (ii = 0; ii < ndata; ii++) {
                if ((ii & 7) == 0) mask = *mask_pointer++;
                data[ii] = (mask & 128) ? data[k++] : UNDEFINED;
                mask <<= 1;
            }
        } else if (bitmap_flag != 255) {
            fatal_error2("unknown bitmap: ", "");
        }
        jas_matrix_destroy(jas_data);
        jas_stream_close(jpcstream);
        jas_image_destroy(image);
        return 0;
    }
#endif
#ifdef USE_PNG
    else if (packing == 41) {  // png
        p = sec[5];
        reference = ieee2flt(p + 11);
        bin_scale = Int_Power(2.0, int2(p + 15));
        dec_scale = Int_Power(10.0, -int2(p + 17));
        nbits = p[19];

        if (nbits == 0) {
            tmp = reference * dec_scale;
            if (bitmap_flag == 255) {
                for (ii = 0; ii < ndata; ii++) {
                    data[ii] = tmp;
                }
                return 0;
            }
            if (bitmap_flag == 0 || bitmap_flag == 254) {
                mask_pointer = sec[6] + 6;
                ieee = sec[7] + 5;
                mask = 0;
                for (ii = 0; ii < ndata; ii++) {
                    if ((ii & 7) == 0) mask = *mask_pointer++;
                    data[ii] = (mask & 128) ? tmp : UNDEFINED;
                    mask <<= 1;
                }
                return 0;
            }
            fatal_error2("unknown bitmap", "");
        }

        if ((c = (unsigned char *)malloc(4 * sizeof(char) * (size_t)ndata)) == NULL)
            fatal_error2("unpk: allocation error", "");

        i = (int)dec_png_clone(sec[7] + 5, &width, &height, (char *)c);
        if (i) fatal_error2("unpk: png decode error ", "");
        mask_pointer = (bitmap_flag == 255) ? NULL : sec[6] + 6;

        //  check sizes

        if (mask_pointer == NULL) {
            if (ndata != width * height) fatal_error2("png size mismatch w*h", "");
        } else {
            if (ndata != width * height + missing_points2(mask_pointer, GB2_Sec3_npts(sec)))
                fatal_error2("png size mismatch", "");
        }

        unpk_0(data, c, mask_pointer, nbits, ndata, reference, bin_scale, dec_scale);
        free(c);
        return 0;
    }
#endif
#ifdef USE_AEC
    else if (packing == 42) {  // aec

        p = sec[5];
        reference = ieee2flt(p + 11);
        bin_scale = Int_Power(2.0, int2(p + 15));
        dec_scale = Int_Power(10.0, -int2(p + 17));
        nbits = p[19];

        if (nbits == 0) {
            tmp = reference * dec_scale;
            if (bitmap_flag == 255) {
                for (ii = 0; ii < ndata; ii++) {
                    data[ii] = tmp;
                }
                return 0;
            }
            if (bitmap_flag == 0 || bitmap_flag == 254) {
                mask_pointer = sec[6] + 6;
                mask = 0;
                for (ii = 0; ii < ndata; ii++) {
                    if ((ii & 7) == 0) mask = *mask_pointer++;
                    data[ii] = (mask & 128) ? tmp : UNDEFINED;
                    mask <<= 1;
                }
                return 0;
            }
            fatal_error2("unknown bitmap", "");
        }

        strm.flags = (int)sec[5][21];
        strm.bits_per_sample = (int)sec[5][19];
        strm.block_size = (int)sec[5][22];
        strm.rsi = uint2(sec[5] + 23);

        strm.next_in = sec[7] + 5;
        strm.avail_in = uint4(sec[7]) - 5;

        numBitsNeeded = (int)sec[5][19];
        size = ((numBitsNeeded + 7) / 8) * (size_t)ndata;

        if ((c = (unsigned char *)malloc(size)) == NULL) fatal_error2("unpk: allocation error", "");

        strm.next_out = c;
        strm.avail_out = size;

        status = aec_buffer_decode(&strm);

        if (status != AEC_OK) fatal_error2("unpk: aec decode error ", "");

        mask_pointer = (bitmap_flag == 255) ? NULL : sec[6] + 6;

        unpk_0(data, c, mask_pointer, ((nbits + 7) / 8) * 8, ndata, reference, bin_scale, dec_scale);

        free(c);
        return 0;
    }
#endif
    fatal_error2("packing type %d not supported", "");
    return 1;
}

/* from missing.c */

static unsigned int bitsum[256] = {
    8, 7, 7, 6, 7, 6, 6, 5, 7, 6, 6, 5, 6, 5, 5, 4, 7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 7, 6, 6, 5, 6,
    5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2, 7, 6, 6, 5, 6, 5, 5, 4, 6, 5,
    5, 4, 5, 4, 4, 3, 6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2, 6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3,
    2, 5, 4, 4, 3, 4, 3, 3, 2, 4, 3, 3, 2, 3, 2, 2, 1, 7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 6, 5, 5, 4,
    5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2, 6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2, 5, 4, 4, 3, 4, 3, 3, 2, 4,
    3, 3, 2, 3, 2, 2, 1, 6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2, 5, 4, 4, 3, 4, 3, 3, 2, 4, 3, 3, 2, 3, 2,
    2, 1, 5, 4, 4, 3, 4, 3, 3, 2, 4, 3, 3, 2, 3, 2, 2, 1, 4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0};

unsigned int missing_points2(unsigned char *bitmap, unsigned int n) {
    unsigned int count, i, j, rem;
    if (bitmap == NULL) return 0;
    /*
       count = 0;
       while (n >= 8) {
       tmp = *bitmap++;
       n -= 8;
       count += bitsum[tmp];
       }
       tmp = *bitmap | ((1 << (8 - n)) - 1);
       count += bitsum[tmp];
       */

    j = n >> 3;
    rem = n & 7;
    count = 0;
#pragma omp parallel for private(i) reduction(+ : count)
    for (i = 0; i < j; i++) {
        count += bitsum[bitmap[i]];
    }
    count += rem ? bitsum[bitmap[j] | ((1 << (8 - rem)) - 1)] : 0;

    return count;
}

/* from addtime.c */

int get_time(unsigned char *p, int *year, int *month, int *day, int *hour, int *minute, int *second) {
    *year = (p[0] << 8) | p[1];
    p += 2;
    *month = (int)*p++;
    *day = (int)*p++;
    *hour = (int)*p++;
    *minute = (int)*p++;
    *second = (int)*p;
    return 0;
}

/* from intpower.c */

double Int_Power(double x, int y) {
    double value;

    if (y < 0) {
        y = -y;
        x = 1.0 / x;
    }
    value = 1.0;

    while (y) {
        if (y & 1) {
            value *= x;
        }
        x = x * x;
        y >>= 1;
    }
    return value;
}

/* from ieee2flt.c */

float ieee2flt(unsigned char *ieee) {
    double fmant;
    int exp;

    if ((ieee[0] & 127) == 0 && ieee[1] == 0 && ieee[2] == 0 && ieee[3] == 0) return (float)0.0;

    exp = ((ieee[0] & 127) << 1) + (ieee[1] >> 7);
    fmant = (double)((int)ieee[3] + (int)(ieee[2] << 8) + (int)((ieee[1] | 128) << 16));
    if (ieee[0] & 128) fmant = -fmant;
    return (float)(ldexp(fmant, (int)(exp - 128 - 22)));
}

float ieee2flt_nan(unsigned char *ieee) {
    double fmant;
    int exp;

    if ((ieee[0] & 127) == 0 && ieee[1] == 0 && ieee[2] == 0 && ieee[3] == 0) return (float)0.0;

    exp = ((ieee[0] & 127) << 1) + (ieee[1] >> 7);

    if (exp == 255) return (float)UNDEFINED;

    fmant = (double)((int)ieee[3] + (int)(ieee[2] << 8) + (int)((ieee[1] | 128) << 16));
    if (ieee[0] & 128) fmant = -fmant;

    return (float)(ldexp(fmant, (int)(exp - 128 - 22)));
}

/* from ffopen.c */

int flush_mode;

struct opened_file {
    char *name;
    FILE *handle;
    int is_read_file;
    int usage_count;
    int do_not_close_flag;
    struct opened_file *next;
};

/* memory buffers */
extern size_t mem_buffer_pos[N_mem_buffers];

/* do not initialize opened_file_start in init_globals */
static struct opened_file *opened_file_start = NULL;

/* do not initialize opened_file_start in init_globals */

FILE *ffopen(const char *filename, const char *mode) {
    struct opened_file *ptr;
    struct stat stat_buf; /* test for pipes */
    int is_read_file, is_write_file, is_append_file, i;
    const char *p;

    /* see if is a read/write file */
    is_read_file = is_write_file = is_append_file = 0;
    p = mode;
    while (*p) {
        if (*p == 'r') is_read_file = 1;
        if (*p == 'w') is_write_file = 1;
        if (*p++ == 'a') is_append_file = 1;
    }
    if (is_read_file + is_write_file + is_append_file != 1) fatal_error2("ffopen: mode is bad %s", mode);

    if (strcmp(filename, "-") == 0) {
        if (is_read_file) return stdin;
        flush_mode = 1;
        return stdout;
    }

    if (strncmp(filename, "@mem:", 5) == 0) {
        fatal_error2("ffopen option does not suport memory files %s", filename);
    }

    /* see if file has already been opened */

    ptr = opened_file_start;
    while (ptr != NULL) {
        if (strcmp(filename, ptr->name) == 0) {
            /* ptr->usage_count   should be zero at start of wgrib2 call
               ptr->is_read_flag, file was opened in read mode
               have to check if the open was for previous wgrib2 call
               */

            /* already open in r/r or (aw)/(aw) modes */
            if (is_read_file == ptr->is_read_file) {
                (ptr->usage_count)++;
                return ptr->handle;
            }

            /* already being used in incompatible mode do not allow r+ w+ or a+ modes */
            if (ptr->usage_count > 0) {
                fatal_error2("ffopen: file %s cannot be used for for read and write at same time", ptr->name);
            }

            /* open a file that has been opened before but not in use */
            /*
               ptr->handle = freopen(NULL,mode, ptr->handle);
               if (!ptr->handle) fatal_error2("ffopen: wgrib2 could not freopen %s", ptr->name);
               */
            /* brute force. close file first. Since not in use, can change ptr->handle */
            i = fclose(ptr->handle);
            if (i)
                fprintf(stderr, "Warning ffopen(%s): closing used file, fclose() failed fclose err=%d\n", ptr->name, i);
            ptr->handle = fopen(ptr->name, mode);
            if (ptr->handle == NULL) fatal_error2("ffopen: failed to reopen %s\n", ptr->name);
            ptr->is_read_file = is_read_file;
            ptr->usage_count = 1;
            return ptr->handle;
        }
        ptr = ptr->next;
    }

    ptr = (struct opened_file *)malloc(sizeof(struct opened_file));
    if (ptr == NULL) fatal_error2("ffopen: memory allocation problem (%s)", filename);
    ptr->name = (char *)malloc(strlen(filename) + 1);
    if (ptr->name == NULL) fatal_error2("ffopen: memory allocation problem (%s)", filename);

    strncpy(ptr->name, filename, strlen(filename) + 1);
    ptr->is_read_file = is_read_file;
    ptr->usage_count = 1;
    ptr->do_not_close_flag = 1;  // default is not to close file file

    /* check if filename is @tmp:XXX or just XXX */
    if (strlen(filename) > 5 && strncmp(filename, "@tmp:", 5) == 0) {
        if (is_read_file) {
            /* remove ptr */
            free(ptr->name);
            free(ptr);
            fatal_error2("ffopen: creating temporary file for read is not a good idea (%s)", filename);
        }
        ptr->handle = tmpfile(); /* ISO C to generate temporary file */
        if (ptr->handle == NULL) {
            /* remove ptr */
            free(ptr->name);
            free(ptr);
            fatal_error2("ffopen: could not open tmpfile %s", filename);
        }
        /* add ptr to linked list */
        ptr->next = opened_file_start;
        opened_file_start = ptr;
        return ptr->handle;
    }

    /* disk file or pipe */
    ptr->handle = fopen(filename, mode);
    if (ptr->handle == NULL) {
        free(ptr->name);
        free(ptr);
        // fatal_error2("ffopen: could not open %s", filename);
        return NULL;
    }

    /* check if output is to a pipe */
    if (is_read_file == 0) {
        if (stat(filename, &stat_buf) == -1) {
            /* remove ptr */
            free(ptr->name);
            free(ptr);
            fatal_error2("ffopen: could not stat file: %s", filename);
        }
        if (S_ISFIFO(stat_buf.st_mode)) {
            flush_mode = 1;
        }
    }

    /* add ptr to linked list */
    ptr->next = opened_file_start;
    opened_file_start = ptr;
    return ptr->handle;
}

int ffclose(FILE *flocal) {
    struct opened_file *ptr, *last_ptr;

    if (flocal == stdin || flocal == stdout) return 0;  // do not close stdin or stdout

    last_ptr = NULL;
    ptr = opened_file_start;
    while (ptr != NULL) {
        if (ptr->handle == flocal) {
            if (ptr->usage_count <= 0) fatal_error2("ffclose: usage_count = %d <= 0, %d", "");
            if (--ptr->usage_count == 0) {  // usage_count is zero, can close file name
                if (ptr->do_not_close_flag == 1) return 0;
                if (last_ptr == NULL)
                    opened_file_start = ptr->next;
                else {
                    last_ptr->next = ptr->next;
                }
                fclose(ptr->handle);
                free(ptr->name);
                free(ptr);
            }
            return 0;
        }
        last_ptr = ptr;
        ptr = ptr->next;
    }
    return 1; /* not found */
}

/* from CodeTable.c */

int code_table_3_1(unsigned char **sec) { return (int)uint2(sec[3] + 12); }

int code_table_4_0(unsigned char **sec) { return GB2_ProdDefTemplateNo(sec); }

int code_table_4_4(unsigned char **sec) {
    unsigned char *p;
    p = code_table_4_4_location(sec);
    if (p == NULL) return -1;
    return (int)*p;
}

unsigned char *code_table_4_4_location(unsigned char **sec) {
    int pdt, center, n;
    pdt = GB2_ProdDefTemplateNo(sec);
    center = GB2_Center(sec);

    switch (pdt) {
        case 0:
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
        case 8:
        case 9:
        case 10:
        case 11:
        case 12:
        case 13:
        case 14:
        case 15:
        case 32:
        case 60:
        case 61:
        case 51:
        case 91:
        case 1000:
        case 1001:
        case 1002:
            return sec[4] + 17;
            break;
        case 40:
        case 41:
        case 42:
        case 43:
            return sec[4] + 19;
            break;
        case 44:
        case 45:
        case 46:
        case 47:
            return sec[4] + 30;
            break;
        case 48:
            return sec[4] + 41;
            break;
        case 52:
            return sec[4] + 20;
            break;
        case 57:
            n = number_of_mode(sec);
            if (n <= 0 || n == 65535) fatal_error2("PDT 4.57 bad number of mode", "");
            return sec[4] + 26 + 5 * n;
            break;
        case 70:
        case 71:
        case 72:
        case 73:
            return sec[4] + 22;
            break;
        case 50008:
            if (center == JMA1 || center == JMA2) return sec[4] + 17;
            return NULL;
            break;
        case 50009:
            if (center == JMA1 || center == JMA2) return sec[4] + 17;
            return NULL;
            break;
        case 50011:
            if (center == JMA1 || center == JMA2) return sec[4] + 17;
            return NULL;
            break;
    }
    return NULL;
}

int code_table_4_5a(unsigned char **sec) {
    unsigned char *p;
    p = code_table_4_5a_location(sec);
    if (p == NULL) return -1;
    return (int)*p;
}

unsigned char *code_table_4_5a_location(unsigned char **sec) {
    int pdt, center, n;
    pdt = GB2_ProdDefTemplateNo(sec);
    center = GB2_Center(sec);

    if (center == JMA1 || center == JMA2) {
        switch (pdt) {
            case 50008:
                return sec[4] + 22;
                break;
            case 50009:
                return sec[4] + 22;
                break;
            case 51020:
            case 51021:
            case 51022:
            case 51122:
            case 52020:
                return NULL;
                break;
            case 50010:
            case 50011:
                return sec[4] + 22;
                break;
            default:
                break;
        }
    }

    switch (pdt) {
        case 0:
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
        case 8:
        case 9:
        case 10:
        case 11:
        case 12:
        case 13:
        case 14:
        case 15:
        case 51:
        case 60:
        case 61:
        case 91:
        case 1100:
        case 1101:
            return sec[4] + 22;
            break;
        case 40:
        case 41:
        case 42:
        case 43:
            return sec[4] + 24;
            break;
        case 44:
            return sec[4] + 33;
            break;
        case 45:
        case 46:
        case 47:
            return sec[4] + 35;
            break;
        case 48:
            return sec[4] + 46;
            break;
        case 52:  // validation
            return sec[4] + 25;
            break;
        case 57:
            n = number_of_mode(sec);
            if (n <= 0 || n == 65535) fatal_error2("PDT 4.57 bad number of mode", "");
            return sec[4] + 31 + 5 * n;
            break;
        case 70:
        case 71:
        case 72:
        case 73:
            return sec[4] + 27;
            break;
        case 20:
        case 30:
        case 31:
        case 32:
        case 1000:
        case 1001:
        case 1002:
        case 254:
            return NULL;
            break;
        default:
            fprintf(stderr, "code_table_4.5a: product definition template #%d not supported\n", pdt);
            return NULL;
            break;
    }
    return NULL;
}

int code_table_4_5b(unsigned char **sec) {
    unsigned char *p;
    p = code_table_4_5b_location(sec);
    if (p == NULL) return -1;
    return *p;
}

unsigned char *code_table_4_5b_location(unsigned char **sec) {
    int pdt, center, n;
    pdt = GB2_ProdDefTemplateNo(sec);
    center = GB2_Center(sec);

    if (center == JMA1 || center == JMA2) {
        switch (pdt) {
            case 50008:
                return sec[4] + 28;
                break;
            case 50009:
                return sec[4] + 28;
                break;
            case 51020:
            case 51021:
            case 51022:
            case 51122:
                return NULL;
                break;
            case 50010:
            case 50011:
                return sec[4] + 28;
                break;
            default:
                break;
        }
    }

    switch (pdt) {
        case 0:
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
        case 8:
        case 9:
        case 10:
        case 11:
        case 12:
        case 13:
        case 14:
        case 15:
        case 51:
        case 60:
        case 61:
        case 91:
        case 1100:
        case 1101:
            return sec[4] + 28;
            break;
        case 40:
        case 41:
        case 42:
        case 43:
            return sec[4] + 30;
            break;
        case 44:
            return sec[4] + 39;
            break;
        case 45:
        case 46:
        case 47:
            return sec[4] + 41;
            break;
        case 48:
            return sec[4] + 52;
            break;
        case 52:
            return NULL;
            break;
        case 57:
            n = number_of_mode(sec);
            if (n <= 0 || n == 65535) fatal_error2("PDT 4.57 bad number of mode", "");
            return sec[4] + 37 + 5 * n;
            break;
        case 70:
        case 71:
        case 72:
        case 73:
            return sec[4] + 33;
            break;
        case 20:
        case 30:
        case 31:
        case 32:
        case 1000:
        case 1001:
        case 1002:
        case 254:
            return NULL;
            break;
        default:
            fprintf(stderr, "code_table_4.5b: product definition template #%d not supported\n", pdt);
            return NULL;
            break;
    }
    return NULL;
}

int code_table_5_0(unsigned char **sec) { return (int)uint2(sec[5] + 9); }

int code_table_5_1(unsigned char **sec) {
    switch (code_table_5_0(sec)) {
        case 0:
        case 1:
        case 2:
        case 3:
        case 40:
        case 41:
        case 42:
            return (int)(sec[5][20]);
            break;
        default:
            return -1;
            break;
    }
    return -1;
}

int code_table_5_4(unsigned char **sec) {
    if (code_table_5_0(sec) < 2 || code_table_5_0(sec) > 3) return -1;
    return (int)(sec[5][21]);
}

int code_table_5_5(unsigned char **sec) {
    if (code_table_5_0(sec) < 2 || code_table_5_0(sec) > 3) return -1;
    return (int)(sec[5][22]);
}

int code_table_5_6(unsigned char **sec) {
    if (code_table_5_0(sec) != 3) return -1;
    return (int)(sec[5][47]);
}

int code_table_6_0(unsigned char **sec) { return (int)(sec[6][5]); }

/* from Code_Values.c */

/*
 * forecast_time_in_units
 *
 * v1.1 4/2015:  allow forecast time to be a signed quantity
 *      old: return unsigned value, ! code_4_4, return 0xffffffff;
 *      new: return signed value    ! code_4_4, return 0
 */

int forecast_time_in_units(unsigned char **sec) {
    unsigned char *p;
    int size;

    if ((p = forecast_time_in_units_location(sec, &size)) == NULL) return 0;
    if (size == 4) {
        return int4(p);
    } else if (size == 2) {
        return int4(p);
    } else
        fatal_error2("forecast_time_in_units size=?", "");
    return 0;
}

/*
 * WMO made a mistake and pdt 4.44 only uses 2 bytes for the forecast length
 */

unsigned char *forecast_time_in_units_location(unsigned char **sec, int *size) {
    unsigned char *code_4_4;

    code_4_4 = code_table_4_4_location(sec);
    if (code_4_4 == NULL) return NULL;
    *size = code_table_4_0(sec) == 44 ? 2 : 4;
    return code_4_4 + 1;
}

unsigned char *stat_proc_verf_time_location(unsigned char **sec) {
    int i, j, center, nb;
    i = code_table_4_0(sec);
    center = GB2_Center(sec);
    j = 0;
    if (i == 8)
        j = 34;
    else if (i == 9)
        j = 47;
    else if (i == 10)
        j = 35;
    else if (i == 11)
        j = 37;
    else if (i == 12)
        j = 36;
    else if (i == 13)
        j = 68;
    else if (i == 14)
        j = 64;
    else if (i == 34) {
        nb = sec[4][22];
        j = 26 + 11 * nb;
    } else if (i == 42)
        j = 36;
    else if (i == 43)
        j = 39;
    else if (i == 46)
        j = 47;
    else if (i == 47)
        j = 50;
    else if (i == 61)
        j = 44;
    else if (i == 50008 && (center == JMA1 || center == JMA2))
        j = 34;
    if (j == 0) return NULL;
    return sec[4] + j;
}

int number_of_mode(unsigned char **sec) {
    int pdt;
    pdt = code_table_4_0(sec);
    if (pdt == 57) return (int)uint2(sec[4] + 13);
    return -1;
}

/* from FlagTable.c */

int flag_table_3_3(unsigned char **sec) {
    unsigned char *p;
    p = flag_table_3_3_location(sec);
    if (p == NULL) return -1;
    return (int)*p;
}

unsigned char *flag_table_3_3_location(unsigned char **sec) {
    int grid_template, center;
    unsigned char *gds;

    grid_template = code_table_3_1(sec);
    center = GB2_Center(sec);
    gds = sec[3];
    switch (grid_template) {
        case 0:
        case 1:
        case 2:
        case 3:
        case 40:
        case 41:
        case 42:
        case 43:
        case 140:
        case 204:
            return gds + 54;
            break;
        case 4:
        case 5:
        case 10:
        case 12:
        case 20:
        case 30:
        case 31:
        case 90:
        case 110:
            return gds + 46;
            break;
#ifdef WMO_VALIDATION
        case 60:
            return gds + 71;
            break;
#endif
        case 32768:
            if (center == NCEP) return gds + 54;
            return NULL;
        case 32769:
            if (center == NCEP) return gds + 54;
            return NULL;
        case 40110:
            if ((center == JMA1) || (center == JMA2)) return gds + 46;
            return NULL;
            break;
        default:
            break;
    }
    return NULL;
}

int flag_table_3_4(unsigned char **sec) {
    unsigned char *p;
    p = flag_table_3_4_location(sec);
    if (p == NULL) return -1;
    return (int)*p;
}

unsigned char *flag_table_3_4_location(unsigned char **sec) {
    int grid_template, center;
    unsigned char *gds;

    gds = sec[3];
    grid_template = code_table_3_1(sec);
    center = GB2_Center(sec);

    switch (grid_template) {
        case 0:
        case 1:
        case 2:
        case 3:
        case 40:
        case 41:
        case 42:
        case 43:
            return gds + 71;
            break;
        case 4:
        case 5:
            return gds + 47;
            break;
        case 10:
        case 12:
            return gds + 59;
            break;
        case 20:
            return gds + 64;
            break;
        case 30:
        case 31:
            return gds + 64;
            break;
        case 50:
        case 51:
        case 52:
        case 53:
            /* spectral modes don't have scan order */
            return NULL;
            break;
        case 90:
        case 140:
            return gds + 63;
            break;
        case 110:
            return gds + 56;
            break;
        case 190:
        case 120:
            return gds + 38;
            break;
        case 204:
            return gds + 71;
            break;
        case 1000:
            return gds + 50;
            break;
#ifdef WMO_VALIDATION
        case 60:
            return gds + 72;
            break;
#endif
        case 32768:
            if (center == NCEP) return gds + 71;
            return NULL;
            break;
        case 32769:
            if (center == NCEP) return gds + 71;
            return NULL;
            break;
        case 40110:
            if ((center == JMA1) || (center == JMA2)) return gds + 56;
            return NULL;
            break;
        default:
            break;
    }
    return NULL;
}

/* from unpk_0.c */
#ifdef USE_OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 1
#endif
static unsigned int mask[] = {0, 1, 3, 7, 15, 31, 63, 127, 255};
static double shift[9] = {1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0};

void unpk_0(float *flt, unsigned char *bits0, unsigned char *bitmap0, int n_bits, unsigned int n, double ref,
            double scale, double dec_scale) {
    unsigned char *bits, *bitmap;

    int c_bits, j_bits, nthreads;
    unsigned int map_mask, bbits, i, j, k, n_missing, ndef, di;
    double jj;

    ref = ref * dec_scale;
    scale = scale * dec_scale;
    bits = bits0;
    bitmap = bitmap0;

    bbits = 0;

    /* assume integer has 32+ bits */
    /* optimized code for n_bits <= 25bits */
    if (n_bits <= 25) {
        n_missing = bitmap ? missing_points2(bitmap0, n) : 0;
        ndef = n - n_missing;

        // 1-cpu: rd_bitstream_flt(bits0, 0, flt+n_missing, n_bits, ndef);
        // 1-cpu: for (j = 0; j < ndef; j++) flt[j+n_missing] = ref + scale*flt[j+n_missing];

#pragma omp parallel private(i, j, k)
        {
#pragma omp single
            {
                nthreads = omp_get_num_threads();
                di = (ndef + nthreads - 1) / nthreads;
                di = ((di + 7) | 7) ^ 7;
            }
#pragma omp for
            for (i = 0; i < ndef; i += di) {
                k = ndef - i;
                if (k > di) k = di;
                rd_bitstream_flt(bits0 + (i / 8) * n_bits, 0, flt + n_missing + i, n_bits, k);
                for (j = i + n_missing; j < i + k + n_missing; j++) {
                    flt[j] = ref + scale * flt[j];
                }
            }
        }
        /*
#pragma omp parallel for private(i,j,k)
for (i = 0; i < ndef; i += CACHE_LINE_BITS) {
k  = ndef - i;
if (k > CACHE_LINE_BITS) k = CACHE_LINE_BITS;
rd_bitstream_flt(bits0 + (i/8)*n_bits, 0, flt+n_missing+i, n_bits, k);
for (j = i+n_missing; j < i+k+n_missing; j++) {
flt[j] = ref + scale*flt[j];
}
}
*/

        if (n_missing != 0) {
            j = n_missing;
            for (i = 0; i < n; i++) {
                /* check bitmap */
                if ((i & 7) == 0) bbits = *bitmap++;
                if (bbits & 128) {
                    flt[i] = flt[j++];
                } else {
                    flt[i] = UNDEFINED;
                }
                bbits = bbits << 1;
            }
        }
    } else {
        /* older unoptimized code, not often used */
        c_bits = 8;
        map_mask = 128;
        while (n-- > 0) {
            if (bitmap) {
                j = (*bitmap & map_mask);
                if ((map_mask >>= 1) == 0) {
                    map_mask = 128;
                    bitmap++;
                }
                if (j == 0) {
                    *flt++ = UNDEFINED;
                    continue;
                }
            }

            jj = 0.0;
            j_bits = n_bits;
            while (c_bits <= j_bits) {
                if (c_bits == 8) {
                    jj = jj * 256.0 + (double)(*bits++);
                    j_bits -= 8;
                } else {
                    jj = (jj * shift[c_bits]) + (double)(*bits & mask[c_bits]);
                    bits++;
                    j_bits -= c_bits;
                    c_bits = 8;
                }
            }
            if (j_bits) {
                c_bits -= j_bits;
                jj = (jj * shift[j_bits]) + (double)((*bits >> c_bits) & mask[j_bits]);
            }
            *flt++ = ref + scale * jj;
        }
    }
    return;
}

/* from unpk_complex.c */
int unpk_complex(unsigned char **sec, float *data, unsigned int ndata) {
    unsigned int i, j, n;
    int k, nbits, ref_group_length;
    unsigned char *p, *d, *mask_pointer;
    double ref_val, factor_10, factor_2, factor;
    float missing1, missing2;
    int n_sub_missing;
    int pack, offset;
    unsigned clocation;
    unsigned int ngroups, ref_group_width, nbit_group_width, len_last, npnts;
    int nbits_group_len, group_length_factor;
    int *group_refs, *group_widths, *group_lengths, *group_offset, *udata;
    unsigned int *group_clocation, *group_location;

    int m1, m2, mask, last, penultimate;
    int extra_vals[2];
    int min_val;
    int ctable_5_4, ctable_5_6, bitmap_flag, extra_octets;

    extra_vals[0] = extra_vals[1] = 0;
    pack = code_table_5_0(sec);
    if (pack != 2 && pack != 3) return 0;

    p = sec[5];
    ref_val = ieee2flt(p + 11);
    factor_2 = Int_Power(2.0, int2(p + 15));
    factor_10 = Int_Power(10.0, -int2(p + 17));
    ref_val *= factor_10;
    factor = factor_2 * factor_10;
    nbits = p[19];
    ngroups = uint4(p + 31);
    bitmap_flag = code_table_6_0(sec);
    ctable_5_6 = code_table_5_6(sec);

    if (pack == 3 && (ctable_5_6 != 1 && ctable_5_6 != 2)) fatal_error2("unsupported: code table 5.6", "");

    extra_octets = (pack == 2) ? 0 : sec[5][48];

    if (ngroups == 0) {
        if (bitmap_flag == 255) {
            for (i = 0; i < ndata; i++) data[i] = ref_val;
            return 0;
        }
        if (bitmap_flag == 0 || bitmap_flag == 254) {
            mask_pointer = sec[6] + 6;
            mask = 0;
            for (i = 0; i < ndata; i++) {
                if ((i & 7) == 0) mask = *mask_pointer++;
                data[i] = (mask & 128) ? ref_val : UNDEFINED;
                mask <<= 1;
            }
            return 0;
        }
        fatal_error2("unknown bitmap", "");
    }

    ctable_5_4 = code_table_5_4(sec);
    ref_group_width = p[35];
    nbit_group_width = p[36];
    ref_group_length = uint4(p + 37);
    group_length_factor = p[41];
    len_last = uint4(p + 42);
    nbits_group_len = p[46];

#ifdef DEBUG
    fprintf(stderr, "ctable 5.4 %d ref_group_width %u nbit_group_width %u ref_group_length %u group_length_factor %d\n",
            ctable_5_4, ref_group_width, nbit_group_width, ref_group_length, group_length_factor);
    fprintf(stderr, "len_last %u nbit_group_len %u\n", len_last, nbits_group_len);
#endif

    npnts = GB2_Sec5_nval(sec);  // number of defined points
    n_sub_missing = sub_missing_values(sec, &missing1, &missing2);

    // allocate group widths and group lengths
    group_refs = (int *)malloc(sizeof(unsigned int) * (size_t)ngroups);
    group_widths = (int *)malloc(sizeof(unsigned int) * (size_t)ngroups);
    group_lengths = (int *)malloc(sizeof(unsigned int) * (size_t)ngroups);
    group_location = (unsigned int *)malloc(sizeof(unsigned int) * (size_t)ngroups);
    group_clocation = (unsigned int *)malloc(sizeof(unsigned int) * (size_t)ngroups);
    group_offset = (int *)malloc(sizeof(unsigned int) * (size_t)ngroups);
    udata = (int *)malloc(sizeof(unsigned int) * (size_t)npnts);
    if (group_refs == NULL || group_widths == NULL || group_lengths == NULL || group_location == NULL ||
        group_clocation == NULL || group_offset == NULL || udata == NULL)
        fatal_error2("unpk_complex: memory allocation", "");

    // read any extra values
    d = sec[7] + 5;
    min_val = 0;
    if (extra_octets) {
        extra_vals[0] = uint_n(d, extra_octets);
        d += extra_octets;
        if (ctable_5_6 == 2) {
            extra_vals[1] = uint_n(d, extra_octets);
            d += extra_octets;
        }
        min_val = int_n(d, extra_octets);
        d += extra_octets;
    }

    if (ctable_5_4 != 1) fatal_error2("internal decode does not support code table 5.4", "");

#pragma omp parallel
    {
#pragma omp sections
        {
#pragma omp section
            {
                // read the group reference values
                rd_bitstream(d, 0, group_refs, nbits, ngroups);
            }

#pragma omp section
            {
                unsigned int i;
                // read the group widths

                rd_bitstream(d + (nbits * ngroups + 7) / 8, 0, group_widths, nbit_group_width, ngroups);
                for (i = 0; i < ngroups; i++) group_widths[i] += ref_group_width;
            }

#pragma omp section
            {
                unsigned int i;
                // read the group lengths

                if (ctable_5_4 == 1) {
                    rd_bitstream(d + (nbits * ngroups + 7) / 8 + (ngroups * nbit_group_width + 7) / 8, 0, group_lengths,
                                 nbits_group_len, ngroups - 1);

                    for (i = 0; i < ngroups - 1; i++) {
                        group_lengths[i] = group_lengths[i] * group_length_factor + ref_group_length;
                    }
                    group_lengths[ngroups - 1] = len_last;
                }
            }
        }

#pragma omp single
        {
            d += (nbits * ngroups + 7) / 8 + (ngroups * nbit_group_width + 7) / 8 + (ngroups * nbits_group_len + 7) / 8;

            // do a check for number of grid points and size
            clocation = offset = n = j = 0;
        }

#pragma omp sections
        {
#pragma omp section
            {
                unsigned int i;
                for (i = 0; i < ngroups; i++) {
                    group_location[i] = j;
                    j += group_lengths[i];
                    n += group_lengths[i] * group_widths[i];
                }
            }

#pragma omp section
            {
                unsigned int i;
                for (i = 0; i < ngroups; i++) {
                    group_clocation[i] = clocation;
                    clocation = clocation + group_lengths[i] * (group_widths[i] / 8) +
                                (group_lengths[i] / 8) * (group_widths[i] % 8);
                }
            }

#pragma omp section
            {
                unsigned int i;
                for (i = 0; i < ngroups; i++) {
                    group_offset[i] = offset;
                    offset += (group_lengths[i] % 8) * (group_widths[i] % 8);
                }
            }
        }
    }

    if (j != npnts) fatal_error2("bad complex packing: n points", "");
    if (d + (n + 7) / 8 - sec[7] != GB2_Sec7_size(sec)) fatal_error2("complex unpacking size mismatch old test", "");

    if (d + clocation + (offset + 7) / 8 - sec[7] != GB2_Sec7_size(sec))
        fatal_error2("complex unpacking size mismatch", "");

#pragma omp parallel for private(i) schedule(static)
    for (i = 0; i < ngroups; i++) {
        group_clocation[i] += (group_offset[i] / 8);
        group_offset[i] = (group_offset[i] % 8);

        rd_bitstream(d + group_clocation[i], group_offset[i], udata + group_location[i], group_widths[i],
                     group_lengths[i]);
    }

    // handle substitute, missing values and reference value
    if (n_sub_missing == 0) {
#pragma omp parallel for private(i, k, j)
        for (i = 0; i < ngroups; i++) {
            j = group_location[i];
            for (k = 0; k < group_lengths[i]; k++) {
                udata[j++] += group_refs[i];
            }
        }
    } else if (n_sub_missing == 1) {
#pragma omp parallel for private(i, m1, k, j)
        for (i = 0; i < ngroups; i++) {
            j = group_location[i];
            if (group_widths[i] == 0) {
                m1 = (1 << nbits) - 1;
                if (m1 == group_refs[i]) {
                    for (k = 0; k < group_lengths[i]; k++) udata[j++] = INT_MAX;
                } else {
                    for (k = 0; k < group_lengths[i]; k++) udata[j++] += group_refs[i];
                }
            } else {
                m1 = (1 << group_widths[i]) - 1;
                for (k = 0; k < group_lengths[i]; k++) {
                    if (udata[j] == m1)
                        udata[j] = INT_MAX;
                    else
                        udata[j] += group_refs[i];
                    j++;
                }
            }
        }
    } else if (n_sub_missing == 2) {
#pragma omp parallel for private(i, j, k, m1, m2)
        for (i = 0; i < ngroups; i++) {
            j = group_location[i];
            if (group_widths[i] == 0) {
                m1 = (1 << nbits) - 1;
                m2 = m1 - 1;
                if (m1 == group_refs[i] || m2 == group_refs[i]) {
                    for (k = 0; k < group_lengths[i]; k++) udata[j++] = INT_MAX;
                } else {
                    for (k = 0; k < group_lengths[i]; k++) udata[j++] += group_refs[i];
                }
            } else {
                m1 = (1 << group_widths[i]) - 1;
                m2 = m1 - 1;
                for (k = 0; k < group_lengths[i]; k++) {
                    if (udata[j] == m1 || udata[j] == m2)
                        udata[j] = INT_MAX;
                    else
                        udata[j] += group_refs[i];
                    j++;
                }
            }
        }
    }

    // post processing

    if (pack == 3) {
        if (ctable_5_6 == 1) {
            last = extra_vals[0];
            i = 0;
            while (i < npnts) {
                if (udata[i] == INT_MAX)
                    i++;
                else {
                    udata[i++] = extra_vals[0];
                    break;
                }
            }
            while (i < npnts) {
                if (udata[i] == INT_MAX)
                    i++;
                else {
                    udata[i] += last + min_val;
                    last = udata[i++];
                }
            }
        } else if (ctable_5_6 == 2) {
            penultimate = extra_vals[0];
            last = extra_vals[1];

            i = 0;
            while (i < npnts) {
                if (udata[i] == INT_MAX)
                    i++;
                else {
                    udata[i++] = extra_vals[0];
                    break;
                }
            }
            while (i < npnts) {
                if (udata[i] == INT_MAX)
                    i++;
                else {
                    udata[i++] = extra_vals[1];
                    break;
                }
            }
            for (; i < npnts; i++) {
                if (udata[i] != INT_MAX) {
                    udata[i] = udata[i] + min_val + last + last - penultimate;
                    penultimate = last;
                    last = udata[i];
                }
            }
        } else
            fatal_error2("Unsupported: code table 5.6", "");
    }

    // convert to float

    if (bitmap_flag == 255) {
#pragma omp parallel for schedule(static) private(i)
        for (i = 0; i < ndata; i++) {
            data[i] = (udata[i] == INT_MAX) ? UNDEFINED : ref_val + udata[i] * factor;
        }
    } else if (bitmap_flag == 0 || bitmap_flag == 254) {
        n = 0;
        mask = 0;
        mask_pointer = sec[6] + 6;
        for (i = 0; i < ndata; i++) {
            if ((i & 7) == 0) mask = *mask_pointer++;
            if (mask & 128) {
                if (udata[n] == INT_MAX)
                    data[i] = UNDEFINED;
                else
                    data[i] = ref_val + udata[n] * factor;
                n++;
            } else
                data[i] = UNDEFINED;
            mask <<= 1;
        }
    } else
        fatal_error2("unknown bitmap", "");

    free(group_refs);
    free(group_widths);
    free(group_lengths);
    free(group_location);
    free(group_clocation);
    free(group_offset);
    free(udata);

    return 0;
}

/* from unpk_run_length.c */

int unpk_run_length(unsigned char **sec, float *data, unsigned int ndata) {
    int i, k, decimal_scale, n_bits;
    int mv, mvl;
    unsigned int j, ncheck;

    double *levels, dec_factor;
    int size_compressed_data, nvals, *vals;
    int v, n, factor, range;
    int bitmap_flag;
    unsigned char *mask_pointer;
    const unsigned int mask[] = {128, 64, 32, 16, 8, 4, 2, 1};

    if (code_table_5_0(sec) != 200) return 0;  // only for DRT 5.200

    if (GB2_Sec3_npts(sec) != ndata) fatal_error2("Run-length decoding: programming error", "");
    n_bits = (int)sec[5][11];
    mv = (int)uint2(sec[5] + 12);
    mvl = (int)uint2(sec[5] + 14);
    decimal_scale = (int)sec[5][16];
    if (decimal_scale > 127) {  // convert to signed negative values
        decimal_scale = -(decimal_scale - 128);
    }
    dec_factor = Int_Power(10.0, -decimal_scale);
    if (mv > mvl) fatal_error2("Run-length decoding: ", "");

#ifdef DEBUG
    printf(" n_bits=%d mv=%d mvl=%d decimal_scale=%d\n", n_bits, mv, mvl, decimal_scale);
#endif

    size_compressed_data = GB2_Sec7_size(sec) - 5;
    nvals = (size_compressed_data * 8) / n_bits;

#ifdef DEBUG
    printf(" size_compressed data %d ndata %u\n", size_compressed_data, ndata);
#endif

    if (nvals == 0 || mv == 0) {
        for (j = 0; j < ndata; j++) {
            data[j] = UNDEFINED;
        }
        return 0;
    }

    levels = (double *)malloc((mvl + 1) * sizeof(double));
    vals = (int *)malloc(nvals * sizeof(int));
    if (levels == NULL || vals == NULL) fatal_error2("Run-length decoding: memory allocation", "");

    /* levels[0..mvl]: levels[0] = UNDEFINED; levels[1..mvl] from table */
    levels[0] = UNDEFINED;
    for (i = 1; i <= mvl; i++) {
        levels[i] = int2(sec[5] + 15 + i * 2) * dec_factor;
    }

#ifdef DEBUG
    for (i = 0; i <= mvl; i++) {
        printf(" lvls[%d] = %lf ", i, levels[i]);
        if (i % 4 == 0) printf("\n");
    }
    printf("\n");
#endif

    rd_bitstream(sec[7] + 5, 0, vals, n_bits, nvals);

    ncheck = i = 0;
    range = (1 << n_bits) - 1 - mv;
    if (range <= 0) fatal_error2("unpk_running_length: range error", "");

    j = 0;
    mask_pointer = sec[6] + 6;
    bitmap_flag = code_table_6_0(sec);
    if (bitmap_flag == 254) bitmap_flag = 0;
    if (bitmap_flag != 0 && bitmap_flag != 255) fatal_error2("unpk_running_length: unsupported bitmap", "");

    while (i < nvals) {
        if (vals[i] > mv) fatal_error2("Run-length decoding: error", "");
        v = vals[i++];

        /* figure out any repeat factor */
        n = 1;
        factor = 1;
        // 12/2014 while (vals[i] > mv && i < nvals) {
        while (i < nvals && vals[i] > mv) {
            n += factor * (vals[i] - mv - 1);
            factor = factor * range;
            i++;
        }

        ncheck += n;
        if (ncheck > ndata) fatal_error2("unpk_run_length: ncheck ,", "");

        if (bitmap_flag != 0) {
            for (k = 0; k < n; k++) {
                data[j++] = levels[v];
            }
        } else {
            for (k = 0; k < n; k++) {
                while (mask_pointer[j >> 3] & mask[j & 7]) {
                    data[j++] = UNDEFINED;
                }
                data[j++] = levels[v];
            }
        }
    }
    if (j != ndata) fatal_error2("Run-length decoding: bitmap problem", "");
    free(levels);
    free(vals);
    return 0;
}

/* from bitstream.c */

static unsigned int ones[] = {0, 1, 3, 7, 15, 31, 63, 127, 255};

void rd_bitstream(unsigned char *p, int offset, int *u, int n_bits, int n) {
    unsigned int tbits;
    int i, t_bits, new_t_bits;

    // not the best of tests

    if (INT_MAX <= 2147483647 && n_bits > 31) fatal_error2("rd_bitstream: n_bits", "");

    if (offset < 0 || offset > 7) fatal_error2("rd_bitstream: illegal offset", "");

    if (n_bits == 0) {
        for (i = 0; i < n; i++) {
            u[i] = 0;
        }
        return;
    }

    t_bits = 8 - offset;
    tbits = (*p++) & ones[t_bits];

    for (i = 0; i < n; i++) {
        while (n_bits - t_bits >= 8) {
            t_bits += 8;
            tbits = (tbits << 8) | *p++;
        }

        if (n_bits > t_bits) {
            new_t_bits = 8 - (n_bits - t_bits);
            u[i] = (int)((tbits << (n_bits - t_bits) | (*p >> new_t_bits)));
            t_bits = new_t_bits;
            tbits = *p++ & ones[t_bits];
        } else if (n_bits == t_bits) {
            u[i] = (int)tbits;
            tbits = t_bits = 0;
        } else {
            t_bits -= n_bits;
            u[i] = (int)(tbits >> t_bits);
            tbits = tbits & ones[t_bits];
        }
    }
}

void rd_bitstream_flt(unsigned char *p, int offset, float *u, int n_bits, int n) {
    unsigned int tbits;
    int i, t_bits, new_t_bits;

    // not the best of tests

    if (INT_MAX <= 2147483647 && n_bits > 31) fatal_error2("rd_bitstream: n_bits", "");

    if (offset < 0 || offset > 7) fatal_error2("rd_bitstream_flt: illegal offset", "");

    if (n_bits == 0) {
        for (i = 0; i < n; i++) {
            u[i] = 0.0;
        }
        return;
    }

    t_bits = 8 - offset;
    tbits = (*p++) & ones[t_bits];

    for (i = 0; i < n; i++) {
        while (n_bits - t_bits >= 8) {
            t_bits += 8;
            tbits = (tbits << 8) | *p++;
        }

        if (n_bits > t_bits) {
            new_t_bits = 8 - (n_bits - t_bits);
            u[i] = (int)((tbits << (n_bits - t_bits) | (*p >> new_t_bits)));
            t_bits = new_t_bits;
            tbits = *p++ & ones[t_bits];
        } else if (n_bits == t_bits) {
            u[i] = (float)tbits;
            tbits = t_bits = 0;
        } else {
            t_bits -= n_bits;
            u[i] = (float)(tbits >> t_bits);
            tbits = tbits & ones[t_bits];
        }
    }
}

#endif
