#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <stdbool.h>

#include "lndsr.h"
#include "keyvalue.h"
#include "const.h"
#include "param.h"
#include "input.h"
#include "prwv_input.h"
#include "lut.h"
#include "output.h"
#include "sr.h"
#include "ar.h"
#include "error.h"
#include "clouds.h"
#include "read_level1_qa.h"
#include "read_grib_tools.h"
#include "sixs_runs.h"

#define AERO_NB_BANDS 3
#define SP_INDEX    0
#define WV_INDEX    1
#define ATEMP_INDEX 2
#define OZ_INDEX    0
#define DEBUG_FLAG  0
/* #define DEBUG_AR 1 */
/* #define DEBUG_CLD 1 */

/* DEM Definition: U_char format, 1 count = 100 meters */
/* 0 = 0 meters */

#define DEMFILE "CMGDEM.hdf"
#define DEM_NBLAT 3600
#define DEM_DLAT 0.05
#define DEM_LATMIN (-90.0)
#define DEM_LATMAX 90.0
#define DEM_NBLON 7200
#define DEM_DLON 0.05
#define DEM_LONMIN (-180.0)
#define DEM_LONMAX 180.0
#define P_DFTVALUE 1013.0

/* Type definitions */

#ifdef DEBUG_AR
FILE *fd_ar_diags = NULL;
int diags_il_ar;
#endif
#ifdef DEBUG_CLD
FILE *fd_cld_diags = NULL;
#endif
/* reading the DEM in hdf */
int32 sds_file_id,sds_id,status;
char sds_name[256];
int32 sds_index;
int32 dim_sizes[2],start[2],stride[2],edges[2];
int32 data_type,n_attrs,rank;
/* Prototypes */
#ifndef  HPUX
#define chand chand_
#define csalbr csalbr_
#endif
void chand(float *phi,float *muv,float *mus,float *tau_ray,float *actual_rho_ray);
void csalbr(float *tau_ray,float *actual_S_r);
static void update_atmos_coefs(atmos_t *atmos_coef, Ar_gridcell_t *ar_gridcell,
                               sixs_tables_t *sixs_tables, int ***line_ar,
                               Lut_t *lut, int nband, int bkgd_aerosol);
void update_gridcell_atmos_coefs(int ipt, atmos_t *atmos_coef,
                                Ar_gridcell_t *ar_gridcell, 
                                sixs_tables_t *sixs_tables, int line_ar,
                                Lut_t *lut, int nband, int bkgd_aerosol);
float calcuoz(short jday,float flat);
float get_dem_spres(short *dem,float lat,float lon);

#ifdef SAVE_6S_RESULTS
#define SIXS_RESULTS_FILENAME "SIXS_RUN_RESULTS.TXT"
int read_6S_results_from_file(char *filename,sixs_tables_t *sixs_tables);
int write_6S_results_to_file(char *filename,sixs_tables_t *sixs_tables);
#endif
void sun_angles (short jday,float gmt,float flat,float flon,float *ts,float *fs);
/* Functions */

int main (int argc, char *argv[]) {
    Param_t *param = NULL;
    Input_t *input = NULL, *input_b6 = NULL;
    InputPrwv_t *prwv_input = NULL;
    InputOzon_t *ozon_input = NULL;
    Lut_t *lut = NULL;
    Output_t *output = NULL;
    int i,j,il, is,ib,ifree;
    int il_start, il_end, il_ar, il_region, is_ar;
    uint16_t *line_out[NBAND_SR_MAX];
    uint16_t *line_out_buf = NULL;
    uint16_t ***line_in = NULL;
    uint16_t **line_in_band_buf = NULL;
    uint16_t *line_in_buf = NULL;
    int ***line_ar = NULL;
    int **line_ar_band_buf = NULL;
    int *line_ar_buf = NULL;
    uint16_t** b6_line = NULL;
    uint16_t* b6_line_buf = NULL;
    float *atemp_line = NULL;
    uint16_t** qa_line = NULL;      /* bqa_pixel */
    uint16_t* qa_line_buf = NULL;
    uint16_t** qa2_line = NULL;     /* bqa_radsat */
    uint16_t* qa2_line_buf = NULL;
    char **ddv_line = NULL;
    char *ddv_line_buf = NULL;
    char **rot_cld[3],**ptr_rot_cld[3],**ptr_tmp_cld;
    char **rot_cld_block_buf = NULL;
    char *rot_cld_buf = NULL;
    char envi_file[STR_SIZE]; /* name of the output ENVI header file */
    char *cptr = NULL;        /* pointer to the file extension */

    Sr_stats_t sr_stats;
    Ar_stats_t ar_stats;
    Ar_gridcell_t ar_gridcell;
    float *prwv_in[NBAND_PRWV_MAX];
    float *prwv_in_buf = NULL;
    int *ozon_in = NULL;
    float corrected_sun_az;   /* (degrees) sun azimuth angle has been corrected
                                  for polar scenes that are ascending or
                                  flipped */
  
    int nbpts;
    int inter_aot;            /* atmospheric opacity */
    float scene_gmt;

    Geoloc_t *space = NULL;
    Space_def_t space_def;
    char *dem_name = NULL;
    Img_coord_float_t img;
    Img_coord_int_t loc;
    Geo_coord_t geo;

    t_ncep_ancillary anc_O3,anc_WV,anc_SP,anc_ATEMP;
    double sum_spres_anc,sum_spres_dem;
    int nb_spres_anc,nb_spres_dem;
    float tmpflt_arr[4];
    double coef;
    int tmpint;
    int osize;
    int debug_flag;

    sixs_tables_t sixs_tables;
    float center_lat,center_lon;
    char tmpfilename[128];
    FILE *fdtmp/*, *fdtmp2 */;
    int tmpid;                  /* file ID for temporary file (ID not used) */
  
    short *dem_array;
    int dem_available;
  
    cld_diags_t cld_diags;

    float flat,flon/*,fts,ffs*/;
    double delta_y,delta_x;
    float adjust_north;
  
    float sum_value,sumsq_value;
    int no_ozone_file;
    short jday;

    Espa_internal_meta_t xml_metadata;  /* XML metadata structure */
    Espa_global_meta_t *gmeta = NULL;   /* pointer to global meta */
    Envi_header_t envi_hdr;             /* output ENVI header information */

    atmos_t atmos_coef;          /* atmospheric coefficients */
    atmos_t atmos_coef_storage;  /* storage space for temporary atmospheric
                                    coefficients */

    debug_flag= DEBUG_FLAG;
    no_ozone_file=0;
  
    /* Read the parameters from the command-line and input parameter file */
    param = GetParam(argc, argv);
    if (param == NULL) EXIT_ERROR("getting runtime parameters", "main");

    printf ("\nRunning lndsr ....\n");
#ifdef _OPENMP
    printf ("lndsr was compiled using multi-threading via -fopenmp\n");
#else
    printf ("lndsr was NOT compiled using multi-threading\n");
#endif

    /* Validate the input metadata file */
    if (validate_xml_file (param->input_xml_file_name) != SUCCESS)
    {  /* Error messages already written */
        EXIT_ERROR("Unable to validate XML file", "main");
    }

    /* Initialize the metadata structure */
    init_metadata_struct (&xml_metadata);

    /* Parse the metadata file into our internal metadata structure; also
       allocates space as needed for various pointers in the global and band
       metadata */
    if (parse_metadata (param->input_xml_file_name, &xml_metadata) != SUCCESS)
    {  /* Error messages already written */
        EXIT_ERROR("parsing XML file", "main");
    }
    gmeta = &xml_metadata.global; /* pointer to global meta */

    /* Open input files; grab QA band for reflectance band */
    input = OpenInput(&xml_metadata, false /* not thermal */);
    if (input == NULL)
        EXIT_ERROR("bad input file", "main");

    input_b6 = OpenInput(&xml_metadata, true /* thermal */);
    if (input_b6 == NULL) {
        param->thermal_band = false;
        printf ("WARNING: no TOA brightness temp band available. "
            "Processing without.");
    }
    else
        param->thermal_band = true;

    if (param->num_prwv_files > 0  && param->num_ncep_files > 0) {
        EXIT_ERROR("both PRWV and PRWV_FIL files specified", "main");
    }

    /* The surface reflectance algorithm cannot be implemented for solar
       zenith angles greater than 76 degrees.  Need to flag if the current
       scene falls into that category. */
    if (input->meta.sun_zen * DEG > 76.0) {
        EXIT_ERROR ("Solar zenith angle is too large to allow for surface "
            "reflectance processing.  Corrections must be limited to top-of-"
            "atmosphere reflectance and brightness temperature corrections. "
            "Use the --process_sr=False command-line argument when running "
            "do_ledaps.py.", "main");
    }

    /* Open prwv input file */
    if (param->num_prwv_files > 0) {
        prwv_input = OpenInputPrwv(param->prwv_file_name);
        if (prwv_input==NULL) EXIT_ERROR("bad input prwv file","main");

        osize= 3 * (prwv_input->size.ntime*prwv_input->size.nlat*
            prwv_input->size.nlon);

        prwv_in_buf = calloc(osize, sizeof(float));
        if (prwv_in_buf == NULL)
            EXIT_ERROR("allocating input prwv buffer", "main");
        prwv_in[0] = prwv_in_buf;
        for (ib = 1; ib < prwv_input->nband; ib++)
            prwv_in[ib] = prwv_in[ib - 1] + (prwv_input->size.ntime*
                prwv_input->size.nlat*prwv_input->size.nlon);

        for (ib = 0; ib < prwv_input->nband; ib++) {
            if (!GetInputPrwv(prwv_input, ib, prwv_in[ib]))
                EXIT_ERROR("reading input prwv data", "main");
        }

        /**** ozone ***/
        if ( param->num_ozon_files<1 )
            no_ozone_file=1;
        else {
            ozon_input = OpenInputOzon(param->ozon_file_name);
            if (ozon_input==NULL) EXIT_ERROR("bad input ozon file", "main");
  
            osize= (ozon_input->size.ntime*ozon_input->size.nlat*
                ozon_input->size.nlon);

            ozon_in = calloc(osize, sizeof(int));
            if (ozon_in == NULL)
                EXIT_ERROR("allocating input ozone buffer", "main");

            if (!GetInputOzon(ozon_input, 0, ozon_in))
                EXIT_ERROR("reading input ozone data", "main");
        }
    }

    /* Get Lookup table, based on reflectance information */
    lut = GetLut(input->nband, &input->meta, &input_b6->meta, &input->size);
    if (lut == NULL) EXIT_ERROR("bad lut file", "main");

    /* Get geolocation space definition */
    if (!get_geoloc_info(&xml_metadata, &space_def))
        EXIT_ERROR("getting space metadata from XML file", "main");
    space = setup_mapping(&space_def);
    if (space == NULL)
        EXIT_ERROR("getting setting up geolocation mapping", "main");

    printf ("Number of input bands: %d\n", input->nband);
    printf ("Number of input lines: %d\n", input->size.l);
    printf ("Number of input samples: %d\n", input->size.s);

    /* If the scene is an ascending polar scene (flipped upside down), then
       the solar azimuth needs to be adjusted by 180 degrees.  The scene in
       this case would be north down and the solar azimuth is based on north
       being up. */
    corrected_sun_az = input->meta.sun_az * DEG;
    if (gmeta->ul_corner[0] < gmeta->lr_corner[0]) {
        corrected_sun_az += 180.0;
        if (corrected_sun_az > 360.0)
            corrected_sun_az -= 360.0;
        printf ("Polar or ascending scene.  Readjusting solar azimuth by "
            "180 degrees.\n  New value: %f radians (%f degrees)\n",
            corrected_sun_az*RAD, corrected_sun_az);
    }

    /* Open the output files and set up the necessary information for appending
       to the XML file */
    output = OpenOutput(&xml_metadata, input, param, lut);
    if (output == NULL) EXIT_ERROR("opening output file", "main");

    /* Open diagnostics files if needed */
#ifdef DEBUG_AR
    strcpy(tmpfilename,param->output_file_name);
    strcat(tmpfilename,".DEBUG_AR");
    fd_ar_diags=fopen(tmpfilename,"w");
    if (fd_ar_diags != NULL) {
        fprintf(fd_ar_diags,"cell_row cell_col total_nb_samples avg_b1 std_b1 avg_b2 std_b2 avg_b3 std_b3 avg_b7 std_b7 szen vzen relaz wv ozone spres fraction_water fraction_clouds fraction_cldshadow fraction_snow spres_ratio tau_ray corrected_T_ray corrected_Sr measured_rho_b1 simulated_b1_01 simulated_b1_02 simulated_b1_03 simulated_b1_04 simulated_b1_05 simulated_b1_06 simulated_b1_07 simulated_b1_08 simulated_b1_09 simulated_b1_10 simulated_b1_11 simulated_b1_12 simulated_b1_13 simulated_b1_14 simulated_b1_15 aot_index coef aot_value new_aot ratio_neg_red\n"); 
    }
#endif
#ifdef DEBUG_CLD
    strcpy(tmpfilename,"tempfile.DEBUG_CLD");
    fd_cld_diags=fopen(tmpfilename,"w");
    if (fd_cld_diags != NULL) {
        fprintf(fd_cld_diags,"cell_row cell_col nb_samples airtemp_2m "
            "avg_t6_clear std_t6_clear avg_b7_clear std_b7_clear\n"); 
    }
#endif

    /* Allocate memory for input lines */
    line_in = calloc(lut->ar_region_size.l, sizeof(uint16_t **));
    if (line_in == NULL) 
        EXIT_ERROR("allocating input line buffer (a)", "main");

    line_in_band_buf = calloc(lut->ar_region_size.l * input->nband,
        sizeof(uint16_t *));
    if (line_in_band_buf == NULL) 
        EXIT_ERROR("allocating input line buffer (b)", "main");

    line_in_buf = calloc(input->size.s * lut->ar_region_size.l * input->nband,
        sizeof(uint16_t));
    if (line_in_buf == NULL) 
        EXIT_ERROR("allocating input line buffer (c)", "main");

    for (il = 0; il < lut->ar_region_size.l; il++) {
        line_in[il] = line_in_band_buf;
        line_in_band_buf += input->nband;
        for (ib = 0; ib < input->nband; ib++) {
            line_in[il][ib] = line_in_buf;
            line_in_buf += input->size.s;
        }
    }

    /* Allocate memory for qa line */
    qa_line = calloc(lut->ar_region_size.l,sizeof(uint16_t *));
    if (qa_line == NULL) EXIT_ERROR("allocating qa line", "main");
    qa_line_buf = calloc(input->size.s * lut->ar_region_size.l,
        sizeof(uint16_t));
    if (qa_line_buf == NULL) EXIT_ERROR("allocating qa line buffer", "main");
    qa2_line = calloc(lut->ar_region_size.l,sizeof(uint16_t *));
    if (qa2_line == NULL) EXIT_ERROR("allocating qa line", "main");
    qa2_line_buf = calloc(input->size.s * lut->ar_region_size.l,
        sizeof(uint16_t));
    if (qa2_line_buf == NULL) EXIT_ERROR("allocating qa line buffer", "main");
    for (il = 0; il < lut->ar_region_size.l; il++) {
        qa_line[il]=qa_line_buf;
        qa_line_buf += input->size.s;
        qa2_line[il]=qa2_line_buf;
        qa2_line_buf += input->size.s;
    }

    /* Allocate memory for one band 6 line */
    if (param->thermal_band) {
        b6_line = calloc(lut->ar_region_size.l,sizeof(uint16_t *));
        if (b6_line == NULL) EXIT_ERROR("allocating b6 line", "main");
        b6_line_buf = calloc(input_b6->size.s * lut->ar_region_size.l,
            sizeof(uint16_t));
        if (b6_line_buf == NULL)
            EXIT_ERROR("allocating b6 line buffer", "main");
        for (il = 0; il < lut->ar_region_size.l; il++) {
            b6_line[il]=b6_line_buf;
            b6_line_buf += input_b6->size.s;
        }
    }

    /* Allocate memory for one air temperature line */
    atemp_line = calloc(input->size.s,sizeof(float));
    if (atemp_line == NULL) EXIT_ERROR("allocating atemp line", "main");

    /* Allocate memory for ddv line */
    ddv_line = calloc(lut->ar_region_size.l,sizeof(char *));
    if (ddv_line == NULL) EXIT_ERROR("allocating ddv line", "main");
    ddv_line_buf = calloc(input->size.s * lut->ar_region_size.l, sizeof(char));
    if (ddv_line_buf == NULL) EXIT_ERROR("allocating ddv line buffer", "main");
    for (il = 0; il < lut->ar_region_size.l; il++) {
        ddv_line[il]=ddv_line_buf;
        ddv_line_buf += input->size.s;
    }

    /* Allocate memory for rotating cloud buffer */
    rot_cld_buf=calloc(input->size.s*lut->ar_region_size.l*3, sizeof(char));
    if (rot_cld_buf == NULL) 
        EXIT_ERROR("allocating roatting cloud buffer (a)", "main");
    rot_cld_block_buf=calloc(lut->ar_region_size.l*3, sizeof(char *));
    if (rot_cld_block_buf == NULL) 
        EXIT_ERROR("allocating rotating cloud buffer (b)", "main");

    for (ib = 0; ib < 3; ib++) {
        rot_cld[ib]=rot_cld_block_buf;
        rot_cld_block_buf += lut->ar_region_size.l;
        for (il = 0; il < lut->ar_region_size.l; il++) {
            rot_cld[ib][il]=rot_cld_buf;
            rot_cld_buf+=input->size.s;
        }
    }

    /* Allocate memory for ar_gridcell */
    ar_gridcell.nbrows=lut->ar_size.l;
    ar_gridcell.nbcols=lut->ar_size.s;
    ar_gridcell.lat=calloc(lut->ar_size.s * lut->ar_size.l, sizeof(float));
    if (ar_gridcell.lat == NULL)
        EXIT_ERROR("allocating ar_gridcell.lat", "main");
    ar_gridcell.lon=calloc(lut->ar_size.s * lut->ar_size.l, sizeof(float));
    if (ar_gridcell.lon == NULL)
        EXIT_ERROR("allocating ar_gridcell.lon", "main");
    ar_gridcell.sun_zen=calloc(lut->ar_size.s * lut->ar_size.l, sizeof(float));
    if (ar_gridcell.sun_zen == NULL)
        EXIT_ERROR("allocating ar_gridcell.sun_zen", "main");
    ar_gridcell.view_zen=calloc(lut->ar_size.s * lut->ar_size.l,sizeof(float));
    if (ar_gridcell.view_zen == NULL)
        EXIT_ERROR("allocating ar_gridcell.view_zen", "main");
    ar_gridcell.rel_az=calloc(lut->ar_size.s * lut->ar_size.l,sizeof(float));
    if (ar_gridcell.rel_az == NULL)
        EXIT_ERROR("allocating ar_gridcell.rel_az", "main");
    ar_gridcell.wv=calloc(lut->ar_size.s * lut->ar_size.l,sizeof(float));
    if (ar_gridcell.wv == NULL)
        EXIT_ERROR("allocating ar_gridcell.wv", "main");
    ar_gridcell.spres=calloc(lut->ar_size.s * lut->ar_size.l,sizeof(float));
    if (ar_gridcell.spres == NULL)
        EXIT_ERROR("allocating ar_gridcell.spres", "main");
    ar_gridcell.ozone=calloc(lut->ar_size.s * lut->ar_size.l,sizeof(float));
    if (ar_gridcell.ozone == NULL)
        EXIT_ERROR("allocating ar_gridcell.ozone", "main");
    ar_gridcell.spres_dem=calloc(lut->ar_size.s * lut->ar_size.l,sizeof(float));
    if (ar_gridcell.spres_dem == NULL)
        EXIT_ERROR("allocating ar_gridcell.spres_dem", "main");

    /* Allocate memory for output lines */
    line_out_buf = calloc(output->size.s * output->nband_out, sizeof(uint16_t));
    if (line_out_buf == NULL) 
        EXIT_ERROR("allocating output line buffer", "main");
    line_out[0] = line_out_buf;
    for (ib = 1; ib < output->nband_out; ib++)
        line_out[ib] = line_out[ib - 1] + output->size.s;

    /* Allocate memory for the aerosol lines */
    line_ar = calloc(lut->ar_size.l, sizeof(int **));
    if (line_ar == NULL) 
        EXIT_ERROR("allocating aerosol line buffer (a)", "main");

    line_ar_band_buf = calloc(lut->ar_size.l * AERO_NB_BANDS, sizeof(int *));
    if (line_ar_band_buf == NULL) 
        EXIT_ERROR("allocating aerosol line buffer (b)", "main");

    line_ar_buf = calloc(lut->ar_size.l * lut->ar_size.s * AERO_NB_BANDS, 
        sizeof(int));
    if (line_ar_buf == NULL) 
        EXIT_ERROR("allocating aerosol line buffer (c)", "main");

    for (il = 0; il < lut->ar_size.l; il++) {
        line_ar[il] = line_ar_band_buf;
        line_ar_band_buf += AERO_NB_BANDS;
        for (ib = 0; ib < AERO_NB_BANDS; ib++) {
            line_ar[il][ib] = line_ar_buf;
            line_ar_buf += lut->ar_size.s;
        }
    }

    /* Initialize the statistics */
    ar_stats.nfill = 0;
    ar_stats.first = true;

    for (ib = 0; ib < output->nband_out; ib++) {
        sr_stats.nfill[ib] = 0;
        sr_stats.nsatu[ib] = 0;
        sr_stats.nout_range[ib] = 0;
        sr_stats.first[ib] = true;
    }

    /****
    Get center lat lon and deviation from true north
    ****/
    img.l=input->size.l/2.; 
    img.s=input->size.s/2.; 
    img.is_fill=false;
    if (!from_space(space, &img, &geo))
        EXIT_ERROR("mapping from space (0)", "main");
    center_lat=geo.lat * DEG;
    center_lon=geo.lon * DEG;

    /* Compute scene gmt time */
    if ( (input->meta.acq_date.hour !=0 ) ||
         (input->meta.acq_date.minute != 0 ) ||
         (input->meta.acq_date.second !=0))
        scene_gmt=input->meta.acq_date.hour + input->meta.acq_date.minute/60. +
          input->meta.acq_date.second/3600.;
    else
        scene_gmt=10.5-center_lon/15.;
    if ( scene_gmt < 0.)
        scene_gmt=scene_gmt+24.;

    printf ("Acquisition Time: %02d:%02d:%fZ\n", input->meta.acq_date.hour,
        input->meta.acq_date.minute, input->meta.acq_date.second);

    /* Read PRWV Data */
    if ( param->num_prwv_files > 0  ) {
        if (!get_prwv_anc(&anc_SP,prwv_input,prwv_in[SP_INDEX],SP_INDEX))
            EXIT_ERROR("Can't get PRWV SP data","main");
        if (!get_prwv_anc(&anc_WV,prwv_input,prwv_in[WV_INDEX],WV_INDEX))
            EXIT_ERROR("Can't get PRWV WV data","main");
        if (!get_prwv_anc(&anc_ATEMP,prwv_input,prwv_in[ATEMP_INDEX],
            ATEMP_INDEX))
            EXIT_ERROR("Can't get PRWV ATEMP data","main");
        if (!no_ozone_file)
            if (!get_ozon_anc(&anc_O3,ozon_input,ozon_in,OZ_INDEX))
                EXIT_ERROR("Can't get OZONE data","main");
    }
    else if ( param->num_ncep_files > 0  ) {
        anc_O3.data[0]=NULL;
        anc_O3.data[1]=NULL;
        anc_O3.data[2]=NULL;
        anc_O3.data[3]=NULL; 
        anc_O3.nblayers=4;
        anc_O3.timeres=6; 
        strcpy (anc_O3.source, "N/A");
        strcpy(anc_O3.filename[0],param->ncep_file_name[0]);
        strcpy(anc_O3.filename[1],param->ncep_file_name[1]);
        strcpy(anc_O3.filename[2],param->ncep_file_name[2]);
        strcpy(anc_O3.filename[3],param->ncep_file_name[3]);

        if (read_grib_anc(&anc_O3,TYPE_OZONE_DATA))
            EXIT_ERROR("Can't read NCEP Ozone data","main");

        anc_WV.data[0]=NULL;
        anc_WV.data[1]=NULL;
        anc_WV.data[2]=NULL;
        anc_WV.data[3]=NULL;
        anc_WV.nblayers=4;
        anc_WV.timeres=6; 
        strcpy (anc_WV.source, "N/A");
        strcpy(anc_WV.filename[0],param->ncep_file_name[0]);
        strcpy(anc_WV.filename[1],param->ncep_file_name[1]);
        strcpy(anc_WV.filename[2],param->ncep_file_name[2]);
        strcpy(anc_WV.filename[3],param->ncep_file_name[3]);
        if (read_grib_anc(&anc_WV,TYPE_WV_DATA))
            EXIT_ERROR("Can't read NCEP WV data","main");

        anc_SP.data[0]=NULL;
        anc_SP.data[1]=NULL;
        anc_SP.data[2]=NULL;
        anc_SP.data[3]=NULL;
        anc_SP.nblayers=4;
        anc_SP.timeres=6; 
        strcpy (anc_SP.source, "N/A");
        strcpy(anc_SP.filename[0],param->ncep_file_name[0]);
        strcpy(anc_SP.filename[1],param->ncep_file_name[1]);
        strcpy(anc_SP.filename[2],param->ncep_file_name[2]);
        strcpy(anc_SP.filename[3],param->ncep_file_name[3]);
        if (read_grib_anc(&anc_SP,TYPE_SP_DATA))
            EXIT_ERROR("Can't read NCEP SP data","main");

        anc_ATEMP.data[0]=NULL;
        anc_ATEMP.data[1]=NULL;
        anc_ATEMP.data[2]=NULL;
        anc_ATEMP.data[3]=NULL;
        anc_ATEMP.nblayers=4;
        anc_ATEMP.timeres=6; 
        strcpy (anc_ATEMP.source, "N/A");
        strcpy(anc_ATEMP.filename[0],param->ncep_file_name[0]);
        strcpy(anc_ATEMP.filename[1],param->ncep_file_name[1]);
        strcpy(anc_ATEMP.filename[2],param->ncep_file_name[2]);
        strcpy(anc_ATEMP.filename[3],param->ncep_file_name[3]);
        if (read_grib_anc(&anc_ATEMP,TYPE_ATEMP_DATA))
            EXIT_ERROR("Can't read NCEP SP data","main");
    }
    else {
        EXIT_ERROR("No input NCEP or PRWV data specified","main");
    }

    /* Convert the units */
    /* convert Pascals into millibars (divide by 100) */
    for (i=0;i<anc_SP.nblayers;i++)
        for (j=0;j<anc_SP.nbrows*anc_SP.nbcols;j++)
            anc_SP.data[i][j] *= 0.01;

    /* convert original PRWV kg/m2 into g/cm2 (divide by 10) */
    for (i=0;i<anc_WV.nblayers;i++)
        for (j=0;j<anc_WV.nbrows*anc_WV.nbcols;j++)
            anc_WV.data[i][j] *= 0.1;

    /* convert O3 to cm-atm (divide by 1000) */
    if (!no_ozone_file) {
        for (i=0;i<anc_O3.nblayers;i++)
            for (j=0;j<anc_O3.nbrows*anc_O3.nbcols;j++)
                anc_O3.data[i][j] *= 0.001;
    }

    /* read DEM file */
    dem_name= (char*)(param->dem_flag ? param->dem_file : DEMFILE );

    /* Open file for SD access */
    sds_file_id = SDstart((char *)dem_name, DFACC_RDONLY);
    if (sds_file_id == HDF_ERROR) {
        EXIT_ERROR("opening dem_file", "OpenDem");
    }
    sds_index=0;
    sds_id= SDselect(sds_file_id,sds_index);
    status=  SDgetinfo(sds_id, sds_name, &rank, dim_sizes, &data_type,&n_attrs);
    start[0]=0;
    start[1]=0;
    edges[0]=3600;   /* number of lines in the DEM data */
    edges[1]=7200;   /* number of samples in the DEM data */
    stride[0]=1;
    stride[1]=1;
    dem_array=(short *)malloc(DEM_NBLAT*DEM_NBLON*sizeof(short));
    status=SDreaddata(sds_id,start, stride, edges,dem_array);
    if (status != 0 ) {
        printf("Fatal error DEM file not read\n");
        exit(EXIT_FAILURE);
    }
    dem_available=1;

    /* Print the ancillary metadata info */
    if ( debug_flag ) {
        print_anc_data(&anc_SP,"SP_DATA");
        print_anc_data(&anc_WV,"WV_DATA");
        print_anc_data(&anc_ATEMP,"ATEMP_DATA");
        if (!no_ozone_file) 
            print_anc_data(&anc_O3,"OZONE_DATA");
    }
    print_anc_data(&anc_O3,"OZONE_DATA");

    /****
    Get center lat lon and deviation from true north
    ****/
    img.l=input->size.l/2.; 
    img.s=input->size.s/2.; 
    img.is_fill=false;
    if (!from_space(space, &img, &geo))
        EXIT_ERROR("mapping from space (0)", "main");
    center_lat=geo.lat * DEG;
    center_lon=geo.lon * DEG;
    printf ("(y0,x0)=(%d,%d)  (lat0,lon0)=(%f,%f)\n",
        (int)img.l,(int)img.s,(float)(geo.lat * DEG),(float)(geo.lon * DEG));

    delta_y=img.l;
    delta_x=img.s;

    img.l=input->size.l/2.-100.; 
    img.s=input->size.s/2.; 
    img.is_fill=false;
    if (!from_space(space, &img, &geo))
        EXIT_ERROR("mapping from space (0)", "main");

    geo.lon=center_lon*RAD;
    geo.is_fill=false;
    if (!to_space(space, &geo, &img))
        EXIT_ERROR("mapping to space (0)", "main");

    delta_y = delta_y - img.l;
    delta_x = img.s - delta_x;
    adjust_north=(float)(atan(delta_x/delta_y)*DEG);

    printf("True North adjustment = %f\n",adjust_north);


#ifdef SAVE_6S_RESULTS
    if (read_6S_results_from_file(SIXS_RESULTS_FILENAME,&sixs_tables)) {
#endif
    /****
    Run 6S and compute atmcor params
    ****/
/*    printf ("DEBUG: Interpolating WV at scene center ...\n"); */
    interpol_spatial_anc(&anc_WV,center_lat,center_lon,tmpflt_arr);
    tmpint=(int)(scene_gmt/anc_WV.timeres);
    if (tmpint>=(anc_WV.nblayers-1))
        tmpint=anc_WV.nblayers-2;
    coef=(double)(scene_gmt-anc_WV.time[tmpint])/anc_WV.timeres;
    sixs_tables.uwv = tmpflt_arr[tmpint]
                    + coef*(tmpflt_arr[tmpint+1] - tmpflt_arr[tmpint]);

    if (!no_ozone_file) {
/*        printf ("DEBUG: Interpolating ozone at scene center ...\n"); */
        interpol_spatial_anc(&anc_O3,center_lat,center_lon,tmpflt_arr);
        tmpint=(int)(scene_gmt/anc_O3.timeres);
        if ( anc_O3.nblayers> 1 ){
            if (tmpint>=(anc_O3.nblayers-1))tmpint=anc_O3.nblayers-2;
            coef=(double)(scene_gmt-anc_O3.time[tmpint])/anc_O3.timeres;
            sixs_tables.uoz = tmpflt_arr[tmpint]
                            + coef*(tmpflt_arr[tmpint+1] - tmpflt_arr[tmpint]);
        }
        else {
            sixs_tables.uoz=tmpflt_arr[tmpint];
        }
    }
    else {
        jday=(short)input->meta.acq_date.doy;
        sixs_tables.uoz=calcuoz(jday,(float)center_lat);
    }

    sixs_tables.target_alt=0.; /* target altitude in km (sea level) */
    sixs_tables.sza=input->meta.sun_zen*DEG;
    sixs_tables.phi=corrected_sun_az;
    sixs_tables.vza=0.;
    sixs_tables.month=9;
    sixs_tables.day=15;
    sixs_tables.srefl=0.14;
/*
    printf ("Center : Lat = %7.2f  Lon = %7.2f \n",center_lat,center_lon);
    printf ("          O3 = %7.2f   SP = %7.2f   WV = %7.2f\n",\
                sixs_tables.uoz,tmpflt,sixs_tables.uwv);
*/
    switch (input->meta.inst) {
        case INST_TM:
            sixs_tables.Inst=SIXS_INST_TM; break;
        case INST_ETM:
            sixs_tables.Inst=SIXS_INST_ETM; break;
        default:
            EXIT_ERROR("Unknown Instrument", "main");
    }
    create_6S_tables(&sixs_tables, &input->meta);
#ifdef SAVE_6S_RESULTS
    write_6S_results_to_file(SIXS_RESULTS_FILENAME,&sixs_tables);
    }
#endif

/***
    interpolate ancillary data for AR grid cells
***/
    img.is_fill=false;
    sum_spres_anc=0.;
    sum_spres_dem=0.;
    nb_spres_anc=0;
    nb_spres_dem=0;
    int index = 0;
    for (il_ar = 0, img.l = 0.5*lut->ar_region_size.l;
         il_ar < lut->ar_size.l; il_ar++, img.l += lut->ar_region_size.l) {
        for (is_ar = 0, img.s = 0.5*lut->ar_region_size.s;
             is_ar < lut->ar_size.s;
             is_ar++, index++, img.s += lut->ar_region_size.s) {
            if (!from_space(space, &img, &geo))
                EXIT_ERROR("mapping from space (1)", "main");

            ar_gridcell.lat[index] = geo.lat*DEG;
            ar_gridcell.lon[index] = geo.lon*DEG;
            ar_gridcell.sun_zen[index]= input->meta.sun_zen*DEG;
            ar_gridcell.view_zen[index] = 3.5;
            ar_gridcell.rel_az[index] = corrected_sun_az;

            interpol_spatial_anc(&anc_WV, ar_gridcell.lat[index],
                                 ar_gridcell.lon[index], tmpflt_arr);
            tmpint=(int)(scene_gmt/anc_WV.timeres);
            if (tmpint>=(anc_WV.nblayers-1))
                tmpint=anc_WV.nblayers-2;
            coef=(double)(scene_gmt-anc_WV.time[tmpint])/anc_WV.timeres;
            ar_gridcell.wv[index] = tmpflt_arr[tmpint]
                           + coef*(tmpflt_arr[tmpint+1] - tmpflt_arr[tmpint]);

            if (!no_ozone_file) {
                interpol_spatial_anc(&anc_O3, ar_gridcell.lat[index],
                                     ar_gridcell.lon[index], tmpflt_arr);
                tmpint=(int)(scene_gmt/anc_O3.timeres);
                if ( anc_O3.nblayers> 1 ){
                    if (tmpint>=(anc_O3.nblayers-1))
                        tmpint=anc_O3.nblayers-2;
                    coef=(double)(scene_gmt-anc_O3.time[tmpint])/anc_O3.timeres;
                    ar_gridcell.ozone[index] = tmpflt_arr[tmpint] +
                        coef*(tmpflt_arr[tmpint+1] - tmpflt_arr[tmpint]);
                   }
                else {
                    ar_gridcell.ozone[index] = tmpflt_arr[tmpint];
                   }
            }
            else {
                jday=(short)input->meta.acq_date.doy;
                ar_gridcell.ozone[index] =
                                 calcuoz(jday, (float)ar_gridcell.lat[index]);
            }

            interpol_spatial_anc(&anc_SP, ar_gridcell.lat[index],
                                 ar_gridcell.lon[index], tmpflt_arr);
            tmpint=(int)(scene_gmt/anc_SP.timeres);
            if (tmpint>=(anc_SP.nblayers-1))
                tmpint=anc_SP.nblayers-2;
            coef=(double)(scene_gmt-anc_SP.time[tmpint])/anc_SP.timeres;
            ar_gridcell.spres[index] = tmpflt_arr[tmpint] +
                coef*(tmpflt_arr[tmpint+1] - tmpflt_arr[tmpint]);
            if (ar_gridcell.spres[index] > 0) {
                sum_spres_anc += ar_gridcell.spres[index];
                nb_spres_anc++;
            }
            if (dem_available) {
                ar_gridcell.spres_dem[index] =
                    get_dem_spres(dem_array, ar_gridcell.lat[index],
                                  ar_gridcell.lon[index]);
                if (ar_gridcell.spres_dem[index] > 0) {
                    sum_spres_dem += ar_gridcell.spres_dem[index];
                    nb_spres_dem++;
                }
            }
        }  /* for is_ar */
    }  /* for il_ar */

    if (dem_available) {
        for (il_ar = 0, index = 0; il_ar < lut->ar_size.l;il_ar++) 
            for (is_ar=0;is_ar < lut->ar_size.s; is_ar++, index++) 
                if (ar_gridcell.spres[index] > 0 &&
                    ar_gridcell.spres_dem[index] > 0)
                    ar_gridcell.spres[index] = ar_gridcell.spres_dem[index]
                                             * ar_gridcell.spres[index]/1013.;
    }

    /* Compute atmospheric coefs for the whole scene with aot550=0.01 for use
       in internal cloud screening : NAZMI */
    nbpts=lut->ar_size.l*lut->ar_size.s;

    /***
    Allocate memory for atmos_coeff
    ***/
    if (allocate_mem_atmos_coeff(nbpts,&atmos_coef))
        EXIT_ERROR("Allocating memory for atmos_coef", "main");

    printf("Compute Atmos Params with aot550 = 0.01\n"); fflush(stdout);
    update_atmos_coefs(&atmos_coef, &ar_gridcell, &sixs_tables, line_ar, lut,
        input->nband, 1);

    /* Read input first time and compute clear pixels stats for internal cloud
       screening */
    /* allocate memory for cld_diags structure and clear sum and nb of obs */
    if (allocate_cld_diags(&cld_diags,CLDDIAGS_CELLHEIGHT_5KM,
        CLDDIAGS_CELLWIDTH_5KM, input->size.l, input->size.s)) {
        EXIT_ERROR("couldn't allocate memory from cld_diags","main");
    }

    /* Allocate memory for the interpolated atmospheric coefficients. */
    if (allocate_mem_atmos_coeff(1, &atmos_coef_storage))
        EXIT_ERROR("Allocating memory for atmos_coef storage", "main");

    /* Screen the clouds */
    for (il = 0; il < input->size.l; il++) {
        if (!(il%100)) 
        {
            printf("First pass cloud screening for line %d\r",il);
            fflush(stdout);
        }

        /* Read each input band */
        for (ib = 0; ib < input->nband; ib++) {
            if (!GetInputLine(input, ib, il, line_in[0][ib]))
                EXIT_ERROR("reading input data for a line (b)", "main");
        }
        if (!GetInputQALine(input, il, qa_line[0], qa2_line[0]))
            EXIT_ERROR("reading input data for qa_line (1)", "main");
        if (param->thermal_band) {
            if (!GetInputLine(input_b6, 0, il, b6_line[0]))
                EXIT_ERROR("reading input data for b6_line (1)", "main");
        }

        tmpint = (int)(scene_gmt / anc_ATEMP.timeres);
        if (tmpint >= anc_ATEMP.nblayers - 1)
            tmpint = anc_ATEMP.nblayers - 2;
        coef = (double)(scene_gmt - anc_ATEMP.time[tmpint]) / anc_ATEMP.timeres;

        img.is_fill = false;
        img.l = il;
#ifdef _OPENMP
        #pragma omp parallel for private (is, geo, flat, flon, tmpflt_arr) \
        firstprivate (img, atemp_line) num_threads(get_num_threads())
#endif
        for (is = 0; is < input->size.s; is++) {
            /* Get the geolocation info for this pixel */
            img.s = is;
            if (!from_space (space, &img, &geo))
                EXIT_ERROR("mapping from space (2)", "main");
            flat = geo.lat * DEG;
            flon = geo.lon * DEG;

            /* Interpolate the anciliary data for this lat/long, then pull the
               information for the scene center time and adjust */
            interpol_spatial_anc (&anc_ATEMP, flat, flon, tmpflt_arr);
            atemp_line[is] = tmpflt_arr[tmpint]
                           + coef*(tmpflt_arr[tmpint+1] - tmpflt_arr[tmpint]);
        }

        /* Run Cld Screening Pass1 and compute stats. This cloud detection
           function contains statistics gathering that needs to be in a
           critical section for multi-threading. */
        if (param->thermal_band)
            if (!cloud_detection_pass1(lut, input->size.s, il, line_in[0],
                               qa_line[0], qa2_line[0], b6_line[0], atemp_line,
                               &atmos_coef, &atmos_coef_storage, &cld_diags))
                EXIT_ERROR("running cloud detection pass 1", "main");
    } /* end for il */
    printf ("\n");

    if (param->thermal_band) {
        int print_step = cld_diags.nbrows/10;
        for (il = 0; il < cld_diags.nbrows; il++) {
            if (!(il%print_step))
            {
                printf("Second pass cloud screening for line %d\r",il);
                fflush(stdout);
            }

            tmpint = (int)(scene_gmt / anc_ATEMP.timeres);
            if (tmpint >= anc_ATEMP.nblayers - 1)
                tmpint = anc_ATEMP.nblayers - 2;
            coef = (double)(scene_gmt - anc_ATEMP.time[tmpint])
                 / anc_ATEMP.timeres;

            img.is_fill = false;
            img.l = (il + 0.5)*cld_diags.cellheight;
            if (img.l >= input->size.l)
                img.l = input->size.l-1;
            for (is = 0; is < cld_diags.nbcols; is++) {
                img.s = (is + 0.5)*cld_diags.cellwidth;
                if (img.s >= input->size.s)
                    img.s = input->size.s-1;
                if (!from_space (space, &img, &geo))
                    EXIT_ERROR("mapping from space (3)", "main");
                flat = geo.lat*DEG;
                flon = geo.lon*DEG;

                interpol_spatial_anc(&anc_ATEMP,flat,flon,tmpflt_arr);
                cld_diags.airtemp_2m[il][is] = tmpflt_arr[tmpint]
                           + coef*(tmpflt_arr[tmpint+1] - tmpflt_arr[tmpint]);

                if (cld_diags.nb_t6_clear[il][is] > 0) {
                    sum_value = cld_diags.avg_t6_clear[il][is];
                    sumsq_value = cld_diags.std_t6_clear[il][is];
                    cld_diags.avg_t6_clear[il][is] = sum_value
                                              / cld_diags.nb_t6_clear[il][is];
                    if (cld_diags.nb_t6_clear[il][is] > 1) {
                        cld_diags.std_t6_clear[il][is] =
                            (sumsq_value -
                             sum_value*sum_value/cld_diags.nb_t6_clear[il][is])
                            /(cld_diags.nb_t6_clear[il][is] - 1);
                        cld_diags.std_t6_clear[il][is] =
                            sqrt(fabs(cld_diags.std_t6_clear[il][is]));
                    }
                    else 
                        cld_diags.std_t6_clear[il][is] = 0;

                    sum_value = cld_diags.avg_b7_clear[il][is];
                    sumsq_value = cld_diags.std_b7_clear[il][is];
                    cld_diags.avg_b7_clear[il][is] = sum_value
                                              / cld_diags.nb_t6_clear[il][is];
                    if (cld_diags.nb_t6_clear[il][is] > 1) {
                        cld_diags.std_b7_clear[il][is] =
                            (sumsq_value -
                             sum_value*sum_value/cld_diags.nb_t6_clear[il][is])
                            /(cld_diags.nb_t6_clear[il][is] - 1);
                        cld_diags.std_b7_clear[il][is] =
                            sqrt(fabs(cld_diags.std_b7_clear[il][is]));
                    }
                    else
                        cld_diags.std_b7_clear[il][is] = 0;
                }
                else {
                    cld_diags.avg_t6_clear[il][is] = -9999.;
                    cld_diags.avg_b7_clear[il][is] = -9999.;
                    cld_diags.std_t6_clear[il][is] = -9999.;
                    cld_diags.std_b7_clear[il][is] = -9999.;
                }
            }  /* end for is */
        }  /* end for il */

        fill_cld_diags(&cld_diags);
#ifdef DEBUG_CLD
        for (il=0;il<cld_diags.nbrows;il++) 
            for (is=0;is<cld_diags.nbcols;is++) 
                if (fd_cld_diags != NULL)
                    fprintf(fd_cld_diags, "%d %d %d %f %f %f %f %f\n",
                            il, is, cld_diags.nb_t6_clear[il][is],
                            cld_diags.airtemp_2m[il][is],
                            cld_diags.avg_t6_clear[il][is],
                            cld_diags.std_t6_clear[il][is],
                            cld_diags.avg_b7_clear[il][is],
                            cld_diags.std_b7_clear[il][is]);
        fclose(fd_cld_diags);
#endif
    }  /* end if thermal band */
    printf ("\n");

/***
    Create dark target temporary file
***/
    strcpy(tmpfilename, "temporary_dark_target_XXXXXX");
    if ((tmpid = mkstemp (tmpfilename)) < 1)
      EXIT_ERROR("creating filename for dark target temporary file", "main");
    close(tmpid);
    if ((fdtmp=fopen(tmpfilename,"w"))==NULL)
      EXIT_ERROR("creating dark target temporary file", "main");

    /* Read input second time and create cloud and cloud shadow masks */
    ptr_rot_cld[0]=rot_cld[0];
    ptr_rot_cld[1]=rot_cld[1];
    ptr_rot_cld[2]=rot_cld[2];

    /* Allocate storage for fill mask. */
    char *fill_mask;
    if ((fill_mask = malloc(lut->ar_region_size.l*input->size.s)) == NULL)
        EXIT_ERROR("allocating space for fill mask", "main");

    for (il_start = 0, il_ar = 0, index = 0; il_start < input->size.l; 
         il_start += lut->ar_region_size.l, il_ar++, index += lut->ar_size.s) {
        ar_gridcell.line_lat = &(ar_gridcell.lat[index]);
        ar_gridcell.line_lon = &(ar_gridcell.lon[index]);
        ar_gridcell.line_sun_zen = &(ar_gridcell.sun_zen[index]);
        ar_gridcell.line_view_zen = &(ar_gridcell.view_zen[index]);
        ar_gridcell.line_rel_az = &(ar_gridcell.rel_az[index]);
        ar_gridcell.line_wv = &(ar_gridcell.wv[index]);
        ar_gridcell.line_spres = &(ar_gridcell.spres[index]);
        ar_gridcell.line_ozone = &(ar_gridcell.ozone[index]);
        ar_gridcell.line_spres_dem = &(ar_gridcell.spres[index]);
        
        il_end = il_start + lut->ar_region_size.l - 1;
        if (il_end >= input->size.l)
            il_end = input->size.l - 1;

        /* Read each input band for each line in region */
        for (il = il_start; il < (il_end + 1); il++) {
            il_region = il - il_start;
            for (ib = 0; ib < input->nband; ib++) {
                if (!GetInputLine(input, ib, il, line_in[il_region][ib]))
                    EXIT_ERROR("reading input data for a line (a)", "main");
            }

            if (!GetInputQALine(input, il, qa_line[il_region],
                                qa2_line[il_region]))
                EXIT_ERROR("reading input data for qa_line (2)", "main");
            if (param->thermal_band) {
                if (!GetInputLine(input_b6, 0, il, b6_line[il_region]))
                    EXIT_ERROR("reading input data for b6_line (2)", "main");

                /* Run Cld Screening Pass2 */
                if (!cloud_detection_pass2(lut, input->size.s, il,
                                           line_in[il_region],
                                           qa_line[il_region],
                                           qa2_line[il_region],
                                           b6_line[il_region], &atmos_coef,
                                           &atmos_coef_storage, &cld_diags,
                                           ptr_rot_cld[1][il_region]))
                    EXIT_ERROR("running cloud detection pass 2", "main");
            }
            else {
                if (!cloud_detection_pass2(lut, input->size.s, il,
                                           line_in[il_region],
                                           qa_line[il_region], 
                                           qa2_line[il_region], 
                                           NULL,
                                           &atmos_coef, &atmos_coef_storage,
                                           &cld_diags,
                                           ptr_rot_cld[1][il_region]))
                    EXIT_ERROR("running cloud detection pass 2", "main");
            }
        }  /* end for il */

        if (param->thermal_band) {
            /* Cloud Mask Dilation : 5 pixels */
            if (!dilate_cloud_mask(lut, input->size.s, ptr_rot_cld, 5))
                EXIT_ERROR("running cloud mask dilation", "main");

            /* Cloud shadow */
            cast_cloud_shadow(lut, input->size.s, il_start, line_in, b6_line,
                &cld_diags,ptr_rot_cld,&ar_gridcell, space_def.pixel_size[0],
                adjust_north);

            /* Dilate Cloud shadow */
            dilate_shadow_mask(lut, input->size.s, fill_mask, ptr_rot_cld, 5);
        }

        /***
        Save cloud and cloud shadow in temporary file
        ***/
        if (il_ar > 0)
            if (fwrite(ptr_rot_cld[0][0],lut->ar_region_size.l*input->size.s,1,
                fdtmp) != 1)
                EXIT_ERROR("writing dark target to temporary file", "main");

        ptr_tmp_cld=ptr_rot_cld[0];
        ptr_rot_cld[0]=ptr_rot_cld[1];
        ptr_rot_cld[1]=ptr_rot_cld[2];
        ptr_rot_cld[2]=ptr_tmp_cld;

        for (i=0;i<lut->ar_region_size.l;i++)
            memset(&ptr_rot_cld[2][i][0],0,input->size.s);
    }  /* end for il_start */

    /** Last Block **/
    dilate_shadow_mask(lut, input->size.s, fill_mask, ptr_rot_cld, 5);
    if (fwrite(ptr_rot_cld[0][0],lut->ar_region_size.l*input->size.s,1,fdtmp)
        != 1)
        EXIT_ERROR("writing dark target to temporary file", "main");

    fclose(fdtmp);

    /* Free temporary storage. */
    free_mem_atmos_coeff(&atmos_coef_storage);
    free(fill_mask);

    /* Done with the cloud diagnostics */
    free_cld_diags (&cld_diags);

    /***
    Open temporary file for read and write
    ***/
    if ((fdtmp=fopen(tmpfilename,"r+"))==NULL)
        EXIT_ERROR("opening dark target temporary file (r+)", "main");

    /* Allocate memory for the Ar() function. */
    if (allocate_mem_atmos_coeff(ar_gridcell.nbrows*ar_gridcell.nbcols,
                                 &atmos_coef_storage))
        EXIT_ERROR("Allocating memory for atmos_coef storage", "main");
    collect_bands_t *cbands;
    if ((cbands = malloc(lut->ar_region_size.s*lut->ar_region_size.l*
                         sizeof(*cbands))) == NULL)
        EXIT_ERROR("Allocating memory for collect bands", "main");

    /* Read input second time and compute the aerosol for each region */
    for (il_start = 0, il_ar = 0, index = 0; il_start < input->size.l; 
         il_start += lut->ar_region_size.l, il_ar++, index += lut->ar_size.s) {
        ar_gridcell.line_lat=&(ar_gridcell.lat[index]);
        ar_gridcell.line_lon=&(ar_gridcell.lon[index]);
        ar_gridcell.line_sun_zen=&(ar_gridcell.sun_zen[index]);
        ar_gridcell.line_view_zen=&(ar_gridcell.view_zen[index]);
        ar_gridcell.line_rel_az=&(ar_gridcell.rel_az[index]);
        ar_gridcell.line_wv=&(ar_gridcell.wv[index]);
        ar_gridcell.line_spres=&(ar_gridcell.spres[index]);
        ar_gridcell.line_ozone=&(ar_gridcell.ozone[index]);
        ar_gridcell.line_spres_dem=&(ar_gridcell.spres[index]);
    
        il_end = il_start + lut->ar_region_size.l - 1;
        if (il_end >= input->size.l)
            il_end = input->size.l - 1;
     
        if (fseek(fdtmp,(long)(il_ar*lut->ar_region_size.l*input->size.s),
            SEEK_SET))
            EXIT_ERROR("seeking in temporary file (r)", "main");
        if (fread(ddv_line[0],lut->ar_region_size.l*input->size.s,1,fdtmp)!=1)
            EXIT_ERROR("reading dark target to temporary file", "main");

        /* Read each input band for each line in region */
        for (il = il_start, il_region = 0; il < (il_end + 1);
             il++, il_region++) {
            for (ib = 0; ib < input->nband; ib++) {
                if (!GetInputLine(input, ib, il, line_in[il_region][ib]))
                    EXIT_ERROR("reading input data for a line (a)", "main");
            }
        }  /* end for il */

        /* Compute the aerosol for the regions */
#ifdef DEBUG_AR
        diags_il_ar=il_ar;
#endif
        if (!Ar(il_ar,lut, &input->size, line_in, ddv_line,
                &atmos_coef_storage, cbands, line_ar[il_ar], &ar_stats,
                &ar_gridcell, &sixs_tables))
            EXIT_ERROR("computing aerosol", "main");

        /***
        Save dark target map in temporary file
        ***/
        if (fseek(fdtmp,il_ar*lut->ar_region_size.l*input->size.s,SEEK_SET))
            EXIT_ERROR("seeking in temporary file (w)", "main");
        if (fwrite(ddv_line[0],lut->ar_region_size.l*input->size.s,1,fdtmp)!=1)
            EXIT_ERROR("writing dark target to temporary file", "main");
    }  /* end for il_start */

    printf("\n");
    fclose(fdtmp);
#ifdef DEBUG_AR
    fclose(fd_ar_diags);
#endif

    /* Free temporary storage. */
    free_mem_atmos_coeff(&atmos_coef_storage);
    free(cbands);

    /***
    Fill Gaps in the coarse resolution aerosol product for bands 1(0), 2(1)
    and 3(2)
    ***/
    Fill_Ar_Gaps(lut, line_ar, 0);

    /* Compute atmospheric coeffs for the whole scene using retrieved aot */
    nbpts=lut->ar_size.l*lut->ar_size.s;

    printf("Compute Atmos Params\n"); fflush(stdout);
#ifdef NO_AEROSOL_CORRECTION
    update_atmos_coefs(&atmos_coef, &ar_gridcell, &sixs_tables, line_ar, lut,
        input->nband, 1);
#else
    update_atmos_coefs(&atmos_coef, &ar_gridcell, &sixs_tables, line_ar, lut,
        input->nband, 0); /*Eric COMMENTED TO PERFORM NO CORRECTION*/
#endif

    /* Re-read input and compute surface reflectance */
    /***
    Open temporary file for read
    ***/
    if ((fdtmp=fopen(tmpfilename,"r"))==NULL)
        EXIT_ERROR("opening dark target temporary file", "main");

    /* Allocate memory for the interpolated atmospheric coefficients. */
    if (allocate_mem_atmos_coeff(1, &atmos_coef_storage))
        EXIT_ERROR("Allocating memory for atmos_coef storage", "main");

    for (il = 0; il < input->size.l; il++) {
        if (!(il%100)) 
        {
            printf("Processing surface reflectance for line %d\r",il);
            fflush(stdout);
        }

        /* Re-read each input band */
        for (ib = 0; ib < input->nband; ib++) {
            if (!GetInputLine(input, ib, il, line_in[0][ib]))
                EXIT_ERROR("reading input data for a line (b)", "main");
        }
    
        if (!GetInputLine(input_b6, 0, il, b6_line[0]))
            EXIT_ERROR("reading input data for b6_line (1)", "main");

        if (!GetInputQALine(input, il, qa_line[0], NULL))
            EXIT_ERROR("reading input data for qa_line (3)", "main");

        /* Compute the surface reflectance */
        if (!Sr(lut, input->size.s, il, &atmos_coef, &atmos_coef_storage,
                line_in[0], qa_line[0], line_out, &sr_stats))
            EXIT_ERROR("computing surface reflectance for a line", "main");

        /***
        Read line from dark target temporary file
        ***/
        if (fread(ddv_line[0],input->size.s,1,fdtmp)!=1)
            EXIT_ERROR("reading line from dark target temporary file", "main");

        loc.l=il;

        for (is=0;is<input->size.s;is++) {
            loc.s=is;

            /* Initialize QA band to off */
            line_out[lut->nband+CLOUD][is] = QA_OFF;

            /* Determine if this is a fill pixel -- mark output band as fill,
               else process SR QA and atmos opacity for each pixel */
            if (!level1_qa_is_fill(qa_line[0][is])) {
                /* AOT / opacity */
                ArInterp(lut, &loc, line_ar, &inter_aot); 
                line_out[lut->nband+ATMOS_OPACITY][is] = inter_aot;

                /* QA is written out in the cloud band as a bit-packed product
                   (16-bit). We will use QA values as-is and no further
                   post-processing QA step will be implemented. We want the QA
                   to reflect the cloud, etc. status that was used in the
                   aerosol and surface reflectance computations. We are not
                   interested in post-processing of the QA information, as
                   there are better QA products available. */
                if (ddv_line[0][is] & AR_DDV)
                    line_out[lut->nband+CLOUD][is] |= (1 << DDV_BIT);

                if (ddv_line[0][is] & AR_ADJ_CLOUD)
                    line_out[lut->nband+CLOUD][is] |= (1 << ADJ_CLOUD_BIT);

                if (!(ddv_line[0][is] & AR_WATER))  /* if water, turn on */
                    line_out[lut->nband+CLOUD][is] |= (1 << LAND_WATER_BIT);

                if (ddv_line[0][is] & AR_CLOUD)
                    line_out[lut->nband+CLOUD][is] |= (1 << CLOUD_BIT);

                if (ddv_line[0][is] & AR_CLOUD_SHADOW)
                    line_out[lut->nband+CLOUD][is] |=
                        (1 << CLOUD_SHADOW_BIT);

                if (ddv_line[0][is] & AR_SNOW)
                    line_out[lut->nband+CLOUD][is] |= (1 << SNOW_BIT);
            }
            else {
                line_out[lut->nband+ATMOS_OPACITY][is]=lut->output_fill_opacity;
            }
        } /* for is */

        /* Write each output band */
        for (ib = 0; ib < output->nband_out; ib++) {
            if (!PutOutputLine(output, ib, il, line_out[ib]))
                EXIT_ERROR("writing output data for a line", "main");
        }
    }  /* for il */
    printf("\n");
    fclose(fdtmp);
    unlink(tmpfilename); 
    
    /* Print the statistics, skip bands that don't exist */
    printf(" total pixels %ld\n", ((long)input->size.l * (long)input->size.s));
    printf(" aerosol coarse  nfill %ld  min  %d  max  %d\n", ar_stats.nfill,
        ar_stats.ar_min, ar_stats.ar_max);

    for (ib = 0; ib < lut->nband; ib++) {
        if (output->metadata.band[ib].name != NULL)
            printf(" sr %s  nfill %ld  nsatu %ld  nout_range %ld  min  %d  "
                "max  %d\n", output->metadata.band[ib].name, 
                sr_stats.nfill[ib], sr_stats.nsatu[ib], sr_stats.nout_range[ib],
                sr_stats.sr_min[ib], sr_stats.sr_max[ib]);
    }

    /* Close input files */
    if (!CloseInput(input)) EXIT_ERROR("closing input file", "main");
    if (!CloseOutput(output)) EXIT_ERROR("closing input file", "main");

    /* Write the ENVI header for reflectance files */
    for (ib = 0; ib < output->nband_out; ib++) {
        /* Create the ENVI header file this band */
        if (create_envi_struct (&output->metadata.band[ib],
            &xml_metadata.global, &envi_hdr) != SUCCESS)
            EXIT_ERROR("Creating the ENVI header structure for this file.",
                "main");

        /* Write the ENVI header */
        strcpy (envi_file, output->metadata.band[ib].file_name);
        cptr = strchr (envi_file, '.');
        strcpy (cptr, ".hdr");
        if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
            EXIT_ERROR("Writing the ENVI header file.", "main");
    }

    /* Append the reflective and thermal bands to the XML file */
    if (append_metadata (output->nband_out, output->metadata.band,
        param->input_xml_file_name) != SUCCESS)
        EXIT_ERROR("appending surfance reflectance and QA bands", "main");

    /* Free the metadata structure */
    free_metadata (&xml_metadata);

    /* Free memory */
    free_mem_atmos_coeff(&atmos_coef);
    free_mem_atmos_coeff(&atmos_coef_storage);

    FreeInput(input);
    FreeInput(input_b6);
    FreeLut(lut);
    if (!FreeOutput(output)) 
        EXIT_ERROR("freeing output file stucture", "main");

    free(space);
    free(line_out[0]);
    free(line_ar[0][0]);
    free(line_ar[0]);
    free(line_ar);
    free(line_in[0][0]);
    free(line_in[0]);
    free(line_in);
    free(qa_line[0]);
    free(qa_line);
    if (param->thermal_band) {
        free(b6_line[0]);
        free(b6_line);
    }
    free(ddv_line[0]);
    free(ddv_line);
    free(rot_cld[0][0]);
    free(rot_cld[0]);
    free(ar_gridcell.lat);
    free(ar_gridcell.lon);
    free(ar_gridcell.sun_zen);
    free(ar_gridcell.view_zen);
    free(ar_gridcell.rel_az);
    free(ar_gridcell.wv);
    free(ar_gridcell.spres);
    free(ar_gridcell.ozone);

    for (ifree=0; ifree<(param->num_ncep_files>0?4:1); ifree++) {
        free(anc_O3.data[ifree]);
        free(anc_WV.data[ifree]);
        free(anc_SP.data[ifree]);
    }
    if (dem_available)
        free(dem_array);
    FreeParam(param);

    /* All done */
    printf ("lndsr complete.\n");
    return (EXIT_SUCCESS);
}

      
int allocate_mem_atmos_coeff(int nbpts,atmos_t *atmos_coef)
{
    int ib;
    if ((atmos_coef->computed=(int *)malloc(nbpts*sizeof(int)))==NULL)
        return -1;
    for (ib=0;ib<7;ib++) {
        if ((atmos_coef->tgOG[ib]=(float *)malloc(nbpts*sizeof(float)))==NULL)
            return -1;
        if ((atmos_coef->tgH2O[ib]=(float *)malloc(nbpts*sizeof(float)))==NULL)
            return -1;
        if ((atmos_coef->td_ra[ib]=(float *)malloc(nbpts*sizeof(float)))==NULL)
            return -1;
        if ((atmos_coef->tu_ra[ib]=(float *)malloc(nbpts*sizeof(float)))==NULL)
            return -1;
        if ((atmos_coef->rho_mol[ib]=(float *)malloc(nbpts*sizeof(float)))==NULL)
            return -1;
        if ((atmos_coef->rho_ra[ib]=(float *)malloc(nbpts*sizeof(float)))==NULL)
            return -1;
        if ((atmos_coef->td_da[ib]=(float *)malloc(nbpts*sizeof(float)))==NULL)
            return -1;
        if ((atmos_coef->tu_da[ib]=(float *)malloc(nbpts*sizeof(float)))==NULL)
            return -1;
        if ((atmos_coef->S_ra[ib]=(float *)malloc(nbpts*sizeof(float)))==NULL)
            return -1;
        if ((atmos_coef->td_r[ib]=(float *)malloc(nbpts*sizeof(float)))==NULL)
            return -1;
        if ((atmos_coef->tu_r[ib]=(float *)malloc(nbpts*sizeof(float)))==NULL)
            return -1;
        if ((atmos_coef->S_r[ib]=(float *)malloc(nbpts*sizeof(float)))==NULL)
            return -1;
        if ((atmos_coef->rho_r[ib]=(float *)malloc(nbpts*sizeof(float)))==NULL)
            return -1;
    }
    return 0;
}

int free_mem_atmos_coeff(atmos_t *atmos_coef)
{
    int ib;
    free(atmos_coef->computed);
    for(ib=0;ib<7;ib++) {
        free(atmos_coef->tgOG[ib]);
        free(atmos_coef->tgH2O[ib]);
        free(atmos_coef->td_ra[ib]);
        free(atmos_coef->tu_ra[ib]);
        free(atmos_coef->rho_mol[ib]);
        free(atmos_coef->rho_ra[ib]);
        free(atmos_coef->td_da[ib]);
        free(atmos_coef->tu_da[ib]);
        free(atmos_coef->S_ra[ib]);
        free(atmos_coef->td_r[ib]);
        free(atmos_coef->tu_r[ib]);
        free(atmos_coef->S_r[ib]);
        free(atmos_coef->rho_r[ib]);
    }
    return 0;
}

/******************************************************************************
!C
!Routine: calcuoz

!Description:  Gets ozone concentration for a particular day and 
               latitude, interpolating if necessary.

!Revision History:
 Original version:    Nazmi Z El Saleous and Eric Vermote

!Input Parameters:
        jday       Julian day
        lat        latitude (in degrees)
        
!Output Parameters:
        uoz        Ozone concentration (in cm-atm) 

!Return value:
        returns 0 (success)
               -1 (latitude is beyond coverage of data, and so
                  the returned ozone value is an approximation).
                
                
!References and Credits:
      Data is from:
    LONDON J.,BOJKOV R.D.,OLTMANS S. AND KELLEY J.I.,1976,
    ATLAS OF THE GLOBAL DISTRIBUTION OF TOTAL OZONE .JULY 1957-JUNE 1967
    NCAR TECHNICAL NOTE, NCAR/TN/113+STR,PP276
      
!Developers:
      Nazmi Z El Saleous
      Eric Vermote
      University of Maryland / Dept. of Geography
      nazmi.elsaleous@gsfc.nasa.gov

!Design Notes:  
      
!END
*******************************************************************************/
float calcuoz(short jday,float flat)
{
    float t,u,tmpf;
    int i1,i2,j1,j2,Minf,Msup,Latinf,Latsup;

/*
CREPARTITION ZONALE PAR BANDE DE 10 DEG DE LATITUDE A PARTIR DE 80 SUD
C            ---   OZONE   ---
*/
    float oz[12][17] = {
    {.315,.320,.315,.305,.300,.280,.260,.240,.240,.240,.250,.280,
     .320,.350,.375,.380,.380},
    {.280,.300,.300,.300,.280,.270,.260,.240,.240,.240,.260,.300,
     .340,.380,.400,.420,.420},
    {.280,.280,.280,.280,.280,.260,.250,.240,.250,.250,.270,.300,
     .340,.400,.420,.440,.440},
    {.280,.280,.280,.280,.280,.260,.250,.250,.250,.260,.280,.300,
     .340,.380,.420,.430,.430},
    {.280,.290,.300,.300,.280,.270,.260,.250,.250,.260,.270,.300,
     .320,.360,.380,.400,.400},
    {.280,.300,.300,.305,.300,.280,.260,.250,.250,.260,.260,.280,
     .310,.330,.360,.370,.370},
    {.290,.300,.315,.320,.305,.280,.260,.250,.240,.240,.260,.270,
     .290,.310,.320,.320,.320},
    {.300,.310,.320,.325,.320,.300,.270,.260,.240,.240,.250,.260,
     .280,.290,.300,.300,.290},
    {.300,.320,.325,.335,.320,.300,.280,.260,.240,.240,.240,.260,
     .270,.280,.280,.280,.280},
    {.320,.340,.350,.345,.330,.300,.280,.260,.240,.240,.240,.260,
     .260,.280,.280,.280,.280},
    {.360,.360,.360,.340,.320,.300,.280,.260,.240,.240,.240,.260,
     .280,.300,.310,.310,.300},
    {.340,.350,.340,.320,.310,.280,.260,.250,.240,.240,.240,.260,
     .300,.320,.330,.340,.330}};
/*     
C      /Begin of interpolation/
C      /Find Neighbours/
C      /First loop for time/
*/
    if (fabs(flat)>=80.) {
        tmpf=.270;
        return -1;
    }
    Minf= (int) ((jday-15.)/30.5);
    if (jday < 15) 
        Minf=Minf-1;
    Latinf=(int) (flat*0.1);
    if (flat < 0.) 
        Latinf=Latinf-1;
    t=((jday-15.)-(30.5*Minf))/30.5;
    u=(flat-10.*Latinf)*0.1;
   
    Minf=(Minf+12)%12;
    Msup=(Minf+13)%12;
    Latsup=Latinf+1;
    i1=Minf;
    j1=Latinf+8;
    i2=Msup;
    j2=Latsup+8;
   
    /* Now Calculate Uo3 at the given point Xlat,Jjulien */
    tmpf = oz[i1][j1]
         + t*(oz[i2][j1] - oz[i1][j1])
         + u*(oz[i1][j2] - oz[i1][j1])
         + t*u*(oz[i1][j1] - oz[i2][j1] - oz[i1][j2] + oz[i2][j2]);
    return tmpf;
}

float get_dem_spres(short *dem,float lat,float lon)
{
    int idem,jdem;
    float dem_spres;
        
    idem=(int)((DEM_LATMAX-lat)/DEM_DLAT+0.5);
    if (idem<0)
        idem=0;
    else if (idem >= DEM_NBLAT)
        idem=DEM_NBLAT-1;
    jdem=(int)((lon-DEM_LONMIN)/DEM_DLON+0.5);
    if (jdem<0)
        jdem=0;
    else if (jdem >= DEM_NBLON)
        jdem=DEM_NBLON-1;
    if (dem[idem*DEM_NBLON+jdem]== -9999)
        dem_spres=1013;
    else
        dem_spres=1013.2*exp(-dem[idem*DEM_NBLON+jdem]/8000.);

    return dem_spres;
}

static void update_atmos_coefs
(
    atmos_t *atmos_coef,
    Ar_gridcell_t *ar_gridcell,
    sixs_tables_t *sixs_tables,
    int ***line_ar,
    Lut_t *lut,
    int nband,
    int bkgd_aerosol
)
{
    int irow, icol, ipt;

    for (irow=0, ipt=0; irow<ar_gridcell->nbrows; irow++)
        for (icol=0; icol<ar_gridcell->nbcols; icol++, ipt++)
            update_gridcell_atmos_coefs(ipt, atmos_coef, ar_gridcell,
                                        sixs_tables, line_ar[irow][0][icol],
                                        lut, nband, bkgd_aerosol);
}

void update_gridcell_atmos_coefs
(
    int ipt,
    atmos_t *atmos_coef,
    Ar_gridcell_t *ar_gridcell,
    sixs_tables_t *sixs_tables,
    int line_ar,
    Lut_t *lut,
    int nband,
    int bkgd_aerosol
)
{
    int ib, k;
    float mus,muv,phi,ratio_spres,tau_ray,aot550;
    double coef;
    float actual_rho_ray,actual_T_ray_up,actual_T_ray_down,actual_S_r;
    float rho_ray_P0,T_ray_up_P0,T_ray_down_P0,S_r_P0;
    float lamda[7]={486.,570.,660.,835.,1669.,0.,2207.};
    float tau_ray_sealevel[7]={0.16511,0.08614,0.04716,0.01835,0.00113,0.00037};
        /* index=5 => band 7 */

    mus=cos(ar_gridcell->sun_zen[ipt]*RAD);
    muv=cos(ar_gridcell->view_zen[ipt]*RAD);
    phi=ar_gridcell->rel_az[ipt];
    ratio_spres=ar_gridcell->spres[ipt]/1013.;
    if (bkgd_aerosol) {
        atmos_coef->computed[ipt]=1;
        aot550=0.01;
    }
    else {
        if (line_ar != lut->aerosol_fill) {
            atmos_coef->computed[ipt]=1;
            aot550=((float)line_ar/1000.)*lamda[0]/550;
        }
        else {
            atmos_coef->computed[ipt]=1;
            aot550=0.01;
        }
    }

    for (k=1;k<SIXS_NB_AOT;k++) {
        if (aot550 < sixs_tables->aot[k])
            break;
    }
    k--;

    if (k>=(SIXS_NB_AOT-1))
        k=SIXS_NB_AOT-2;
    coef=(aot550-sixs_tables->aot[k])/(sixs_tables->aot[k+1]-
        sixs_tables->aot[k]);

    for (ib=0;ib < nband; ib++) {
        atmos_coef->tgOG[ib][ipt]=sixs_tables->T_g_og[ib];                
        atmos_coef->tgH2O[ib][ipt]=sixs_tables->T_g_wv[ib];                
        atmos_coef->td_ra[ib][ipt] = sixs_tables->T_ra_down[ib][k]
                                   + coef*(sixs_tables->T_ra_down[ib][k+1] -
                                           sixs_tables->T_ra_down[ib][k]);
        atmos_coef->tu_ra[ib][ipt] = sixs_tables->T_ra_up[ib][k]
                                   + coef*(sixs_tables->T_ra_up[ib][k+1] -
                                           sixs_tables->T_ra_up[ib][k]);
        atmos_coef->rho_mol[ib][ipt]=sixs_tables->rho_r[ib];                
        atmos_coef->rho_ra[ib][ipt] = sixs_tables->rho_ra[ib][k]
                                    + coef*(sixs_tables->rho_ra[ib][k+1] -
                                            sixs_tables->rho_ra[ib][k]);
        atmos_coef->td_da[ib][ipt] = sixs_tables->T_a_down[ib][k]
                                   + coef*(sixs_tables->T_a_down[ib][k+1] -
                                           sixs_tables->T_a_down[ib][k]);
        atmos_coef->tu_da[ib][ipt] = sixs_tables->T_a_up[ib][k]
                                   + coef*(sixs_tables->T_a_up[ib][k+1] -
                                           sixs_tables->T_a_up[ib][k]);
        atmos_coef->S_ra[ib][ipt] = sixs_tables->S_ra[ib][k]
                                  + coef*(sixs_tables->S_ra[ib][k+1] -
                                          sixs_tables->S_ra[ib][k]);

        /**
        compute DEM-based pressure correction for each grid point
        **/
        tau_ray=tau_ray_sealevel[ib]*ratio_spres;
        chand(&phi,&muv,&mus,&tau_ray,&actual_rho_ray);
        actual_T_ray_down = (2./3. + mus + (2./3. - mus)*exp(-tau_ray/mus))
                          / (4./3. + tau_ray); /* downward */
        actual_T_ray_up = (2./3. + muv + (2./3. - muv)*exp(-tau_ray/muv))
                        / (4./3. + tau_ray);   /* upward */

        csalbr(&tau_ray,&actual_S_r);
                        
        rho_ray_P0=sixs_tables->rho_r[ib];
        T_ray_down_P0=sixs_tables->T_r_down[ib];
        T_ray_up_P0=sixs_tables->T_r_up[ib];
        S_r_P0=sixs_tables->S_r[ib];

        atmos_coef->rho_ra[ib][ipt] += actual_rho_ray - rho_ray_P0;
                                          /* will need to correct for uwv/2 */
        atmos_coef->td_ra[ib][ipt] *= actual_T_ray_down/T_ray_down_P0;
        atmos_coef->tu_ra[ib][ipt] *= actual_T_ray_up/T_ray_up_P0;
        atmos_coef->S_ra[ib][ipt] += actual_S_r - S_r_P0;
        atmos_coef->td_r[ib][ipt] = actual_T_ray_down;
        atmos_coef->tu_r[ib][ipt] = actual_T_ray_up;
        atmos_coef->S_r[ib][ipt] = actual_S_r;
        atmos_coef->rho_r[ib][ipt] = actual_rho_ray;
    }  /* for ib */
}


void sun_angles
(
    short jday,
    float gmt,
    float flat,
    float flon,
    float *ts,
    float *fs
)
{
    double mst,tst,tet,et,ha,delta;
    double dlat,amuzero,elev,az,caz,azim;
    double A1=.000075,A2=.001868,A3=.032077,A4=.014615,A5=.040849;
    double B1=.006918,B2=.399912,B3=.070257,B4=.006758;
    double B5=.000907,B6=.002697,B7=.001480;
      
    dlat=(double)flat*M_PI/180.;
                          
    /*      
    SOLAR POSITION (ZENITHAL ANGLE ThetaS,AZIMUTHAL ANGLE PhiS IN DEGREES)
    J IS THE DAY NUMBER IN THE YEAR
    MEAN SOLAR TIME (HEURE DECIMALE)
    */
    mst=gmt+(flon)/15.;
    tet=2.*M_PI*(double)jday/365.;

    /* TIME EQUATION (IN MN.DEC) */
    et=A1+A2*cos(tet)-A3*sin(tet)-A4*cos(2.*tet)-A5*sin(2.*tet);
    et=et*12.*60./M_PI;

    /* TRUE SOLAR TIME */
    tst=mst+et/60.;
    tst=(tst-12.);

    /* HOUR ANGLE */
    ha=tst*15.*M_PI/180.;

    /* SOLAR DECLINATION   (IN RADIAN) */
    delta=B1-B2*cos(tet)+B3*sin(tet)-B4*cos(2.*tet)+B5*sin(2.*tet)-
        B6*cos(3.*tet)+B7*sin(3.*tet);

    /* ELEVATION,AZIMUTH */
    amuzero=sin(dlat)*sin(delta)+cos(dlat)*cos(delta)*cos(ha);
    elev=asin(amuzero);
    az=cos(delta)*sin(ha)/cos(elev);
    if (az<-1.)
        az=-1;
    else if (az>1.)
        az=1.;
    caz=(-cos(dlat)*sin(delta)+sin(dlat)*cos(delta)*cos(ha))/cos(elev);
    azim=asin(az);
    if (caz < 0.) 
        azim=M_PI-azim;
    if ((caz > 0.) && (az < 0.)) 
        azim=2*M_PI+azim;
    azim=azim+M_PI;
    if (azim > (2.*M_PI)) 
        azim=azim-2.*M_PI;
    elev=elev*180./M_PI;

    /* CONVERSION IN DEGREES */
    *ts=90.-elev;
    *fs=azim*180./M_PI;
}
