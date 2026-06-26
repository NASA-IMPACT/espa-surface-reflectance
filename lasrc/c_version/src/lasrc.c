#include <sys/stat.h>
#include <unistd.h>
#include "lasrc.h"

/******************************************************************************
MODULE:  lasrc (Land Surface Reflectance Code - LaSRC)

PURPOSE:  Computes the surface reflectance values for the Landsat 8/9 and
Sentinel 2 products.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           An error occurred during processing of the surface reflectance
SUCCESS         Processing was successful

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
**Landsat 8/9:
1. Bands 1-7 are corrected to surface reflectance.  Band 8 (pand band) is not
   processed.  Band 9 (cirrus band) is corrected to TOA reflectance.  Bands
   10 and 11 are corrected to brightness temperature.
2. SDstart and SDreaddata have minor memory leaks.  Ultimately both call
   HAregister_atom which makes a malloc call and the memory is never freed.
3. Conversion algorithms for TOA reflectance and TOA brightness temperature are
   available from
   http://landsat.usgs.gov/Landsat8_Using_Product.php
4. The TOA and SR corrections utilize solar and view angles.  Previous versions
   simply corrected every pixel for the solar angles at the scene center.
   However as of version 1.0.0, the solar and view angles are read from the
   per-pixel angle bands and better reflect the angle at that location within
   the scene.  These solar/view angles are scaled int16s and therefore need to
   be unscaled (multiply by 0.01) before using.  They are in units of degrees.
5. The solar angles from band to band are fairly stable, thus a single array of
   per-pixel angles will be used.  The array can be from a "representative
   band" for the reflectance bands or it can be an average of the reflectance
   bands. We will use band 4 as the representative band for these per-pixel
   angle values.
6. The view/observation angles from band to band are a bit more unstable
   (particularly the view azimuth) near nadir.  A single array of per-pixel
   angles will still be used, however the information near nadir may be slightly
   affected.  We will use band 4 as the representative band for these per-pixel
   angle values.

**Sentinel 2:
1. Bands 1-12, including band 8a, are corrected to surface reflectance.
2. All bands are converted to 10m band resolution. Thus 20m and 60m bands are
   read and then converted to 10m.
******************************************************************************/
int main (int argc, char *argv[])
{
    bool verbose;            /* verbose flag for printing messages */
    char FUNC_NAME[] = "main"; /* function name */
    char errmsg[STR_SIZE];   /* error message */
    char envi_file[STR_SIZE];/* ENVI filename */
    char *aux_path = NULL;   /* path for Landsat auxiliary data */
    char *xml_infile = NULL; /* input XML filename */
    char *aux_infile = NULL; /* input auxiliary filename for water vapor
                                and ozone*/
    char *aux2_infile = NULL; /* second input auxiliary filename for water vapor
                                and ozone -- 23-FEB-26, JPR */
    char *cptr = NULL;       /* pointer to the file extension */
    char aux_year[5];        /* string to contain the year of auxiliary file */
    aux_src_t aux_src;       /* identifies the source of atmospheric aux data */

    int retval;              /* return status */
    int ib;                  /* looping variable for input bands */
    int i;                   /* looping variables */
    Sat_t sat;               /* satellite */
    Input_t *input = NULL;       /* input structure for the Landsat product */
    Output_t *toa_output = NULL; /* output structure and metadata for the TOA
                                    product */
    Espa_internal_meta_t xml_metadata;  /* XML metadata structure */
    Espa_global_meta_t *gmeta = NULL;   /* pointer to global meta */
    Envi_header_t envi_hdr;      /* output ENVI header information */
    struct stat statbuf;     /* buffer for the file stat function */

    int16 *sza = NULL;       /* Landsat per-pixel solar zenith angles,
                                nlines x nsamps */
    int16 *saa = NULL;       /* Landsat per-pixel solar azimuth angles,
                                nlines x nsamps */
    int16 *vza = NULL;       /* Landsat per-pixel view zenith angles,
                                nlines x nsamps */
    int16 *vaa = NULL;       /* Landsat per-pixel view azimuth angles,
                                nlines x nsamps */
    float **toaband = NULL;  /* Sentinel TOA unscaled reflectance bands,
                                nlines x nsamps */
    float **sband = NULL;    /* output surface reflectance and brightness
                                temp bands, qa band is separate as a uint16 */
    uint16 *out_band = NULL; /* scaled output */
    uint16 *qaband = NULL;   /* Landsat input QA band, Sentinel generated QA
                                band, nlines x nsamps */
    float xts;               /* scene center solar zenith angle (deg) */
    float xmus;              /* cosine of solar zenith angle */
    bool process_sr = true;  /* this is set to false if the solar zenith
                                is too large and the surface reflectance
                                cannot be calculated or if the user specifies
                                that surface reflectance processing will not
                                be completed and only TOA processing will be
                                done */
    bool write_toa = false;  /* this is set to true if the user specifies
                                TOA products should be output for delivery */
    bool use_orig_aero=false;/* this is set to true if the user specifies
                                the original aerosol algorithm should be used,
                                o/w the semi-empirical approach is used for
                                faster processing */
    float pixsize;           /* pixel size for the reflectance bands */
    int nlines, nsamps;      /* number of lines/samples in the reflectance and
                                thermal (Landsat) bands */

    /* The following arguments are all names of the LUTs. The first five are
       all tables of coefficients generated by the 6S software and provided
       as input to this application. */
    char anglehdf[STR_SIZE];  /* angle HDF filename */
    char intrefnm[STR_SIZE];  /* intrinsic reflectance filename */
    char transmnm[STR_SIZE];  /* transmission filename */
    char spheranm[STR_SIZE];  /* spherical albedo filename */
    char cmgdemnm[STR_SIZE];  /* climate modeling grid DEM filename */
    char rationm[STR_SIZE];   /* ratio averages filename ("ratio map" used by
                                 the aerosol retrieval algorithm) */
    char auxnm[STR_SIZE];     /* auxiliary filename for ozone and water vapor*/
    char aux2nm[STR_SIZE];    /* second auxiliary filename for ozone and water vapor*/
    
    /* GRIB-handling code -- 24-FEB-26, JPR) */
    int isGrb2, isAlsoGrb2, Grb2_file(char *filename); 
    int read_grib_date(FILE *input, char *what, char *where, char *date, int *frcst_hrs);  
    int read_grib2_date(char *filename, char *what, char *where, char *date, int *frcst_hrs);
    void get_new_where(char *where);
    FILE *fd;
    char where[50];
    char GRBdate[256];
    int frcst_hrs;

      
    /* Read the command-line arguments */
    retval = get_args (argc, argv, &xml_infile, &aux_infile, &aux2_infile, &process_sr,
        &write_toa, &use_orig_aero, &verbose);
    if (retval != SUCCESS)
    {   /* get_args already printed the error message */
        exit (ERROR);
    }

#ifdef _OPENMP
    printf ("LaSRC was compiled using multi-threading via -fopenmp\n");
#else 
    printf ("LaSRC was NOT compiled using multi-threading\n");
#endif
    printf ("Starting TOA and surface reflectance processing ...\n");
    if (use_orig_aero)
        printf ("Using the original aerosol inversion algorithm\n");
    else
        printf ("Using the semi-empirical approach for aerosol inversion\n");

    /* Provide user information if verbose is turned on */
    if (verbose)
    {
        printf ("  XML input file: %s\n", xml_infile);
        printf ("  AUX input file: %s\n", aux_infile);
        if (aux2_infile != NULL) printf ("  AUX2 input file: %s\n", aux2_infile);
        if (!process_sr)
        {
            printf ("    **Surface reflectance corrections will not be "
                "completed.  Only top of atmosphere corrections will be "
                "completed.\n");
        }
    }

    /* Validate the input metadata file */
    if (validate_xml_file (xml_infile) != SUCCESS)
    {  /* Error messages already written */
        exit (ERROR);
    }

    /* Initialize the metadata structure */
    init_metadata_struct (&xml_metadata);

    /* Parse the metadata file into our internal metadata structure; also
       allocates space as needed for various pointers in the global and band
       metadata */
    if (parse_metadata (xml_infile, &xml_metadata) != SUCCESS)
    {  /* Error messages already written */
        exit (ERROR);
    }

    /* Open the reflectance product, set up the input data structure, and
       allocate memory for the data buffers */
    input = open_input (&xml_metadata, use_orig_aero, process_sr);
    if (input == (Input_t *) NULL)
    {
        sprintf (errmsg, "Error opening/reading the input DN data: %s",
            xml_infile);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }
    gmeta = &xml_metadata.global;

    /* Verify that we have either Landsat 8/9 or Sentine-2 products */
    sat = input->meta.sat;
    if (sat != SAT_LANDSAT_8 && sat != SAT_LANDSAT_9 && sat != SAT_SENTINEL_2)
    {
        sprintf (errmsg, "The satellite identified from the input metadata "
            "is not Landsat 8/9 or Sentinel-2.  This application only supports "
            "Landsat 8/9 or Sentinel-2 products.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Output some information from the input files if verbose */
    if (verbose)
    {
        printf ("  Nband: %d\n", input->nband);
        printf ("  Number of lines/samples: %d/%d\n", input->size.nlines,
            input->size.nsamps);
        printf ("  Pixsize: %f,%f\n\n", input->size.pixsize[0],
            input->size.pixsize[1]);

        printf ("  Nband thermal: %d\n", input->nband_th);
        printf ("  Number of thermal lines/samples: %d/%d\n",
            input->size_th.nlines, input->size_th.nsamps);
        printf ("  Pixsize: %f,%f\n\n", input->size_th.pixsize[0],
            input->size_th.pixsize[1]);

        printf ("  Nband QA: %d\n", input->nband_qa);
        printf ("  Number of qa lines/samples: %d/%d\n",
            input->size_qa.nlines, input->size_qa.nsamps);
        printf ("  Pixsize: %f,%f\n\n", input->size_qa.pixsize[0],
            input->size_qa.pixsize[1]);

        printf ("  Fill value: %d\n", input->meta.fill);
        printf ("  Solar zenith: %f\n", xml_metadata.global.solar_zenith);
        printf ("  Solar azimuth: %f\n", xml_metadata.global.solar_azimuth);
    }

    /* Pull the needed metadata from the XML file and input structure */
    xts = gmeta->solar_zenith;
    pixsize = (float) input->size.pixsize[0];
    xmus = cos (xts * DEG2RAD);
    nlines = input->size.nlines;
    nsamps = input->size.nsamps;

    /* If this is a Sentinel product and TOA reflectance was requested, then
       it will not be written */
    if (sat == SAT_SENTINEL_2 && write_toa)
    {
        sprintf (errmsg, "This is a Sentinel-2 product and therefore TOA "
            "reflectance will not be written as an output.  Sentinel-2 L1C "
            "products are already TOA reflectance.  Surface reflectance "
            "corrections will be applied.");
        error_handler (false, FUNC_NAME, errmsg);
    }

    /* The surface reflectance algorithm cannot be implemented for solar
       zenith angles greater than 76 degrees.  Need to flag if the current
       scene falls into that category. */
    if (xts > 76.0 && process_sr)
    {
        sprintf (errmsg, "Solar zenith angle is too large to allow for surface "
            "reflectance processing.  Corrections must be limited to top-of-"
            "atmosphere reflectance and brightness temperature corrections. "
            "Use the --process_sr=false command-line argument. "
            "(solar zenith angle out of range)");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Allocate memory for all the data arrays. Note: per-pixel solar/view
       angles are only used for Landsat, toaband is for Sentinel only. */
    if (verbose)
        printf ("Allocating memory for the data arrays ...\n");
    retval = memory_allocation_main (sat, nlines, nsamps, use_orig_aero, &sza,
        &saa, &vza, &vaa, &qaband, &out_band, &sband, &toaband);
    if (retval != SUCCESS)
    {   /* get_args already printed the error message */
        sprintf (errmsg, "Error allocating memory for the data arrays from "
            "the main application.");
        error_handler (false, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Read the QA band for Landsat */
    if (sat == SAT_LANDSAT_8 || sat == SAT_LANDSAT_9)
    {
        if (get_input_qa_lines (input, 0, 0, nlines, qaband) != SUCCESS)
        {
            sprintf (errmsg, "Reading QA band");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }

    /* Read the scaled solar zenith per pixel angle bands, in degrees, for
       Landsat */
    if (sat == SAT_LANDSAT_8 || sat == SAT_LANDSAT_9)
    {
        if (get_input_ppa_lines (input, 0, nlines, use_orig_aero, sza, saa,
            vza, vaa) != SUCCESS)
        {
            sprintf (errmsg, "Reading per-pixel solar angle bands");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }

    /* Is the atmospheric aux file VIIRS (VJ104ANC/VNP04ANC) or MODIS (L8ANC)  */
    
    aux_src = AUXNULL;
    if (strstr(aux_infile, "04ANC") != NULL)
        aux_src = VIIRS;
    if (strstr(aux_infile, "L8ANC") != NULL)
        aux_src = MODIS;
	
    /* Get the auxiliary directory and the full pathname of the auxiliary
       files to be read if processing surface reflectance */
    if (process_sr)
    {
        /* Get the path for the auxiliary products from the LASRC_AUX_DIR
           environment variable.  If it isn't defined, then assume the products
           are in the local directory. */
        aux_path = getenv ("LASRC_AUX_DIR");
        if (aux_path == NULL)
        {
            aux_path = ".";
            sprintf (errmsg, "LASRC_AUX_DIR environment variable isn't "
                "defined. It is assumed the auxiliary products will be "
                "available from the local directory.");
            error_handler (false, FUNC_NAME, errmsg);
        }
	
        /* Is the atmospheric aux file VIIRS (VJ104ANC/VNP04ANC) or MODIS (L8ANC) ... or GRIB (23-FEB-26, JPR) */
	
        /* Note: GRIB data has filenames with considerable variation, and you can't 
      	 * know what you have just from the filename.  Need to run Grb2_file(char *filename) on the
     	 * filename, and for that the files have to be at LASRC_AUX_DIR/LADS/, NOT LASRC_AUX_DIR/LADS/<year>,
     	 * because <year> also cannot be reliably derived from the filename.
     	 */
	if ((aux_src != VIIRS) && (aux_src != MODIS)) {
	 
            sprintf (auxnm, "%s/LADS/%s", aux_path, aux_infile);
            if ((isGrb2 = Grb2_file(auxnm)) > -1) 
     	       aux_src = GRIB;
             }
	     
        /* If we have a second file and the first was GRIB, make sure the second one is the same type of GRIB -- 23-FEB-26, JPR */
        if ((aux2_infile != NULL) && (aux_src == GRIB)) {  /* if first isn't GRIB, ignore the second -- in other words we'll only use the second for GRIB files.  */
	       
              sprintf (aux2nm, "%s/LADS/%s", aux_path, aux2_infile);
              if ((isAlsoGrb2 = Grb2_file(aux2nm)) > -1) {
     	          if (isGrb2 != isAlsoGrb2) {
     	            sprintf (errmsg, "both GRIB auxiliary products must match.");
    	            error_handler (false, FUNC_NAME, errmsg);
		     }
                  }
              else {
     	         sprintf (errmsg, "both auxiliary products must match.");
    	         error_handler (false, FUNC_NAME, errmsg);
    	      }
    	   }
	   

        /* Grab the year of the auxiliary input file to be used for the correct
           location of the auxiliary file in the auxiliary directory */
        if (aux_src == VIIRS)
            strncpy (aux_year, &aux_infile[10], 4);
        else if (aux_src == MODIS)
            strncpy (aux_year, &aux_infile[5], 4);
        else {  /* perhaps not necessary for GRIB files, but just in case... */
	    strcpy(where, "atmos col");
	    if (isGrb2 == 1) {
	       fd = fopen(auxnm, "rb");
	       read_grib_date(fd, "PWAT", where, GRBdate, &frcst_hrs);  
	       //printf("  xxxxxxxx GRBdate is %s\n", GRBdate);	       
	       fclose(fd);
	       }
	    else if (isGrb2 == 2) {
	       get_new_where(where);
	       read_grib2_date(auxnm, "PWAT", where, GRBdate, &frcst_hrs);  
	       //printf("  rrxxxxxxxx GRBdate is %s\n", GRBdate);	       
	       }	       
	    strncpy (aux_year, &GRBdate[0], 4);
	    }
        aux_year[4] = '\0';

        /* Set up the look-up table files and make sure they exist */
        if (sat == SAT_LANDSAT_8 || sat == SAT_LANDSAT_9)
        {
            sprintf (anglehdf, "%s/LDCMLUT/ANGLE_NEW.hdf", aux_path);
            sprintf (intrefnm, "%s/LDCMLUT/RES_LUT_V3.0-URBANCLEAN-V2.0.hdf",
                aux_path);
            sprintf (transmnm,"%s/LDCMLUT/TRANS_LUT_V3.0-URBANCLEAN-V2.0.ASCII",
                aux_path);
            sprintf (spheranm, "%s/LDCMLUT/AERO_LUT_V3.0-URBANCLEAN-V2.0.ASCII",
                aux_path);
        }
        else if (sat == SAT_SENTINEL_2)
        {
            sprintf (anglehdf, "%s/MSILUT/ANGLE_NEW.hdf", aux_path);
            sprintf (intrefnm, "%s/MSILUT/RES_LUT_V3.0-URBANCLEAN-V3.0.hdf",
                aux_path);
            sprintf (transmnm,"%s/MSILUT/TRANS_LUT_V3.0-URBANCLEAN-V3.0.ASCII",
                aux_path);
            sprintf (spheranm, "%s/MSILUT/AERO_LUT_V3.0-URBANCLEAN-V3.0.ASCII",
                aux_path);
        }

        sprintf (cmgdemnm, "%s/CMGDEM.hdf", aux_path);
        sprintf (rationm, "%s/ratiomapndwiexp.hdf", aux_path);
        if (aux_src != GRIB ) {  /* if it's GRIB, they've already been set */
	    sprintf (auxnm, "%s/LADS/%s/%s", aux_path, aux_year, aux_infile);   
            sprintf (aux2nm, "%s/LADS/%s/%s", aux_path, aux_year, aux2_infile);
	   }

        if (stat (anglehdf, &statbuf) == -1)
        {
            sprintf (errmsg, "Could not find anglehdf data file: %s\n  Check "
                "LASRC_AUX_DIR environment variable.", anglehdf);
            error_handler (false, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        if (stat (intrefnm, &statbuf) == -1)
        {
            sprintf (errmsg, "Could not find intrefnm data file: %s\n  Check "
                "LASRC_AUX_DIR environment variable.", intrefnm);
            error_handler (false, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        if (stat (transmnm, &statbuf) == -1)
        {
            sprintf (errmsg, "Could not find transmnm data file: %s\n  Check "
                "LASRC_AUX_DIR environment variable.", transmnm);
            error_handler (false, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        if (stat (spheranm, &statbuf) == -1)
        {
            sprintf (errmsg, "Could not find spheranm data file: %s\n  Check "
                "LASRC_AUX_DIR environment variable.", spheranm);
            error_handler (false, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        if (stat (cmgdemnm, &statbuf) == -1)
        {
            sprintf (errmsg, "Could not find cmgdemnm data file: %s\n  Check "
                "LASRC_AUX_DIR environment variable.", cmgdemnm);
            error_handler (false, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        if (stat (rationm, &statbuf) == -1)
        {
            sprintf (errmsg, "Could not find rationm data file: %s\n  Check "
                "LASRC_AUX_DIR environment variable.", rationm);
            error_handler (false, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        if (stat (auxnm, &statbuf) == -1)
        {
            sprintf (errmsg, "Could not find auxnm data file: %s\n  Check "
                "LASRC_AUX_DIR environment variable.", auxnm);
            error_handler (false, FUNC_NAME, errmsg);
            exit (ERROR);
        }
        if ((aux2_infile != NULL) && (stat (aux2nm, &statbuf) == -1) && (aux_src == GRIB))
        {
            sprintf (errmsg, "Could not find aux2nm data file: %s\n  Check "
                "LASRC_AUX_DIR environment variable.", aux2nm);
            error_handler (false, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }

    /* Landsat needs the TOA reflectance and brightness temp computed from the
       input Level-1 data product.  Sentinel is already TOA reflectance so just
       read it and convert all the bands to the same resolution. */
    if (sat == SAT_LANDSAT_8 || sat == SAT_LANDSAT_9)
    {
        /* Compute the TOA reflectance and TOA brightness temp */
        printf ("Calculating Landsat TOA reflectance and TOA brightness "
                "temps...\n");
        retval = compute_landsat_toa_refl (input, &xml_metadata, qaband, nlines,
            nsamps, gmeta->instrument, sza, sband);
        if (retval != SUCCESS)
        {
            sprintf (errmsg, "Error computing Landsat TOA reflectance and TOA "
                "brightness temperatures.");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }
    else if (sat == SAT_SENTINEL_2)
    {
        /* Read the TOA reflectance bands and return unscaled values */
        printf ("Reading Sentinel-2 TOA reflectance bands...\n");
        retval = read_sentinel_toa_refl (input, &xml_metadata, toaband);
        if (retval != SUCCESS)
        {
            sprintf (errmsg, "Error reading Sentinel-2 TOA reflectance bands.");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }

    /* Write the output TOA and BT bands for Landsat. Nothing is done for
       Sentinel. */
    if (sat == SAT_LANDSAT_8 ||sat == SAT_LANDSAT_9)
    {
        /* Open the TOA output file, and set up the bands according to whether
           the TOA reflectance bands will be written. */
        toa_output = open_output (&xml_metadata, input, OUTPUT_TOA, NULL);
        if (toa_output == NULL)
        {   /* error message already printed */
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
        printf ("Writing TOA reflectance corrected data to the output "
                "files ...\n");

        /* If we are writing the TOA data, do so now for bands 1-7.  This will
           occur if the user specified TOA to be written or if the surface
           reflectance processing will not be completed. */
        if (write_toa || !process_sr)
        {
            for (ib = SRL_BAND1; ib <= SRL_BAND7; ib++)
            {
                /* Convert and scale the floating point values */
                printf ("  Band %d: %s\n", ib+1,
                    toa_output->metadata.band[ib].file_name);
                convert_output (sband, ib, nlines, nsamps, false, out_band);

                /* Write the scaled integers */
                if (put_output_lines (toa_output, out_band, ib, 0, nlines,
                    sizeof (uint16)) != SUCCESS)
                {
                    sprintf (errmsg, "Writing output TOA data for band %d",
                        ib+1);
                    error_handler (true, FUNC_NAME, errmsg);
                    exit (ERROR);
                }

                /* Create the ENVI header file this band */
                if (create_envi_struct (&toa_output->metadata.band[ib],
                    &xml_metadata.global, &envi_hdr) != SUCCESS)
                {
                    sprintf (errmsg, "Creating ENVI header structure.");
                    error_handler (true, FUNC_NAME, errmsg);
                    exit (ERROR);
                }
          
                /* Write the ENVI header */
                strcpy (envi_file, toa_output->metadata.band[ib].file_name);
                cptr = strchr (envi_file, '.');
                strcpy (cptr, ".hdr");
                if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
                {
                    sprintf (errmsg, "Writing ENVI header file.");
                    error_handler (true, FUNC_NAME, errmsg);
                    exit (ERROR);
                }
            }

            /* Append the TOA reflectance bands, bands 1-7, to the XML file */
            if (append_metadata (7, toa_output->metadata.band, xml_infile) !=
                SUCCESS)
            {
                sprintf (errmsg, "Appending TOA reflectance bands to XML "
                    "file.");
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }
        }

        /* Write bands 9-11 (cirrus and thermals), which don't get any further
           processing. */
        for (ib = SRL_BAND9; ib <= SRL_BAND11; ib++)
        {
            /* If processing OLI-only, then bands 10 and 11 don't exist */
            if (!strcmp (gmeta->instrument, "OLI") &&
                (ib == SRL_BAND10 || ib == SRL_BAND11))
                continue;
            
            /* Convert and scale the floating point values. Band 9 TOA has a
               different scale factor than bands 10 and 11 BT. */
            printf ("  Band %d: %s\n", ib+2,
                toa_output->metadata.band[ib].file_name);
            if (ib == SRL_BAND9)
                convert_output (sband, ib, nlines, nsamps, false, out_band);
            else
                convert_output (sband, ib, nlines, nsamps, true, out_band);

            /* Write the scaled integers */
            if (put_output_lines (toa_output, out_band, ib, 0, nlines,
                sizeof (uint16)) != SUCCESS)
            {
                sprintf (errmsg, "Writing output TOA data for band %d", ib+2);
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }

            /* Create the ENVI header file this band */
            if (create_envi_struct (&toa_output->metadata.band[ib],
                &xml_metadata.global, &envi_hdr) != SUCCESS)
            {
                sprintf (errmsg, "Creating ENVI header structure.");
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }
  
            /* Write the ENVI header */
            strcpy (envi_file, toa_output->metadata.band[ib].file_name);
            cptr = strchr (envi_file, '.');
            strcpy (cptr, ".hdr");
            if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
            {
                sprintf (errmsg, "Writing ENVI header file.");
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }

            /* Append the TOA cirrus/thermal band to the XML file */
            if (append_metadata (1, &toa_output->metadata.band[ib], xml_infile)
                != SUCCESS)
            {
                sprintf (errmsg, "Appending TOA cirrus/thermal band to XML "
                    "file.");
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }
        }

        /* Close the output TOA products, cleanup bands, and free the memory */
        close_output (sat, toa_output, OUTPUT_TOA);
        if (process_sr && !write_toa)
        {
            /* Remove the TOA bands 1-7 that were created by the open routine,
               since they aren't actually used */
            for (ib = SRL_BAND1; ib <= SRL_BAND7; ib++)
                unlink (toa_output->metadata.band[ib].file_name);
        }
        free_output (toa_output, OUTPUT_TOA);
    }

    /* Only continue with the surface reflectance corrections if SR processing
       has been requested and is possible due to the solar zenith angle */
    if (process_sr)
    {
        /* Perform atmospheric correction for the reflectance bands and write
           the data to the SR output file */
        printf ("Performing atmospheric corrections for each reflectance "
            "band ...\n");
        if (sat == SAT_LANDSAT_8 || sat == SAT_LANDSAT_9)
        {
            retval = compute_landsat_sr_refl (input, &xml_metadata, xml_infile,
                qaband, out_band, nlines, nsamps, pixsize, sband, sza, saa,
                vza, vaa, xts, xmus, use_orig_aero, anglehdf, intrefnm,
                transmnm, spheranm, cmgdemnm, rationm, auxnm, aux2nm, aux_src);
            if (retval != SUCCESS)
            {
                sprintf (errmsg, "Error computing Landsat surface reflectance");
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }
        }
        else if (sat == SAT_SENTINEL_2)
        {
            retval = compute_sentinel_sr_refl (input, &xml_metadata, xml_infile,
                qaband, nlines, nsamps, pixsize, toaband, sband, out_band, xts,
                xmus, use_orig_aero, anglehdf, intrefnm, transmnm, spheranm,
                cmgdemnm, rationm, auxnm, aux2nm, aux_src);
            if (retval != SUCCESS)
            {
                sprintf (errmsg, "Error computing Sentinel-2 surface "
                    "reflectance");
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }
        }
    }  /* end if process_sr */
  
    /* Free the metadata structure */
    free_metadata (&xml_metadata);

    /* Close the input product */
    printf ("Closing input/output and freeing pointers ...\n");
    close_input (input, use_orig_aero);
    free_input (input, use_orig_aero);

    /* Free the filename pointers */
    free (xml_infile);
    free (aux_infile);
    if (aux2_infile != NULL) free (aux2_infile);

    /* Free memory for band data */
    free (qaband);
    free (out_band);
    if (sat == SAT_LANDSAT_8 || sat == SAT_LANDSAT_9)
    {
        free (sza);
        if (use_orig_aero)
        {
            free (saa);
            free (vza);
            free (vaa);
        }
        for (i = 0; i < NBANDL_TTL_OUT-1; i++)
            free (sband[i]);
    }
    else if (sat == SAT_SENTINEL_2)
    {
        for (i = 0; i < NBANDS_TTL_OUT-1; i++)
        {
            free (toaband[i]);
            free (sband[i]);
        }
        free (toaband);
    }
    free (sband);

    /* Indicate successful completion of processing */
    printf ("Surface reflectance processing complete!\n");
    exit (SUCCESS);
}


/******************************************************************************
MODULE:  usage

PURPOSE:  Prints the usage information for this application.

RETURN VALUE:
Type = None

NOTES:
******************************************************************************/
void usage ()
{
    printf ("LaSRC (Land Surface Reflectance Code) computes the surface "
            "reflectance values for the input Landsat 8/9 DN and Sentinel-2 "
            "L1C TOA reflectance products.\n"
            "Landsat 8/9: Surface reflectance correction and/or top of "
            "atmosphere correction is applied and written for bands 1-7.  "
            "Top of atmosphere and corrections are applied and written for "
            "bands 9 (cirrus), 10 (thermal), and 11 (thermal).  Surface "
            "reflectance corrections are available for OLI_TIRS products.  "
            "OLI-only scenes are corrected up through TOA and not surface "
            "reflectance.\n"
            "Sentinel 2: Surface reflectance correction is applied and "
            "written for bands 1-12.\n");
    printf ("version: %s\n\n", SR_VERSION);

    printf ("usage: lasrc "
            "--xml=input_xml_filename "
            "--aux=input_auxiliary_filename "
            "[--aux2=input_auxiliary_filename] "
            "--process_sr=true:false --write_toa [--use_orig_aero_alg] "
            "[--verbose] [--version]\n");

    printf ("\nwhere the following parameters are required:\n");
    printf ("    -xml: name of the input XML file to be processed\n");
    printf ("    -aux: name of the input auxiliary file containing ozone "
            "and water vapor for the scene date.  The file is expected to "
            "live in the $LASRC_AUX_DIR/LADS directory or in the local "
            "directory.\n");

    printf ("\nwhere the following parameters are optional:\n");
    printf ("    -process_sr: the default is to process surface reflectance, "
            "however if this flag is set to false then only the TOA "
            "reflectance processing (Landsat) and brightness temperature will "
            "be done. If set to false for Sentinel products, then nothing is "
            "done.\n");
    printf ("    -aux2: name of a second input auxiliary file containing ozone "
            "and water vapor for the scene date.  The file is expected to "
            "live in the $LASRC_AUX_DIR/LADS directory or in the local "
            "directory.  If both aux files are given, data is temporally interpolated "
	    "between them.  When two aux files are to be processed, both must be "
	    "in the GRIB or GRIB2 format\n");
    printf ("    -write_toa: the intermediate Landsat TOA reflectance products "
            "for bands 1-7 as well as the brightness temperature, are written "
            "to the output file. This argument has no relevance for Sentinel-2 "
            "products, since they are input as TOA reflectance, and therefore "
            "is ignored.\n");
    printf ("    -use_orig_aero_alg: specifies if the original FORTRAN-based "
            "algorithm should be used for the aerosol inversion. (default is "
            "to use the semi-empirical approach for accessing the look-up "
            "table which provides a speed-up in processing)\n");
    printf ("    -verbose: should intermediate messages be printed? (default "
            "is false)\n");
    printf ("    -version: print the LaSRC version. When this parameter is "
            "used, none of the other parameters are used or required.\n");

    printf ("\nlasrc --help will print the usage statement\n");
    printf ("\nExample: lasrc "
            "--xml=LC08_L1TP_041027_20130630_20140312_02_T1.xml "
            "--aux=VJ104ANC.A2013181.002.2022179071944.h5 --verbose\n");
    printf ("   ==> Writes bands 9-11 as TOA reflectance and brightness "
            "temperature.  Writes bands 1-7 as surface reflectance.\n\n");


    printf ("\nExample: lasrc "
            "--xml=LC08_L1TP_078014_20220306_20220314_02_T1.xml "
            "--aux=gdas1.PGrbF00.220306.18z --aux2=gdas1.PGrbF00.220307.00z --verbose\n");
    printf ("   ==> Writes bands 9-11 as TOA reflectance and brightness "
            "temperature.  Writes bands 1-7 as surface reflectance.  Temporally "
	    "interpolates water vapor and ozone data from two GRIB format files "
	    "bracketting the scene acquisition time.\n\n");

    printf ("\nExample: lasrc "
            "--xml=LC08_L1TP_041027_20130630_20140312_02_T1.xml "
            "--aux=VJ104ANC.A2013181.002.2022179071944.h5 --write_toa "
            "--verbose\n");
    printf ("   ==> Writes bands 1-11 as TOA reflectance and brightness "
            "temperature.  Writes bands 1-7 as surface reflectance.\n");

    printf ("\nExample: lasrc "
            "--xml=LC08_L1TP_041027_20130630_20140312_02_T1.xml "
            "--aux=VJ104ANC.A2013181.002.2022179071944.h5 --process_sr=false "
            "--verbose\n");
    printf ("   ==> Writes bands 1-11 as TOA reflectance and brightness "
            "temperature.  Surface reflectance corrections are not applied.\n");

    printf ("\nExample: lasrc "
            "--xml=L1C_T10TFR_A016455_20180816T190038.xml "
            "--aux=VJ104ANC.A2018228.002.2022179071944.h5 --verbose\n");
    printf ("   ==> Writes bands 1-12 as surface reflectance.\n\n");

    printf ("\nExample: lasrc "
            "--xml=L1C_T10TFR_A016455_20180816T190038.xml "
            "--aux=VJ104ANC.A2018228.002.2022179071944.h5 --use_orig_aero_alg "
            "--verbose\n");
    printf ("   ==> Writes bands 1-12 as surface reflectance but uses the "
            "original FORTRAN aerosol inversion algorithm.\n\n");

    printf ("\nExample: lasrc "
            "--xml=L1C_T10TFR_A016455_20180816T190038.xml "
            "--aux=VJ104ANC.A2018228.002.2022179071944.h5 --write_toa "
            "--verbose\n");
    printf ("   ==> Writes bands 1-12 as surface reflectance.\n\n");

    printf ("\nExample: lasrc "
            "--xml=L1C_T10TFR_A016455_20180816T190038.xml "
            "--aux=VJ104ANC.A2018228.002.2022179071944.h5 --process_sr=false "
            "--verbose\n");
    printf ("   ==> Does nothing.\n\n");
}


/******************************************************************************
MODULE:  btest

PURPOSE:  Tests to see if bit n is set in the byte_val variable.

RETURN VALUE:
Type = bool
Value      Description
-----      -----------
false      bit n is not set in byte_val
true       bit n is set in byte_val

NOTES:
******************************************************************************/
bool btest
(
    uint8 byte_val,   /* I: byte value to be tested with the bit n */
    byte n            /* I: bit number to be tested (0 is rightmost bit) */
)
{
    /* Take 2 ** n, then AND that result with the byte value */
    return (byte_val & (1 << n));
}

