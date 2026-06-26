/******************************************************************************
FILE: ncep_reanalysis_surface_repackage

PURPOSE: Defines the functions used to repackage the NCEP Reanalysis datasets
into products used as input by LEDAPS.

PROJECT: Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE: NASA Open Source Agreement Version 1.3

HISTORY:
Date        Programmer       Reason
----------  ---------------  -------------------------------------
10/27/2014  Gail Schmidt     Modified original application to read the new
                             netCDF4 files (vs. netCDF3 files) from NOAA/NCEP.
                             Modifications also included the need to handle
                             new data types for the NCEP variables.

NOTES:
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "netcdf.h"
#include "hdf.h"
#include "mfhdf.h"

/* Defines for processing */
#define CLIMATE_NDIMS 3
#define CLIMATE_NVARS 4
#define CLIMATE_XDIM_NAME "lon"
#define CLIMATE_YDIM_NAME "lat"
#define CLIMATE_TDIM_NAME "time"

/* Prototypes for accessory functions */
int copy_sds (int ncid, int nvars, size_t dimsizes[], char *var_name,
    int32 first_time_index, int32 sdout_id, int verbose);
int get_size (nc_type data_type);
char *get_dt_string (nc_type data_type);
int32 get_hdf_dt (nc_type data_type);


/******************************************************************************
METHOD: main

PURPOSE: Handles the end to end processing of the annual NCEP input files and
repackages them to output daily HDF files (for the user-specified DOY),
containing multiple NCEP variables (air temp, precipitable water,
surface pressure).  Global attributes are copied to the new HDF file.  The
lat/long dimensions are also copied, along with their attributes.

RETURN VALUE:
Type = int
Value  Description
-----  -----------
-1     Error processing
-99    Error reading the input NCEP file
0      Successful processing

HISTORY:
Date        Programmer       Reason
----------  ---------------  -------------------------------------
10/27/2014  Gail Schmidt     Modified original application to read the new
                             netCDF4 files (vs. netCDF3 files) from NOAA/NCEP.
                             Modifications also included the need to handle
                             new data types for the NCEP variables.

NOTES:
******************************************************************************/
int main(int argc,char **argv) {
    int32 sdout_id;          /* HDF ID for output file */
    int verbose;             /* should status messages be written? */
    int index;               /* index for variables and attributes */
    size_t start[MAX_VAR_DIMS], cnt[MAX_VAR_DIMS];
    int32 first_time_index;  /* location in the time dimension for the 3D
                                variable to start reading (-1 if this isn't
                                one of the 3D variables) */
    char name[MAX_NC_NAME];  /* attribute name */
    void *buffer = NULL;     /* buffer for reading attribute values */
    double *timebuff = NULL; /* array of time values */
    int write_metadata;      /* should the global metadata and lat/long
                                dimensions be written to the output HDF file? */
    int doy;                 /* DOY for the current run */
    int base_date[3];        /* base date for NCEP file (year, month, day) */

    /* netCDF input vars */
    char varname[MAX_NC_NAME+1]; /* var names as read from netCDF file */
    char dimname[MAX_NC_NAME+1]; /* dim names as read from netCDF file */
    int ncid;                /* netCDF ID */
    int ndims;               /* number of input dimensions in netCDF file */
    int nvars;               /* number of variables in the netCDF file */
    int primary_index;       /* index of the primary variable (time) */
    int nb_globattrs;        /* number of global attributes in netCDF file */
    int var_ndims;           /* num dims for each var */
    int var_dimids[MAX_VAR_DIMS]; /* array for the dimension IDs */
    int var_natts;           /* number of var attributes */
    nc_type data_type;       /* data type for each variable */
    size_t dimsizes[MAX_VAR_DIMS]; /* dimension sizes */
    size_t count;            /* count of the attributes */

    if (argc != 4) {
        fprintf (stderr, "usage: %s <input> <output> <doy>\n", argv[0]);
        exit(-1);
    }
    verbose = 1;

/****
    open input netCDF4 file
****/
    if (nc_open (argv[1], NC_NOWRITE, &ncid)) {
        fprintf (stderr, "Error opening netCDF file: %s\n", argv[1]);
        exit(-99);
    }

/****
    check and open output, if the file doesn't exist then create it and
    set the write metadata flag so that the file attributes and lat/long
    dimensions are written
****/
    write_metadata = 0;
    if ((sdout_id = SDstart(argv[2], DFACC_RDONLY)) < 0) {
        if ((sdout_id = SDstart(argv[2], DFACC_CREATE)) < 0) {
           fprintf (stderr, "can't create output %s\n", argv[2]);
           exit(-1);
        }
        write_metadata = 1;
    } else {
        SDend (sdout_id);
        if ((sdout_id = SDstart (argv[2], DFACC_WRITE)) < 0) {
           fprintf (stderr, "can't open output %s\n", argv[2]);
           exit(-1);
        }
    }
    doy = (int16) atoi (argv[3]);        
    
/****
    determine how many netCDF variables, dimensions, and global attributes
    are in the file; also the dimension id of the unlimited dimension, if
    there is one.
****/
    if (nc_inq (ncid, &ndims, &nvars, &nb_globattrs, NULL)) {
        fprintf (stderr, "Error inquiring about the variables, dimensions, "
            "global attributes, etc. for the netCDF file: %s\n", argv[1]);
        exit(-1);
    }

    if (verbose) {
        printf ("Number of global attributes = %d\n",nb_globattrs);
        printf ("Number of variables = %d\n", nvars);
        printf ("Number of dimensions = %d\n", ndims);
    }

/****
    Get information about each of the dimensions
****/
    for (index = 0; index < ndims; index++) {
        if (nc_inq_dim (ncid, index, dimname, &dimsizes[index])) {
            fprintf (stderr, "Error inquiring about dimension %d", index);
            exit(-1);
        }
        printf ("Dimension %d: %s - %d\n", index, dimname,
            (int) dimsizes[index]);
    }

/****
    Copy global attributes and lat/long dimensions to output if creating
    a new file
****/
    if (write_metadata) {
        for (index = 0; index < nb_globattrs; index++) {
            if (nc_inq_attname (ncid, NC_GLOBAL, index, name)) {
                fprintf (stderr, "Error inquiring about global attribute "
                    "%d (0-based)\n", index);
                exit(-1);
            }

            if (nc_inq_att (ncid, NC_GLOBAL, name, &data_type, &count) == -1) {
                fprintf (stderr, "Error inquiring about global attribute "
                    "%s\n", name);
                exit(-1);
            }

#ifdef DEBUG
            if (verbose)
                printf ("attribute #%d = %s (%s - %d bytes) size = %d\n",
                    index, name, get_dt_string (data_type),
                    get_size (data_type), (int) count);
#endif

            /* Allocate a buffer to hold the attribute data. */
            buffer = malloc ((int) count * get_size(data_type));
            if (buffer == NULL) {
                fprintf (stderr, "Error allocating memory for global attr "
                    "%s\n", name);
                exit(-1);
            }

            /* Read the attribute data. */
            if (nc_get_att (ncid, NC_GLOBAL, name, buffer) != NC_NOERR) {
                fprintf (stderr, "Error getting global attribute %s\n", name);
                exit(-1);
            }

            /* Write the attribute data. */
            if (SDsetattr (sdout_id, name, get_hdf_dt (data_type), (int) count,
                buffer) < 0) {
                fprintf (stderr, "Error writing global attribute %s\n", name);
                exit(-1);
            }

            /* Free buffer */
            free (buffer);
        } /* end for */

        /* Get information about the variables in the file */
        primary_index = -1;     /* initialize index to invalid value */
        for (index = 0; index < nvars; index++) {
            if (nc_inq_var (ncid, index, varname, &data_type, &var_ndims,
                var_dimids, &var_natts)) {
                fprintf (stderr, "Error inquiring about variable %d\n", index);
                exit(-1);
            }

            if (verbose) {
                printf ("Variable %d: %s\n", index, varname);
                printf ("  ndims - %d\n", var_ndims); 
                printf ("  natts - %d\n", var_natts); 
            }

            /* If this variable is the time variable, then store the index */
            if (!strcmp (varname, CLIMATE_TDIM_NAME))
                primary_index = index;
        }

        /* Make sure the "time" variable was found */
        if (primary_index == -1) {
            fprintf (stderr, "%s variable was not found in the netCDF "
                "dataset.\n", CLIMATE_TDIM_NAME);
            exit(-1);
        }

/****
        Determine the base_date of this file by reading the time variable
        this is a double variable.
****/
        if (nc_inq_var (ncid, primary_index, varname, &data_type, &var_ndims,
            var_dimids, &var_natts)) {
            fprintf (stderr, "Error inquiring about variable %d\n",
                primary_index);
            exit(-1);
        }
        if (data_type != NC_DOUBLE) {
            fprintf (stderr, "Error: the %s variable is expected to be "
                "a double, but it's %s\n", CLIMATE_TDIM_NAME,
                get_dt_string (data_type));
            exit(-1);
        }

        timebuff = malloc (dimsizes[var_dimids[0]] * sizeof (double));
        if (timebuff == NULL) {
            fprintf (stderr, "Error allocating memory for timebuff\n");
            exit(-1);
        }

        start[0] = 0;
        cnt[0] = var_dimids[0];
        if (nc_get_vara_double (ncid, primary_index, start, cnt, timebuff)) {
            fprintf (stderr, "Error reading data from %s variable\n",
                CLIMATE_TDIM_NAME);
            exit(-1);
        }

        /* compute the year using the first time value in the file; these
           values represent hours since 1800-01-01 00:00:0.0 */
        base_date[0] = (int16) (timebuff[0] / 8765.81277) + 1800;
        printf ("year %d\n", base_date[0]);
        base_date[1] = 1;
        base_date[2] = 1;

/****
        Write base_date and Day Of Year to output
****/
        if (SDsetattr (sdout_id, "base_date", DFNT_INT16, 3, &base_date) < 0) {
            fprintf (stderr, "Error writing global attribute base_date\n");
            exit(-1);
        }
        if (SDsetattr (sdout_id, "Day Of Year", DFNT_INT16, 1, &doy) < 0) {
            fprintf (stderr, "Error writing global attribute Day of Year\n");
            exit(-1);
        }
        
/****
        Copy lat/lon SDS
****/
        if (copy_sds (ncid, nvars, dimsizes, CLIMATE_YDIM_NAME, -1, sdout_id,
            verbose)) {
            fprintf (stderr, "ERROR: couldn't copy SDS %s ... ABORT\n",
                CLIMATE_YDIM_NAME);
            exit(-1);
        }

        if (copy_sds (ncid, nvars, dimsizes, CLIMATE_XDIM_NAME, -1, sdout_id,
            verbose)) {
            fprintf (stderr, "ERROR: couldn't copy SDS %s ... ABORT\n",
                CLIMATE_XDIM_NAME);
            exit(-1);
        }        
    }  /* end write metadata */

    if (verbose)      
        printf ("Number of variables = %d\n", nvars);

/****
    Loop over variables, get name and number of attributes for each one
****/
    for (index = 0; index < nvars; index++) {
        if (nc_inq_var (ncid, index, varname, &data_type, &var_ndims,
            var_dimids, &var_natts)) {
            fprintf (stderr, "Error inquiring about variable %d\n", index);
            exit(-1);
        }

        if (!strcmp (varname, "pres") || !strcmp (varname, "pr_wtr") ||
            !strcmp (varname, "slp") || !strcmp (varname, "air")) {
            if (verbose) {
                printf ("SDS #%d = %s (%s - %d bytes) Rank = %d "
                    "Nb Attrs = %d\n", index, varname,
                    get_dt_string (data_type), get_size (data_type),
                    var_ndims, var_natts);

                printf("\tDims = ");
                for (index = 0; index < var_ndims; index++) {
                    if (index == 0)
                        printf ("%d", (int) dimsizes[var_dimids[index]]);
                    else
                        printf ("x%d", (int) dimsizes[var_dimids[index]]);
                }
                printf ("\n");
            }

            /* Copy this SDS, starting at the specified index position within
               the netCDF file for the time dimension */
            first_time_index = (doy - 1) * 4;
            if (copy_sds (ncid, nvars, dimsizes, varname, first_time_index,
                sdout_id, verbose)) {
                fprintf (stderr, "ERROR: couldn't copy SDS %s ... ABORT\n",
                    varname);
                exit(-1);
            }        
        }
    }  /* end for */

/****
    Close input & output
****/
    nc_close (ncid);
    SDend (sdout_id);
    return 0;
}   


/******************************************************************************
METHOD: copy_sds

PURPOSE: Copies the specified netCDF variable to the output HDF file as a
new SDS.  Copies the attributes to the SDS as well as sets the dimension
names.

RETURN VALUE:
Type = int
Value  Description
-----  -----------
-1     Error processing
0      Successful processing

HISTORY:
Date        Programmer       Reason
----------  ---------------  -------------------------------------
10/27/2014  Gail Schmidt     Modified original application to read the new
                             netCDF4 files (vs. netCDF3 files) from NOAA/NCEP.
                             Modifications also included the need to handle
                             new data types for the NCEP variables.
11/13/2014  Gail Schmidt     Added a check to the variable data type to make
                             sure it is a floating point, just to make sure the
                             old NCEP int16 products aren't being processed.

NOTES:
1. The NCEP variables are global datasets.  The starting point for the
   longitudinal values is 0 degrees.  The output HDF files need these variables
   to start at -180.  So, some minor data manipulation is done after reading
   the NCEP variables and before writing them to the output HDF file.
******************************************************************************/
int copy_sds
(
    int ncid,                 /* I: netCDF file ID for input */
    int nvars,                /* I: number of variables in netCDF file */
    size_t dimsizes[],        /* I: dimension sizes for netCDF file */
    char *var_name,           /* I: name of the variable to be copied to the
                                    output file */
    int32 first_time_index,   /* I: location in the time dimension for the 3D
                                    variable to start reading */
    int32 sdout_id,           /* I: HDF file ID for output */
    int verbose               /* I: should status messages be written? */
)
{
    int32 sdsout_id, start_hdf[MAX_VAR_DIMS], edge_hdf[MAX_VAR_DIMS];
    int32 dimout_id;
    int32 hdf_dim_sizes[MAX_VAR_DIMS];
    int buf_size;
    int il;

    int index;               /* index for variables and attributes */
    size_t start[MAX_VAR_DIMS], cnt[MAX_VAR_DIMS];
    char name[MAX_NC_NAME];  /* attribute name */
    char *buffer = NULL;     /* buffer for reading attribute values */
    char *oneline = NULL;    /* one line of input data */

    /* netCDF input vars */
    char varname[MAX_NC_NAME+1]; /* var names as read from netCDF file */
    int primary_index;       /* index of the primary variable (time) */
    int var_ndims;           /* num dims for each var */
    int var_dimids[MAX_VAR_DIMS]; /* array for the dimension IDs */
    int var_natts;           /* number of var attributes */
    nc_type data_type;       /* data type for each variable */
    size_t count;            /* count of the attributes */

    /* Get information about the variables in the file */
    primary_index = -1;     /* initialize index to invalid value */
    for (index = 0; index < nvars; index++) {
        if (nc_inq_var (ncid, index, varname, &data_type, &var_ndims,
            var_dimids, &var_natts)) {
            fprintf (stderr, "Error inquiring about variable %d\n", index);
            return (-1);
        }

        /* If this variable is the specified variable, then store the index */
        if (!strcmp (varname, var_name)) {
            primary_index = index;
            break;
        }
    }

    /* Make sure the variable was found */
    if (primary_index == -1) {
        fprintf (stderr, "%s variable was not found in the netCDF "
            "dataset.\n", var_name);
        return (-1);
    }

    /* If this is not one of the lat/long dimensions, then it's our NCEP
       variable.  In that case, we need to make sure this is a floating point
       data type and not int16 as the previous NCEP products were. */
    if (strcmp (var_name, CLIMATE_XDIM_NAME) &&
        strcmp (var_name, CLIMATE_YDIM_NAME) && 
        strcmp (var_name, CLIMATE_TDIM_NAME)) {
        if (data_type != NC_FLOAT) {
            fprintf (stderr, "Error: Non-dimensional variable (%s) should be "
                "floating point.\n", var_name);
            return -1;
        }
    }

    /* Set up the dimensions of the output HDF SDS */
    for (index = 0; index < var_ndims; index++)
        hdf_dim_sizes[index] = dimsizes[var_dimids[index]];

    /* If this is not one of the lat/long dimensions, then it's our NCEP
       variable.  In that case, we are only pulling data for the current DOY
       which is 4 values which represent data once every 6 hours. */
    if (strcmp (var_name, CLIMATE_XDIM_NAME) &&
        strcmp (var_name, CLIMATE_YDIM_NAME) && 
        strcmp (var_name, CLIMATE_TDIM_NAME))
        hdf_dim_sizes[0] = 4; 

    /* Create an SDS in the output file for this variable */
    if ((sdsout_id = SDcreate (sdout_id, var_name, get_hdf_dt (data_type),
        var_ndims, hdf_dim_sizes)) < 0) {
        fprintf (stderr, "Error creating SDS in output HDF file for %s\n",
            var_name);
        return -1;
    }
    
/****
    Read and write the data.  Dimensions get read and written in whole.
    The NCEP variables only read and write the four hours pertaining to the
    specified DOY.
****/
    if (first_time_index < 0) {
        /* Determine the size of the buffer needed for all the dimensions
           and data type */
        buf_size = 1;
        for (index = 0; index < var_ndims; index++) {
            buf_size *= hdf_dim_sizes[index];
            start[index] = 0;
            start_hdf[index] = 0;
            cnt[index] = hdf_dim_sizes[index];
            edge_hdf[index] = hdf_dim_sizes[index];
        }
        buf_size *= get_size (data_type);

        /* Allocate memory */
        buffer = malloc (buf_size);
        if (buffer == NULL) {
            fprintf (stderr, "Error allocating memory for %s\n", var_name);
            return (-1);
        }

        /* Read the data from the netCDF file */
        if (nc_get_vara (ncid, primary_index, start, cnt, buffer)) {
            fprintf (stderr, "Error reading data from %s variable\n", var_name);
            return (-1);
        }

        /* Write the data to the HDF file */
        if (SDwritedata (sdsout_id, start_hdf, NULL, edge_hdf, buffer) < 0) {
            fprintf (stderr, "Error writing %s data to SDS\n", var_name);
            return (-1);
        }

        /* Free the memory */
        free(buffer);
    } else {
        /* Determine the size of the buffer needed for the lat/long dimensions
           and data type */
        buf_size = hdf_dim_sizes[1] * hdf_dim_sizes[2] * get_size (data_type);

        /* Allocate memory */
        buffer = malloc (buf_size);
        if (buffer == NULL) {
            fprintf (stderr, "Error allocating memory for %s\n", var_name);
            return (-1);
        }

        /* Allocate memory for one line */
        oneline = malloc (hdf_dim_sizes[2] * get_size (data_type));
        if (oneline == NULL) {
            fprintf (stderr, "Error allocating memory for %s\n", var_name);
            return (-1);
        }

        /* Set up the start and edge arrays for reading and writing one line
           of time data at a time */
        start[0] = first_time_index;
        start_hdf[0] = 0;
        cnt[0] = 1;
        edge_hdf[0] = 1;
        start[1] = 0;
        start_hdf[1] = 0;
        cnt[1] = hdf_dim_sizes[1];
        edge_hdf[1] = hdf_dim_sizes[1];
        start[2] = 0;
        start_hdf[2] = 0;
        cnt[2] = hdf_dim_sizes[2];
        edge_hdf[2] = hdf_dim_sizes[2];

        /* Loop through the four time datasets for the specified DOY, using
           the first_time_index as the start of the current DOY */
        for (index = 0; index < 4; index++) {
            /* Read the data from the netCDF file */
            if (nc_get_vara (ncid, primary_index, start, cnt, buffer)) {
                fprintf (stderr, "Error reading data from %s variable for time "
                    "%d\n", var_name, index);
                return (-1);
            }

            /* Rearrange the NCEP variable data values since they start at
               0 degrees and we want to write the global values starting at
               -180 degrees */
            for (il = 0; il < hdf_dim_sizes[1]; il++) {
                memcpy (oneline,
                   &buffer[(il*hdf_dim_sizes[2]+hdf_dim_sizes[2]/2)*get_size(data_type)],
                   (hdf_dim_sizes[2]/2)*get_size(data_type));

                memcpy(&oneline[(hdf_dim_sizes[2]/2)*get_size(data_type)],
                    &buffer[il*hdf_dim_sizes[2]*get_size(data_type)],
                    (hdf_dim_sizes[2]-hdf_dim_sizes[2]/2)*get_size(data_type));

                memcpy(&buffer[il*hdf_dim_sizes[2]*get_size(data_type)],
                    oneline, hdf_dim_sizes[2]*get_size(data_type));
            }

            /* Write the data to the HDF file */
            if (SDwritedata (sdsout_id, start_hdf, NULL, edge_hdf, buffer)
                < 0) {
                fprintf (stderr, "Error writing %s data to SDS\n", var_name);
                return (-1);
            }

            /* Increment the start pointers for time dimension */
            start[0]++;
            start_hdf[0]++;
        }

        /* Free the memory */
        free (buffer);
        free (oneline);
    }

/****
    Loop over SDS attributes and copy
****/
    for (index = 0; index < var_natts; index++) {
        if (nc_inq_attname (ncid, primary_index, index, name)) {
            fprintf (stderr, "Error inquiring about variable attribute %d "
                "(0-based)\n", index);
            return (-1);
        }

        if (nc_inq_att (ncid, primary_index, name, &data_type, &count) == -1) {
            fprintf (stderr, "Error inquiring about variable attribute %s\n",
                name);
            return (-1);
        }

#ifdef DEBUG
        if (verbose)
            printf ("attribute #%d = %s (%s - %d bytes) size = %d\n",
                index, name, get_dt_string (data_type),
                get_size (data_type), (int) count);
#endif

        /* Allocate a buffer to hold the attribute data. */
        buffer = malloc ((int) count * get_size(data_type));
        if (buffer == NULL) {
            fprintf (stderr, "Error allocating memory for attr %s\n", name);
            return (-1);
        }

        /* Read the attribute data. */
        if (nc_get_att (ncid, primary_index, name, buffer) != NC_NOERR) {
            fprintf (stderr, "Error getting variable attribute %s\n", name);
            return (-1);
        }

        /* Write the attribute data. */
        if (SDsetattr (sdsout_id, name, get_hdf_dt (data_type), (int) count,
            buffer) < 0) {
            fprintf (stderr, "Error writing SDS attribute %s\n", name);
            return (-1);
        }

        /* Free buffer */
        free (buffer);
    } /* end for */

/****
    Set the SDS dimensions.  The XDIM, YDIM, and TDIM variables are 1D.
    The others are 3D.
****/
    if (!strcmp (var_name, CLIMATE_XDIM_NAME)) {
        dimout_id = SDgetdimid (sdsout_id, 0);
        SDsetdimname (dimout_id, CLIMATE_XDIM_NAME); 
    }
    else if (!strcmp (var_name, CLIMATE_YDIM_NAME)) {
        dimout_id = SDgetdimid (sdsout_id, 0);
        SDsetdimname (dimout_id, CLIMATE_YDIM_NAME); 
    }
    else if (!strcmp (var_name, CLIMATE_TDIM_NAME)) {
        dimout_id = SDgetdimid (sdsout_id, 0);
        SDsetdimname (dimout_id, CLIMATE_TDIM_NAME); 
    }
    else {
        dimout_id = SDgetdimid (sdsout_id, 0);
        SDsetdimname (dimout_id, CLIMATE_TDIM_NAME); 
        dimout_id = SDgetdimid (sdsout_id, 1);
        SDsetdimname (dimout_id, CLIMATE_YDIM_NAME); 
        dimout_id = SDgetdimid (sdsout_id, 2);
        SDsetdimname (dimout_id, CLIMATE_XDIM_NAME); 
    }

/****
    Close SDS
****/
    SDendaccess(sdsout_id);
    return 0;
}


/******************************************************************************
METHOD: get_size

PURPOSE: Determine the size (in bytes) of the netCDF data type specified

RETURN VALUE:
Type = int
Value  Description
-----  -----------
-99    Data type isn't supported or it was NAT (Not a Type)
Else   Size of the dataype (in bytes)

HISTORY:
Date        Programmer       Reason
----------  ---------------  -------------------------------------
10/27/2014  Gail Schmidt     Original development

NOTES:
******************************************************************************/
int get_size
(
    nc_type data_type       /* I: netCDF data type for each variable */
)
{
    switch (data_type)
    {
        case NC_BYTE:
        case NC_UBYTE:
        case NC_CHAR:
            return 1;

        case NC_SHORT:
        case NC_USHORT:
            return 2;

        case NC_INT:   /* covers NC_LONG which is deprecated */
        case NC_UINT:
            return 4;

        case NC_INT64:
        case NC_UINT64:
            return 8;

        case NC_FLOAT:
            return 4;

        case NC_DOUBLE:
            return 8;

        default:
            return -99;
    }
}


/******************************************************************************
METHOD: get_size

PURPOSE: Determine the data type string of the netCDF data type specified

RETURN VALUE:
Type = int
Value                    Description
-----                    -----------
"Unsupported data type"  Data type not supported
Else                     String for the input data type

HISTORY:
Date        Programmer       Reason
----------  ---------------  -------------------------------------
10/27/2014  Gail Schmidt     Original development

NOTES:
******************************************************************************/
char *get_dt_string
(
    nc_type data_type       /* I: netCDF data type for each variable */
)
{
    switch (data_type)
    {
        case NC_BYTE:
            return "NC_BYTE";
        case NC_UBYTE:
            return "NC_UBYTE";
        case NC_CHAR:
            return "NC_CHAR";
        case NC_SHORT:
            return "NC_SHORT";
        case NC_USHORT:
            return "NC_USHORT";
        case NC_INT:   /* covers NC_LONG which is deprecated */
            return "NC_INT";
        case NC_UINT:
            return "NC_UINT";
        case NC_INT64:
            return "NC_INT64";
        case NC_UINT64:
            return "NC_UINT64";
        case NC_FLOAT:
            return "NC_FLOAT";
        case NC_DOUBLE:
            return "NC_DOUBLE";
        case NC_NAT:
            return "NC_NAT";
        default:
            return "Unsupported data type";
    }
}


/******************************************************************************
METHOD: get_size

PURPOSE: Determine the HDF data type for the netCDF data type specified

RETURN VALUE:
Type = int
Value  Description
-----  -----------
-99    Data type isn't supported or it was NAT (Not a Type)
Else   HDF data type

HISTORY:
Date        Programmer       Reason
----------  ---------------  -------------------------------------
10/27/2014  Gail Schmidt     Original development

NOTES:
******************************************************************************/
int32 get_hdf_dt
(
    nc_type data_type       /* I: netCDF data type for each variable */
)
{
    switch (data_type)
    {
        case NC_BYTE:
            return DFNT_CHAR;
        case NC_UBYTE:
            return DFNT_UCHAR;
        case NC_CHAR:
            return DFNT_CHAR;
        case NC_SHORT:
            return DFNT_INT16;
        case NC_USHORT:
            return DFNT_UINT16;
        case NC_INT:   /* covers NC_LONG which is deprecated */
            return DFNT_INT32;
        case NC_UINT:
            return DFNT_UINT32;
        case NC_INT64:
            return DFNT_INT64;
        case NC_UINT64:
            return DFNT_UINT64;
        case NC_FLOAT:
            return DFNT_FLOAT;
        case NC_DOUBLE:
            return DFNT_FLOAT64;
        default:
            return -99;
    }
}
