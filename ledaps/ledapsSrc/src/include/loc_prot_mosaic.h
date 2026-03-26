/******************************************************************************

FILE:  loc_prot_mosaic.h

PURPOSE:  Function prototypes

HISTORY:
Version  Date   Programmer       Code  Reason
-------  -----  ---------------  ----  -------------------------------------
         03/02  Gail Schmidt           Initial development            
 
HARDWARE AND/OR SOFTWARE LIMITATIONS:  
  None

PROJECT:    MODIS Reprojection Tool

NOTES:

******************************************************************************/
#ifndef _LOC_PROT_MOSAIC_H_
#define _LOC_PROT_MOSAIC_H_

/*** Prototypes for Mosaicking ***/
intn AppendMetadataMosaic
(
    int numh_tiles,             /* I: number of horiz tiles in the mosaic */
    int numv_tiles,             /* I: number of vert tiles in the mosaic */
    int **tile_array,           /* I: 2D array of size [numv_tiles][numh_tiles]
                                      specifying which input file represents
                                      that tile location */
    MosaicDescriptor infiles[], /* I: file descriptor array for the input
                                      files */
    MosaicDescriptor *mosaic    /* I: mosaic info */
);

int CheckFieldsMosaic
(
    int nstrings,
    char *HeaderStrings[],
    int StringsPresent[],
    MosaicDescriptor *P
);

int CheckMosaicArgs
(
    int argc,              /* I: number of arguments */
    char *argv[],          /* I: argument strings */
    char ifilenames[][LARGE_STRING],
                           /* O: input filenames specified in the arguments */
    int *num_files,        /* O: number of input filenames specified */
    char ofilename[],      /* O: output filename specified in the arguments */
    char bandstr[],        /* O: string of spectral subset bands */
    int *determine_tiles,  /* O: was -t switch specified to determine the
                                 tiles for each input filename? */
    int *write_tmphdr,     /* O: was -h switch specified to write the
                                 raw binary header info for the mosaic? */
    int *spectral_subset   /* O: was -s switch specified for spectral
                                 subsetting? */
);

int CompareProducts
(
    int num_infiles,            /* number of input files */
    MosaicDescriptor infiles[]  /* file descriptor array for the input files */
);

void CopyMosaicDescriptor
(
    MosaicDescriptor *in,    /* I:  the source descriptor to be copied */
    char output_filename[],  /* I:  name of the output file */
    MosaicDescriptor *out    /* O:  the descriptor to be copied to */
);

FileDescriptor *CreateFileDescriptorMosaic
(
    MosaicDescriptor *mosaic,   /* I:  session info */
    int in_bandnum,             /* I:  number of the input band */
    int out_bandnum,            /* I:  number of the output band */
    FileOpenType fileopentype,  /* I:  reading or writing */
    char *filename              /* I:  the file name */
);

int CreateHdfEosFieldMosaic
(
    MosaicDescriptor *mosaic,   /* I: session info */
    int curr_band,              /* I: current output band to be processed */
    FileDescriptor *output      /* O: output file descriptor */
);

int CreateHdfEosGridMosaic
(
    MosaicDescriptor *mosaic,   /* I: session info */
    int curr_band,              /* I: current output band to be processed */
    FileDescriptor *output      /* O: output file descriptor */
);

void FillBufferBackground
(
    BandType *bandinfo,     /* I: band information for current band */
    int numh_tiles,         /* I: number of horizontal tiles */
    double *buffer          /* O: output buffer */
);

int GetHdfEosFieldMosaic
(
    HdfEosFD *hdfptr,           /* I/O: file to get field number */
    int fieldnum                /* I:   field to find */
);

int GetInputFilenames
(
    char *str,             /* I: argument string to be parsed */
    char ifilenames[][LARGE_STRING]
                           /* O: input filenames specified in the arguments
                                 (can't be more than NUM_ISIN_TILES strings) */
);

ProjInfo *GetInputProjectionMosaic
(
    MosaicDescriptor *mosaic     /* I:  session info */
);

int GetSpectralSubsetMosaic
(
    MosaicDescriptor *Mosaic,   /* O:  session info */
    char *str,                  /* I:  string to be parsed */
    int spectral_subset         /* I:  was spectral subsetting specified? */
);

void InitializeMosaicDescriptor
(
    MosaicDescriptor *P    /* I/O:  the descriptor to initialize */
);

FileDescriptor *MakeHdfEosFDMosaic
(
    MosaicDescriptor *mosaic,   /* I:  session info */
    HdfEosFD *hdfptr,           /* I:  HDF file descriptor */
    FileOpenType mode,          /* I:  read or write */
    int in_bandnum,             /* I:  number of the input band */
    int out_bandnum,            /* I:  number of the output band */
    int *status                 /* I/O:  error status code */
);

int MosaicTiles
(
    int numh_tiles,      /* I: number of horiz tiles in the mosaic */
    int numv_tiles,      /* I: number of vert tiles in the mosaic */
    int **tile_array,    /* I: 2D array of size [numv_tiles][numh_tiles]
                               specifying which input file represents that
                               tile location */
    int num_infiles,     /* I: number of input files */
    MosaicDescriptor infiles[],
                         /* I: file descriptor array for the input files */
    MosaicDescriptor *mosaicfile
                         /* I: file descriptor for the output mosaic file */
);

void MosaicUsage
(
    void
);

FileDescriptor *OpenInImageMosaic
(
    MosaicDescriptor *mosaic,   /* I:  session info */
    int in_bandnum,             /* I:  number of the input band */
    int out_bandnum,            /* I:  number of the output band */
    int *status                 /* O:  Error status */
);

FileDescriptor *OpenMultiFileMosaic
(
    MosaicDescriptor *mosaic,   /* I:  session info */
    FileOpenType mode,          /* I:  reading or writing */
    int in_bandnum,             /* I:  number of the input band */
    int out_bandnum,            /* I:  number of the output band */
    int *status                 /* O:  error status */
);

FileDescriptor *OpenOutImageMosaic
(
    MosaicDescriptor *mosaic,   /* I:  session info */
    int in_bandnum,             /* I:  number of the input band */
    int out_bandnum,            /* I:  number of the output band */
    int *status                 /* O:  Error status */
);

int OutputHdrMosaic
(
    MosaicDescriptor *mosaic,
    char *output_filename
);

void PrintOutputFileInfoMosaic
(
    MosaicDescriptor *mosaic
);

int read_tile_number_hdf
(
    char filename[],          /* I: the filename of the input files */
    int *horiz,               /* O: horizontal tile number */
    int *vert                 /* O: vertical tile number */
);

int read_tile_number_rb
(
    char filename[],          /* I: the filename of the input files */
    int *horiz,               /* O: horizontal tile number */
    int *vert                 /* O: vertical tile number */
);

int ReadHDFHeaderMosaic
(
    MosaicDescriptor *mosaic     /* I/O:  session info */
);

int ReadHeaderFileMosaic
(
    MosaicDescriptor *mosaic     /* I/O:  session info */
);

int ReadMetaMosaic
(
    MosaicDescriptor *P     /* I/O:  session info */
);

int ReadTiles
(
    char filenames[][LARGE_STRING], /* I: files to be mosaicked */
    int num_files                   /* I: number of files to be mosaicked */
);

int SortProducts
(
    int num_infiles,            /* I: number of input files */
    MosaicDescriptor infiles[], /* I: file descriptor array for the input
                                      files */
    char output_filename[],     /* I: name of the output file */
    MosaicDescriptor *mosaicfile, /* O: file descriptor for the output mosaic
                                        file */
    int *numh_tiles,            /* O: number of horiz tiles in the mosaic */
    int *numv_tiles,            /* O: number of vert tiles in the mosaic */
    int ***tile_array           /* O: 2D array of size [numv_tiles][numh_tiles]
                                      specifying which input file represents
                                      that tile location (space is allocated
                                      for this array) */
);

#endif
