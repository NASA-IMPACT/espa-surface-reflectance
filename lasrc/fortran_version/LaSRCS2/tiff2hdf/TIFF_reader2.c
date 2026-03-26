/*
 * TIFF_reader.c -- TIFF-format file input module.
 *
 *  After: Niles D. Ritter
 *
 *
 *  Jim Ray, SSAI
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "geotiff.h"
#include "geotiffio.h"
#include "xtiffio.h"
#include "geo_normalize.h"
#include "geovalues.h"
#include "tiffio.h"
#include "geo_tiffp.h"
#include "geo_keyp.h"
extern int makelatlon;

enum {VERSION=0,MAJOR,MINOR};
int tiff_read(char *filename, int *rwidth, int *rheight, int *rnbands, short **data[], unsigned char **data8[], 
              float **latbuffer, float **lonbuffer, int *latindex, char *TOC)
	{
		TIFF *tif=(TIFF*)0;  /* TIFF-level descriptor */
		GTIF *gtif=(GTIF*)0; /* GeoKey-level descriptor */
		int versions[3];
		int cit_length;
		geocode_t model;    /* all key-codes are of this type */
		char *citation;
		int i, size, width, height, orientation, nbands, bitspersample, sampleformat, retval, stlen, config;
		uint16 xuint16;
	        double tiepoints[24], pixscale[3];
                tagtype_t type;
                void ReadImage16(TIFF *tif, int width, int height, int nbands, short **image[], int config);
                void ReadImage8(TIFF *tif, int width, int height, int nbands, unsigned char **image[], int config);
		int make_lat_lon(GTIF *gtif, double tiepoints[24], int width, int height, float **latbuffer, float **lonbuffer);

	        char *localTOC;
		
		
		/* Open TIFF descriptor to read GeoTIFF tags */
		tif=XTIFFOpen(filename,"r");  
		if (!tif) {
		    printf("Error in XTIFFOpen, exiting\n");
		    exit(-1);
		 }
		
		/* Open GTIF Key parser; keys will be read at this time. */
		gtif = GTIFNew(tif);
		if (!gtif) {
		    printf("Error in GTIFNew, exiting\n");
		    exit(-1);
		 }

		/* Get the GeoTIFF directory info */
		GTIFDirectoryInfo(gtif,versions,0);
		if (versions[MAJOR] > 1)
		{
		    printf("Error in GTIFDirectoryInfo, exiting\n");
		    exit(-1);		 
		}
		/* Some TIF images give an error here: and since 'model' isn't 
		   used here any further, why bother?
		   This may perhaps be necessary for GEOTIFF files, but some TIFs
		   don't like it.
		
		if (!GTIFKeyGet(gtif, GTModelTypeGeoKey, &model, 0, 1))
		{
		    printf("Error in GTIFKeyGet, exiting\n");
		    exit(-1);
		}*/
		
		/* ASCII keys are variable-length; compute size */
		cit_length = GTIFKeyInfo(gtif,GTCitationGeoKey,&size,&type);
		if (cit_length > 0)
		{
		
			citation = (char *)malloc(size*cit_length);
			if (!citation) {printf("malloc error,exiting\n");exit(-2);}
			GTIFKeyGet(gtif, GTCitationGeoKey, citation, 0, cit_length);
			/*printf("Citation:%s\n",citation);*/
		}

		/* Get some TIFF info on this image */
		/*TIFFGetField(tif,TIFFTAG_IMAGEWIDTH,    &width);*/
	 	
	        TIFFGetField(tif,TIFFTAG_IMAGEWIDTH,    &width);
	        TIFFGetField(tif,TIFFTAG_IMAGELENGTH,   &height);
		
		/* Apparently, you don't need to allocate memory for 
		   whatever is read with TIFFTAG_IMAGEDESCRIPTION.
		   But as of 27-MAY-10 (and probably much earlier),
		   you do...  
		 * /
		localTOC = (char *) malloc(1000000);
	        TIFFGetField(tif,TIFFTAG_IMAGEDESCRIPTION,   &localTOC); 
 
                stlen = strlen(localTOC);
	        for(i=0;i<stlen;i++) TOC[i] = localTOC[i];
	        TOC[i] = '\0';
	        free(localTOC);
		Forget it -- this doesn't work consistently */

	        TIFFGetField(tif,TIFFTAG_ORIENTATION,   &xuint16);
		orientation = (int) xuint16;
		
	        if (TIFFGetField(tif,TIFFTAG_SAMPLESPERPIXEL, &xuint16) != 1) nbands = 1;
		else nbands = (int) xuint16;
		
                TIFFGetField(tif,TIFFTAG_BITSPERSAMPLE, &xuint16);
		bitspersample = (int) xuint16;
		
                TIFFGetField(tif,TIFFTAG_SAMPLEFORMAT,  &sampleformat);

                TIFFGetField(tif, TIFFTAG_PLANARCONFIG,  &config);

		
	        /* As of 27-MAY-10 (and probably much earlier),
		   the TIFF library can't handle reading tiepoints.
		   Unknown why.  		
		TIFFGetField(tif,TIFFTAG_GEOTIEPOINTS,  tiepoints);*/
		
		
	        /* I don't think it ever had much success reading pixscale values;
		   at least not in *this* program.		   
		TIFFGetField(tif,TIFFTAG_GEOPIXELSCALE, pixscale);*/
		
		/*printf("width = %d height = %d orientation = %d nbands = %d bitspersample = %d sampleformat = %d\n",
		    width, height, orientation, nbands, bitspersample, sampleformat);*/
	        /*for(i=0;i<24;i++) printf("tiepoint %d = %lf\n", i, tiepoints[i]);*/
	        /*for(i=0;i<3;i++) printf("pixscale %d = %f\n", i, pixscale[i]);*/

 		 /*  printf("width = %d height = %d orientation = %d nbands = %d bitspersample = %d sampleformat = %d config = %d\n",
		           width, height, orientation, nbands, bitspersample, sampleformat, config);*/
			   
               if (bitspersample == 8) {
                   *data8 = (unsigned char **)malloc(nbands*sizeof(unsigned char *));
		   for(i=0;i<nbands;i++) (*data8)[i] = (unsigned char *)malloc(width*height*sizeof(unsigned char));
                   ReadImage8(tif, width, height, nbands, data8, config);
		    }
                else if (bitspersample == 16) {
                   *data = (short **)malloc(nbands*sizeof(short *));
		   for(i=0;i<nbands;i++) (*data)[i] = (short *)malloc(width*height*sizeof(short));
                   ReadImage16(tif, width, height, nbands, data, config);
		    }
                else {
		   printf("Can't handle file of that format: \n");
		   printf("width = %d height = %d orientation = %d nbands = %d bitspersample = %d sampleformat = %d\n",
		           width, height, orientation, nbands, bitspersample, sampleformat);
		   printf("exiting...\n");
		   return(1);
	           }
		 
		(*latbuffer) = (float *)malloc(width*height*sizeof(float));
		(*lonbuffer) = (float *)malloc(width*height*sizeof(float));
		if (makelatlon == 1)
		   retval = make_lat_lon(gtif, tiepoints, width, height, latbuffer, lonbuffer);
		else 
		    retval = 1;
		if (retval == 0) {   /* success */
		   *latindex = nbands;
		   nbands += 2;
		   }
                else {
		   free(*latbuffer);
		   free(*lonbuffer);
		   *latindex = -1;
		   }

      		/* get rid of the key parser */
		GTIFFree(gtif);
		
		/* close the TIFF file descriptor */
		XTIFFClose(tif);
		
		*rwidth = width;
		*rheight = height;
		*rnbands = nbands;
		
		return (0);
	}


void get_next_tiffname(char name[200], char *tifftableofcontents, int length, int k)
{  /* assumes tifftableofcontents is a quoted, comma-separated list of SDS names, e. g., 
      "1km Band 1","1km Band 2","1km Band 3",
      It probably can work without the quotes, provided the names don't have any
      internal commas.
    */      
int i, n, nc, inside;
char x;
/*printf("%s\n", tifftableofcontents);*/

nc = n = inside = 0;
for(i=0;i<length;i++) {
    x = tifftableofcontents[i];
    if ( x == '"') {
       if (inside == 0) inside = 1;
       else if (inside == 1) inside = 0;
       continue;
       }
    else if (( x == ',') && (inside == 0)) nc++;	
    else {
        if (nc == k) name[n++] = x;
       }
   }
name[n] = '\0';   
return;
}


int is_a_tiffile(char *filename)
{
TIFF *tif=(TIFF*)0;  /* TIFF-level descriptor */
uint16 xuint16;
int ret, nbands;

tif=XTIFFOpen(filename,"r");  
if (!tif) return(-1);
else {
   ret = TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &xuint16);
   nbands = (int) xuint16;
   XTIFFClose(tif);
   if (ret != 1) nbands = 1;
  }
return(nbands);
}


void ReadImage16(TIFF *tif, int width, int height, int nbands, short **image[], int config)
{
	int i,j,k,kk,ii;
	unsigned char *buffer, bf1, bf2;
	short data;
	int band;
	
   if (config == 2) {
	buffer = (unsigned char *)malloc(width*2*sizeof(unsigned char));

 	for(band=0;band<nbands;band++) {
	   ii = 0;
	   
 	   for (i=0;i<height;i++) {
	
              if (!TIFFReadScanline(tif, buffer, i, band))
			TIFFError("ReadImage","failure in ReadScanline\n");
			
	      for (k=0;k<width;k++) {
	   
	         kk = k*2;
		 
		 bf1 = buffer[kk++];
		 bf2 = buffer[kk++];
	         data = (bf2<<8) + bf1;
		 (*image)[band][ii] = data;
		 /*printf("row %d col %d band %d val1 %d val2 %d value %d\n",
		     i, k, band, bf1, bf2, data);*/
	   
	         ii++;
	   
	          }
	       }  
	   }
  
      }
   else {   

	buffer = (unsigned char *)malloc(width*nbands*2*sizeof(unsigned char));
	ii = 0;
	for (i=0;i<height;i++) {
	
           if (!TIFFReadScanline(tif, buffer, i, 0))
			TIFFError("ReadImage","failure in ReadScanline\n");

	   for (k=0;k<width;k++) {
	   
	      kk = k*nbands*2;
	      
	      for(band=0;band<nbands;band++) {
	      		 
		 bf1 = buffer[kk++];
		 bf2 = buffer[kk++];
	         data = (bf2<<8) + bf1;
		 (*image)[band][ii] = data;
		 /*printf("row %d col %d band %d val1 %d val2 %d value %d\n",
		     i, k, band, bf1, bf2, data);*/
		 
	         }

	      ii++;

	      }
	   
         }
      }
         free(buffer);
return;
}

void ReadImage8(TIFF *tif, int width, int height, int nbands, unsigned char **image[], int config)
{
	int i,j,k,kk,ii;
	unsigned char *buffer;
	char data;
	int band;
	
    if (config == 2) {
	buffer = (unsigned char *)malloc(width*2*sizeof(unsigned char));

 	for(band=0;band<nbands;band++) {
	   ii = 0;
	   
 	   for (i=0;i<height;i++) {
	
              if (!TIFFReadScanline(tif, buffer, i, band))
			TIFFError("ReadImage","failure in ReadScanline\n");
			
	      for (k=0;k<width;k++) {
	   		 
	         data = buffer[k];
		 (*image)[band][ii] = data;
	   
	         ii++;
	   
	          }
	       }  
	   }
  
      }
    else {
  	buffer = (unsigned char *)malloc(width*nbands*2*sizeof(unsigned char));
	
	ii = 0;
	for (i=0;i<height;i++) {

           if (!TIFFReadScanline(tif, buffer, i, 0))
			TIFFError("ReadImage","failure in ReadScanline\n");


	   for (k=0;k<width;k++) {
	   
	      kk = k*nbands;
	      
	      for(band=0;band<nbands;band++) {
	      		 
		 data = buffer[kk++];
		 (*image)[band][ii] = data;
		 
	         }

	      ii++;

	      }
	   
         }
     }
         free(buffer);
return;
}

int make_lat_lon(GTIF *gtif, double tiepoints[24], int width, int height, float **latbuffer, float **lonbuffer)
{
int i, k, ii;
GTIFDefn defn;
double dlat, dlon;
int toPCSret, toLatLonret, ret = 0;
extern int VERBOSE, VVERBOSE;

/*for (i=0;i<24;i++) printf("%lf\n", tiepoints[i]);*/


GTIFGetDefn( gtif, &defn );
if (VERBOSE || VVERBOSE) printf("Calculating lat, lon...\n");
	ii = 0;
	for (i=0;i<height;i++) {
	
	   for (k=0;k<width;k++) {

               dlat = (double)i;
               dlon = (double)k;
               toPCSret =  GTIFImageToPCS( gtif, &dlon, &dlat );
               toLatLonret = GTIFProj4ToLatLong( &defn, 1, &dlon, &dlat );
	       
               ret += toLatLonret;
	       
	       (*latbuffer)[ii]   = dlat;
	       (*lonbuffer)[ii++] = dlon;
	   }
	}
	
if (VERBOSE || VVERBOSE) printf(" ...finished\n");
	
if (ret == height*width) 	
   return(0);
else
   return(1);   
}

