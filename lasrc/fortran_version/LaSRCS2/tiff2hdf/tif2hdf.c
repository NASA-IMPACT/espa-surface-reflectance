#include <stdio.h>
#include <ctype.h>
#include <strings.h>
#include <string.h>
#include <stdlib.h>  
#include "mfhdf.h"
#include "math.h"

/* Compile with:

   gcc -O0 -o tiff2hdf *.c -I/usr/local/TOOLKIT/SDPToolkit_gfortran44-5.2.17/hdf/linux64/hdf-4.2.5/include/ /usr/lib64/libtiff.so \
      /usr/local/geotiff/lib/libgeotiff.a -L/usr/local/TOOLKIT/SDPToolkit_gfortran44-5.2.17/hdf/linux64/hdf-4.2.5/lib/ -lmfhdf \\
      -ldf -ljpeg -lz -lsz -lm -lproj
*/

int VERBOSE = 0;
int VVERBOSE = 0;
int makelatlon;
int main(int argc, char **argv)
{
FILE *fd;
int32 i, ii, j, k, n, sd, ret, n_sets, n_gattr, n_val, tif_rgb_val, allzero, n_OK, n_rez, dimid;
int latindex; 
int is_a_tiffile(char *filename);
int xdimone, xdimtwo, xn_sets, tiffret;
int tiff_read(char *filename, int *rwidth, int *rheight, int *rnbands, short **data[], unsigned char **data8[], 
              float **latbuffer, float **lonbuffer, int *latindex, char *tifftableofcontents); /* make lat lon construction the default */
int16 **TIFdata;
unsigned char **other_TIFdata;
char *tifftableofcontents;
float *lat_buffer, *lon_buffer;
int tiffTOCln;
void get_next_tiffname(char name[200], char *tifftableofcontents, int tiffTOCln, int k);

int32 sds, start[2], dims[2];
char name[200];

makelatlon = 0;

if (argc < 3) {
   printf("\n%s reads TIFF files, writes HDF files\n", argv[0]);
   printf("usage: %s <TIFF file> <HDF file> {-l}\n", argv[0]);
   printf("where {-l} (optional) is for calculating lat/lon for each pixel (very slow)\n");
   printf("Jim Ray, SSAI, %s\n\n", __DATE__);
   exit(0);
   }
if (argc > 3) {
   if ((argv[3][0] == '-') && (argv[3][1] == 'l')) makelatlon = 1;
   }
/*printf("%s %s %d\n", argv[1], argv[2], makelatlon);*/

ret = is_a_tiffile(argv[1]);
if (ret < 0) {
      printf("Error opening file %s: check file\n", argv[1]); 
      exit(-1);
  }
else if (ret == 0) {
      printf("File %s apparently has no bands: check file\n", argv[1]); 
      exit(-2);
  }
else {
      n_sets = ret;
  }
	     
sd = SDstart(argv[2], DFACC_CREATE);
	     
   xdimone = xdimtwo = 0;
   TIFdata = (int16 **)NULL;  
   other_TIFdata = (unsigned char **)NULL;  
   lat_buffer = (float *)NULL;  
   lon_buffer = (float *)NULL;  
   tifftableofcontents = (char *)malloc(256*60);  /* "256*60" comes from make_GEOTIFF_output(), in program hdf_to_TIFF() (13-MAY-08) */
   tifftableofcontents[0] = '\0';
   tiffret = tiff_read(argv[1], &xdimtwo, &xdimone, &xn_sets, &TIFdata, &other_TIFdata, &lat_buffer, &lon_buffer, &latindex, tifftableofcontents);
   /*printf("tifftableofcontents = %s\n", tifftableofcontents);*/
   if (tiffret != 0) exit(-2);
   
   if (makelatlon == 1) xn_sets -= 2;  /* xn_sets includes lat/lon, but we don't need that here */
   
   tiffTOCln = strlen(tifftableofcontents);
   /* Note; the TIFF library read of TIFFTAG_IMAGEDESCRIPTION can report 
      some pixels read when nothing is there... */
   if (tiffTOCln < 10) {
      tiffTOCln = 0; 
      for(k=0;k<xn_sets;k++) {
         
         sprintf(name, "BAND %d\0", k);
	 if (xn_sets == 3) {  /* assume R, G, B */
             switch (k) {
	         case (0) : strcat(name, " red\0"); break;
	         case (1) : strcat(name, " green\0"); break;
	         case (2) : strcat(name, " blue\0"); break;
	        }
	    }
	 
	 dims[0] = xdimone;
	 dims[1] = xdimtwo;
	 start[0] = start[1] = 0;
	 if (TIFdata != (int16 **)NULL)  {
	    sds = SDcreate(sd, name, DFNT_INT16, 2, dims);
	    if((dimid=SDgetdimid(sds, 0))!=-1) SDsetdimname(dimid, "YDim\0");
	    if((dimid=SDgetdimid(sds, 1))!=-1) SDsetdimname(dimid, "XDim\0");
	    SDwritedata(sds, start, NULL, dims, TIFdata[k]);
	    }
	 else {
	    sds = SDcreate(sd, name, DFNT_UINT8, 2, dims);
	    if((dimid=SDgetdimid(sds, 0))!=-1) SDsetdimname(dimid, "YDim\0");
	    if((dimid=SDgetdimid(sds, 1))!=-1) SDsetdimname(dimid, "XDim\0");
	    SDwritedata(sds, start, NULL, dims, other_TIFdata[k]);
	  }
	 SDendaccess(sds);
	 
         }
          
      }
   else {
      /*printf("tifftableofcontents = %d\n", tiffTOCln);*/
   
      /* v4.1: check all the names in the tifftableofcontents */
      i = 0;
      for(k=0;k<xn_sets;k++) {
         get_next_tiffname(name, tifftableofcontents, tiffTOCln, k);
	 
	 dims[0] = xdimone;
	 dims[1] = xdimtwo;
	 start[0] = start[1] = 0;
	 if (TIFdata != (int16 **)NULL)  {
	    sds = SDcreate(sd, name, DFNT_INT16, 2, dims);
	    if((dimid=SDgetdimid(sds, 0))!=-1) SDsetdimname(dimid, "YDim\0");
	    if((dimid=SDgetdimid(sds, 1))!=-1) SDsetdimname(dimid, "XDim\0");
	    SDwritedata(sds, start, NULL, dims, TIFdata[k]);
	    }
	 else {
	    sds = SDcreate(sd, name, DFNT_UINT8, 2, dims);
	    if((dimid=SDgetdimid(sds, 0))!=-1) SDsetdimname(dimid, "YDim\0");
	    if((dimid=SDgetdimid(sds, 1))!=-1) SDsetdimname(dimid, "XDim\0");
	    SDwritedata(sds, start, NULL, dims, other_TIFdata[k]);
	  }
	 SDendaccess(sds);
	 
         if (strlen(name) > 0) i++;
         }
      }
      if (makelatlon == 1) {
	    sds = SDcreate(sd, "latitude", DFNT_FLOAT, 2, dims);
	    if((dimid=SDgetdimid(sds, 0))!=-1) SDsetdimname(dimid, "YDim\0");
	    if((dimid=SDgetdimid(sds, 1))!=-1) SDsetdimname(dimid, "XDim\0");
	    SDwritedata(sds, start, NULL, dims, lat_buffer);
	    SDendaccess(sds);
	    sds = SDcreate(sd, "longitude", DFNT_FLOAT, 2, dims);
	    if((dimid=SDgetdimid(sds, 0))!=-1) SDsetdimname(dimid, "YDim\0");
	    if((dimid=SDgetdimid(sds, 1))!=-1) SDsetdimname(dimid, "XDim\0");
	    SDwritedata(sds, start, NULL, dims, lon_buffer);
	    SDendaccess(sds);
         
         }
  
SDend(sd);
  
  
} 
