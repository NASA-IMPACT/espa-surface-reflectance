#include <stdio.h>
#include <string.h>
#include "mfhdf.h"

/* Compile with:

   cc -O0 -c setprojinfo.c

 */
int set_proj_(int32 *sds, int *xd, int *yd, float *uplx, float *uply, float *lorx, float *lory, int *zone,int *sphere, float *wbc, float *ebc,
float *nbc , float *sbc ) {
char structstring[1000] =
 {"GROUP=SwathStructure\n"
    "END_GROUP=SwathStructure\n"
    "GROUP=GridStructure\n"
    "\tGROUP=GRID_1\n"
    "\t\tGridName=\"Grid\"\n"
    "\t\tXDim=%d\n"
    "\t\tYDim=%d\n"
    "\t\tUpperLeftPointMtrs=(%f,%f)\n"
    "\t\tLowerRightMtrs=(%f,%f)\n"
    "\t\tProjection=GCTP_UTM\n"
    "\t\tZoneCode=%d\n"
    "\t\tSphereCode=%d\n"
    "\tEND_GROUP=GRID_1\n"
    "END_GROUP=GridStructure\n"
    "GROUP=PointStructure\n"
    "END_GROUP=PointStructure\n"
    "END\n\0"};
char complete_structstring[1500];

float  westboundcoord= *wbc ;
float  eastboundcoord= *ebc ;
float  northboundcoord= *nbc ;
float southboundcoord = *sbc ;

sprintf(complete_structstring, structstring, *xd, *yd, *uplx, *uply, *lorx, *lory, *zone, *sphere);
/*                                                                xdim, ydim, UpperLeftPointMtrs, et cetera                               */

SDsetattr(*sds, "StructMetadata.0",  DFNT_CHAR, strlen(complete_structstring), (VOIDP)complete_structstring);
SDsetattr(*sds, "WestBoundingCoordinate", DFNT_FLOAT, 1, (VOIDP)&westboundcoord); /* if they're floats, else DFNT_DOUBLE */
SDsetattr(*sds, "EastBoundingCoordinate", DFNT_FLOAT, 1, (VOIDP)&eastboundcoord);
SDsetattr(*sds, "NorthBoundingCoordinate", DFNT_FLOAT, 1, (VOIDP)&northboundcoord);
SDsetattr(*sds, "SouthBoundingCoordinate", DFNT_FLOAT, 1, (VOIDP)&southboundcoord);

sprintf(complete_structstring, "UTMSR.gatesproject.hdf" ); /* whatever */
SDsetattr(*sds, "LocalGranuleID",  DFNT_CHAR, strlen(complete_structstring), (VOIDP)complete_structstring);
return(0);
}
