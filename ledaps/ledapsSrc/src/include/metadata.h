/*****************************************************************************
NAME:                         METADATA.H 

PURPOSE:  Include file containing number of fields in standard metadata
	  format and an enumeration for the field names.

PROGRAM HISTORY:
VERSION         DATE    AUTHOR          DESCRIPTION
-------         ----    ------          -----------
  1.0           6/96    C. Engebretson  Original Development

******************************************************************************/


/* Number of fields in the standard metadata file.
--------------------------------------------------*/

#define FIELDS 68


/* Mnemonic names for each field.  An expanded discussion of the meanings of
   these fields can be found in the online LAS documentation.
----------------------------------------------------------------------------*/

enum MetaField { SCENE_ID,		/* New scene ID			*/
		 OLD_SCENE,		/* Old scene ID			*/
		 ORIG_SOURCE,		/* Recieving station		*/
		 PROJECT,		/* Project this scene is for	*/
		 OFFSET,		/* WRS offset from scene center */
		 USE_RESTR,		/* Flag for client restriction	*/
		 PATH_NBR,		/* WRS Path Number		*/
		 ROW_NBR,		/* WRS Row Number		*/
		 WRS_F_ROW,		/* First row (for multi-scene)	*/
		 WRS_L_ROW,		/* Last row (for multi-scene)	*/
		 WRS_TYPE,		/* WRS type			*/	
		 ACQU_DATE,		/* Date of acquisition		*/
		 CLOUD_COVER,		/* Percentage of cloud cover	*/
		 RECD_TECH,		/* Recording technique		*/
		 SUN_AZIMUTH,		/* Sun azimuth angle		*/
		 SUN_ELEVATION,		/* Sun elevation angle		*/
		 CENTER_LAT,		/* Scene center latitude	*/
		 CENTER_LON,		/* Scene center longitude	*/	
		 QUALITY,		/* Scene quality		*/
		 SATELLITE,		/* Satellite number		*/
		 MICRO_FRAME,
		 CONDITION,		/* Special image conditions	*/
		 RECORD_STAT,		
		 RESAMP_TECH,		/* Resampling technique		*/
		 USER_ORIENT,		/* User Orientation		*/
		 ORIEN_ANGLE,		/* Orientation Angle		*/
		 HORIZ_DATUM,		/* Horizontal Datum		*/
		 EDC_PROCESS,		/* EDC Processing Level		*/
		 BROW_AVAIL,		/* Browse image availability	*/
		 PROD_MEDIA,		/* Product Media		*/
		 MGMT_CODE,		/* Management Code		*/
		 MAP_PROJEC,		/* Map Projection		*/
		 IMS_VISIBL,		/* IMS visible flag		*/
		 INTERLVNG,		/* Band leaving on tape(s)	*/
		 PROD_SIZE,		/* Product Size			*/
		 PIXEL_WIDTH,
		 PIXEL_HGHT,
		 DEM_AVAIL,		/* DEM Availability Flag	*/
		 SPECT_BAND,		/* Spectral bands available	*/
		 RAD_E,			/* Radiometric Enhancement	*/
		 RAD_C,			/* Radiometric Calibration 	*/
		 UPP_L_LAT,		/* Corner latitude/longitude	*/
		 UPP_L_LON,
		 UPP_R_LAT,
		 UPP_R_LON,
		 LOW_R_LAT,
		 LOW_R_LON,
		 LOW_L_LAT,
		 LOW_L_LON,
		 MAP_ZONE,		/* USGS Map Zone		*/
		 PROJ_PAR1,		/* Projection Parameters	*/
		 PROJ_PAR2,
		 PROJ_PAR3,
		 PROJ_PAR4,
		 PROJ_PAR5,
		 PROJ_PAR6,
		 PROJ_PAR7,
		 PROJ_PAR8,
		 PROJ_PAR9,
		 PROJ_PAR10,
		 PROJ_PAR11,
		 PROJ_PAR12,
		 PROJ_PAR13,
		 PROJ_PAR14,
		 PROJ_PAR15,
		 ROOT_ERROR,		/* Root Mean Square Error	*/	
		 ELLIP_MAJ,		/* Semi-major axis		*/
		 ELLIP_MIN };		/* Semi-minor axis		*/
