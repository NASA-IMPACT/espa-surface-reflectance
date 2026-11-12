#ifndef LOCAL_PROTOS_H
#define LOCAL_PROTOS_H

int read_tps( long fp,		/* File pointer */
	      struct TPSDATA *pt_data	/* Tie point data */
     );

int read_mtp( long fp,		/* File pointer */
	      struct MTPDATA *pt_data	/* Tie point data */
     );

int read_tpl( long fp,		/* File pointer */
	      struct TPLDATA *pt_data	/* Tie point data */
     );

int write_tps( long fp,		/* File pointer */
	       struct TPSDATA *pt_data	/* Tie point data */
     );

int write_mtp( long fp,		/* File pointer */
	       struct MTPDATA *pt_data	/* Tie point data */
     );

int write_tpl( long fp,		/* File pointer */
	       struct TPLDATA *pt_data	/* Tie point data */
     );

int ndcon( double *incoor,	/* input coordinates                    */
	   long *conv,		/* 1 = 27 - 83, 2 = 27 - 83             */
	   double *outcoor,	/* output coordinates                   */
	   char *errkey2,	/* error key                            */
	   char *message2	/* error message                        */
     );

int proj_err( long err );

void proj_type( char text[75]	/* Projection name */
     );

void proj_gen( double val, char *text );

void proj_dms( double angle, char *text );

void proj_spheroid( double major, double e2 );

void proj_rad( double radius );

void false_ne( double east, double north );

void proj_zone( long zone );

void proj_nad( double val );

void report_type( char text[75]	/* Projection name */
     );

void report_gen( double val, char *text );

void report_deg( double angle, char *text );

void report_spheroid( double major, double e2 );

void report_rad( double radius );

void report_zone( long zone );

void report_lf( void );

void spcs_zone( long zone	/* State Plane Coordinate System zone id */
     );

#endif
