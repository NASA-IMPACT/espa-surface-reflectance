#ifndef LTWG_H
#define LTWG_H

/*****************************************************************************
NAME:                           LTWG.H

PURPOSE:  Include file defining structure for LTWG files.

PROGRAM HISTORY:
VERSION         DATE    AUTHOR          DESCRIPTION
-------         ----    ------          -----------
  6.0            3/94   D. Bridges      Original Development

******************************************************************************/
#include <sysdef.h>
#include <typlim.h>
#include <worgen.h>

struct LTWGVDR
  {
  /* VOLUME DESCRIPTOR RECORD 
  ---------------------------*/
  long   rec_number;		/* Record Sequence Number 		*/
  long   frst_rec_code;		/* First record sub-type code  		*/
  long   rec_code;		/* Record type code 			*/
  long   scnd_rec_code;		/* Second record sub-type code		*/
  long   thrd_rec_code;		/* Third record sub-type code		*/
  long 	 rec_length;		/* Length of this record		*/
  char   a_e_flag[3];		/* ASCII/EBCDIC flag			*/
  char   control_num[13];	/* Superstructure control doc number 	*/
  char   control_rev_num[3];	/* Superstruct contrl doc rev number 	*/
  char   format_rev[3];		/* Superstruct record format rev letter	*/
  char   release_num[13];	/* Software release number		*/
  char   tape_id[17];		/* Tape ID for physical volume		*/
  char   logical_vol_id[17];	/* Logical volume ID			*/
  char   vol_set_id[17];	/* Volume set ID			*/
  long   num_phy_vol;		/* Number of pyhsical volumes in set	*/
  long   start_log_vol;		/* Phys vol num, start of logical vol	*/
  long   end_log_vol;		/* Phys vol num, end of logical vol	*/
  long   volume_num;		/* Phys vol num containing vol discript */
  long   frst_ref_file;		/* First referenced file num in phy vol */
  long   log_vol_num_set;	/* Logical vol num within volume set 	*/
  long   log_vol_num_phys;	/* Logical vol num within Physical vol	*/
  char   log_creat_date[9];	/* Logical volume creation date		*/
  char   log_creat_time[9];	/* Logical volume creation time 	*/
  char 	 log_gen_country[13];	/* Logical volume generating country	*/
  char   log_gen_agency[9];	/* Logical volume generating agency	*/
  char   log_gen_facility[13];	/* Logical volume generating facility	*/
  long   num_point_rec;		/* number pointer records in vol direct */
  long   num_records;		/* number of records in volume directory*/
  };

struct LTWGFPR
  {
  /* FILE POINTER RECORD 
  ----------------------*/
  long   rec_number;		/* Record Sequence Number 		*/
  long   frst_rec_code;		/* First record sub-type code  		*/
  long   rec_code;		/* Record type code 			*/
  long   scnd_rec_code;		/* Second record sub-type code		*/
  long   thrd_rec_code;		/* Third record sub-type code		*/
  long	 rec_length;		/* Length of this record		*/
  char   a_e_flag[3];		/* ASCII/EBCDIC flag			*/
  long   num;			/* Leader file number			*/
  char   name[17];		/* Leader file name			*/
  char   class[29]; 		/* Leader file class			*/
  char   cl_code[5];		/* Leader file class code		*/
  char   datatyp[29];		/* Leader file data type		*/
  char   dt_code[5];		/* Leader file data type code		*/
  long   num_rec;		/* Number of files in leader file 	*/
  long   des_rec_lth; 		/* Leader file descriptor record length */
  long   max_lth;		/* Leader file maximum record length	*/
  char   rec_lth_typ[13];	/* Leader file record length type	*/
  char   rec_lth_code[15];	/* Leader file record length type code	*/
  long   vol_num_start;		/* Leader file physical vol num, start	*/
  long   vol_num_end;		/* Leader file physical vol num, end	*/
  long   file_por_frst;		/* Leader file portion, first rec num	*/
  long   file_por_last;		/* Leader file portion, last rec num	*/
  };

struct LTWGTRIV
  {
  /* VOLUME DIRECTORY TEXT RECORD - IMAGERY VOLUME 
  ------------------------------------------------*/
  long   rec_number;		/* Record Sequence Number 		*/
  long   frst_rec_code;		/* First record sub-type code  		*/
  long   rec_code;		/* Record type code 			*/
  long   scnd_rec_code;		/* Second record sub-type code		*/
  long   thrd_rec_code;		/* Third record sub-type code		*/
  long	 rec_length;		/* Length of this record		*/
  char   a_e_flag[3];		/* ASCII/EBCDIC flag			*/
  char   cont_flag[3];		/* Continuation flag			*/
  char   product_type[51];	/* Product type				*/
  char   loc_and_date[59];	/* Location, date &time of prod creation*/
  char   in_scene_id[50];	/* Input scene identification		*/
  char   tape_id[44];		/* Physical tape identification		*/
  char   wrs_id[29];		/* WRS & product identification		*/
  char   level_cor[24];		/* level of processing corrections appld*/
  };

struct LTWGLFD
  {
  /*  LEADER FILE - FILE DESCRIPTOR RECORD
  ----------------------------------------*/

  long   rec_number;		/* Record Sequence Number 		*/
  long   frst_rec_code;		/* First record sub-type code  		*/
  long   rec_code;		/* Record type code 			*/
  long   scnd_rec_code;		/* Second record sub-type code		*/
  long   thrd_rec_code; 	/* Third record sub-type code		*/
  long	 rec_length;		/* Length of this record		*/
  char   a_e_flag[3];		/* ASCII/EBCDIC flag			*/
  char   control_doc_num[13];	/* Control document number for data file*/
  char   control_doc_rev[3];	/* Control document revision number	*/
  char   file_des_rev[3];	/* File design descriptor revision num	*/
  char   rel_num[13];		/* Software release Number		*/
  long   file_num;		/* File number				*/
  char   file_name[16];		/* File name				*/
  char   record_seq[5];		/* Record sequence and location type flg*/
  long   seq_num_loc;   	/* Sequence number location		*/
  long	 seq_num_lgt;		/* Sequence number field length		*/
  char   rec_loc_flg[5];	/* Record code and location type flag	*/
  long   rec_loc;   		/* Record code location			*/
  long 	 rec_lth;		/* Record code field length		*/
  char   rec_lth_flg[5];	/* Record length and location type flag	*/
  long   rec_lth_loc;   	/* Record length location		*/
  long	 rec_fld_lth;		/* Record length field length		*/
  char   data_interp_des[2];	/* Flag if data interpretation included */
				/* In descriptor record			*/
  char   data_interp_oth[2];	/* Flag if data interp inc in other file*/
  char   data_disp_des[2];	/* Flag if data display information inc	*/
				/* In file discriptor record		*/
  char   data_disp_oth[2];	/* Flag if data disp inc in other file  */
  long   num_header;   		/* Number of header records		*/
  long   lth_header;		/* Header record length			*/
  long   num_map_proj;   	/* Number of map projection ancillary	*/
  long   lth_map_proj;   	/* Map project ancillary record length	*/
  long   num_ancill_rec;	/* Number of radiometric capibrat anc rec*/
  long   lth_ancill_rec;	/* Radiometric ancillary record length	*/
  char   scene_id_loc[17];	/* Scene identification field locator	*/
  char   wrs_id_loc[17];	/* WRS identification field locator	*/
  char   miss_id_loc[17];	/* mission identification field locator */
  char   sensor_id_loc[17];	/* Sensor identification field locator 	*/
  char   exp_date_loc[17];	/* Exposre date-time field locator	*/
  char   geog_ref_loc[17];	/* Geographic reference field locator	*/
  char   imag_proc_loc[17];	/* Image processing performed fld locat */
  char   imag_for_loc[17];	/* Imagery format indicator locator	*/
  char   band_ind_loc[17];	/* Band indicator locator		*/
  char   quad_ind_loc[17];	/* quadrant indicator locator		*/
  char   inter_scale_loc[17];	/* Inter-pixel &inter-line scale locator*/
  char   quad_ver_ovrlp[17];	/* Quadrant product vert overlap indicat*/
  char   quad_hor_ovrlp[17];    /* Quadrand product horz overlap indicat*/
  };

struct LTWGSHR
  {
  /*  LEADER FILE - SCENE HEADER RECORD  
  -------------------------------------*/
      /* SCENE PARAMETERS */
  long   rec_number;		/* Record Sequence Number 		*/
  long   frst_rec_code;		/* First record sub-type code  		*/
  long   rec_code;		/* Record type code 			*/
  long   scnd_rec_code;		/* Second record sub-type code		*/
  long   thrd_rec_code; 	/* Third record sub-type code		*/
  long	 rec_length;		/* Length of this record		*/
  long   header_seq_num;	/* Header record sequence number	*/
  char   prod_type[17];		/* Product type identification		*/
  char   in_scene_id[17];	/* Input scene identification		*/
  double in_scene_centre_lat;	/* Input scene centre latitude in degrees*/
  double in_scene_centre_lon;	/* Input scene centre longitude in degrees*/
  double in_line_num_cent;	/* Line number at input scene centre	*/
  double in_pixel_num_cent;	/* Pixel number at input scene centre	*/
  char   time_centre[33];	/* Input scene centre time		*/
  float  time_offset;		/* Time offset form WRS frame (msecs)	*/
  char   wrs_design[17];	/* WRS designator (path and row)	*/
  char   wrs_cycle[17];		/* WRS cycle				*/
  char   pro_scene_id[17];	/* Processed scene identification	*/
  double pro_scene_centre_lat;	/* Proc scene centre latitude in degrees*/
  double pro_scene_centre_lon;	/* Proc scene centre longitude in degrees*/
  double pro_line_num_cent;	/* Line number at procsd scene centre	*/
  double pro_pixel_num_cent;	/* Pixel number at procsd scene centre	*/
  long   ovrlap_lines;		/* Count of overlap lines		*/
  long   ovrlap_pixels;		/* Count of overlap pixels		*/
      /* MISSION PARAMETERS */
  char   miss_id[17];		/* Mission identification		*/
  char   sensor_id[17];		/* Sensor identification		*/
  long   orbit_num;		/* Orbit number				*/
  char   asc_des_flag[17];	/* Ascending/descending flag		*/
      /* SENSOR PARAMETERS */
  long   wavelgth_range[64][2]; /* Upper and lower wavelgth range (nm)	*/
  long   num_bands;		/* Number of active bands in procsd image*/
  long   pixels_per_line;	/* Number of pixels per line in procsd im*/
  long   num_lines;		/* Number of scene lines in image	*/
      /* PROCESSING PARAMETERS */
  char   rad_cal_des[17];	/* Radiometric calibration designator 	*/
  long   rad_res_des;            /* Radiometric resoultion designator    */
  char   scenic_rad_des[17];	/* Scenic radiometric correction desgntr*/
  char   geom_corr_des[17];	/* geometric correction designator	*/
  char   resamp_des[17];	/* Resampling designator		*/
  char   map_proj_id[17];	/* Map projection identifier		*/
  long   process_level;		/* Product processing level		*/
  long   num_map_proj_anc;	/* Number of map projection ancillary rec*/
  char   replac_detect[17];	/* Technique to replace failed detectors */
  char   kernal_for_replac[17];	/* Kernal used to replace failed detector*/
  long   num_rad_ancill;	/* Number of radiometric ancillary record*/
  char   active_bands[65];	/* Number of active bands		*/
  char   interleave_ind[17];	/* interleaving idicator		*/
  long   detector_subst[100];	/* Detector substitution array		*/
  char   detector_smooth[101];	/* Detector smoothing array		*/
  double mirror_scan_vel[12];	/* Mirror scan velocity profile coeff	*/
  long   dect_adjust[4][16];	/* Detector adjustment array		*/
  };

struct LTWGMPR
  {
  /*  LEADER FILE - MAP PROJECTION ANCILLARY RECORD 
  -------------------------------------------------*/
  long   rec_number;		/* Record Sequence Number 		*/
  long   frst_rec_code;		/* First record sub-type code  		*/
  long   rec_code;		/* Record type code 			*/
  long   scnd_rec_code;		/* Second record sub-type code		*/
  long   thrd_rec_code; 	/* Third record sub-type code		*/
  long	 rec_length;		/* Length of this record		*/
      /*  Inut scene related data  */
  long   nominal_num_pixels;	/* Nominal number of scene data pixels	*/
  long   nominal_num_lines;	/* Nominal number of scene data lines	*/
  double nominal_scale_pixel;	/* Nominal scale of inter-pixel dist (m)*/
  double nominal_scale_lines;	/* Nominal scale of inter-line dist (m) */
  double image_skew;		/* Image skew at centre from rotation(deg)*/
      /*  UTM/PS Related data for input scene  */
  char   utm_datum[7]; 		/* UTM datum for input image WRS centre	*/
  long   utm_zone;		/* UTM zone number for input image 	*/
  double wrs_northing;		/* Northing of WRS centre in meters 	*/
  double wrs_easting;		/* Easting of WRS centre in meters	*/
  double input_northing;	/* Northing of input image centre-meters*/
  double input_easting;		/* Easting of input image centre-meters */
  double vert_offset;		/* Verticlal offset of centre to WRS nominal*/
  double horz_offset;		/* Horizontal offset of center to WRS 	*/
  double orientation;		/* Orient of input image center in degrees*/
				/* relative to the UTM grid		*/
      /*  Processed scene related data */
  double prod_pixel_line;	/* Product pixels per line		*/
  double prod_num_lines;	/* Number of lines per processed image	*/
  double inter_pixel_scale;	/* Scale of processed inter-pixel dist(m)*/
  double inter_line_scale;	/* Scale of processed inter-line dist(m)*/
  char   proc_utm_datum[7];	/* UTM datum for processed image	*/
  long	 proc_utm_zome;		/* UTM zone number for processed image	*/
  double proc_center_line;	/* Line num of WRS centre in processed im*/
  double proc_center_pix;	/* Pixel num of WRS cent. in processed im*/
  double proc_orient;		/* Orientation of proc. image centre (deg)*/
      /*  Sensor data  */
  double inclination;		/* Sat. orbital inclin from polar orbit */
  double ascend_node;		/* Spacecraft ascending node		*/
  double altitude;		/* Satellite altit. at WRS centre (m)	*/
  double ground_speed;		/* Sat grnd speed at WRS centre (m/s)	*/
  double sat_heading;		/* Satellite heading (degrees)		*/
  double cross_track_fov;	/* Cross-track field of view (degrees)	*/
  double scan_rate;		/* Sensor scan rate (scans/sec)		*/
  double samp_rate;		/* Sensor sampling rate (samples/sec)	*/
  double sun_elev;		/* Sun elevation angle at WRS centre(deg)*/
  double sun_azimuth;		/* Sun azimuth angle at WRS centre(deg)	*/
  double utm_corner_1[2];	/* UTM coord for corner 1 (meters)	*/
  double utm_corner_2[2];	/* UTM coord for corner 2 (meters)	*/
  double utm_corner_3[2];	/* UTM coord for corner 3 (meters)	*/
  double utm_corner_4[2];	/* UTM coord for corner 4 (meters)	*/
  double lat_long_corner1[2];	/* Lat Long of corner 1 (degrees)	*/
  double lat_long_corner2[2];	/* Lat Long of corner 2 (degrees)	*/
  double lat_long_corner3[2];	/* Lat Long of corner 3 (degrees)	*/
  double lat_long_corner4[2];	/* Lat Long of corner 4 (degrees)	*/
  double pix_line_corner1[2];   /* Coord of line expressed in input pixels */
				/* and lines of corner 1 		*/
  double pix_line_corner2[2];   /* Coord of line expressed in input pixels */
				/* and lines of corner 2		*/
  double pix_line_corner3[2];   /* Coord of line expressed in input pixels */
				/* and lines of corner 3		*/
  double pix_line_corner4[2];   /* Coord of line expressed in input pixels */
				/* and lines of corner 4		*/
  };

struct LTWGRCR
  {
  /*  RADIOMETRIC CALIBRATION ANCILLARY RECORD 
  --------------------------------------------*/
  long   rec_number;		/* Record Sequence Number 		*/
  long   frst_rec_code;		/* First record sub-type code  		*/
  long   rec_code;		/* Record type code 			*/
  long   scnd_rec_code;		/* Second record sub-type code		*/
  long   thrd_rec_code; 	/* Third record sub-type code		*/
  long	 rec_length;		/* Length of this record		*/
  long   band_number;		/* Band number				*/
  long   reflectance_lower;	/* Lower reflectance limit in percentage*/
  long   reflectance_upper;	/* upper reflectance limit in percentage*/
  double a0;			/* Offset coefficient			*/
  double a1;			/* Gain coefficient			*/
  long   detector_data[16][256]; /* data for detectors			*/
  };

struct LTWGIFD
  {
  /*   IMAGERY FILE - FILE DESCRIPTOR RECORD
  ------------------------------------------*/
  long   rec_number;		/* Record Sequence Number 		*/
  long   frst_rec_code;		/* First record sub-type code  		*/
  long   rec_code;		/* Record type code 			*/
  long   scnd_rec_code;		/* Second record sub-type code		*/
  long   thrd_rec_code; 	/* Third record sub-type code		*/
  long	 rec_length;		/* Length of this record		*/
  char   a_e_flag[3];		/* ASCII/EBCDIC flag			*/
  char   control_doc_num[13];	/* Control document number for data file*/
  char   control_doc_rev[3];	/* Control document revision number	*/
  char   file_des_rev[3];	/* File design descriptor revision num	*/
  char   rel_num[13];		/* Software release Number		*/
  long   file_num;		/* File number				*/
  char   file_name[16];		/* File name				*/
  char   record_seq[5];		/* Record sequence and location type flg*/
  long   seq_num_loc;   	/* Sequence number location		*/
  long	 seq_num_lgt;		/* Sequence number field length		*/
  char   rec_loc_flg[5];	/* Record code and location type flag	*/
  long   rec_loc;   		/* Record code location			*/
  long 	 rec_lth;		/* Record code field length		*/
  char   rec_lth_flg[5];	/* Record length and location type flag	*/
  long   rec_lth_loc;   	/* Record length location		*/
  long	 rec_fld_lth;		/* Record length field length		*/
  char   data_interp_des[2];	/* Flag if data interpretation included */
				/* In descriptor record			*/
  char   data_interp_oth[2];	/* Flag if data interp inc in other file*/
  char   data_disp_des[2];	/* Flag if data display information inc	*/
				/* In file discriptor record		*/
  char   data_disp_oth[2];	/* Flag if data disp inc in other file  */
  long   number_img_rec;	/* Number of image records		*/
  long   record_lth;		/* Record length			*/
  long   bits_per_pixel;	/* Number of bits per pixel		*/
  long   pixel_per_group;	/* Number of pixel/data group		*/
  long   bytes_per_group;	/* Number of bytes/data group		*/
  char   pix_just[5];		/* Pixel justification - RJLR = pixels are*/
				/* right justified with pirst pixel leftmost*/
  long   num_bands;		/* number of images(bands) this file	*/
  long   lines_band;		/* lines per band			*/
  long   left_pix_line;		/* Left border pixel/line		*/
  long   img_pix_line;		/* Image pixels per line		*/
  long   right_pix_line;	/* Right border pixel/line		*/
  long   top_lines;		/* Number of top border lines		*/
  long   bot_lines;		/* Number of bottom border lines	*/
  char   file_type[5];		/* Type of file BSQ or BIL		*/
  long   records_line;		/* Number of physical records per line	*/
  long   recs_multispec_line;	/* Physical records/multispectrial line	*/
  long   size_prefix;		/* Bytes of prefix data per logical rec	*/
  long   size_image;		/* Bytes of image data per logical rec	*/
  long   size_suffix;		/* Bytes of suffix data per logical rec	*/
  char   line_num_loc[9];	/* Scan line number locator		*/
  char   band_num_loc[9];	/* Banc number locator			*/
  char   location_time[9];	/* Scan location time			*/
  char   left_fill_loc[9];	/* Left fill count locator		*/
  char   right_fill_loc[9];	/* Right fill count locator		*/
  char   unknown1[9];		/* Do not know what field represents	*/
  char   unknown2[9];		/* Do not know what field represents	*/
  long   left_fill;             /* Left fill/pixel			*/
  long   right_fill;		/* Right fill/pixel			*/
  char   unknown3[9];		/* Do not know what field represents	*/
  }; 

struct LTWGIDR
  {
  /*  IMAGERY FILE _ IMAGE DATA RECORD
  ------------------------------------*/
     /* PREFIX DATA */
  long   scan_line;		/* Scan line number			*/
  long   band_num;		/* Logical band number			*/
  long   start_time;		/* Time in GMT at start of scan (millsec)*/
  long   left_fill;		/* Count of left fill pixels		*/
  long   right_fill;		/* Count of right fill pixels		*/
     /* SUFFIX DATA */
  long   sync_loss;		/* sync loss indicator 1-yes 0-no	*/
  long   local_qual_flag1[2];	/* Local use quality code		*/
  long   local_qual_flag2[4];	/* Local use quality code		*/
  long   detect_sub;		/* Detector substitution ind 1-yes 0-no */
  long   cal_pulse_width;	/* Calibration pulse width		*/
  long   line_lgth;		/* Counted full-scan line length	*/
  long   embedded_lgth;		/* Embedded line length			*/
  long   time_err_start;	/* Time error - line start to midscan	*/
  long   time_err_end;		/* Time error - midscan to line end	*/
  long   scan_direction;	/* scan line direction 0-forward 1 reverse*/
  long   half_line_lngth;	/* current half-scan line length	*/
  char   sat_time_code[17];	/* Satellite time code, start of scan	*/
     /* CALIBRATION INFORMATION */
  long   detect_id;		/* Detector identification		*/
  long   calib_lamp_qual;	/* Calibration lamp value quality	*/
  long   calib_lamp_state;	/* Calibration lamp state		*/
  long   calib_state_seq;	/* Calibration state sequence number	*/
  long   low_calib_value_bef;	/* Low level calib val (thousandths- levels)*/
				/* before dc restore			*/
  long   high_calib_value;	/* High level calib val (millionths  units)*/ 
  long   comput_gain;		/* Computed gain val (millionths of units)*/
  long   comput_bais;		/* Computed bais val (millionths of units)*/
  long   app_gain;		/* Applied gain val (millionths of units)*/
  long   app_bais;		/* Applied bais val (millionths of units)*/
  long   low_calib_value_aft;	/* Low level calib val (thousandths- levels)*/
				/* after dc restore			*/
  
  };

struct LTWGTFD
  {
  /*   TRAILER FILE - FILE DESCRIPTOR 
  -----------------------------------*/
  long   rec_number;		/* Record Sequence Number 		*/
  long   frst_rec_code;		/* First record sub-type code  		*/
  long   rec_code;		/* Record type code 			*/
  long   scnd_rec_code;		/* Second record sub-type code		*/
  long   thrd_rec_code; 	/* Third record sub-type code		*/
  long	 rec_length;		/* Length of this record		*/
  char   a_e_flag[3];		/* ASCII/EBCDIC flag			*/
  char   control_doc_num[13];	/* Control document number for data file*/
  char   control_doc_rev[3];	/* Control document revision number	*/
  char   file_des_rev[3];	/* File design descriptor revision num	*/
  char   rel_num[13];		/* Software release Number		*/
  long   file_num;		/* File number				*/
  char   file_name[16];		/* File name				*/
  char   record_seq[5];		/* Record sequence and location type flg*/
  long   seq_num_loc;   	/* Sequence number location		*/
  long	 seq_num_lgt;		/* Sequence number field length		*/
  char   rec_loc_flg[5];	/* Record code and location type flag	*/
  long   rec_loc;   		/* Record code location			*/
  long 	 rec_lth;		/* Record code field length		*/
  char   rec_lth_flg[5];	/* Record length and location type flag	*/
  long   rec_lth_loc;   	/* Record length location		*/
  long	 rec_fld_lth;		/* Record length field length		*/
  char   data_interp_des[2];	/* Flag if data interpretation included */
				/* In descriptor record			*/
  char   data_interp_oth[2];	/* Flag if data interp inc in other file*/
  char   data_disp_des[2];	/* Flag if data display information inc	*/
				/* In file discriptor record		*/
  char   data_disp_oth[2];	/* Flag if data disp inc in other file  */
  long   num_trailer;		/* Number of trailer records		*/
  long	 lth_trailer;		/* Trailer record length		*/
  char   parity_count_loc[17];	/* Parity error count field locator	*/
  char   quality_code_loc[17];	/* Quality code summary map field locator*/
  };

struct LTWGTR
  {
  /*  TRAILER RECORDS
  -------------------*/
  long   rec_number;            /* Record Sequence Number               */
  long   frst_rec_code;         /* First record sub-type code           */
  long   rec_code;              /* Record type code                     */
  long   scnd_rec_code;         /* Second record sub-type code          */
  long   thrd_rec_code;         /* Third record sub-type code           */
  long   rec_length;            /* Length of this record                */
  long   record_seq_num;	/* Trailer record sequence number	*/
  long   seq_in_band;		/* Sequence number of trailer rec w/in bnd*/
  long   detect_histo[4][256];	/* Histogram of RAW data for 4 detectors*/
				/* for the band				*/
  long   parity_count;		/* Parity error count			*/
  char   quality_summary[201];	/* Quality summary 			*/
  };

struct LTWGSTR
  {     
  /*  SUPPLEMENTAL VOLUME TEXT RECORD 
  -----------------------------------*/ 
  long   rec_number;		/* Record Sequence Number 		*/
  long   frst_rec_code;		/* First record sub-type code  		*/
  long   rec_code;		/* Record type code 			*/
  long   scnd_rec_code;		/* Second record sub-type code		*/
  long   thrd_rec_code;		/* Third record sub-type code		*/
  long	 rec_length;		/* Length of this record		*/
  char   a_e_flag[3];		/* ASCII/EBCDIC flag			*/
  char   continuation_flg[3];	/* Continuation flag			*/
  char   supp_file_descr[41];	/* Supplemental file description	*/
  char   location_time[59];	/* Location and date/time of prod creation*/
  char   orbit_id[50];		/* Orbit identification			*/
  char   tape_id[44];		/* Physical tape identification		*/
  };

struct LTWGSFD
  {
  /*   SUPPLEMENTAL FILE - FILE DESCRIPTOR RECORD 
  -----------------------------------------------*/
  long   rec_number;		/* Record Sequence Number 		*/
  long   frst_rec_code;		/* First record sub-type code  		*/
  long   rec_code;		/* Record type code 			*/
  long   scnd_rec_code;		/* Second record sub-type code		*/
  long   thrd_rec_code; 	/* Third record sub-type code		*/
  long	 rec_length;		/* Length of this record		*/
  char   a_e_flag[3];		/* ASCII/EBCDIC flag			*/
  char   control_doc_num[13];	/* Control document number for data file*/
  char   control_doc_rev[3];	/* Control document revision number	*/
  char   file_des_rev[3];	/* File design descriptor revision num	*/
  char   rel_num[13];		/* Software release Number		*/
  long   file_num;		/* File number				*/
  char   file_name[16];		/* File name				*/
  char   record_seq[5];		/* Record sequence and location type flg*/
  long   seq_num_loc;   	/* Sequence number location		*/
  long	 seq_num_lgt;		/* Sequence number field length		*/
  char   rec_loc_flg[5];	/* Record code and location type flag	*/
  long   rec_loc;   		/* Record code location			*/
  long 	 rec_lth;		/* Record code field length		*/
  char   rec_lth_flg[5];	/* Record length and location type flag	*/
  long   rec_lth_loc;   	/* Record length location		*/
  long	 rec_fld_lth;		/* Record length field length		*/
  char   data_interp_des[2];	/* Flag if data interpretation included */
				/* In descriptor record			*/
  char   data_interp_oth[2];	/* Flag if data interp inc in other file*/
  char   data_disp_des[2];	/* Flag if data display information inc	*/
				/* In file discriptor record		*/
  char   data_disp_oth[2];	/* Flag if data disp inc in other file  */
  long   num_header_rec;	/* Number of interval-related header rec*/
  long   lth_header_rec;	/* Interval-releated header length	*/
  long   num_tm_house;		/* Number of TM housekeeping records	*/
  long   lth_tm_house;		/* Length of TM housekeeping records	*/
  long   num_ephemeris;		/* Number of ephemeris data records	*/
  long	 lth_ephemeris;		/* Length of ephemeris data records 	*/
  long   num_scene_def;		/* Number of scene definition		*/
  long   lth_scene_def;		/* Length of scene definition		*/
  long   num_scene_qual;	/* Number of scene quality data records */
  long   lth_scene_qual;	/* Length of scene quality data records	*/
  long   num_geometric_mod;	/* Number of geometric modelling records*/
  long   lth_geometric_mod;	/* Length of geometric modelling records*/
  long   num_sparse_matrix;	/* Number of sparse matrices records	*/
  long   lth_sparse_matrix;	/* Length of sparse matrices records	*/
  long   num_gcd_start_time;	/* Number of GCD mirror scan start time rec*/
  long   lth_gcd_start_time;	/* Length of GCD mirror scan start time rec*/
  long   num_across_scan;	/* Number of high freq across-scan lines*/
	 			/* matrix records(raw fitter meas)	*/
  long   lth_across_scan;	/* Length of high freq across-scan lines*/
				/* matrix records(raw fitter meas)	*/
  long   num_cross_scan;	/* Number of high freq cross-scan lines */
	 			/* matrix records                       */
  long   lth_cross_scan;        /* Length of high freq cross-scan lines */
				/* matrix records                       */
  long   num_mission_rec;	/* Number of mission telemetry ancillary rec*/
  long   lth_mission_rec;	/* Length of mission telemetry ancillary rec*/
  long   num_grnd_cont_rec;	/* Number of ground control point records*/
  long   lth_grnd_cont_rec;	/* Length of ground control point records*/
  char   start_time_loc[17];	/* Interval data start time locator	*/
  char   stop_time_loc[17];	/* Interval data stop time locator	*/
  char   orbit_fld_loc[17];	/* Orbit field locator			*/
  };

struct LTWGSIH
  {
  /*   SUPPLEMENTAL FILE - INTERVAL HEADER FILE
  ---------------------------------------------*/
  long   rec_number;		/* Record Sequence Number 		*/
  long   frst_rec_code;		/* First record sub-type code  		*/
  long   rec_code;		/* Record type code 			*/
  long   scnd_rec_code;		/* Second record sub-type code		*/
  long   thrd_rec_code;		/* Third record sub-type code		*/
  long	 rec_length;		/* Length of this record		*/
  long   header_num;		/* Interval header record sequence num	*/
  char   input_start_time[17];	/* Input scene start time		*/
  char   input_stop_time[17];	/* Input scene stop time		*/
  char   pcd_start_time[17];	/* PCD telemetry start time		*/
  char   pcd_stop_time[17];	/* PCD telemetry stop time		*/
  long   pcd_frames;		/* Number of PCD major frames		*/
  long	 orbit_num;		/* Orbit number				*/
  };

struct LTWGTMHK
  {
  /*   SUPPLEMENTAL FILE - TM HOUSEKEEPING DATA ANCILLARY RECORDS 
  ---------------------------------------------------------------*/
  long   rec_number;		/* Record Sequence Number 		*/
  long   frst_rec_code;		/* First record sub-type code  		*/
  long   rec_code;		/* Record type code 			*/
  long   scnd_rec_code;		/* Second record sub-type code		*/
  long   thrd_rec_code;		/* Third record sub-type code		*/
  long	 rec_length;		/* Length of this record		*/
  long   rec_num;		/* TM housekeeping data record sequence num*/
  long   pcd_ident;		/* PCD major frame identifier		*/
      /*   PCD MAJOR FRAME START TIME    */
  long   pcd_start_day;		/* Day of year at start of PCD frame	*/
  long   pcd_start_msec;	/* Milliseconds of day at start of PCD frame*/
  long   pcd_start_usec;	/* Microseconds at start of PCD frame   */
  long   pcd_start_tenths;	/* Tenths of a microsecond at start of PCD*/
      /*   TM HOUSEKEEPING DATA    */
  /* Blackbody temperature sample and conversion to Degrees C */
  long   bb_time_offset;	/* Blackbody samp time offset (microsec)*/
  long	 bb_hk_sample;		/* Blackbody TM housekeeping sample	*/
  double bb_A0;			/* Blackbody A0 coefficient		*/
  double bb_A1;			/* Blackbody A1 coefficient		*/
  double bb_A2;			/* Blackbody A2 coefficient		*/
  double bb_A3;			/* Blackbody A3 coefficient		*/
  double bb_A4;			/* Blackbody A4 coefficient		*/
  double bb_A5;			/* Blackbody A5 coefficient		*/	
  /* Focal plane temperature sample and conversion to Degrees C */
  long   ff_time_offset;	/* Focal plane samp time offset (microsec)*/
  long	 ff_hk_sample;		/* Focal plane TM housekeeping sample	*/
  double ff_A0;			/* Focal plane A0 coefficient		*/
  double ff_A1;			/* Focal plane A1 coefficient		*/
  double ff_A2;			/* Focal plane A2 coefficient		*/
  double ff_A3;			/* Focal plane A3 coefficient		*/
  double ff_A4;			/* Focal plane A4 coefficient		*/
  double ff_A5;			/* Focal plane A5 coefficient		*/	
  /* Calib shutter temperature sample and conversion to Degrees C */
  long   st_time_offset;	/* Calib shutter samp time offset (microsec)*/
  long	 st_hk_sample;		/* Calib shutter TM housekeeping sample	*/
  double st_A0;			/* Calib shutter A0 coefficient		*/
  double st_A1;			/* Calib shutter A1 coefficient		*/
  double st_A2;			/* Calib shutter A2 coefficient		*/
  double st_A3;			/* Calib shutter A3 coefficient		*/
  double st_A4;			/* Calib shutter A4 coefficient		*/
  double st_A5;			/* Calib shutter A5 coefficient		*/	
  /* Baffle temperature sample and conversion to Degrees C */
  long   bt_time_offset;	/* Baffle samp time offset (microsec)*/
  long	 bt_hk_sample;		/* Baffle TM housekeeping sample	*/
  double bt_A0;			/* Baffle A0 coefficient		*/
  double bt_A1;			/* Baffle A1 coefficient		*/
  double bt_A2;			/* Baffle A2 coefficient		*/
  double bt_A3;			/* Baffle A3 coefficient		*/
  double bt_A4;			/* Baffle A4 coefficient		*/
  double bt_A5;			/* Baffle A5 coefficient		*/	
  /* Cold stage temperature sample and conversion to Degrees C */
  long   cs_time_offset;	/* Cold stage samp time offset (microsec)*/
  long	 cs_hk_sample;		/* Cold stage TM housekeeping sample	*/
  double cs_A0;			/* Cold stage A0 coefficient		*/
  double cs_A1;			/* Cold stage A1 coefficient		*/
  double cs_A2;			/* Cold stage A2 coefficient		*/
  double cs_A3;			/* Cold stage A3 coefficient		*/
  double cs_A4;			/* Cold stage A4 coefficient		*/
  double cs_A5;			/* Cold stage A5 coefficient		*/	
  /* Scan line corrector temperature sample and conversion to Degrees C */
  long   lc_time_offset;	/* Line corrector samp time offset (microsec)*/
  long	 lc_hk_sample;		/* Line corrector TM housekeeping sample*/
  double lc_A0;			/* Line corrector A0 coefficient	*/
  double lc_A1;			/* Line corrector A1 coefficient	*/
  double lc_A2;			/* Line corrector A2 coefficient	*/
  double lc_A3;			/* Line corrector A3 coefficient	*/
  double lc_A4;			/* Line corrector A4 coefficient	*/
  double lc_A5;			/* Line corrector A5 coefficient	*/	
  /* Shutter hub temperature sample and conversion to Degrees C */
  long   sh_time_offset;	/* Shutter hub samp time offset (microsec)*/
  long	 sh_hk_sample;		/* Shutter hub TM housekeeping sample	*/
  double sh_A0;			/* Shutter hub A0 coefficient		*/
  double sh_A1;			/* Shutter hub A1 coefficient		*/
  double sh_A2;			/* Shutter hub A2 coefficient		*/
  double sh_A3;			/* Shutter hub A3 coefficient		*/
  double sh_A4;			/* Shutter hub A4 coefficient		*/
  double sh_A5;			/* Shutter hub A5 coefficient		*/	
  /* Relay optics temperature sample and conversion to Degrees C */
  long   ro_time_offset;	/* Relay optics samp time offset (microsec)*/
  long	 ro_hk_sample;		/* Relay optics TM housekeeping sample	*/
  double ro_A0;			/* Relay optics A0 coefficient		*/
  double ro_A1;			/* Relay optics A1 coefficient		*/
  double ro_A2;			/* Relay optics A2 coefficient		*/
  double ro_A3;			/* Relay optics A3 coefficient		*/
  double ro_A4;			/* Relay optics A4 coefficient		*/
  double ro_A5;			/* Relay optics A5 coefficient		*/	
  /* Unpacked as serial word B  */
  long   wb_time_offset;	/* Serial word B samp time offset (microsec)*/
  long	 wb_hk_sample;		/* Serial word B TM housekeeping sample	*/
  double wb_A0;			/* Serial word B A0 coefficient		*/
  double wb_A1;			/* Serial word B A1 coefficient		*/
  double wb_A2;			/* Serial word B A2 coefficient		*/
  double wb_A3;			/* Serial word B A3 coefficient		*/
  double wb_A4;			/* Serial word B A4 coefficient		*/
  double wb_A5;			/* Serial word B A5 coefficient		*/	
  /* Unpacked as serial word D  */
  long   wd_time_offset;	/* Serial word D samp time offset (microsec)*/
  long	 wd_hk_sample;		/* Serial word D TM housekeeping sample	*/
  double wd_A0;			/* Serial word D A0 coefficient		*/
  double wd_A1;			/* Serial word D A1 coefficient		*/
  double wd_A2;			/* Serial word D A2 coefficient		*/
  double wd_A3;			/* Serial word D A3 coefficient		*/
  double wd_A4;			/* Serial word D A4 coefficient		*/
  double wd_A5;			/* Serial word D A5 coefficient		*/	
  /* Unpacked as serial word E  */
  long   we_time_offset;	/* Serial word E samp time offset (microsec)*/
  long	 we_hk_sample;		/* Serial word E TM housekeeping sample	*/
  double we_A0;			/* Serial word E A0 coefficient		*/
  double we_A1;			/* Serial word E A1 coefficient		*/
  double we_A2;			/* Serial word E A2 coefficient		*/
  double we_A3;			/* Serial word E A3 coefficient		*/
  double we_A4;			/* Serial word E A4 coefficient		*/
  double we_A5;			/* Serial word E A5 coefficient		*/	
  /* Unpacked as serial word F  */
  long   wf_time_offset;	/* Serial word F samp time offset (microsec)*/
  long	 wf_hk_sample;		/* Serial word F TM housekeeping sample	*/
  double wf_A0;			/* Serial word F A0 coefficient		*/
  double wf_A1;			/* Serial word F A1 coefficient		*/
  double wf_A2;			/* Serial word F A2 coefficient		*/
  double wf_A3;			/* Serial word F A3 coefficient		*/
  double wf_A4;			/* Serial word F A4 coefficient		*/ 
  double wf_A5;			/* Serial word F A5 coefficient		*/	
  /* Unpacked as serial word G  */
  long   wg_time_offset;	/* Serial word G samp time offset (microsec)*/
  long	 wg_hk_sample;		/* Serial word G TM housekeeping sample	*/
  double wg_A0;			/* Serial word G A0 coefficient		*/
  double wg_A1;			/* Serial word G A1 coefficient		*/
  double wg_A2;			/* Serial word G A2 coefficient		*/ 
  double wg_A3;			/* Serial word G A3 coefficient		*/ 
  double wg_A4;			/* Serial word G A4 coefficient		*/
  double wg_A5;			/* Serial word G A5 coefficient		*/	
  /* Unpacked as serial word L  */
  long   wl_time_offset;	/* Serial word L samp time offset (microsec)*/
  long	 wl_hk_sample;		/* Serial word L TM housekeeping sample	*/
  double wl_A0;			/* Serial word L A0 coefficient		*/
  double wl_A1;			/* Serial word L A1 coefficient		*/
  double wl_A2;			/* Serial word L A2 coefficient		*/
  double wl_A3;			/* Serial word L A3 coefficient		*/
  double wl_A4;			/* Serial word L A4 coefficient		*/
  double wl_A5;			/* Serial word L A5 coefficient		*/	
  /* Primary mirror temperature sample and conversion to Degrees C */
  long   pm_time_offset;	/* Primary mirror samp time offset (microsec)*/
  long	 pm_hk_sample;		/* Primary mirror TM housekeeping sample*/
  double pm_A0;			/* Primary mirror A0 coefficient	*/
  double pm_A1;			/* Primary mirror A1 coefficient	*/
  double pm_A2;			/* Primary mirror A2 coefficient	*/
  double pm_A3;			/* Primary mirror A3 coefficient	*/
  double pm_A4;			/* Primary mirror A4 coefficient	*/
  double pm_A5;			/* Primary mirror A5 coefficient	*/	
  /* Secondary mirror temperature sample and conversion to Degrees C */
  long   sm_time_offset;	/* Second mirror samp time offset (microsec)*/
  long	 sm_hk_sample;		/* Secondary mirror TM housekeeping sample*/
  double sm_A0;			/* Secondary mirror A0 coefficient	*/
  double sm_A1;			/* Secondary mirror A1 coefficient	*/
  double sm_A2;			/* Secondary mirror A2 coefficient	*/
  double sm_A3;			/* Secondary mirror A3 coefficient	*/
  double sm_A4;			/* Secondary mirror A4 coefficient	*/
  double sm_A5;			/* Secondary mirror A5 coefficient	*/	
  char   ser_wd_b[9];		/* Serial word B			*/
  char   ser_wd_d[9];		/* Serial word D			*/
  char   ser_wd_e[9];		/* Serial word E			*/
  char   ser_wd_f[9];		/* Serial word F			*/
  char   ser_wd_g[9];		/* Serial word G			*/
  char   ser_wd_l[9];		/* Serial word L			*/
  };

struct LTWGEAA 
  {
  /* EPHOMERIS & ATTITUDE ANCILLARY RECORD 
  ----------------------------------------*/
  long   rec_number;		/* Record Sequence Number 		*/
  long   frst_rec_code;		/* First record sub-type code  		*/
  long   rec_code;		/* Record type code 			*/
  long   scnd_rec_code;		/* Second record sub-type code		*/
  long   thrd_rec_code;		/* Third record sub-type code		*/
  long	 rec_length;		/* Length of this record		*/
  long   rec_seq_num;		/* Ephemeris & attitude data record seq num*/
  long   ephemeris_flag;	/* Ephemeris updated (0-NO, 1-YES)	*/
  long   attitude_flag;		/* Attitude updated (0-NO, 1-YES)	*/
  long   gyro_flag;		/* Gyro updated (0-NO, 1-YES)		*/
  long   gyro_drift_flag;	/* Gyro drift updated (0-NO, 1-YES)	*/
      /*   PCD MAJOR FRAME START TIME    */
  long   pcd_start_day;		/* Day of year at start of PCD frame	*/
  long   pcd_start_msec;	/* Milliseconds of day at start of PCD frame*/
  long   pcd_start_usec;	/* Microseconds at start of PCD frame   */
  long   pcd_start_tenths;	/* Tenths of a microsecond at start of PCD*/
     /* EPHEMERIS DATA */
  long   eph_time_offset;	/* Ephemeris measure time offset (microsec)*/
  double craft_pos[3];		/* Spacecraft position components x,y,z	*/
  double craft_vel[3];		/* Spacecraft velosity components x,y,z	*/
     /* ATTITUDE DATA */
  long   att_time_offset;	/* Attitude measure time offset (micorsec)*/
  double epa1;			/* EPA1					*/
  double epa2;			/* EPA2					*/
  double epa3;			/* EPA3					*/
  double epa4;			/* EPA4					*/
     /* GYRO DATA  */
  long   gyro_time_offset;	/* Gyro measure time offset (microsec)	*/
  double gyro_meas[3][64];	/* sets of gyro measurments components	*/
     /* GYRO DRIFT DATA */
  long   gyro_drift_time_offset;/* Gyro drift measure time offset (microsec)*/
  double gyro_drift[3];		/* Gyro drift components		*/
  };

struct LTWGRJR
  {
  /* RAW JITTER MEASUREMENTS ANCILLARY RECORD	
  -------------------------------------------*/
  long   rec_number;		/* Record Sequence Number 		*/
  long   frst_rec_code;		/* First record sub-type code  		*/
  long   rec_code;		/* Record type code 			*/
  long   scnd_rec_code;		/* Second record sub-type code		*/
  long   thrd_rec_code;		/* Third record sub-type code		*/
  long	 rec_length;		/* Length of this record		*/
  long   rec_seq_num;		/* Raw jitter measurements seq num	*/ 
  long   seq_num;		/* Sequence number within PCD major frame*/
  long   day_start;		/* Day of year at start of PCD frame STC*/
  long   millisec_start;	/* Millisec of day at start of frame	*/
  long   microsec_start;	/* Microsec of day at start of frame	*/
  long   tenths_microsec_start;	/* Tenths of a microsec of day at start */
    /* ADS TEMPERATURE DATA */
  long   ads_temp1_time_offset;	/* ADS temp (1) time offset (microsec)	*/
  double ads_temp1;		/* ADS temperature 1 (degrees C)	*/
  long   ads_temp2_time_offset;	/* ADS temp (2) time offset (microsec)	*/
  double ads_temp2;		/* ADS temperature 2 (degrees C)	*/
  long   ads_temp3_time_offset;	/* ADS temp (3) time offset (microsec)	*/
  double ads_temp3;		/* ADS temperature 3 (degrees C)	*/
  long   ads_temp4_time_offset;	/* ADS temp (4) time offset (microsec)	*/
  double ads_temp4;		/* ADS temperature 4 (degrees C)	*/
    /* ADS MEASUREMENTS */
  long   ads_measur1[3];	/* ADS measurment (1) (3 axes) (microsec)*/
  long   ads_measur[3][1024];	/* Array of ADS measurments for 3 axes	*/
  };

struct LTWGMTR
  {
  /* MISSION TELEMETRY ANCILLARY RECORD
  -------------------------------------*/
  long   rec_number;            /* Record Sequence Number               */
  long   frst_rec_code;         /* First record sub-type code           */
  long   rec_code;              /* Record type code                     */
  long   scnd_rec_code;         /* Second record sub-type code          */
  long   thrd_rec_code;         /* Third record sub-type code           */
  long   rec_length;            /* Length of this record                */
  long   rec_seq_num;           /* Mission telemetry data record seq num*/
    /* TELEMETRY MAJOR FRAME START TIME */
  long   day_start;             /* Day of year at start of PCD frame STC*/
  long   millisec_start;        /* Millisec of day at start of frame    */
  long   microsec_start;        /* Microsec of day at start of frame    */
  long   tenths_microsec_start; /* Tenths of a microsec of day at start */
    /* EPHEMERIS DATA */
  long   eph_time_offset;       /* flight software time offset          */
  double craft_pos[3];          /* spacecraft position x,y, and z comp  */
  double craft_vel[3];          /* spacecraft velocity x,y, and z comp  */
    /* ATTIDTUDE DATA */
  long   att_time_offset;       /* flight software time offset          */
  double epa1;                  /* EPA1                                 */
  double epa2;                  /* EPA2                                 */
  double epa3;                  /* EPA3                                 */
  double epa4;                  /* EPA4                                 */
    /* ROLL, PITCH, YAW AXIS ANGULAR INCREMENTS */
  long   inc_time_offset;	/* flight software time offset		*/
  double roll_inc;		/* role axis angular increment		*/
  double pitch_inc;		/* pitch axis angular increment		*/
  double yaw_inc;		/* Yaw axis angular increment		*/
    /* ROLL, PITCH, YAW ANGULAR RATE */
  long   rate_time_offset;	/* flight software time offset		*/ 
  double roll_rate;		/* role axis angular rate		*/
  double pitch_rate;		/* Pitch axis angular rate		*/
  double yaw_rate;		/* Yaw axis angular rate		*/
    /* ROLL, PITCH, YAW ATTITUDE ERROR */
  long   err_time_offset;	/* flight software time offset		*/
  double roll_err;              /* role axis attitude error             */
  double pitch_err;             /* Pitch axis attitude error            */
  double yaw_err;               /* Yaw axis attitude error              */
    /* HOUSEKEEPING DATA SAMPLES */
  double item1_coef[6];		/* A0 - A5 coefficients for item 1	*/
  double item2_coef[6];         /* A0 - A5 coefficients for item 2      */
  double item3_coef[6];         /* A0 - A5 coefficients for item 3      */
  double item4_coef[6];         /* A0 - A5 coefficients for item 4      */
  double item5_coef[6];         /* A0 - A5 coefficients for item 5      */
  double item6_coef[6];         /* A0 - A5 coefficients for item 6      */
  double item7_coef[6];         /* A0 - A5 coefficients for item 7      */
  double item8_coef[6];         /* A0 - A5 coefficients for item 8      */
  double item9_coef[6];         /* A0 - A5 coefficients for item 9      */
  double item10_coef[6];        /* A0 - A5 coefficients for item 10     */
  double item11_coef[6];        /* A0 - A5 coefficients for item 11     */
  double item12_coef[6];        /* A0 - A5 coefficients for item 12     */
  long   samples[12][128];	/* 128 samples, S, of each of a maximum */
				/* 12 data items, D, - stored in the order */
				/* S1D1, S1D2, ..., S1D12,..., S128,D12 */
  };

#endif
