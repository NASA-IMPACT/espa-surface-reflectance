#ifndef EDIPS_H
#define EDIPS_H

/*****************************************************************************
NAME:				EDIPS.H

PURPOSE:  Include file defining structure for EDIPS files.

PROGRAM HISTORY:
VERSION		DATE	AUTHOR		DESCRIPTION
-------		----	------		-----------
  1.0		12/92	S. Nelson	Original Development

******************************************************************************/
#include "sysdef.h"
#include "typlim.h"
#include "worgen.h"

struct EDIPSH
  {
  /*  A. Scene identification
      ----------------------*/
  long	record_number;		/* Record number for each file present	*/
  long	record_type_code;	/* Record type code			*/
  long	record_length;		/* Record length			*/
  char	image_id[13];		/* image ID of the form ADDDDHHMMSB
				   where  A    = Landsat mission
					  DDDD = day number, relative to 
						 launch
					  HH   = Hour of observation
					  MM   = Minute of observation
					  S    = Tens of seconds of observation
					  B    = IPF ID or Band ident. code */
  char	wrs_designator[9];	/* Terrestrial image identifier of the form
				   MPPPRR where
					  M   = A (ascending) or D (decending)
					  PPP = WRS path number (001-251)
					  RRR = WRS row number (001-248)*/
  long	tape_generation_day;	/* date of tape generation		*/
  long  tape_generation_mon;	/* month of tape generation		*/
  long	tape_generation_yr;	/* year of tape generation		*/
  /*  B. Spacecraft Description
      ------------------------*/
  char	sensor_id[4];		/* Sensor ID (either MSS or RBV)	*/
  long  mission_number;		/* Mission number			*/
  long  orbit_number;		/* Orbit number of the spacecraft	*/
  long  status_bnd_detect[5][6];/* Active detector status one detector for each
				   band in the form [detector][band]:
					0 = inactive	1 = active	*/
  long  active_detector_cnt;	/* Active detector/camera count		*/
  long  pix_original_image;	/* Number of pixels in original image	*/
  /*  C. Time of Exposure/WRS designator
      ---------------------------------*/
  long	wrs_scan_line;		/* Scan line number of WRS		*/
  long  wrs_pixel_number;	/* Pixel number of WRS			*/
  char  center_year[3];		/* Center picture exposure time year	*/
  char  center_julian_day[4];	/* Center picture julian day		*/
  char  center_hour[3];		/* Center picture hour			*/
  char  center_min[3];		/* Center picture minutes		*/
  char  center_sec[3];		/* Center picture seconds		*/
  char  center_milliseconds[4];	/* Center picture milliseconds		*/
  /*  D. Data Identification and Characteristics
      -----------------------------------------*/
  long  header_length;		/* Header record length			*/
  long  number_header_records;	/* Number of header records		*/
  long  header_bytes;		/* Number of bytes of header data	*/
      /*  Annotation Data Characteristics
          ------------------------------*/
  long	annotat_length;		/* Annotation record length		*/
  long	number_annotat_records;	/* Number of annotation records		*/
      /*  Ancillary Data Characteristics
          -----------------------------*/
  long  ancill_length;		/* Ancillary record length		*/
  long  number_ancill_records;	/* Number of ancillary records		*/
  char  geom_correct[4];	/* Geometric corrections applied (yes/no) */
  char	geom_data[4];		/* Geometric data present (yes/no)	*/
  char  radiometric_correct[4];	/* Radiometric corrections applied (yes/no) */
  char	radiometric_data[4];	/* Radiometric data present (yes/no)	*/
	/*  Image Data Characteristics
   	    -------------------------*/
  long  image_record_length;	/* Image record length			*/
  long  calibration_words_line;	/* Number of calibration/quality data words
				   per scan line			*/
  char  image_data_format[17];	/* (Unframed/Framed rect/Framed square)	*/
  char  interleav_type[4];	/* Interleaving type indicator (BSQ/BIL) */
  long	bil_interleaving_count; /* BIL line interleaving count (0 for N/A) */
  long  bits_per_pixel;		/* Number of bits per pixel		*/
  long	resampling;		/* Resampling applied where:
					192 = None
					  9 = Cubic convolution
					 18 = Nearest neighbor		*/
  long  map_projection;		/* Map projection applied to first ancillary
				   and annotation data sections where:
					192 = None
					  9 = Universal Transverse Mercator
					 18 = Polar Stereographic
					219 = Hotline Oblique Mercator
					237 = Space Oblique Mercator	*/
  long	wrs_center_offset;	/* WRS offset from fully processed image
				   center - right(positive) left (negative) */
  char	data_justification[6];	/* Image data justification (left/right)*/
  char  most_signifcnt_bit[6];	/* Locatn of most significant bit (left/right)*/
  long  pixels_per_line;	/* Number of pixels per scan line	*/
  long  images_per_scene;	/* Number of usable images per scene(1,4 or 5)*/
  char  mss_band_number[2];	/* MSS band number or RBV camera number	*/
  long	support_data;		/* Number of bits of support data or "one"
				   filler in the end of each image data record
				   where:  252 = PM
					   378 = PR
					   252 = AM (landsat 1, 2, and 3)
					     7 = AR
				 	     0 = AM (landsat 4)		*/
	/* Trailer Data Characteristics
	   ---------------------------*/
  long	trailer_length;		/* Trailer record length		*/
  long	number_trailer_record;	/* Number of trailer records (1 per image)*/ 
	/* Special Purpose Fields
	   ---------------------*/
  char	orbit_direction[6];	/* Orbital direction where (day/night)	*/
  double orientation_angle;	/* Image orientation angle in		*/
  long	sensor_mode;		/* Sensor mode where:
				      7 = low gain lanear transmission
				     56 = low gain compressed transmission
				     63 = high gain linear transmission
				    192 = high gain compressed transmission */
	/* Temporal Registration Data--NO data available for Imagery without
	   Geometric Corrections
	   ----------------------------------------------------------------*/
  char	temporal_scene_id[13];	/* Scene ID used for temporal registration
				   processing of the form:  ADDDDHHMMSB where
					A    = Landsat mission
					DDDD = Day number, relative to launch,
					       at time of observation
					HH   = Hour at time of observation
					MM   = Minute at time of observation
					S    = Tens of seconds
					B    = IPF ID code or Band	*/
  char	temp_wrs_designator[9];	/* WRS designator - terrestrial image
				   identifier of the form:  MPPPRRR where
					M   = A (ascending) or D (descending)
					PPP = nominal WRS path number
					RRR = nominal WRS row number	*/
  long	cur_imag_line[4];	/* Current image scan line for temporal
				   registration marks (for four points)	*/
  long  cur_imag_pixel[4];	/* Current image pixel number for temporal
				   registration marks (for four points) */
  long 	ref_imag_line[4];	/* Reference image scan line for temporal
				   registration marks (for four points) */
  long	ref_imag_pixel[4];	/* Reference image pixel number for temporal
				   registration marks (for four points) */
  long	first_overlap_line;	/* Scan line of 1st overlap mark(upper left)*/
  long	first_overlap_pixel;	/* Pixel number of 1st overlap mark	*/
  long	second_overlap_line;	/* Scan line of 2nd overlap mark(upper right)*/
  long	second_overlap_pixel;	/* Pixel number of 2nd overlap mark	*/
  long	third_overlap_line;	/* Scan line of 3rd overlap mark(lower left)*/
  long	third_overlap_pixel;	/* Pixel number of 3rd overlap mark	*/
  long	fourth_overlap_line;	/* Scan line of 4th overlap mark(lower right)*/
  long	fourth_overlap_pixel;	/* Pixel number of 4th overlap mark	*/
  long	overlap_pixel_offset;	/* Nominal overlap mark pixel offset	*/
  char	quality_code[2];	/* Geometric correction quality code	*/
  long	tick_marks_top;		/* Actual number of tick marks for top	*/
  long	tick_marks_left;	/* Actual number of tick marks for left	*/
  long	tick_marks_right;	/* Actual number of tick marks for right */
  long	tick_marks_bottom;	/* Actual number of tick marks for bottom */
	/* data available for both corrected and uncorrected data types
	   -----------------------------------------------------------*/
  char	band_quality[2];	/* Overall band quality indicator	*/
  long	radiometric_calibratn;	/* Radiometric calibration Method where:
					 0 = No corrections applied
					 9 = Histogram method
					27 = Cal wedge values only (no hist.)
					45 = Non-standard corrections applied*/
  float calibration_accuracy;	/* Relative Calibration Accuracy (RCA), maximum
				   difference between detector means for the
				   image, where:
					  0 <= RCA <= 1.0  Excellent
					1.0 <  RCA <= 2.0  Good
					2.0 <  RCA         Acceptable	*/
	/* Input Data Quality Indicators - data generating partially
	   processed image
	   --------------------------------------------------------*/
 	   /* Telemetry:
	      ----------*/ 
  long tel_ephemeris_points;	/* Number of ephemeris data points in interval*/
  long tel_rejected_ephemeris;	/* Rejected ephemeris data points in interval */
  long tel_attitude_points;	/* Number of attitude data points in interval */
  long tel_rejected_attitude;	/* Rejected attitude data points in interval  */
  float tel_interval_length;	/* Length of telemetry interval in seconds    */

 		/* Accuracy of ephemeris and attitude fit, RMS difference in
		   meters between fit and data points for ephemeris--altitude,
		   along-track position, and across-track position and 
		   attitude--pitch, roll, and yaw.
		   ----------------------------------------------------------*/
  float tel_ephm_accuracy_alt;	 /* Accuracy of ephemeris fit, altitude	    */
  float tel_ephm_acc_along_trck; /* Accuracy of ephemeris fit, along-track  */
  float tel_ephm_acc_acros_trck; /* Accuracy of ephemeris fit, across-track */

  float tel_att_accuracy_pitch;	/* Accuracy of attitude fit, pitch	*/
  float tel_att_accuracy_roll;	/* Accuracy of attitude fit, roll	*/
  float tel_att_accuracy_yaw;	/* Accuracy of attitude fit, yaw	*/

	    /* Control Points:
	       --------------*/
  char	cp_band_quality[4][2];	/* Overall Band qualities of scene from which
				   control points were extracted (4 bands)  */
  long	cp_geodetic_points;	/* Number of geodetic points used to correct
				   reference image			    */
  long	cp_registratn_success;	/* Average previous registration success    */
  float cp_autocorrelatn_init;	/* Average initial autocorrelation peak value
				   in control point generation process	    */
  float cp_err_ellipse_along;	/* 90% error ellipse of geodetic c.p. location
				   in corrected image--along-track.	    */
  float cp_err_ellipse_across;	/* 90% error ellipse of geodetic c.p. location
				   in corrected image--across-track.	    */
  float cp_autocorrelatn_ref;	/* Average autocorrelation peak values of 
				   control points--reference image.	    */
  float cp_suitability;		/* Average c.p. suitability measure	    */
  char	data_source[2];		/* Data source where:
				    W = TDRSS/White Sands	S = Simulator
				    U = Alaska			N = NTTF
				    F = Foreign  		G = Goldstone
				    T = Transportable Ground Station	    */
  long	uncorrect_ecc_count;	/* Uncorrectable ECC count for the scene    */
  long	bit_error_rate;		/* Indication of bit error rate		    */
  long	cal_wedge_value_use;	/* Use of Nominal Calibraton Wedge Values(CWVs)
					 0 = Not used
					 7 = Used for comparison only
					56 = Used to replace CWVs outside window
					     but not in radiometric calibration
					63 = Used to replace CWVs outside window
					     and in radiometric calibration */
  long	window_size;		/* Neighborhood which actual CWVs are compared*/
  long	cal_wedge_values[36];	/* Nominal calibration wedge values	    */
  long	wed_qual_sens_samp[6][6];/* Calibration wedge quality for each
				   sample and sensor.			    */
  float wrs_center_lat;		/* WRS center latitude in radians	    */
  float wrs_center_long;	/* WRS center longitude in redians	    */
  long	contrast_enhancement;	/* EDIPS performed contrast enhancement
					0 = false	255 = true	    */
  long	atmospheric_scatter;	/* EDIPS performed atmospheric scatter
					0 = false	255 = true	    */
  long	edge_enhancement;	/* EDIPS performed edge enhancement
					0 = false	255 = true	    */
  char	sensor_gain_opt[4][2];	/* Sensor gain options for each band where:
					H = high gain	L = low gain	    */
  char	transmission_type[4][2];/* Type of MSS transmission for each band
					1 = linear mode	2 = compressed mode */

  };

	/* Trailer record information
	----------------------------*/
struct EDIPST
  {
  long	record_number;		/* Record number			*/
  long	record_type;		/* Record type code			*/
  long	last_scene_data_acqu;	/* Flag for last scene in a data acquisition
				   pass or swath:
					0 = No	255 = Yes		*/
  long	last_scene_hdt;		/* Flag for last scene on the HDT
					0 = No	255 = Yes		*/
  long	destriped_image;	/* Destriped image			*/
  char	contrast_stretch[2];	/* Units of contrast stretch values,
					grey level			*/
  long	minimum_contrast_enh;	/* Units of contrast enhancement minimum*/
  long	maximum_contrast_enh;	/* Units of contrast enhancement maximum*/
  long	radiance_value_scatter;	/* Radiance value used for atmospheric
				   scatter compensation			*/
  long	edge_enh_kernal_size_x;	/* Edge enhancement kernal size in x direction*/
  long	edge_enh_kernal_size_y;	/* Edge enhancement kernal size in y direction*/
  };
	/* Annotation record information
	-------------------------------*/
struct EDIPSA
  {
  long	record_number;		/* Record number			*/
  long	record_type;		/* Record type code			*/
  char	exposure_day_mon_yr[9];	/* Day, month, year of exposure		*/
  char	format_ctr_lat_long[17];/* Format center lat/long in deg and minutes*/
  char	path_row_identifier[9];	/* Path row identifier			*/
  char	wrs_ctr_lat_long[17];	/* WRS center lat and long		*/
  char	sensor_band_id_code[8];	/* Sensor spectral-band identification code */
  char	transmission[2];	/* Transmission: D = Direct transmission
				      R = Stored data from WBUT recorder    */
  char	sun_angles[14];		/* Sun-elevation and -azimuth angle measured
				   clockwise from true North at m.p. of MSS
				   frame to nearest degree		*/
  char	correction_applied[2];	/* Type of correction applied
					U = uncorrected
					S = system level corrected
					G = geometrically using geodetic info
					    (no temporal registration)
					T = temporal reg. using geodetic
					    inform. from single reference scene
					R = temporal reg. to single ref. scene
					    (no geodetic info)		*/
	/* Landsat 1-3 only */
  char	scale_of_image[2];	/* Scale of image
					1 = 185 km x 185 km
					2 = 99 km x 99 km
					3 = 185 km x 170 km		*/
  char	projection[2];		/* Projection
					L = Lambert projection
					P = Polar stereographic
					S = Space oblique mercator
					U = Universal transverse mercator
					H = Hotline oblique mercator	*/
  char	resample_algorithm[2];	/* Resampling algorithm
					C = cubic convolution
					N = Nearest neighbor		*/
  char	type_ephemeris_data[2];	/* Type of ephemeris data to compute image cntr
					P = predictive
					D = definitive
					G = global positioning system	*/
  char	process_procedure[2];	/* Processing procedure
					A = abnormal
					N = normal			*/
  char	sensor_gain_opt[2];	/* Sensor gain options
					H = high gain
					L = low gain			*/
  char	mss_transmission[2];	/* Type of MSS transmission
					1 = linear mode
					2 = compressed mode		*/
  char  agency_and_project[13];	/* Agency and project			*/
  char	frame_id_number[16];	/* Frame identification number in the form
				   E-ADDDD-HHMMS-B
					E    = Project identifier
					A    = Mission
					DDDD = Day number, relative to launch
					HH   = Minute
					S    = Tens of seconds
					B    = Image code		*/

	/* Tick Mark data
	   -------------*/
  long  tm_top_location[16];	/* Location of Tick Mark--pixel number	*/
  char  tm_top_loctn_flg[16][2];/* Flag of location with respect to tick mark
				    1 = above or to left of tick mark with
					leading blanks
				    0 = below or to right of tick mark with
					trailing blanks			*/
  char  tm_top_coord[16][8];	/* Coordinate Data			*/
  long  tm_left_location[25];	/* Location of Tick Mark--line number	*/
  char  tm_left_loctn_flg[25][2];/* Flag of location with respect to tick mark
				    1 = above or to left of tick mark with
					leading blanks
				    0 = below or to right of tick mark with
					trailing blanks			*/
  char  tm_left_coord[25][8];	/* Coordinate Data			*/
  long  tm_rgt_location[25];	/* Location of Tick Mark--line number	*/
  char  tm_rgt_loctn_flg[25][2];/* Flag of location with respect to tick mark
				    1 = above or to left of tick mark with
					leading blanks
				    0 = below or to right of tick mark with
					trailing blanks			*/
  char  tm_rgt_coord[25][8];	/* Coordinate Data			*/
  long  tm_bot_location[16];	/* Location of Tick Mark--pixel number	*/
  char  tm_bot_loctn_flg[16][2];/* Flag of location with respect to tick mark
				    1 = above or to left of tick mark with
					leading blanks
				    0 = below or to right of tick mark with
					trailing blanks			*/
  char  tm_bot_coord[16][8];	/* Coordinate Data			*/
  };
struct EDPSAC
  {
	/* RECORD 1
	   -------*/
  long	r1_record_number;	/* Record number			*/
  long	r1_record_type_code;	/* Record type code			*/
  long  r1_record_length;	/* Record length			*/
  long	input_pixels_per_line;	/* number of pixels per input line	*/
  long	input_lines_per_image;	/* number of input lines in image	*/
  double input_inter_pixel_dis;	/* scale of input inter-pixel distance
				   in meters/pixel 			*/
  double input_inter_line_dis;	/* scale of input inter-line distance
				   in meters/pixel 			*/
  long	output_pixels_per_line;	/* pixels per output line of processed image*/
  long	output_lines_per_image;	/* lines per output image of processed image*/
  double output_inter_pixl_dis;	/* scale of output inter-pixel distance
				   in meters/pixel 			*/
  double output_inter_line_dis;	/* scale of output inter-line distance
				   in meters/pixel 			*/
  double spacecraft_altitude;	/* spacecraft altitude in meters	*/
  double input_image_width;	/* input image width in meters		*/
  double mss_mirror_coef[4];	/* MSS mirror model coefficients	*/
  double mss_max_mirror_ang;	/* MSS maximum mirror angle in radians	*/
  double scan_skew_constant;	/* Scan skew constant			*/
  double tim_bet_mirror_sweeps;	/* time between successive MSS mirror sweeps
				   in seconds				*/
  double tim_active_mirror_swp;	/* time for the active portion of an MSS
				   mirror sweep in seconds		*/
  double semi_major_axis;	/* semi-major axis of Earth ellipsoid	*/
  double semi_minor_axis;	/* semi-major axis of Earth ellipsoid	*/
  double earth_curvature_const;	/* Earth curvature constant		*/
  float	mss_samp_delay_const[24];/* MSS sampling delay constants--input image
		 		   along-scan pixels			*/
  float mss_band_offsets[3];	/* MSS band-to-band offsets with respect to
				   band 1 (one for each of bands 2,3,4)	*/
	/* RECORD 2 
	   -------*/
  long	r2_record_number;	/* Record number			*/
  long	r2_record_type_code;	/* Record type code			*/
  long  r2_record_length;	/* Record length			*/
  char	wrs_frame_number[5];	/* WRS frame number			*/
  char	wrs_orbit_number[5];	/* WRS orbit number			*/
  double wrs_ctr_lat;		/* WRS format center latitude in radians */
  double wrs_ctr_long;		/* WRS format center latitude in radians */
  char  center_time[15];	/* Spacecraft time of frame center
				  Yr,Yr,D,D,D,Hr,Hr,Min,Min,Sec,Sec,ms,ms,ms */
	/* Scene Center information
  	   -----------------------*/
  double scene_ctr_lat;		/* Scene center latitude in radians	*/
  double scene_ctr_long;	/* Scene center longitude in radians	*/
  double sc_earth_cntrd_coords[3];/* Scene center in Earth-centered Earth-fixed
				   coordinates in meters		*/
  double sc_heading_angle;	/* Spacecraft heading angle at scene
				   center in radians			*/
  double sc_scan_line_coord;	/* Scan line coordinate of scene center in 
				   partially-processed image		*/
  double sc_pixel_coord;	/* Pixel coordinate of scene center in
				   partially-processed image		*/
  double sc_norm_velocity_err;	/* Normalize spacecraft velocity error from
				   nominal at scene center		*/
  double sc_earth_rotat_velcty;	/* Earth rotation velocity at scene center
				   in meters per second			*/
  float earth_rotation_param;	/* Earth rotation parameter (image skew) in
				   radians				*/
	/* Spacecraft state vector at scene center
	   --------------------------------------*/
  double sc_state_vector_pitch;	/* Pitch in radians			*/
  double sc_state_vector_roll;	/* Roll in radians			*/
  double sc_state_vector_yaw;	/* Yaw in radians			*/
  double sc_state_vector_x;	/* X in Km				*/
  double sc_state_vector_y;	/* Y in Km				*/
  double sc_state_vector_z;	/* Z in Km				*/
  double sc_s_vctr_delta_pitch;	/* Delta pitch in radians/sec		*/
  double sc_s_vctr_delta_roll;	/* Delta roll in radians/sec		*/
  double sc_s_vctr_delta_yaw;	/* Delta yaw in radians/sec		*/
  double sc_s_vctr_delta_x;	/* Delta X in Km/sec			*/
  double sc_s_vctr_delta_y;	/* Delta Y in Km/sec			*/
  double sc_s_vctr_delta_z;	/* Delta Z in Km/sec			*/
  long	cps_in_attitude;	/* Number of CPs in attitude/ephemeris fit */
  long	gcps_in_attitude;	/* Number of GCPs in attitude/ephemeris fit */
  long	cp_correlatns_attempt;	/* Number of CP correlations attempted	*/
  long	cps_reject_correlatns;	/* Number of correlated CPs rejected	*/
  float	rms_along_track_error;	/* RMS along-track geom modeling error/meters*/
  float	rms_across_track_error;	/* RMS across-track geom modeling error/meters*/
  long	cps_per_zone[25];	/* Distribution of CPs used of the WRS frame */
  char	cp_identificatn[25][9];	/* Identification of CPs used of the form
				    BTXXYYY where:
					B = band number	T = type
				       XX = Zone      YYY = Sequence	*/

  };
	/* Geometric Correction Parameters Ephemeris Data
	   ---------------------------------------------*/
struct EDPSAC_EPH
  {
  char	time_of_first_set[15];	/* Time of the first set of ephemeris entries
			          Yr,Yr,D,D,D,Hr,Hr,Min,Min,Sec,Sec,ms,ms,ms */
  float	interval_between_sets;	/* Time interval between successive sets/sec */
  long number_of_sets;	/* Number of sets of ephemeris entries	*/
  float spccrft_locatn_x[16];  	/* Spacecraft location - x		*/
  float spccrft_locatn_y[16];	/* Spacecraft location - y		*/
  float spccrft_locatn_z[16];	/* Spacecraft location - z		*/
  float spccrft_velocity_x[16];	/* Spacecraft velocity - Vx		*/
  float spccrft_velocity_y[16];	/* Spacecraft velocity - Vy		*/
  float spccrft_velocity_z[16];	/* Spacecraft velocity - Vz		*/
  long	data_quality[16];	/* Data quality indicator where:
				     0 = Corresponding input data/valid
				     1 = No corresponding input data
				     2 = Corresponding input data/not valid */
  };
	/* Attitude Data
	   ------------*/
struct EDPSAC_ATT
  {
  char	time_of_first_set[15];	/* Time of first set of attitude entries
			          Yr,Yr,D,D,D,Hr,Hr,Min,Min,Sec,Sec,ms,ms,ms */
  float	interval_between_sets;	/* Time interval between successive sets/sec */
  long 	number_of_sets;		/* Number of sets of attitude entries	*/
  float	pitch_angle[60];	/* pitch angle in radians		*/
  float	roll_angle[60];		/* roll angle in radians		*/
  float	yaw_angle[60];		/* yaw angle in radians			*/
  long	data_quality[60];	/* Data quality indicator where:
				     0 = Valid data
				     1 = Angular increment data not valid,
					 replaced by last good value
				     2 = Angular increment data not valid,
					 replaced by 0
				     3 = No valid drift bias, replaced by 0
				     4 = Angular increment data not available
				 	 initial attitude information used
				     5 = initial attitude information not 
					 available, replaced by 0	*/
	/* Partial Derivatives for SOM Projection
	   --------------------------------------*/
  float pd_along_trk_x[3][5];	/* Partial derivative along-track 3x5 matrix */
  float pd_along_trk_y[3][5];	/* Partial derivative along-track 3x5 matrix */
  float pd_acros_trk_x[3][5];	/* Partial derivative across-track 3x5 matrix */
  float pd_acros_trk_y[3][5];	/* Partial derivative across-track 3x5 matrix */
  float pd_altitude_x[3][5];	/* Partial derivative altitude 3x5 matrix */
  float pd_altitude_y[3][5];	/* Partial derivative altitude 3x5 matrix */
  float pd_pitch_x[3][5];	/* Partial derivative pitch 3x5 matrix	*/
  float pd_pitch_y[3][5];	/* Partial derivative pitch 3x5 matrix	*/
  float pd_roll_x[3][5];	/* Partial derivative roll 3x5 matrix	*/
  float pd_roll_y[3][5];	/* Partial derivative roll 3x5 matrix	*/
  float pd_yaw_x[3][5];		/* Partial derivative yaw 3x5 matrix	*/
  float pd_yaw_y[3][5];		/* Partial derivative yaw 3x5 matrix	*/
	      /* Multiplicative and additive radiometric correction constants
	         -----------------------------------------------------------*/
  float rc_mult_band_det[4][6];	/* Multiplicative radiometric correction
				   constants of the form [band][detector] */
  float rc_add_band_det[4][6];	/* Additive radiometric correction
				   constants of the form [band][detector] */
  };
struct EDPSAC_RS
  {
  /* Data for HRS (Horizontal Resampling) and VRS (Vertical Resampling) 
     ----------------------------------------------------------------*/
		/* Record type code for each record number
			byte 5 = 36
			byte 6 = 36
			byte 7 = 219 for UTM/PS
				 237 for SOM/HOM
			byte 8 = 18					*/
  long	rn4_rec_type_code;	/* Record type code for record num 4/12	*/
  long	rn5_rec_type_code;	/* Record type code for record num 5/13	*/
  long	rn6_rec_type_code;	/* Record type code for record num 6/14	*/
  long	rn7_rec_type_code;	/* Record type code for record num 7/15	*/
  long	rn8_rec_type_code;	/* Record type code for record num 8/16	*/
  long	rn9_rec_type_code;	/* Record type code for record num 9/17	*/
  long	rn10_rec_type_code;	/* Record type code for record num 10/18*/
  long	rn11_rec_type_code;	/* Record type code for record num 11/19*/
  double hrs_row_coord[51][61]; /* HRS Pixel Coordinates for in the form
				    of [row][coord]             	*/
  double hrs_ln_fill_left[51]; 	/* HRS line fill left count [row]	*/
  double hrs_ln_fill_rgt[51]; 	/* HRS line fill right count [row]	*/
  double vrs_row_coord[44][61];	/* VRS Pixel Coordinates for in the form
				    of [row][coord]			*/
  long	pixel_wrs_cntr;		/* Pixel number of WRS center in fully
				   processed image			*/
  long	offset_wrs_cntr;	/* Offset of WRS center from fully processed
				   image center pixel			*/
  char	temp_reg_scene_id[21];	/* Temporal Registration Scene ID form
					ADDDDHHMMSB where:
					A = Landsat Mission *
					DDDD = Day number
					HH = Hour at time of observation
					MM = Minute at time of observation
					S = Tens of seconds at time of obs.
					B = Band ID code		*/
  long	tr_cur_imag_line[4];	/* Current image scan line for temporal
				   registration marks (for four points)	*/
  long  tr_cur_imag_pixel[4];	/* Current image pixel number for temporal
				   registration marks (for four points) */
  long 	tr_ref_imag_line[4];	/* Reference image scan line for temporal
				   registration marks (for four points) */
  long	tr_ref_imag_pixel[4];	/* Reference image pixel number for temporal
				   registration marks (for four points) */
  long	first_overlap_line;	/* Scan line of 1st overlap mark(upper left)*/
  long	first_overlap_pixel;	/* Pixel number of 1st overlap mark	*/
  long	second_overlap_line;	/* Scan line of 2nd overlap mark(upper right)*/
  long	second_overlap_pixel;	/* Pixel number of 2nd overlap mark	*/
  long	third_overlap_line;	/* Scan line of 3rd overlap mark(lower left)*/
  long	third_overlap_pixel;	/* Pixel number of 3rd overlap mark	*/
  long	fourth_overlap_line;	/* Scan line of 4th overlap mark(lower right)*/
  long	fourth_overlap_pixel;	/* Pixel number of 4th overlap mark	*/
  long	tick_marks_top;		/* Actual number of tick marks for top	*/
  long	tick_marks_left;	/* Actual number of tick marks for left	*/
  long	tick_marks_right;	/* Actual number of tick marks for right */
  long	tick_marks_bottom;	/* Actual number of tick marks for bottom */
  long	samp_corner_pts[4];	/* Input sample value of four corner points
				   in output image (band independent)	*/
  double image_orientatn_ang;	/* Image Orientation Angle with respect
				   to center line--in radians		*/
  long	num_sweeps;		/* Number of sweeps prior to scene center
				   at which the HRS, VRS grid points begin */
  };

#endif
