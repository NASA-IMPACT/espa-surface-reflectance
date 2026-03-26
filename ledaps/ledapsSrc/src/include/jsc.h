#ifndef JSC_H
#define JSC_H

/*****************************************************************************
NAME:				EDIPS.H

PURPOSE:  Include file defining structure for EDIPS files.

PROGRAM HISTORY:
VERSION		DATE	AUTHOR		DESCRIPTION
-------		----	------		-----------
  1.0		01/93	S. Nelson	Original Development

******************************************************************************/
#include "sysdef.h"
#include "typlim.h"
#include "worgen.h"

struct JSCUH
  {
  /*  Universal Header Record
      -----------------------*/
  char	computing_sys_id[33];	/* Computing system identification	*/
  char	tape_library_id[21];	/* Tape library identification of master tape */
  char	sensor_id[9];		/* Sensor identification		*/
  long	master_tape_gen_day;	/* Day of master tape generation	*/
  long	master_tape_gen_mon;	/* Month of master tape generation	*/
  long	master_tape_gen_yr;	/* Year of master tape generation	*/
  long	tape_seq_number;	/* Tape sequence number			*/
  long	mission_number;		/* Mission number			*/
  long	row_number;		/* Site number--Row number from WRS	*/
  long	path_number;		/* Line number--Path number from WRS	*/
  long  run_number;		/* Run number--Cycle number for LANDSAT	*/
  long	orbit_number;		/* Orbit number				*/
  long	first_scan_millisec;	/* Time of first scan--tenths of millisecond*/
  long	first_scan_seconds;	/* Time of first scan--seconds		*/
  long	first_scan_minutes;	/* Time of first scan--minutes		*/
  long	first_scan_hours;	/* Time of first scan--hours		*/
  long	first_scan_days;	/* Time of first scan--days		*/
  long	first_scan_month_num;	/* Time of first scan--month number	*/
  long	first_scan_year;	/* Time of first scan--year (last two digits) */
  long	channels_active[64];	/* Channels active in this job--1 if active
				   0 if otherwise			*/
  long	processing_flag;	/* Processing flag--0 if raw data, 1 if radio-
				   metrically or geometrically processed */
  long	num_channels_on_tape;	/* Number of channels on this tape	*/
  long	num_bits_in_picture;	/* Number of bits in a picture element	*/
  long	byt_loctn_strt_of_video_data;	/* Byte location of start of effective
				   video data within a scan line	*/
  long	byt_loctn_strt_of_calibratn;	/* Byte location of start of first
				   calibration area within scan line.	*/
  long  num_video_elemnts_per_scan;	/* Number of video elements per scan
				   within a single channel		*/
  long	num_calibratn_elemnts_scan;	/* Number of calibration elements in
				   the first calibration area within the scan
				   in a single channel			*/
  long	record_size_in_bytes;	/* Physical record size in bytes	*/
  long	channels_per_record;	/* Number of channels per physical record */
  long	num_records_per_scan;	/* Number of physical records per scan per
				   channel				*/
  long	num_records_data_set;	/* Number of physical records to make a
				   complete data set			*/
  long	byts_ancill_data_per_data_set;	/* Number of bytes of ancillary
				   data included in each data set	*/
  long	data_order_indicator;	/* Data order indicator where:
				      0 = video ordered by channel in a record
				      1 = video ordered by pixel in a record */
  long	start_pixel_number;	/* Number of the first pixel per scan on this
				   tape referenced to original image	*/
  long	stop_pixel_number;	/* Number of the last pixel per scan on this
				   tape referenced to original image	*/

  long	a0_coefficient[64]; 	/* Ao coefficient			*/
  long	e0_exponent[64];	/* Eo exponent				*/
  long	a1_coefficient[64];	/* A1 coefficient			*/
  long	e1_exponent[64];	/* E1 coefficient			*/

  long	color_code_data_chanl[64];	/* Colour code data per channel where
						0 = not active
						1 = red
						2 = green
						3 = blue		*/
  long	picture_field_width;	/* Picture element field width--smallest
				   multiple of 8 which equals or is greater
				   than the number of bits in a pixel	*/
  long	pixel_registration;	/* Number of bits that LSB of pixel is
				   displaced from LSB of pixel field 	*/
  long	computer_word_size;	/* Word size of generating computer or smallest
				   quanitity, in bits, that machine can write
				   on tape				*/
  char	lowr_lmt_wavelength[64][9];	/* Lower limits of wavelength range
				   per channel in nanometres            */
  char	uppr_lmt_wavelength[64][9];	/* Upper limits of wavelength range
				   per channel in nanometres		*/
  long	data_sets_per_record;	/* Number of data sets per physical record */
  long	byte_address_start_second_calibratn;	/* Byte address of start of
				   second calibration within scan	*/
  long	calibration_elements_second_area;	/* Number of calibration
				   elements in the second calibration area
				   within the scan in a single channel	*/
  long	calibration_source_indicator;	/* Calibration source indicator	*/
  long	file_skip_flag;		/* File skip flag.  Skip the next "n" end of
				   file marks, before the video data	*/
  long	bands_in_first_record;	/* Number of bands in the first physical
				   record of the data set		*/
  long	bytes_per_scan;		/* Total number of bytes per scan per channel */
  long	pixel_skip_factor;	/* The number to be added to the number of the
				   last pixel processed to yield the number of
				   the next pixel processed		*/
  long	scan_skip_factor;	/* The quantity to be added to the number of
				   the last scan processed to yeild the number
				   of the next scan to be processed	*/
  char	general_information[4][9];	/* General information		*/
  char	sun_elevation_angle_centre_deg[9];	/* Sun elevation angle in 
				   degrees at centre			*/
  char	sun_azimuth_angle_centre_deg[9];	/* Sun azimuth angle in 
				   degrees at centre			*/
  long	number_auxiliary_files;		/* Number of Auxiliary files between the
				   Universal Header file and the data file */
  long	start_scan_line;	/* The number of the first scan line on this 
				   tape referenced to the original image */
  long	stop_scan_line;		/* The number of the last scan line on this 
				   tape referenced to the original image */
  long	scan_lins_per_frame_per_thous;	/* "n" Thousnad scan lines per frame */
/* Image Annotation Data
   ---------------------*/
  char	ann_satellite_station_id[2];	/* Satellite station identification
				            P = PASS	S = SCSS	*/
  char	ann_path[4];		/* Path					*/
  char  ann_row[4];		/* Row: Canadian day rows 1 - 40
					Canadian night rows 202 - 242 	*/
  char	ann_day[3];		/* Day					*/
  char	ann_month[3];		/* Month				*/
  char	ann_year[3];		/* Year					*/
  char	ann_mission[2];		/* Mission				*/
  char	ann_days_since_launch[5];	/* Days since launch		*/
  char	ann_hours[3];		/* Hours				*/
  char	ann_minutes[3];		/* Minutes				*/
  char	ann_tens_of_seconds[2];	/* Tens of seconds			*/

  long	nominal_altitude_met;	/* Nominal altitude in metres		*/
  long	nominal_ground_speed_met_sec;	/* Nominal ground speed in met/sec */
  long	scan_type;		/* Scan type where 0 = linear
					   	   1 = smoothed
						0 = uncorrected for LLC
						1 = corrected for LLC	*/
  long	arc_angle_deg;		/* Angle of arc in degrees		*/
  long	camera;			/* Camera: 0 = 70mm
					   1 = 127mm
					   2 = 240mm
					 255 = no camera		*/
  long	input_device;		/* Input device: 0 = CCT
						 1 = High density 	*/
  long	truncation;		/* Truncation: 0 = 2 low order bits
					       1 = 2 high order bits
					       2 = no truncation 	*/
  long	channels_requested[64];	/* Channels requested where:
				   1 = requested       0 = otherwise 	*/
  long	processing_mode;	/* Processing Mode where:
				   0 = serially	   1 = concurrently 	*/
  long	color_select;		/* Color select where:  0 = no colour
			   	   1 = assigned colour  2 = false colour */
  long  image_format;		/* Image format where:
					0 = single image 1 = overlay image
					2 = abut image 3 = offset image	*/
  long	repeat_pixels_per_scan;	/* Repeat of pixels per scan		*/
  long	repeat_of_scan;		/* Repeat of scan			*/
  long	partial_scan_strt_pixel;  /* Partial scan--starting pixel number*/
  long	partial_scan_end_pixel;	/* Partial scan--ending pixel number	*/
  long	sensor_scan_rate;	/* Sensor scan rate in scans per second	*/
  long	pixel_size;		/* Pixel size				*/
  long	angle_of_drift_int;	/* Angle of drift in degrees--integer	*/
  long	angle_of_drift_fract;	/* Angle of drift in degrees--fraction	*/

  };
struct JSCLH
  {
  char  productn_syst_id_code[11]; /* Production system identification code.
				   This is of the form (100 x generating
				   centre code) + (duplicating centre code)
				   1 = NASA, GSFC	2 = NASA, GSFC
				   3 = EDC, U.S. Dept. of the Int.
				   4 = CCRS, Ottawa	5 = CCRS, PASS
				   6 = CCRS, SCSS  */
  char	id_code_description[31]; /* Identification code description	*/
  char	mission_id[11];		/* Mission identification		*/
  char	title1[31];		/* LANDSAT				*/
  char	day_num_since_launch[11]; /* Day number since launch--2nd - 5th
				     digits of frame id		 	*/
  char	title2[31];		/* DAY NUMBER SINCE LAUNCH		*/
  char	orbit_number[11];	/* Orbit number				*/
  char	title3[31];		/* ORBIT NUMBER				*/
  char	frame_centre_time[11];	/* Frame centre time--final 5 digits of frame
				   id corresponding to hours, min, & tens of
				   seconds GMT where centre point was imaged */
  char	title4[31];		/* FRAME ID (HHMMT)			*/
  char	frame_centre_lat[11];	/* Frame centre latitude in degrees and
				   minutes (DDDMM)			*/
  char	title5[31];		/* CENTRE LATITUDE (DDDMM)		*/
  char	frame_centre_long[11];	/* Frame centre longitude in deg and min */
  char 	title6[31];		/* CENTRE LONGITUDE (DDDMM)		*/
  char	utm_zone_number[11];	/* UTM zone number			*/
  char	title7[31];		/* UTM ZONE NUMBER			*/
  char	path_number[11];	/* Path number				*/
  char	title8[31];		/* PATH NUMBER				*/
  char	row_number[11];		/* Row number				*/
  char	title9[31];		/* ROW NUMBER				*/
  char	cycle_number[11];	/* Cycle number--18-day orbital cycle	*/
  char	title10[31];		/* CYCLE NUMBER				*/

  char	data_acquisition_ddmmyy[11];	/* Date of data acquisition in the
				   form:  Day, Month number, Year (DDMMYY) */
  char	title11[16];		/* FRAME IMAGED				*/
  char	data_acquisition_ddmmmyy[16];	/* Date of data acquisition in the
				   form:  Day, Month name, Year (DDMMMYY) */

  char	master_cct_producd_ddmmyy[11];	/* Date master CCT produced in the
				   form:  Day, Month number, Year (DDMMMYY) */
  char	title12[21];		/* MASTER GENERATED			*/
  char	master_cct_producd_ddmmmyy[11]; /* Date master CCT produced in the
				   form:  Day, Month name, Year (DDMMMYY) */

  char	date_copy_producd_ddmmyy[11];	/* Date copy produced (DDMMYY)	*/
  char	title13[16];		/* COPY PRODUCED			*/
  char	date_copy_producd_ddmmmyy[16];	/* Date copy produced (DDMMMYY)	*/

  char	heading[11];		/* Heading				*/
  char	title14[31];		/* SATELLITE HEADING (DDDMM)		*/
  char	seq_number_of_tape[11];	/* Sequential number of tape		*/
  char	title15[31];		/* SEQUENTIAL NUMBER OF TAPE		*/
  char	tape_start_time[11];	/* Time of tape start since beginning of 
				   standard frame--in microseconds	*/
  char	title16[31];		/* TAPE START TIME			*/
  char	title17[31];		/* PROCESS FLAGS			*/
  char	radiometric_calibration[2];  /* 0 = no radiometric cal
					1 = radiometric striping removal only
					2 = CAL2 radiometric calibration
					3 = CAL3 radiometric calibration  */
  char	radiometric_levels[2];	/* 9 = 64 levels
				   0 = 128 levels
				   1 = 256 levels			*/
  char	earth_rotatn_correctn[2];    /* 0 = no earth rot corr
					1 = earth rot corr		*/
  char	mirror_scan_velocity_corr[2];/* 0 = no mirror scan vel corr
					1 = mirror scan velocity corr
					2 = panoramic distortion and earth
					    curvature correction
					3 = mirror scan velocity, panoramic
					    distortion, and earth curv corr */
  char	data_representation[2];	/* 0 = linear (decompressed) representation
				   1 = logarithmic representation	*/
  char	line_length_correction[2];  /*	0 = no line length correction
				 	1 = line length correction	*/
  char	character_code_for_lhdr[2];	/* Character code of LANDSAT Header file
					   0 = EBCDIC
					   1 = ASCII			*/
	/* Geometric Transformation Record
   	   ------------------------------*/
  char	gt_utm_zone_number[21];	/* UTM zone number of greatest coverage	*/
  char	gt_utm_coordinate_northing[21];	/* UTM coordinate (northing) in metres*/
  char	gt_utm_coordinate_easting[21];	/* UTM coordinate (easting) in metres*/
  char	gt_orientation[21];	/* Orientation				*/
  };
	/* Radiometric Transformation Record
	   --------------------------------*/
struct JSCRT
  {
  char	radiometric_transformation_detector_1[64][5];
				/* records for detector	1	*/
  char	radiometric_transformation_detector_2[64][5];
				/* records for detector	2	*/
  char	radiometric_transformation_detector_3[64][5];
				/* records for detector	3	*/
  char	radiometric_transformation_detector_4[64][5];
				/* records for detector	4	*/
  char	radiometric_transformation_detector_5[64][5];
				/* records for detector	5	*/
  char	radiometric_transformation_detector_6[64][5];
				/* records for detector	6	*/
  };
   struct JSCIN 
	{
  	short anc_start;
  	short anc_stop;
  	long  anc_line_length;
  	long  anc_days;
  	long  anc_hours;
  	long  anc_minutes;
  	long  anc_seconds;
  	long  anc_hundreth_seconds;
  	long  anc_sensor;
  	long  anc_cal_band_4[6];
  	long  anc_cal_band_5[6];
  	long  anc_cal_band_6[6];
  	long  anc_cal_band_7[6];
  	long  anc_cal_thermal[6];
	};

#endif
