#ifndef SIAT_H
#define SIAT_H

/**************************************************************************** 
*   stucture for CCT-AT and CCT-PT tapes which are read in by ccttipsp      * 
*									    *
*   See "Generation and Physical Characteristics of the Landsat -1, -2  and *
*   -3 Mss Computer Compatible Tapes" for more information on specific      *
*   members of the structure		    				    *
*									    *
****************************************************************************/

struct SIAT 
   {
   /* tape related information
   -------------------------*/
   char  tape_number[13];		/* tape number			*/ 
   char  tape_sequence[3];		/* logical sequence number	*/
   char  tape_date[11];			/* date of tape preparation	*/
   char  tape_siat_number[15];		/* siat number and seq number	*/
   char  tape_rbv_number[9];		/* rbv tape number		*/
   char  tape_mss_number[9];		/* mss tape number		*/
   short tape_number_data_files;	/* number of files on logical siat*/
   short tape_number_rbv_vtc;		/* number of rbv/vtc		*/	
   short tape_number_mss_vtc;		/* number of mss/vtc		*/
   short tape_number_rbv_tfc;		/* number of rbv/tfc		*/
   short tape_number_mss_tfc;		/* number of mss/tfc		*/
   char  tape_rvb_scene_id[64][13];	/* 1st-64th rbv scene id's	*/
   char  tape_mss_scene_id[64][13];	/* 1st-64th mss scene id's	*/

   /* calibration modifiers record
   -----------------------------*/
   short calibration_id;		/* satellite id			*/
   short calibration_days;		/* days since launch		*/
   float calibration_band4_low_compressed_multiplicative_sensor[6];
   float calibration_band4_low_compressed_additive_sensor[6];
   float calibration_band4_low_linear_multiplicative_sensor[6];
   float calibration_band4_low_linear_additive_sensor[6];
   float calibration_band4_high_compressed_multiplicative_sensor[6];
   float calibration_band4_high_compressed_additive_sensor[6];
   float calibration_band4_high_linear_multiplicative_sensor[6];
   float calibration_band4_high_linear_additive_sensor[6];
   float calibration_band5_low_compressed_multiplicative_sensor[6];
   float calibration_band5_low_compressed_additive_sensor[6];
   float calibration_band5_low_linear_multiplicative_sensor[6];
   float calibration_band5_low_linear_additive_sensor[6];
   float calibration_band5_high_compressed_multiplicative_sensor[6];
   float calibration_band5_high_compressed_additive_sensor[6];
   float calibration_band5_high_linear_multiplicative_sensor[6];
   float calibration_band5_high_linear_additive_sensor[6];
   float calibration_band6_low_compressed_multiplicative_sensor[6];
   float calibration_band6_low_compressed_additive_sensor[6];
   float calibration_band6_low_linear_multiplicative_sensor[6];
   float calibration_band6_low_linear_additive_sensor[6];
   float calibration_band7_low_linear_multiplicative_sensor[6];
   float calibration_band7_low_linear_additive_sensor[6];


   /* processing instruction record
   ------------------------------*/
   short processing_scenes_remaining_rbv_vfc;	/* number of scenes remaining */
   short processing_scenes_remaining_mss_vfc;	/* number of scenes remaining */
   short processing_scenes_remaining_rbv_vtc;	/* number of scenes remaining */
   short processing_scenes_remaining_mss_vtc;	/* number of scenes remaining */
   char  processing_scene_id[13];		/* scene id		      */
   long  processing_mission_number;		/* mission number	      */
   long  processing_day_number_from_launch;	/* days from launch	      */
   long  processing_hours;			/* hours of the day	      */
   long  processing_minutes;			/* minute		      */
   long  processing_seconds;			/* second		      */
   long  processing_mission_number2;		/* mission number	      */
   long  processing_day_number_from_launch2;	/* days from launch	      */
   long  processing_hours2;			/* hours		      */
   long  processing_minutes2;			/* minutes		      */
   long  processing_seconds2;			/* seconds		      */
   short processing_code_slat;			/* processing code slat	      */
   short processing_code_mss;			/* processing code mss	      */
   short processing_code_psp;			/* processing code PS project */
   char  processing_preceding_id[11];		/* preceding closeset RCI ID  */
   char  processing_succeedingid[11];		/* succeding closeset RCI ID  */
   char  processing_band_information[8][11];	/* array of band information  */
   char  processing_instructions[73];		/* special prcessing instruct */



   /* use data from spacecraft performance record
   --------------------------------------------*/
   char  spacecraft_rbv1_transmission[9];	/* mode of transmission	*/
   char  spacecraft_rbv1_duration[3];		/* exposure duration	*/
   char  spacecraft_rbv1_aperture[3];		/* aperture correction	*/
   char  spacecraft_rbv2_transmission[9];	/* mode of transmission	*/
   char  spacecraft_rbv2_duration[3];		/* exposure duration	*/
   char  spacecraft_rbv2_aperture[3];		/* aperture correction	*/
   char  spacecraft_rbv3_transmission[9];	/* mode of transmission	*/
   char  spacecraft_rbv3_duration[3];		/* exposure duration	*/
   char  spacecraft_rbv3_aperture[3];		/* aperture correction	*/
   char  spacecraft_rbv4_transmission[13];	/* mode of transmission	*/
   char  spacecraft_rbv5_transmission[13];	/* mode of transmission	*/
   char  spacecraft_rbv6_transmission[13];	/* mode of transmission	*/
   char  spacecraft_rbv7_transmission[13];	/* mode of transmission	*/
   char  spacecraft_rbv8_transmission[13];	/* mode of transmission	*/
   long  spacecraft_sensor_gain_band5;		/* mss sensor gain	*/
   long  spacecraft_sensor_gain_band4;		/* mss sensor gain	*/
   long  spacecraft_sensor_encoding_band6;	/* mss sensor encodeing	*/
   long  spacecraft_sensor_encoding_band5;	/* mss sensor encodeing	*/
   long  spacecraft_sensor_encoding_band4;	/* mss sensor encodeing	*/
   char  spacecraft_tape_id[9];			/* spdt tape id		*/
   char  spacecraft_sun_day[5];			/* mss sun cal day	*/
   short processing_mss_sun_sensor[24];		/* mss sun cal sensor	*/
   char  spacecraft_sun_day2[5];		/* mss sun cal day desired */
   char  spacecraft_sun_flag[5];		/* mss sun cal flag	*/
 
 
   /* use data from annotation record
   --------------------------------*/
   char  annotation_day[3];			/* day of exposure	*/
   char  annotation_month[4];			/* month of exposure	*/
   char  annotation_year[4];			/* year of exposure	*/
   char  annotation_latitude[7];		/* lat of format center	*/
   char  annotation_longitude[8];		/* long of format center*/
   char  annotation_direction[2];		/* descending or ascending*/
   char  annotation_path_row[8];		/* path/row number	*/ 
   char  annotation_nadir_latitude[7];		/* nadir latitude	*/
   char  annotation_nadir_longitude[9];		/* nadir longitude	*/
   char  annotation_sensor_code[2];		/* sensor code		*/
   char  annotation_designated_band[8];		/* designated band	*/
   char  annotation_transmission_mode[2];	/* transmission mode	*/
   char  annotation_sun_elevation[10];		/* sun elev of nadir	*/
   char  annotation_sun_azimuth[7];		/* sun azimuth of nadir	*/
   char  annotation_image_prcessing_codes[13];	/* image processing codes*/
   char  annotation_scene_id[14];		/* scene id		*/
   char  annotation_rbv[25];			/* ???????????????????? */
   char  annotation_mss_mode_site[5];		/* mode and acquistion site */
   char  annotation_heading[5];			/* heading of satellite */
   char  annotation_rev[7];			/* revision number	*/
   char  annotation_data_acquistion[5];		/* rbv data acquisition	*/
   char  annotation_orbit_type[3];		/* type of orbit data	*/
   char  annotation_rci_calibration[2];		/* calibration level	*/
   char  annotation_rbv_mode1[6];		/* transmition mode	*/
   char  annotation_rbv_shutter1[3];		/* shutter setting	*/
   char  annotation_rbv_mode2[7];		/* transmition mode	*/
   char  annotation_rbv_shutter2[3];		/* shutter setting	*/
   char  annotation_rbv_mode3[7];		/* transmition mode	*/
   char  annotation_rbv_shutter3[3];		/* shutter setting	*/
   char  annotation_mss_mode[6];		/* transmition mode	*/

 
   /* use data from rbv record
   -------------------------*/
   short rbv_spacecraft_ddd;		/* spacecraft day		*/
   short rbv_spacecraft_hh;		/* spacecraft hour 		*/
   short rbv_spacecraft_mm;		/* spacecraft minute		*/
   short rbv_spacecraft_ss;		/* spacecraft second		*/
   short rbv_spacecraft_cc;		/* spacecraft decimal seconds	*/
   short rbv_greenwich_ddd;		/* greenwich day		*/
   short rbv_greenwich_hh;		/* greenwich hour		*/
   short rbv_greenwich_mm;		/* greenwich minute		*/
   short rbv_greenwich_ss;		/* greenwich second		*/
   short rbv_greenwich_mmm;		/* greenwich decimal seconds	*/
   float rbv_altitude_change;		/* normalized altitude change	*/
   char  rbv_gmt_date[11];		/* date of exposure		*/
   char  rbv_gmt_time[9];		/* time of exposure		*/
   float rbv_latitude_format_center;	/* latitude of format center	*/
   float rbv_longitude_format_center;	/* longitude of format center	*/
   float rbv_latitude_nadir;		/* latitude of nadir		*/
   float rbv_longitude_nadir;		/* longitude of nadir		*/
   float rbv_altitude;			/* sapcecraft altitude		*/
   float rbv_gmt_exposure;		/* GMT of exposure		*/
   float rbv_flight_heading;		/* s/c flight path heading	*/
   float rbv_pitch;			/* pitch			*/
   float rbv_roll;			/* roll				*/
   float rbv_yaw;			/* yaw				*/


   /* use data from mss record
   -------------------------*/
   short mss_timesc_ddd;		/* spacecraft day		*/
   short mss_timesc_hh;			/* spacecraft hour		*/
   short mss_timesc_mm;			/* spacecraft minute		*/
   short mss_timesc_ss;			/* spacecraft second		*/
   short mss_timesc_cc;			/* spacecraft decimal second	*/
   short mss_gmtsc_ddd;			/* greenwich day		*/
   short mss_gmtsc_hh;			/* greenwich houre		*/
   short mss_gmtsc_mm;			/* greenwich minute		*/
   short mss_gmtsc_ss;			/* greenwich second		*/
   short mss_gmtsc_mmm;			/* greenwich decimal second	*/
   float mss_norm_alt_change_ic[9];	/* normalized altitude change	*/
   float mss_alt[9];			/* array of altitude values	*/
   float mss_vehicle_roll_ic;		/* roll at image center		*/
   float mss_vehicle_pitch_ic;		/* pitch at image center	*/
   float mss_vehicle_yaw_ic;		/* yaw at image center		*/
   float mss_roll[11];			/* array of roll values		*/
   float mss_pitch[11];			/* array of pitch values	*/
   float mss_yaw[11];			/* array of yaw values		*/
   float mss_skew;			/* image skew			*/
   float mss_norm_vel_change;		/* nomralized velocity change	*/
   float mss_mean_pitch;		/* mean pitch			*/
   float mss_mean_roll;			/* mean roll			*/
   float mss_mean_yaw;			/* mean yaw			*/
   float mss_mean_pitch_rate;		/* mean pitch rate		*/
   float mss_mean_roll_rate;		/* mean roll rate		*/
   float mss_mean_yaw_rate;		/* mean yaw rate		*/
   float mss_mean_alt;			/* mean altitude		*/
   float mss_mean_alt_rate;		/* mean altitude rate		*/
   long  mss_gmt_ms[11];		/* array of GMT time		*/
   float mss_latitude[11];		/* array of latitudes		*/ 
   float mss_longitude[11];		/* array of longitudes		*/
   float  mss_latlon_alt[11];		/* array of altitude		*/

   /* use data from image location record
   ------------------------------------*/
   float image_rbv_top_point[6];	     /* array of top tick marks      */
   char  image_rbv_top_point_ascii[6][9];    /* array of top ascii lat/long  */
   float image_rbv_left_point[6];	     /* array of left tick marks     */
   char  image_rbv_left_point_ascii[6][9];   /* array of left ascii lat/long */
   float image_rbv_right_point[6];	     /* array of right tick marks    */
   char  image_rbv_right_point_ascii[6][9];  /* array of right ascii lat/long*/
   float image_rbv_bottom_point[6];	     /* array of bottom tick marks   */
   char  image_rbv_bottom_point_ascii[6][9]; /* array of bottom ascii lat/long*/
   float image_mss_top_point[6];             /* array of top tick marks      */
   char  image_mss_top_point_ascii[6][9];    /* array of top ascii lat/long  */
   float image_mss_left_point[6];            /* array of left tick marks     */
   char  image_mss_left_point_ascii[6][9];   /* array of left ascii lat/long */
   float image_mss_right_point[6];           /* array of right tick marks    */
   char  image_mss_right_point_ascii[6][9];  /* array of right ascii lat/long*/
   float image_mss_bottom_point[6];          /* array of bottom tick marks   */
   char  image_mss_bottom_point_ascii[6][9]; /* array of bottom ascii lat/long*/

   };

#endif
