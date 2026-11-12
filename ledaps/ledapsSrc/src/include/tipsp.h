#ifndef TIPSP_H
#define TIPSP_H

/**************************************************************************** 
*   stucture for CCT-AT and CCT-PT tapes which are read in by ccttipsp      * 
*									    *
*   See "User's Guide for Landsat Thematic Mapper Computer-Compatible Tapes"*
*   for more information on specific members of the structure		    *
*									    *
****************************************************************************/

struct TIPSP
   {
   /* header related information
   ---------------------------*/
   float head_center_latitude;		/* full scene center latitude	*/
   float head_center_longitude;		/* full scene center longitude	*/
   float head_center_line;		/* full scene center line number*/
   float head_center_sample;		/* full scene center sample # 	*/
   float head_pquadrant_lat;		/* processed quadrant center lat*/
   float head_pquadrant_long;		/* processed quadrant center long*/
   float head_quadrant_line;		/* quadrant center line		*/
   float head_quadrant_sample;		/* quadrant center sample	*/

   long  head_average_samples;		/* average number of data pixels*/
   long  head_bands;			/* number of active bands	*/
   long  head_bands_present[7];		/* bands present		*/
   long  head_center_time_offset;	/* scene center time offset	*/
   long  head_detector_map[100];	/* detector substitution map	*/
   long  head_hquadrant_overlap;	/* horizontal quadrant overlap	*/
   long  head_lines;			/* number of scene lines in quadrant*/
   long  head_orbit_number;		/* orbit number			*/
   long  head_vquadrant_overlap;	/* vertical quadrant overlap	*/
   long  head_wrs_cycle;		/* WRS cycle			*/

   char  head_center_time[33];		/* full scene center time	*/
   char  head_dectector_smooth_array[101]; /*dectector smoothing array	*/
   char  head_geometric_correction[17];	/* correction designator    	*/
   char  head_interleave[17];		/* interleaving indicator   	*/
   char  head_map_ancillary_records[17];/* number of map ancillary 
					   records		    	*/
   char  head_map_projection[17];	/* projection code	    	*/
   char  head_mission_id[17];		/* mission id			*/
   char  head_orbit_direction[17];	/* orbit direction		*/
   char  head_parent_id[17];		/* parent scene id		*/
   char  head_product_class[17];	/* product class id	    	*/
   char  head_product_id[17];   	/* product id			*/
   char  head_radiometric_ancillary_records[17];/* number of radiometric 
						   ancillary records	    */
   char  head_radiometric_correction[17];	/* radiometric correction   */
   char  head_radiometric_resolution[17];	/* radiometric resolution   */
   char  head_resampling_designator[17];	/* resampling designator    */	
   char  head_scene_id[17];		/* scene identification		*/
   char  head_sensor_id[17];		/* sensor id			*/
   char  head_wrs_designator[17];	/* WRS designator		*/


   /* annotation related information
   -------------------------------*/
   char  annot_acquistion_date[9];		/* acquisition date	      */
   char  annot_agency[14];			/* agency, project id	      */
   char  annot_band_id_code[11];		/* sensor band id code	      */
   char  annot_center_latitude_longitude[18];	/* center latitude, longitude */
   char  annot_ephemeris[2];			/* type of emphemeris used    */
   char  annot_geometric_correction_header[2];  /* geometric correction	      */
   char  annot_processing_procedure[2];		/* processing procedure	      */
   char  annot_projection[2];			/* map projection code	      */
   char  annot_resampling_algorithm[2];		/* resampling algorithm type  */
   char  annot_scene_identification[26];	/* scene id number	      */
   char  annot_sun_angle[15];			/* sun angle	  	      */
   char  annot_wrs_latitude_longitude[18];	/* wrs center lat, long	      */
   char  annot_wrs_path[13];			/* WRS path and row identifier*/
 
   /* tick mark related information
   ------------------------------*/
   short annot_top_tick[8];		/* x, or y pixel location	*/
   short annot_left_tick[8];		/* x, or y pixel location	*/
   short annot_right_tick[8];		/* x, or y pixel location	*/
   short annot_bottom_tick[8];		/* x, or y pixel location	*/

   char  annot_bottom_tick_ascii[8][8];	/* coord in meters or deg - min */
   char  annot_left_tick_ascii[8][8];	/* coord in meters or deg - min */
   char  annot_right_tick_ascii[8][8];	/* coord in meters or deg - min */
   char  annot_top_tick_ascii[8][8];	/* coord in meters or deg - min */

/* Map projection related information
-----------------------------------*/
   /* Input scene related information
   --------------------------------*/
   long  map_nominal_pixels;		/* nominal pixels per line	*/
   long  map_nominal_lines;		/* nominal number lines per scene*/
   float map_nominal_interpixel;	/* nominal interpixel distance	*/
   float map_nominal_interline;		/* nominal interline distance	*/
   float map_skew;			/* image skew at scene center	*/

   /* UTM/PS related information
   ---------------------------*/
   char  map_utm_datum[7];		/* zone number			*/
   long  map_utm_zone;			/* datum number			*/
   float map_utm_wrs_northing;		/* northing of wrs center	*/
   float map_utm_wrs_easting;		/* easting of wrs center	*/
   float map_utm_center_northing;	/* northing of full scene center*/
   float map_utm_center_easting;	/* easting of full scene center	*/
   float map_utm_vertical_offset;	/* verticle offset of wrs center*/
   float map_utm_horizontal_offset;	/* horizontaloffset of wrs center*/
   float map_utm_rotation;		/* display rotation angle	*/


   /* SOM related information
   ------------------------*/
   float map_som_wrs_horizontal;	/* horizontal coordinate of WRS center*/
   float map_som_wrs_vertical;		/* vertical coordinate of WRS center  */
   float map_som_center_horizontal;	/* horizontal coord of full scene cent*/
   float map_som_center_vertical;	/* vertical coord of full scene center*/
   float map_som_vertical_offset;	/* verticle offset of wrs center      */
   float map_som_horizontal_offset;	/* horizontal offset of wrs center    */
   float map_som_rotation;		/* display rotation angle	      */

   /* Processed Quadrant related information
   ---------------------------------------*/
   float map_pixels;			/* pixels per line		*/
   float map_lines;			/* lines per quadrant		*/
   float map_corrected_interpixel;	/* corrected interpixel distance*/
   float map_corrected_interline;	/* corrected intersample distance*/
   float map_processed_zone;		/* zone number			*/
   long  map_center_line_number;	/* line number at wrs center	*/
   long  map_center_sample_number;	/* sample number at wrs center	*/
   float map_rotation;			/* display rotation angle	*/
   float map_orbit_inclination;		/* satellite orbital inclination*/
   float map_ascending_longitude;	/* ascending node longitude	*/
   float map_altitude;			/* satellite altitude		*/
   float map_ground_speed;		/* ground speed			*/
   float map_satellite_heading;		/* satellite heading at center  */
   float map_cross_track_field_of_view;	/* cross track field of view	*/
   float map_scan_rate;			/* sensor scan rate 		*/
   float map_sampling_rate;		/* sensor sampling rate		*/
   float map_sun_elevation;		/* sun elevation		*/
   float map_sun_azimuth;		/* sun azimuth  		*/

/* Radiometric calibration related information
--------------------------------------------*/
   /* Read data from the Radiometric calibration record which is in the header
   -------------------------------------------------------------------------*/
   char   radio_band[5];		/* band number			*/
   char   radio_lower_reflectance[5];	/* lower reflectance limit	*/
   char   radio_upper_reflectance[5];	/* upper reflectance limit	*/
   double radio_offset;			/* offset coefficient		*/
   double radio_gain;			/* gain coefficient		*/


/* Band quality related information
---------------------------------*/
   char band_quality[2];		/* overall band quality code	*/

   /* radiometric correction data
   ----------------------------*/
   char  band_calibration[2];			/* calibration mode	*/
   char  band_calibration_lamp_use[3];		/* calibration lamp use */
   char  band_lamp_mode[2];			/* internal lamp mode	*/
   char  band_selected_histogram_processing[2];	/* histogram flag	*/

   long  band_number_of_scenes;		     /* number of scenes on the HDT-AT*/
   long  band_scene_sequence;		     /* scene sequence number 	*/
   long  band_number_radiometric_corrections;/* # radiometric corrections (RC)*/
   long  band_number_scans_per_rc;	     /* # scans per RC sections	*/
   long  band_number_subsegments_per_rc;	/* segments per RC	*/

   short band_number_hist_iterations;		/* # itterations allowed*/

   float band_minimum_radiance;			/* min radiance value	*/
   float band_maximum_radiance;			/* max radiance value	*/

   /* Mirror scan correction information
   -----------------------------------*/
   float band_nominal_line_forward;		/* Nominal line length	*/
   float band_nominal_line_backward;		/* Nominal line length	*/
   float band_maximum_line_error_forward;	/* maximum line error	*/
   float band_maximum_line_error_backward;	/* maximum line error	*/
   float band_minimum_line_error_forward;	/* minimum line error	*/
   float band_minimum_line_error_backward;	/* minimum line error	*/
   float band_mean_line_error_forward;		/* mean line error forward */
   float band_mean_line_error_backward;		/* mean line error backward */
   float band_sigma_line_error_forward;		/* sigam line error forward */
   float band_sigma_line_error_backward;	/* sigam line error backward */
   float band_nominal_first_line_error_forward; /* nominal first 1/2 length */ 
   float band_nominal_first_line_error_backward;/* nominal first 1/2 length */
   float band_maximum_first_line_error_forward; /* maximum first 1/2 length */  
   float band_maximum_first_line_error_backward;/* maximum first 1/2 length */
   float band_minimum_first_line_error_forward; /* minimum first 1/2 length */ 
   float band_minimum_first_line_error_backward;/* minimum first 1/2 length */ 
   float band_mean_first_line_error_forward;	/* mean first 1/2 length    */
   float band_mean_first_line_error_backward;	/* mean first 1/2 length    */
   float band_sigma_first_line_error_forward;	/* sigma first 1/2 length   */
   float band_sigma_first_line_error_backward;	/* sigma first 1/2 length   */
   float band_nom_second_line_error_forward;	/* nominal second 1/2 length*/
   float band_nom_second_line_error_backward;	/* nominal second 1/2 length*/
   float band_max_second_line_error_forward;	/* maximum second 1/2 length*/
   float band_max_second_line_error_backward;	/* maximum second 1/2 length*/
   float band_min_second_line_error_forward;	/* minimum second 1/2 length*/
   float band_min_second_line_error_backward;	/* minimum second 1/2 length*/
   float band_mean_second_line_error_forward;	/* mean second 1/2 length   */
   float band_mean_second_line_error_backward;	/* mean second 1/2 length   */
   float band_sigma_second_line_error_forward;	/* sigma second 1/2 length  */
   float band_sigma_second_line_error_backward;	/* sigma second 1/2 length  */
   float band_max_nonlin_line_error_forward;	/* maximum nonlinearity	    */
   float band_max_nonlin_line_error_backward;	/* maximum nonlinearity	    */
   float band_min_nonlin_line_error_forward;	/* minimum nonlinearity	    */
   float band_min_nonlin_line_error_backward;	/* minimum nonlinearity	    */
   float band_mean_nonlin_line_error_forward;	/* mean nonlinearity	    */
   float band_mean_nonlin_line_error_backward;	/* mean nonlinearity	    */
   float band_sigma_nonlin_line_error_forward;	/* sigma nonlinearity	    */
   float band_sigma_nonlin_line_error_backward;	/* sigma nonlinearity	    */

   /* Control point historical information
   -------------------------------------*/
   char  band_quality_codes[15];		/* quality codes for refernce*/ 
   char  band_scene_id[13];			/* reference scene id	    */

   float band_number_scences_in_cp;		/* number of scenes in cp   */
   float band_sequence_number_in_cp;		/* sequence number of cp    */
   float band_number_geodetic_points_in_cp;	/* number of geodetic points*/
   float band_number_geodetic_points_in_reference;/* number of geodetic points*/
   float band_average_autocorrelation;		/* average initial autocor  */
   float band_average_initial_peak;		/* curvature from reference */
   float band_elipse_along_track_interval_error; /* ellipse error	    */ 
   float band_elipse_across_track_interval_error; /* ellipse error	    */ 
   float band_elipse_along_track_ref_error;	/* ellipse error	    */
   float band_elipse_across_track_ref_error;	/* ellipse error	    */
   float band_elipse_average_success;		/* previous registration suc*/

   /* Current control point information
   ----------------------------------*/
   long  band_number_scenes_interval;		/* scenes in interval	      */
   long  band_number_cp_interval;		/* # cp in interval	      */
   long  band_number_cp_scene;			/* # cp in scene	      */
   long  band_number_cp_previous_scenes;	/* # cp in previous scene     */
   long  band_number_geodetic_cp;		/* # geodetic cp	      */
   long  band_number_cp_correlations_interval;	/* # geodetic geometric 
						     corrections per cp       */
   long  band_number_cp_rejected_interval;	/* # cp rejected per int      */
   long  band_number_cp_rejected_prefiltering_interval;/* # cp rejected in
						prefiltering operations	      */
   long  band_number_cp_rejected_filtering_interval;/* # cp rejected in
						   filtering operations	      */
   long  band_number_cp_correlations_scene;	/* # cp correlations attempted*/
   long  band_number_cp_rejected_scene;		/* # cp rejected scene	      */
   long  band_number_cp_rejected_prefiltering_scene;/* # cp rejected in
						prefiltering operations	      */
   long  band_number_cp_rejected_filtering_scene;/* # cp rejected in
						   filtering operations	      */
   long  band_sequence_number_in_interval;	/* sequence in interval	      */

   char  band_interval_correction[3];		/* inveterval correction type */
   char  band_scene_correction[3];		/* scene correction type      */
   char  band_cp_band_mode[3];			/* cp band mode		      */
   char  band_cp_id[25][16];			/* cp identifications	      */
   char  band_cp_interval_mode[3];		/* cp interval mode	      */

   float band_average_cp_curve;			/* ave cp correlation curve   */
   float band_average_cp_peak;			/* ave cp correlation peak    */
   float band_cp_line[25];			/* cp line number	      */
   float band_cp_sample[25];			/* cp sample number	      */

   /* Geometric correction information
   ---------------------------------*/
   char  band_geometric_quality[2];		/* overall geometric quality  */

   float band_error_matrix_along_track_sigma;	/* measurement error	      */
   float band_error_matrix_across_track_sigma;	/* measurement error	      */
   float band_geometric_along_track_interval;	/* RMS geometric modeling err */
   float band_geometric_across_track_interval;	/* RMS geometric modeling err */
   float band_geometric_along_track_scene;	/* RMS geometric modeling err */
   float band_geometric_across_track_scene;	/* RMS geometric modeling err */
   float band_horizontal_strip_maximum;		/* horizontal strip informatio*/
   float band_horizontal_strip_minimum;		/* horizontal strip informatio*/
   float band_horizontal_strip_mean;		/* horizontal strip informatio*/
   float band_horizontal_strip_maximum_rotation;/* horizontal strip informatio*/
   float band_scan_gap_maximum;			/* scan gap size	      */
   float band_scan_gap_minimum;			/* scan gap size	      */
   float band_scan_gap_mean;			/* scan gap size	      */
   float band_scan_skew_maximum;		/* scan gap skew	      */
   float band_scan_skew_minimum;		/* scan gap skew	      */
   float band_scan_skew_mean;			/* scan gap skew	      */
   float band_state_vector_radial_position;	/* state vector information   */
   float band_state_vector_along_track_position;/* state vector information   */
   float band_state_vector_across_track_position;/* state vector information  */
   float band_state_vector_radial_rate;	        /* state vector information   */
   float band_state_vector_along_track_rate;	/* state vector information   */
   float band_state_vector_across_track_rate;	/* state vector information   */
   float band_state_vector_angle_roll;	        /* state vector information   */
   float band_state_vector_angle_pitch;	        /* state vector information   */
   float band_state_vector_angle_yaw;	        /* state vector information   */
   float band_state_vector_roll_angular_rate;	/* state vector information   */
   float band_state_vector_pitch_angular_rate;	/* state vector information   */
   float band_state_vector_yaw_angular_rate;	/* state vector information   */
   float band_state_error_cov_matrix[12][12];	/* covariance matrix	      */
   float band_value_for_e;			/* value for E		      */

   long  band_horizontal_strip_exceeding_maximum;/*horizontal strip informatio*/
   long  band_horizontal_strip_exceeding_minimum;/*horizontal strip informatio*/
   long  band_scan_gap_exceeding_maximum;	/* scan gap size	      */
   long  band_scan_gap_exceeding_minimum;	/* scan gap size	      */
   long  band_scan_skew_exceeding_maximum;	/* scan gap skew	      */
   long  band_scan_skew_exceeding_minimum;	/* scan gap skew	      */

   /* trailer information
   --------------------*/
   char  trail_number_scan_lines[5];		/* number scan lines in HDT-AT*/
   char  trail_quality_indicator_summary[21];	/* summary counts 	      */
   char  trail_number_major_sync_losses[5];	/* major sync losses	      */
   char  trail_number_minor_sync_losses[5];	/* minor sync losses	      */
   char  trail_number_minor_sync_errors[5];	/* minor sync errors	      */
   char  trail_number_bit_slips[5];		/* number of bit slips	      */
   char  trail_number_time_code_substitutions[5];/* number of substitutions   */
   char  trail_quality_summary[157];		/* quality summary	      */
   };

#endif
