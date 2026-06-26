#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "lndsr.h"
#include "sixs_runs.h"

struct etm_spectral_function_t {
	int nbvals[SIXS_NB_BANDS];
	float wlinf[SIXS_NB_BANDS];
	float wlsup[SIXS_NB_BANDS];
	float response[SIXS_NB_BANDS][155];
} etm_spectral_function_t;

#define SIXS_STR_LENGTH 27  /* length of 6s report strings */

int create_6S_tables(sixs_tables_t *sixs_tables, Input_meta_t *meta) {
	char cmd[128],sixs_cmd_filename[1024],sixs_out_filename[1024],line_in[256];
    /* char tmp_file[1024], cmd_string[1024]; */
	int i,j,k;
	FILE *fd;
	float tgoz,tgco2,tgo2,tgno2,tgch4,tgco;
	int tm_band[SIXS_NB_BANDS]={25,26,27,28,29,30};
    char short_name[1024];
    char local_granule_id[1024];
    char acq_date_string[MAX_DATE_LEN + 1];
    const char *sat_names[SAT_MAX] = {"1", "2", "3", "4", "5", "7"};
    const char *inst_names[INST_MAX] = {"M", "T", "E"};
    const char *wrs_names[WRS_MAX] = {"1", "2"};
	
	struct etm_spectral_function_t etm_spectral_function = {
		{54,61,65,81,131,155},
		{0.420,0.500,0.580,0.730,1.501,2.0},
		{0.550,0.650,0.740,0.930,1.825,2.386},
		{
			{0.000,0.000,0.000,0.000,0.000,0.000,0.016,0.071,0.287,0.666,
             0.792,0.857,0.839,0.806,0.779,0.846,0.901,0.900,0.890,0.851,
             0.875,0.893,0.884,0.930,0.958,0.954,0.980,0.975,0.965,0.962,
             0.995,0.990,0.990,0.979,0.983,0.969,0.960,0.768,0.293,0.054,
             0.009,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
             0.000,0.000,0.000,0.000},
			{0.001,0.002,0.003,0.012,0.026,0.074,0.174,0.348,0.552,0.696,
             0.759,0.785,0.822,0.870,0.905,0.929,0.947,0.952,0.952,0.951,
             0.953,0.950,0.954,0.967,0.959,0.941,0.933,0.938,0.951,0.956,
             0.955,0.956,0.973,0.992,1.000,0.976,0.942,0.930,0.912,0.799,
             0.574,0.340,0.185,0.105,0.062,0.038,0.021,0.011,0.005,0.002,
             0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
             0.000},
			{0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
             0.000,0.000,0.000,0.000,0.000,0.001,0.002,0.010,0.047,0.174,
             0.419,0.731,0.921,0.942,0.937,0.937,0.949,0.965,0.973,0.970,
             0.958,0.955,0.962,0.980,0.993,0.998,1.000,0.995,0.992,0.988,
             0.977,0.954,0.932,0.880,0.729,0.444,0.183,0.066,0.025,0.012,
             0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
             0.000,0.000,0.000,0.000,0.000},
			{0.000,0.000,0.000,0.000,0.000,0.002,0.004,0.002,0.001,0.020,
             0.032,0.052,0.069,0.110,0.175,0.271,0.402,0.556,0.705,0.812,
             0.871,0.896,0.908,0.918,0.926,0.928,0.930,0.926,0.925,0.928,
             0.923,0.916,0.908,0.903,0.909,0.924,0.946,0.954,0.971,0.969,
             0.967,0.965,0.967,0.961,0.949,0.931,0.925,0.929,0.943,0.961,
             0.985,0.992,0.998,0.992,0.994,0.997,0.998,1.000,0.991,0.988,
             0.969,0.926,0.868,0.817,0.819,0.880,0.854,0.572,0.256,0.104,
             0.044,0.022,0.011,0.007,0.000,0.000,0.000,0.000,0.000,0.000,
             0.000},
			{0.000,0.003,0.000,0.001,0.007,0.008,0.008,0.012,0.012,0.028,
             0.041,0.062,0.087,0.114,0.176,0.230,0.306,0.410,0.481,0.543,
             0.598,0.642,0.686,0.719,0.750,0.785,0.817,0.845,0.867,0.881,
             0.902,0.900,0.896,0.892,0.899,0.882,0.872,0.872,0.872,0.878,
             0.868,0.860,0.877,0.884,0.897,0.895,0.898,0.912,0.921,0.927,
             0.937,0.947,0.948,0.954,0.961,0.962,0.962,0.964,0.969,0.956,
             0.952,0.951,0.952,0.953,0.939,0.934,0.928,0.943,0.945,0.935,
             0.944,0.947,0.944,0.949,0.960,0.966,0.971,0.978,0.993,0.998,
             0.996,0.996,0.997,0.986,0.990,0.988,0.992,0.985,0.982,0.978,
             0.970,0.966,0.952,0.927,0.883,0.832,0.751,0.656,0.577,0.483,
             0.393,0.310,0.239,0.184,0.142,0.104,0.080,0.063,0.049,0.041,
             0.036,0.023,0.021,0.019,0.012,0.006,0.008,0.000,0.000,0.000,
             0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
             0.000},
			{0.004,0.001,0.003,0.000,0.002,0.001,0.002,0.002,0.012,0.008,
             0.009,0.018,0.017,0.031,0.037,0.046,0.058,0.076,0.088,0.110,
             0.149,0.196,0.242,0.303,0.367,0.437,0.519,0.610,0.677,0.718,
             0.756,0.774,0.784,0.775,0.789,0.782,0.778,0.766,0.762,0.768,
             0.775,0.769,0.788,0.808,0.794,0.823,0.811,0.819,0.836,0.837,
             0.836,0.851,0.859,0.855,0.871,0.873,0.875,0.859,0.872,0.859,
             0.872,0.863,0.865,0.868,0.877,0.873,0.869,0.876,0.868,0.879,
             0.873,0.876,0.880,0.874,0.870,0.858,0.863,0.859,0.844,0.859,
             0.854,0.863,0.868,0.856,0.847,0.861,0.851,0.852,0.838,0.847,
             0.840,0.831,0.836,0.838,0.822,0.838,0.839,0.842,0.854,0.862,
             0.873,0.868,0.879,0.891,0.898,0.919,0.920,0.926,0.928,0.934,
             0.936,0.953,0.954,0.952,0.960,0.973,0.985,0.972,0.970,0.994,
             0.989,0.975,1.000,0.991,0.968,0.966,0.956,0.929,0.929,0.926,
             0.903,0.924,0.929,0.928,0.920,0.853,0.775,0.659,0.531,0.403,
             0.275,0.218,0.131,0.104,0.075,0.052,0.029,0.028,0.014,0.019,
             0.013,0.007,0.015,0.000,0.004}
		}
	};
	sixs_tables->aot[0]=0.01;
	sixs_tables->aot[1]=0.05;
	sixs_tables->aot[2]=0.10;
	sixs_tables->aot[3]=0.15;
	sixs_tables->aot[4]=0.20;
	sixs_tables->aot[5]=0.30;
	sixs_tables->aot[6]=0.40;
	sixs_tables->aot[7]=0.60;
	sixs_tables->aot[8]=0.80;
	sixs_tables->aot[9]=1.00;
	sixs_tables->aot[10]=1.20;
	sixs_tables->aot[11]=1.40;
	sixs_tables->aot[12]=1.60;
	sixs_tables->aot[13]=1.80;
	sixs_tables->aot[14]=2.00;
	
    /* Determine the 6s command and output filenames */
    if (sprintf(short_name, "L%s%s%s", sat_names[meta->sat], 
        inst_names[meta->inst], "SR") < 0) {
        fprintf(stderr, "ERROR:creating short name\n");
        exit(-1);
    }

    if (!FormatDate(&meta->acq_date, DATE_FORMAT_DATEB, acq_date_string)) {
        fprintf(stderr, "ERROR:formatting acquisition date\n");
        exit(-1);
    }
    acq_date_string[4] = '\0';

    sprintf(local_granule_id, "%s.a%4s%3s.w%1sp%03dr%03d",
        short_name, acq_date_string, &acq_date_string[5],
        wrs_names[meta->wrs_sys], meta->ipath, meta->irow);
	
	/* Run 6s */
#ifdef _OPENMP
        #pragma omp parallel for private (i, j, k, sixs_cmd_filename, \
                                          sixs_out_filename, fd, cmd, line_in, \
                                          tgoz, tgco2, tgo2, tgno2, tgch4, \
                                          tgco) num_threads (get_num_threads())
#endif
	for (i=0;i<SIXS_NB_BANDS;i++) {
        int tg_read = 0;  /* number of "tg" components read for this band */
		for (j=0;j<SIXS_NB_AOT;j++) {
            sprintf (sixs_cmd_filename, "sixs_cmd_%s_%d_%d",
                     local_granule_id, i+1, j+1);
            sprintf (sixs_out_filename, "sixs_output_%s_%d_%d",
                     local_granule_id, i+1, j+1);
			printf("Processing 6S for band %d  AOT %2d\r",i+1,j+1);
            fflush(stdout);
			if ((fd=fopen(sixs_cmd_filename,"w"))==NULL) {
				fprintf(stderr,"ERROR: creating temporary file %s\n",
                        sixs_cmd_filename);
				exit(-1);
			}

			fprintf(fd,"%s <<+ >%s\n",SIXS_APP,sixs_out_filename);
			fprintf(fd,"0 (user defined)\n");
			fprintf(fd,"%.2f %.2f %.2f %.2f %d %d (geometrical conditions "
                    "sza saz vza vaz month day)\n", sixs_tables->sza,
                    sixs_tables->phi, sixs_tables->vza, 0.,
                    sixs_tables->month, sixs_tables->day);
			fprintf(fd,"8 (option for water vapor and ozone)\n");
			fprintf(fd,"%.2f %.2f (water vapor and ozone)\n",
                    sixs_tables->uwv, sixs_tables->uoz);
			fprintf(fd,"1 (continental model)\n"
                       "0 (option for optical thickness at 550 nm)\n");
			fprintf(fd,"%.3f (value of aot550\n",sixs_tables->aot[j]);
			fprintf(fd,"%f (target level)\n",sixs_tables->target_alt);
			fprintf(fd,"-1000 (sensor level : -1000=satellite level)\n");
			switch (sixs_tables->Inst) {
				case SIXS_INST_TM:
					fprintf(fd,"%d (predefined band)\n",tm_band[i]);
				break;
				case SIXS_INST_ETM:
					fprintf(fd,"1 (user defined filter function)\n");
					fprintf(fd,"%05.3f %05.3f (wlinf wlsup)\n",
                            etm_spectral_function.wlinf[i],
                            etm_spectral_function.wlsup[i]);
					for (k=0;k<etm_spectral_function.nbvals[i];k++) {
						fprintf(fd,"%05.3f ",
                                etm_spectral_function.response[i][k]);
						if (!((k+1)%10))
							fprintf(fd,"\n");
					}
					if (k%10)
						fprintf(fd,"\n");
				break;
				default:
					fprintf(stderr,"ERROR: Unknown Instrument in six_run "
                            "parameters\n");
					exit(-1);
			}
			fprintf(fd,"0 (homogeneous surface)\n"
			           "0 (no directional effects)\n"
                       "0 (constant value for rho)\n");
            fprintf(fd,"%.3f (value of rho)\n",sixs_tables->srefl);
			fprintf(fd,"-1 (no atmospheric correction)\n"
			           "0\n"
			           "+\n");
			fclose(fd);

            /* Modified 9/26/2014 to run bash shell vs. sh */
            sprintf(cmd,"bash %s",sixs_cmd_filename);
            if (system(cmd)) {
                fprintf(stderr,"ERROR: Can't run 6S \n");
                exit(-1);
            }

			if ((fd=fopen(sixs_out_filename,"r"))==NULL) {
				fprintf(stderr,"ERROR: reading temporary file %s\n",
                        sixs_out_filename);
				exit(-1);
			}

            while (fgets(line_in,256,fd)) {
                line_in[strlen(line_in)-1]='\0';
                if (j==0) {
                    if (!strncmp(line_in,"*      rayl.  sca. trans. :",
                                 SIXS_STR_LENGTH)) {
                        sscanf(&line_in[SIXS_STR_LENGTH], "%f %f %f",
                               &sixs_tables->T_r_down[i],
                               &sixs_tables->T_r_up[i], &sixs_tables->T_r[i]);
                    }
                    else if (!strncmp(line_in,"*      water   \"     \"    :",
                                      SIXS_STR_LENGTH)) {
                        sscanf(&line_in[SIXS_STR_LENGTH], "%*f %*f %f",
                               &sixs_tables->T_g_wv[i]);
                    }
                    else if (!strncmp(line_in,"*      ozone   \"     \"    :",
                                      SIXS_STR_LENGTH)) {
                        sscanf(&line_in[SIXS_STR_LENGTH], "%*f %*f %f", &tgoz);
                        tg_read++;
                    }
                    else if (!strncmp(line_in,"*      co2     \"     \"    :",
                                      SIXS_STR_LENGTH)) {
                        sscanf(&line_in[SIXS_STR_LENGTH], "%*f %*f %f", &tgco2);
                        tg_read++;
                    }
                    else if (!strncmp(line_in,"*      oxyg    \"     \"    :",
                                      SIXS_STR_LENGTH)) {
                        sscanf(&line_in[SIXS_STR_LENGTH], "%*f %*f %f", &tgo2);
                        tg_read++;
                    }
                    else if (!strncmp(line_in,"*      no2     \"     \"    :",
                                      SIXS_STR_LENGTH)) {
                        sscanf(&line_in[SIXS_STR_LENGTH], "%*f %*f %f", &tgno2);
                        tg_read++;
                    }
                    else if (!strncmp(line_in,"*      ch4     \"     \"    :",
                                      SIXS_STR_LENGTH)) {
                        sscanf(&line_in[SIXS_STR_LENGTH], "%*f %*f %f", &tgch4);
                        tg_read++;
                    }
                    else if (!strncmp(line_in,"*      co      \"     \"    :",
                                      SIXS_STR_LENGTH)) {
                        sscanf(&line_in[SIXS_STR_LENGTH], "%*f %*f %f", &tgco);
                        tg_read++;
                    }
                    if (tg_read == 6)
                    {
                        sixs_tables->T_g_og[i] = tgoz*tgco2*tgo2*tgno2*tgno2
                                               * tgch4*tgco;
                        tg_read++; /* Increment so this computation isn't
                                      repeated. */
                    }
                } /* j = 0 */
                if (!strncmp(line_in,"*      spherical albedo   :",
                             SIXS_STR_LENGTH)) {
                    if (j==0)
                        sscanf(&line_in[SIXS_STR_LENGTH], "%f %*f %f",
                               &sixs_tables->S_r[i], &sixs_tables->S_ra[i][j]);
                    else
                        sscanf(&line_in[SIXS_STR_LENGTH], "%*f %*f %f",
                               &sixs_tables->S_ra[i][j]);
                }
                else if (!strncmp(line_in,"*      optical depth total:",
                                  SIXS_STR_LENGTH)) {
                    sscanf(&line_in[SIXS_STR_LENGTH], "%*f %f",
                           &sixs_tables->aot_wavelength[i][j]);
                }
                else if (!strncmp(line_in,"*      aeros. sca.   \"    :",
                                  SIXS_STR_LENGTH)) {
                    sscanf(&line_in[SIXS_STR_LENGTH], "%f %f %f",
                           &sixs_tables->T_a_down[i][j],
                           &sixs_tables->T_a_up[i][j],
                           &sixs_tables->T_a[i][j]);
                }
                else if (!strncmp(line_in,"*      total  sca.   \"    :",
                                  SIXS_STR_LENGTH)) {
                    sscanf(&line_in[SIXS_STR_LENGTH], "%f %f %f",
                           &sixs_tables->T_ra_down[i][j],
                           &sixs_tables->T_ra_up[i][j],
                           &sixs_tables->T_ra[i][j]);
                }
                else if (!strncmp(line_in,"*      reflectance I      :",
                                  SIXS_STR_LENGTH)) {
                    if (j==0)
                        sscanf(&line_in[SIXS_STR_LENGTH], "%f %f %f",
                               &sixs_tables->rho_r[i],
                               &sixs_tables->rho_a[i][j],
                               &sixs_tables->rho_ra[i][j]);
                    else
                        sscanf(&line_in[SIXS_STR_LENGTH], "%*f %f %f",
                               &sixs_tables->rho_a[i][j],
                               &sixs_tables->rho_ra[i][j]);
                }
            } /* fgets() loop */
			fclose(fd);
/* For OZONE debugging:
	        sprintf (tmp_file, "%s_b%d_aot%02d", sixs_cmd_filename, i, j);
            sprintf (cmd_string, "cp %s %s", sixs_cmd_filename, tmp_file);
            system (cmd_string);
	        sprintf (tmp_file, "%s_b%d_aot%02d", sixs_out_filename, i, j);
            sprintf (cmd_string, "cp %s %s", sixs_out_filename, tmp_file);
            system (cmd_string);
*/
	        unlink(sixs_cmd_filename);
	        unlink(sixs_out_filename);
		}  /* for j */
	}  /* for i */
	printf ("\n");
	return 0;
}


#ifdef SAVE_6S_RESULTS
int read_6S_results_from_file(char *filename,sixs_tables_t *sixs_tables) {
	FILE *fd;
	int run_sixs,i,j;

	if ((fd=fopen(filename,"r"))==NULL)
		return 1;
	run_sixs=0;
	if (fscanf(fd,"%d %d",&sixs_tables->month,&sixs_tables->day) != 2)
		run_sixs=1;
	if (fscanf(fd,"%f",&sixs_tables->srefl)!=1)
		run_sixs=1;
	if (fscanf(fd, "%f %f %f", &sixs_tables->sza, &sixs_tables->vza,
               &sixs_tables->phi) != 3)
		run_sixs=1;
	if (fscanf(fd, "%f %f %f", &sixs_tables->uwv, &sixs_tables->uoz,
               &sixs_tables->target_alt) != 3)
		run_sixs=1;
	for (i=0;i<SIXS_NB_AOT;i++)
		if (fscanf(fd,"%f ",&sixs_tables->aot[i])!=1)
			run_sixs=1;
	for (i=0;i<SIXS_NB_BANDS;i++) {
		if (fscanf(fd, "%f %f %f %f %f %f %f", &sixs_tables->S_r[i],
                   &sixs_tables->T_r_up[i], &sixs_tables->T_r_down[i],
                   &sixs_tables->T_r[i], &sixs_tables->T_g_wv[i],
                   &sixs_tables->T_g_og[i], &sixs_tables->rho_r[i]) != 7)
			run_sixs=1;
		for (j=0;j<SIXS_NB_AOT;j++)
			if (fscanf(fd, "%f %f %f %f %f %f %f %f %f %f %f",
                       &sixs_tables->aot_wavelength[i][j],
                       &sixs_tables->T_a_up[i][j],
                       &sixs_tables->T_a_down[i][j], &sixs_tables->T_a[i][j],
                       &sixs_tables->rho_ra[i][j], &sixs_tables->rho_a[i][j],
                       &sixs_tables->S_ra[i][j], &sixs_tables->T_ra_up[i][j],
                       &sixs_tables->T_ra_down[i][j],
                       &sixs_tables->T_ra[i][j],
                       &sixs_tables->rho_toa[i][j]) != 11)
				run_sixs=1;
	}
	fclose(fd);
	return run_sixs;
	
}

int write_6S_results_to_file(char *filename,sixs_tables_t *sixs_tables) {
	FILE *fd;
	int i,j;

	if ((fd=fopen(filename,"w"))==NULL)
		return -1;
	fprintf(fd,"%02d %03d\n",sixs_tables->month,sixs_tables->day);
	fprintf(fd,"%010.6f\n",sixs_tables->srefl);
	fprintf(fd, "%010.6f %010.6f %010.6f\n", sixs_tables->sza,
            sixs_tables->vza, sixs_tables->phi);
	fprintf(fd, "%010.6f %010.6f %010.2f\n", sixs_tables->uwv,
            sixs_tables->uoz, sixs_tables->target_alt);
	for (i=0;i<SIXS_NB_AOT;i++)
		fprintf(fd,"%07.4f ",sixs_tables->aot[i]);
	fprintf(fd,"\n");
	for (i=0;i<SIXS_NB_BANDS;i++) {
		fprintf(fd, "%010.6f %010.6f %010.6f %010.6f %010.6f %010.6f "
                "%010.6f\n", sixs_tables->S_r[i], sixs_tables->T_r_up[i],
                sixs_tables->T_r_down[i], sixs_tables->T_r[i],
                sixs_tables->T_g_wv[i], sixs_tables->T_g_og[i],
                sixs_tables->rho_r[i]);
			for (j=0;j<SIXS_NB_AOT;j++)
				fprintf(fd, "%07.4f %010.6f %010.6f %010.6f %010.6f %010.6f "
                        "%010.6f %010.6f %010.6f %010.6f %010.6f\n",
                        sixs_tables->aot_wavelength[i][j],
                        sixs_tables->T_a_up[i][j],
                        sixs_tables->T_a_down[i][j], sixs_tables->T_a[i][j],
                        sixs_tables->rho_ra[i][j], sixs_tables->rho_a[i][j],
                        sixs_tables->S_ra[i][j], sixs_tables->T_ra_up[i][j],
                        sixs_tables->T_ra_down[i][j], sixs_tables->T_ra[i][j],
                        sixs_tables->rho_toa[i][j]);
	}
	fclose(fd);
	return 0;
}
#endif
