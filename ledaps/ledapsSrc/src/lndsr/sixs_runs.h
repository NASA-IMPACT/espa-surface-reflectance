#ifndef SIXS_H
#define SIXS_H
#include "input.h"

#define SIXS_NB_AOT 15
#define SIXS_NB_BANDS 6
#define SIXS_APP "sixsV1.0B"

typedef enum {
  SIXS_INST_NULL = -1,
  SIXS_INST_MSS = 0, 
  SIXS_INST_TM,
  SIXS_INST_ETM, 
  SIXS_INST_MAX
} Sixs_Inst_t;

typedef struct {
	Sixs_Inst_t Inst;
	int month,day;
	float sza,vza,phi,uwv,uoz,srefl;
	float aot_wavelength[SIXS_NB_BANDS][SIXS_NB_AOT];
	float aot[SIXS_NB_AOT];
	float S_r[SIXS_NB_BANDS];	/* Rayleigh spherical albedo */
	float T_a_up[SIXS_NB_BANDS][SIXS_NB_AOT];  /* aerosol transmittance Up */
	float T_a_down[SIXS_NB_BANDS][SIXS_NB_AOT];  /* aerosol transmittance Down */
	float T_a[SIXS_NB_BANDS][SIXS_NB_AOT];  /* aerosol transmittance */
	float rho_ra[SIXS_NB_BANDS][SIXS_NB_AOT];  /* rayleigh + aerosol reflectance */
	float rho_a[SIXS_NB_BANDS][SIXS_NB_AOT];  /* aerosol reflectance */

	float S_ra[SIXS_NB_BANDS][SIXS_NB_AOT];	/* Aerosol+Rayleigh spherical albedo */
	float T_ra_up[SIXS_NB_BANDS][SIXS_NB_AOT];  /* Rayleigh+aerosol transmittance Up */
	float T_ra_down[SIXS_NB_BANDS][SIXS_NB_AOT];  /* Rayleigh+aerosol transmittance Down */
	float T_ra[SIXS_NB_BANDS][SIXS_NB_AOT];  /* Rayleigh+aerosol transmittance */

	float T_r_up[SIXS_NB_BANDS];  /* Rayleigh transmittance Up */
	float T_r_down[SIXS_NB_BANDS];  /* Rayleigh transmittance Down */
	float T_r[SIXS_NB_BANDS];  /* Rayleigh transmittance */
	float T_g_wv[SIXS_NB_BANDS];  /* water vapor transmittance */
	float T_g_og[SIXS_NB_BANDS];  /* Other gases transmittance */
	float rho_r[SIXS_NB_BANDS];  /* rayleigh reflectance */
	float rho_toa[SIXS_NB_BANDS][SIXS_NB_AOT];  /* Integrated apparent reflectance */
	float target_alt;					/* target altitude in km */
} sixs_tables_t;

typedef struct {
	int band,month,day;
	float sza,vza,phi,uwv,uoz,aot,srefl;
	float S_r;	/* Rayleigh spherical albedo */
	float T_r_up,T_r_down;  /* Rayleigh transmittance */
	float T_a_up,T_a_down;  /* Aerosol transmittance */
	float T_g_wv;  /* water vapor transmittance */
	float T_g_og;  /* Other gases transmittance */
	float rho_r;  /* rayleigh reflectance */
	float rho_a;  /* aerosol reflectance */
} sixs_atmos_params_t;

int create_6S_tables(sixs_tables_t *sixs_tables, Input_meta_t *meta);

#endif
