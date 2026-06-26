#ifndef _LUT_SUBR_H_
#define _LUT_SUBR_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "common.h"
#include "espa_metadata.h"
#include "espa_geoloc.h"
#include "error_handler.h"
#include "date.h"

/* Prototypes */
void atmcorlamb2_new
(
    Sat_t sat,                /* I: satellite */
    float tgo,                /* I: other gaseous transmittance  */
    float roatm_upper,        /* I: roatm upper bound poly_fit, given band */
    float roatm_coef[NCOEF],  /* I: poly_fit coefficients for roatm  */
    float ttatmg_coef[NCOEF], /* I: poly_fit coefficients for ttatmg */
    float satm_coef[NCOEF],   /* I: poly_fit coefficients for satm */
    float raot550nm,          /* I: nearest value of AOT */
    int iband,                /* I: band index (0-based) */
    float normext_ib_0_3,     /* I: normext[iband][0][3] */
    float rotoa,              /* I: top of atmosphere reflectance */
    float *roslamb,           /* O: lambertian surface reflectance */
    float eps                 /* I: angstroem coefficient; spectral dependency
                                    of the AOT */
);

int atmcorlamb2
(
    Sat_t sat,                   /* I: satellite */
    float xts,                   /* I: solar zenith angle (deg) */
    float xtv,                   /* I: observation zenith angle (deg) */
    float xmus,                  /* I: cosine of solar zenith angle */
    float xmuv,                  /* I: cosine of observation zenith angle */
    float xfi,                   /* I: azimuthal difference between sun and
                                       observation (deg) */
    float cosxfi,                /* I: cosine of azimuthal difference */
    float raot550nm,             /* I: nearest value of AOT */
    int iband,                   /* I: band index (0-based) */
    float pres,                  /* I: surface pressure */
    float tpres[NPRES_VALS],     /* I: surface pressure table */
    float aot550nm[NAOT_VALS],   /* I: AOT look-up table */
    float *rolutt,               /* I: intrinsic reflectance table
                                       [NSR_BANDS x NPRES_VALS x NAOT_VALS x
                                        NSOLAR_VALS] */
    float *transt,               /* I: transmission table
                                       [NSR_BANDS x NPRES_VALS x NAOT_VALS x
                                        NSUNANGLE_VALS] */
    float xtsstep,               /* I: solar zenith step value */
    float xtsmin,                /* I: minimum solar zenith value */
    float xtvstep,               /* I: observation step value */
    float xtvmin,                /* I: minimum observation value */
    float *sphalbt,              /* I: spherical albedo table
                                       [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *normext,              /* I: aerosol extinction coefficient at
                                       the current wavelength (normalized
                                       at 550nm)
                                       [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *tsmax,                /* I: maximum scattering angle table
                                       [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *tsmin,                /* I: minimum scattering angle table
                                       [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *nbfic,                /* I: communitive number of azimuth angles
                                       [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *nbfi,                 /* I: number of azimuth angles
                                       [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float tts[NSOLAR_ZEN_VALS],  /* I: sun angle table */
    int32 indts[NSUNANGLE_VALS], /* I: index for the sun angle table */
    float *ttv,                  /* I: view angle table
                                       [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float uoz,                   /* I: total column ozone */
    float uwv,                   /* I: total column water vapor (precipital
                                       water vapor) */
    float tauray[NSR_BANDS],     /* I: molecular optical thickness coeff */
    double ogtransa1[NSR_BANDS], /* I: other gases transmission coeff */
    double ogtransb0[NSR_BANDS], /* I: other gases transmission coeff */
    double ogtransb1[NSR_BANDS], /* I: other gases transmission coeff */
    double wvtransa[NSR_BANDS],  /* I: water vapor transmission coeff */
    double wvtransb[NSR_BANDS],  /* I: water vapor transmission coeff */
    double oztransa[NSR_BANDS],  /* I: ozone transmission coeff */
    float rotoa,                 /* I: top of atmosphere reflectance */
    float *roslamb,              /* O: lambertian surface reflectance */
    float *tgo,                  /* O: other gaseous transmittance */
    float *roatm,                /* O: intrinsic atmospheric reflectance */
    float *ttatmg,               /* O: total atmospheric transmission */
    float *satm,                 /* O: spherical albedo */
    float *xrorayp,              /* O: reflectance of the atmosphere due to
                                       molecular (Rayleigh) scattering */
    float eps                    /* I: angstroem coefficient; spectral
                                       dependency of the AOT */
);

void subaeroret_new
(
    Sat_t sat,                             /* I: satellite */
    bool water,                            /* I: water pixel flag */
    int iband1,                            /* I: band 1 index (0-based) */
    float erelc[NSR_BANDS],                /* I: band ratio variable */
    float troatm[NSR_BANDS],               /* I: toa reflectance */
    float tgo_arr[NREFL_BANDS],            /* I: per-band other gaseous
                                                 transmittance */
    int roatm_iaMax[NREFL_BANDS],          /* I: roatm_iaMax */
    float roatm_coef[NREFL_BANDS][NCOEF],  /* I: per band polynomial
                                                 coefficients for roatm */
    float ttatmg_coef[NREFL_BANDS][NCOEF], /* I: per band polynomial
                                                 coefficients for ttatmg */
    float satm_coef[NREFL_BANDS][NCOEF],   /* I: per band polynomial
                                                 coefficients for satm */
    float normext_p0a3_arr[NREFL_BANDS],   /* I: normext[iband][0][3] */
    float *raot,     /* O: AOT reflectance */
    float *residual, /* O: model residual */
    int *iaots,      /* I/O: AOT index that is passed in and out for multiple
                             calls (0-based) */
    float eps        /* I: angstroem coefficient; spectral dependency of AOT */
);

int subaeroret
(
    Sat_t sat,                   /* I: satellite */
    bool water,                  /* I: water pixel flag */
    int iband1,                  /* I: band 1 index (0-based) */
    float xts,                   /* I: solar zenith angle (deg) */
    float xtv,                   /* I: observation zenith angle (deg) */
    float xmus,                  /* I: cosine of solar zenith angle */
    float xmuv,                  /* I: cosine of observation zenith angle */
    float xfi,                   /* I: azimuthal difference between sun and
                                       observation (deg) */
    float cosxfi,                /* I: cosine of azimuthal difference */
    float pres,                  /* I: surface pressure */
    float uoz,                   /* I: total column ozone */
    float uwv,                   /* I: total column water vapor (precipital
                                       water vapor) */
    float erelc[NSR_BANDS],      /* I: band ratio variable */
    float troatm[NSR_BANDS],     /* I: atmospheric reflectance table */
    float tpres[NPRES_VALS],     /* I: surface pressure table */
    float *rolutt,               /* I: intrinsic reflectance table
                          [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSOLAR_VALS] */
    float *transt,               /* I: transmission table
                       [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSUNANGLE_VALS] */
    float xtsstep,               /* I: solar zenith step value */
    float xtsmin,                /* I: minimum solar zenith value */
    float xtvstep,               /* I: observation step value */
    float xtvmin,                /* I: minimum observation value */
    float *sphalbt,              /* I: spherical albedo table
                                       [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *normext,              /* I: aerosol extinction coefficient at the
                                       current wavelength (normalized at 550nm)
                                       [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *tsmax,                /* I: maximum scattering angle table
                                       [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *tsmin,                /* I: minimum scattering angle table
                                       [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *nbfic,                /* I: communitive number of azimuth angles
                                       [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *nbfi,                 /* I: number of azimuth angles
                                       [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float tts[NSOLAR_ZEN_VALS],  /* I: sun angle table */
    int32 indts[NSUNANGLE_VALS], /* I: index for the sun angle table */
    float *ttv,                  /* I: view angle table
                                       [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float tauray[NSR_BANDS],     /* I: molecular optical thickness coeff */
    double ogtransa1[NSR_BANDS], /* I: other gases transmission coeff */
    double ogtransb0[NSR_BANDS], /* I: other gases transmission coeff */
    double ogtransb1[NSR_BANDS], /* I: other gases transmission coeff */
    double wvtransa[NSR_BANDS],  /* I: water vapor transmission coeff */
    double wvtransb[NSR_BANDS],  /* I: water vapor transmission coeff */
    double oztransa[NSR_BANDS],  /* I: ozone transmission coeff */
    float *raot,                 /* O: AOT reflectance */
    float *residual,             /* O: model residual */
    int *iaots,                  /* I/O: AOT index that is passed in and out
                                         for multiple calls (0-based) */
    float eps                    /* I: angstroem coefficient; spectral
                                       dependency of the AOT */
);

void local_chand
(
    float xphi,    /* I: azimuthal difference between sun and observation
                         (deg) */
    float xmuv,    /* I: cosine of observation zenith angle */
    float xmus,    /* I: cosine of solar zenith angle */
    float xtau,    /* I: molecular optical depth */
    float *xrray   /* O: molecular reflectance, 0.0 to 1.0 */
);

void comptg
(
    int iband,                   /* I: band index (0-based) */
    float xmus,                  /* I: cosine of solar zenith angle */
    float xmuv,                  /* I: cosine of observation zenith angle */
    float uoz,                   /* I: total column ozone */
    float uwv,                   /* I: total column water vapor (precipital
                                       water vapor) */
    float atm_pres,              /* I: pressure at sea level */
    double ogtransa1[NSR_BANDS], /* I: other gases transmission coeff */
    double ogtransb0[NSR_BANDS], /* I: other gases transmission coeff */
    double ogtransb1[NSR_BANDS], /* I: other gases transmission coeff */
    double wvtransa[NSR_BANDS],  /* I: water vapor transmission coeff */
    double wvtransb[NSR_BANDS],  /* I: water vapor transmission coeff */
    double oztransa[NSR_BANDS],  /* I: ozone transmission coeff */
    float *tgoz,                 /* O: ozone transmission */
    float *tgwv,                 /* O: water vapor transmission */
    float *tgwvhalf,             /* O: water vapor transmission, half content */
    float *tgog                  /* O: other gases transmission */
);

void compsalb
(
    int ip1,            /* I: index variable for surface pressure */
    int ip2,            /* I: index variable for surface pressure */
    int iaot1,          /* I: index variable for AOT */
    int iaot2,          /* I: index variable for AOT */
    float raot550nm,    /* I: nearest value of AOT */
    int iband,          /* I: band index (0-based) */
    float pres,         /* I: surface pressure */
    float tpres[NPRES_VALS],   /* I: surface pressure table */
    float aot550nm[NAOT_VALS], /* I: AOT look-up table */
    float *sphalbt,     /* I: spherical albedo table
                              [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *normext,     /* I: aerosol extinction coefficient at the current
                              wavelength (normalized at 550nm)
                              [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *satm         /* O: spherical albedo */
);

void comptrans
(
    int ip1,            /* I: index variable for surface pressure */
    int ip2,            /* I: index variable for surface pressure */
    int iaot1,          /* I: index variable for AOT */
    int iaot2,          /* I: index variable for AOT */
    float xts,          /* I: zenith angle */
    float raot550nm,    /* I: nearest value of AOT */
    int iband,          /* I: band index (0-based) */
    float pres,         /* I: surface pressure */
    float tpres[NPRES_VALS],   /* I: surface pressure table */
    float aot550nm[NAOT_VALS], /* I: AOT look-up table */
    float *transt,      /* I: transmission table
                       [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSUNANGLE_VALS] */
    float xtsstep,      /* I: zenith angle step value */
    float xtsmin,       /* I: minimum zenith angle value */
    float tts[NSOLAR_ZEN_VALS], /* I: sun angle table */
    float *xtts         /* O: downward transmittance */
);

void comproatm
(
    int ip1,            /* I: index variable for surface pressure */
    int ip2,            /* I: index variable for surface pressure */
    int iaot1,          /* I: index variable for AOT */
    int iaot2,          /* I: index variable for AOT */
    float xts,          /* I: solar zenith angle (deg) */
    float xtv,          /* I: observation zenith angle (deg) */
    float xmus,         /* I: cosine of solar zenith angle */
    float xmuv,         /* I: cosine of observation zenith angle */
    float cosxfi,       /* I: cosine of azimuthal difference */
    float raot550nm,    /* I: nearest value of AOT */
    int iband,          /* I: band index (0-based) */
    float pres,         /* I: surface pressure */
    float tpres[NPRES_VALS],   /* I: surface pressure table */
    float aot550nm[NAOT_VALS], /* I: AOT look-up table */
    float *rolutt,      /* I: intrinsic reflectance table
                           [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSOLAR_VALS] */
    float *tsmax,       /* I: maximum scattering angle table
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *tsmin,       /* I: minimum scattering angle table
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *nbfic,       /* I: communitive number of azimuth angles
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *nbfi,        /* I: number of azimuth angles
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float tts[NSOLAR_ZEN_VALS],  /* I: sun angle table */
    int32 indts[NSUNANGLE_VALS], /* I: index for the sun angle table */
    float *ttv,         /* I: view angle table
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float xtsstep,      /* I: solar zenith step value */
    float xtvstep,      /* I: observation step value */
    float xtvmin,       /* I: minimum observation value */
    int its,            /* I: index for the sun angle table */
    int itv,            /* I: index for the view angle table */
    float *roatm        /* O: atmospheric intrinsic reflectance */
);

int readluts
(
    Sat_t sat,                  /* I: satellite */
    float *tsmax,               /* O: maximum scattering angle table
                                      [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *tsmin,               /* O: minimum scattering angle table
                                      [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *ttv,                 /* O: view angle table
                                      [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float tts[NSOLAR_ZEN_VALS], /* O: sun angle table */
    float *nbfic,               /* O: communitive number of azimuth angles
                                      [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *nbfi,                /* O: number of azimuth angles
                                      [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    int32 indts[NSUNANGLE_VALS],/* O: index for the sun angle table */
    float *rolutt,              /* O: intrinsic reflectance table
                           [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSOLAR_VALS] */
    float *transt,              /* O: transmission table
                        [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSUNANGLE_VALS] */
    float *sphalbt,             /* O: spherical albedo table
                                      [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *normext,             /* O: aerosol extinction coefficient at the
                                      current wavelength (normalized at 550nm)
                                      [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float xtsstep,              /* I: solar zenith step value */
    float xtsmin,               /* I: minimum solar zenith value */
    char anglehdf[STR_SIZE],    /* I: angle HDF filename */
    char intrefnm[STR_SIZE],    /* I: intrinsic reflectance filename */
    char transmnm[STR_SIZE],    /* I: transmission filename */
    char spheranm[STR_SIZE]     /* I: spherical albedo filename */
);

int memory_allocation_main
(
    Sat_t sat,           /* I: satellite */
    int nlines,          /* I: number of lines in the scene */
    int nsamps,          /* I: number of samples in the scene */
    bool use_orig_aero,  /* I: use the original aerosol handling if specified,
                               o/w use the semi-empirical approach */
    int16 **sza,         /* O: solar zenith angle, nlines x nsamps  */
    int16 **saa,         /* O: solar azimuth angle table, nlines x nsamps */
    int16 **vza,         /* O: view zenith angle, nlines x nsamps  */
    int16 **vaa,         /* O: view azimuth angle table, nlines x nsamps */
    uint16 **qaband,     /* O: QA band for the input image, nlines x nsamps */
    uint16 **out_band,   /* O: scaled output, nlines x nsamps */
    float ***sband,      /* O: unscaled surface reflectance and brightness temp
                               bands, nlines x nsamps */
    float ***toaband     /* O: Sentinel unscaled TOA reflectance bands,
                               nlines x nsamps */
);

int landsat_memory_allocation_sr
(
    int nlines,          /* I: number of lines in the scene */
    int nsamps,          /* I: number of samples in the scene */
    float **aerob1,      /* O: atmospherically corrected band 1 data
                              (unscaled TOA refl), nlines x nsamps */
    float **aerob2,      /* O: atmospherically corrected band 2 data
                              (unscaled TOA refl), nlines x nsamps */
    float **aerob4,      /* O: atmospherically corrected band 4 data
                              (unscaled TOA refl), nlines x nsamps */
    float **aerob5,      /* O: atmospherically corrected band 5 data
                              (unscaled TOA refl), nlines x nsamps */
    float **aerob7,      /* O: atmospherically corrected band 7 data
                              (unscaled TOA refl), nlines x nsamps */
    uint8 **ipflag,      /* O: QA flag to assist with aerosol interpolation,
                               nlines x nsamps */
    float **twvi,        /* O: interpolated water vapor value,
                               nlines x nsamps */
    float **tozi,        /* O: interpolated ozone value, nlines x nsamps */
    float **tp,          /* O: interpolated pressure value, nlines x nsamps */
    float **taero,       /* O: aerosol values for each pixel, nlines x nsamps */
    float **teps,        /* O: eps (angstrom coefficient) for each pixel,
                               nlines x nsamps*/
    int16 **dem,         /* O: CMG DEM data array [DEM_NBLAT x DEM_NBLON] */
    int16 **andwi,       /* O: avg NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 **sndwi,       /* O: standard NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 **ratiob1,     /* O: mean band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **ratiob2,     /* O: mean band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **ratiob7,     /* O: mean band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **intratiob1,  /* O: band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **intratiob2,  /* O: band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **intratiob7,  /* O: band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **slpratiob1,  /* O: slope band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **slpratiob2,  /* O: slope band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **slpratiob7,  /* O: slope band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    uint16 **wv,         /* O: water vapor values [CMG_NBLAT x CMG_NBLON] */
    uint8 **oz,          /* O: ozone values [CMG_NBLAT x CMG_NBLON] */
    float **rolutt,      /* O: intrinsic reflectance table
                         [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSOLAR_VALS] */
    float **transt,      /* O: transmission table
                        [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSUNANGLE_VALS] */
    float **sphalbt,     /* O: spherical albedo table
                               [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float **normext,     /* O: aerosol extinction coefficient at the current
                               wavelength (normalized at 550nm)
                               [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float **tsmax,       /* O: maximum scattering angle table
                               [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float **tsmin,       /* O: minimum scattering angle table
                               [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float **nbfic,       /* O: communitive number of azimuth angles
                               [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float **nbfi,        /* O: number of azimuth angles
                               [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float **ttv          /* O: view angle table
                               [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
);

int sentinel_memory_allocation_sr
(
    int nlines,          /* I: number of lines in the scene */
    int nsamps,          /* I: number of samples in the scene */
    uint8 **ipflag,      /* O: QA flag to assist with aerosol interpolation,
                               nlines x nsamps */
    float **twvi,        /* O: interpolated water vapor value,
                               nlines x nsamps */
    float **tozi,        /* O: interpolated ozone value, nlines x nsamps */
    float **tp,          /* O: interpolated pressure value, nlines x nsamps */
    float **taero,       /* O: aerosol values for each pixel, nlines x nsamps */
    float **teps,        /* O: eps (angstrom coefficient) for each pixel,
                               nlines x nsamps*/
    int16 **dem,         /* O: CMG DEM data array [DEM_NBLAT x DEM_NBLON] */
    int16 **andwi,       /* O: avg NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 **sndwi,       /* O: standard NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 **ratiob1,     /* O: mean band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **ratiob2,     /* O: mean band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **ratiob7,     /* O: mean band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **intratiob1,  /* O: band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **intratiob2,  /* O: band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **intratiob7,  /* O: band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **slpratiob1,  /* O: slope band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **slpratiob2,  /* O: slope band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **slpratiob7,  /* O: slope band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    uint16 **wv,         /* O: water vapor values [CMG_NBLAT x CMG_NBLON] */
    uint8 **oz,          /* O: ozone values [CMG_NBLAT x CMG_NBLON] */
    float **rolutt,      /* O: intrinsic reflectance table
                         [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSOLAR_VALS] */
    float **transt,      /* O: transmission table
                        [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSUNANGLE_VALS] */
    float **sphalbt,     /* O: spherical albedo table
                               [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float **normext,     /* O: aerosol extinction coefficient at the current
                               wavelength (normalized at 550nm)
                               [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float **tsmax,       /* O: maximum scattering angle table
                               [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float **tsmin,       /* O: minimum scattering angle table
                               [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float **nbfic,       /* O: communitive number of azimuth angles
                               [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float **nbfi,        /* O: number of azimuth angles
                               [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float **ttv          /* O: view angle table
                               [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
);

int read_auxiliary_files
(
    char *cmgdemnm,     /* I: climate modeling grid DEM filename */
    char *rationm,      /* I: ratio averages filename */
    char *auxnm,        /* I: auxiliary filename for ozone and water vapor */
    char *aux2nm,       /* I: 2nd auxiliary filename for ozone and water vapor (23-FEB-26, JPR) */
    aux_src_t aux_src,  /* I: indentifies the source of atmospheric aux data */
    int16 *dem,         /* O: CMG DEM data array [DEM_NBLAT x DEM_NBLON] */
    int16 *andwi,       /* O: avg NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 *sndwi,       /* O: standard NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob1,     /* O: mean band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob2,     /* O: mean band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob7,     /* O: mean band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *intratiob1,  /* O: band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *intratiob2,  /* O: band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *intratiob7,  /* O: band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *slpratiob1,  /* O: slope band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *slpratiob2,  /* O: slope band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *slpratiob7,  /* O: slope band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    uint16 *wv,         /* O: water vapor values [CMG_NBLAT x CMG_NBLON] */
    uint8 *oz,           /* O: ozone values [CMG_NBLAT x CMG_NBLON] */
    Date_t date
);

#endif
