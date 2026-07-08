/* cref-owned C shims: adapt production C functions whose signatures take heavy
 * ESPA structs into flat, FFI-friendly entry points. No production C is edited.
 */
#include <string.h>
#include "aero_interp.h"   /* pulls lasrc.h + espa_metadata.h (Espa_internal_meta_t) */

/* aerosol_interp_landsat takes an Espa_internal_meta_t only to compute
 * refl_indx, which it then never uses. Pass a zeroed metadata (nbands = 0) so
 * the band-search loop is skipped, and forward the real array arguments. */
void cref_aerosol_interp_landsat(
    int aero_window,
    int half_aero_window,
    uint16 *qaband,
    uint8 *ipflag,
    float *taero,
    int nlines,
    int nsamps)
{
    Espa_internal_meta_t meta;
    memset(&meta, 0, sizeof(meta));  /* nbands = 0 -> refl_indx defaults, unused */
    aerosol_interp_landsat(&meta, aero_window, half_aero_window,
                           qaband, ipflag, taero, nlines, nsamps);
}
