/*****************************************************************************
FILE: read_level2_qa.h
  
PURPOSE: Contains function prototypes for the LaSRC Level-2 QA band
manipulation for collection products.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
*****************************************************************************/

#ifndef READ_LEVEL2_QA_H
#define READ_LEVEL2_QA_H

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "error_handler.h"
#include "espa_metadata.h"
#include "parse_metadata.h"
#include "raw_binary_io.h"

/* Defines */
/* Define the constants used for shifting bits and ANDing with the bits to
   get to the desire quality bits */
#define ESPA_L2_SINGLE_BIT 0x01            /* 00000001 */
#define ESPA_L2_DOUBLE_BIT 0x03            /* 00000011 */

/* LaSRC QA bits - aerosol */
#define LASRC_FILL_VAL 0x01                /* bit 0 */
#define LASRC_VALID_AEROSOL_RET_VAL 0x02   /* bit 1 */
#define LASRC_WATER_VAL 0x04               /* bit 2 */
#define LASRC_AEROSOL_INTERP_VAL 0x20      /* bit 5 -- 32 decimal = 02 hex */
#define LASRC_AEROSOL_LEVEL_BIT 6          /* bits 6&7 */


/* Inline Function Prototypes */

/*** LaSRC AEROSOL ***/
/******************************************************************************
MODULE:  lasrc_qa_is_fill

PURPOSE: Determines if the LaSRC aerosol QA pixel is fill

RETURN VALUE:
Type = boolean
Value           Description
-----           -----------
true            Pixel is fill
false           Pixel is not fill

NOTES:
1. This is an inline function so it should be fast as the function call overhead
   is eliminated by dropping the code inline with the original application.
******************************************************************************/
static inline bool lasrc_qa_is_fill
(
    uint8_t l2_qa_pix      /* I: Level-2 QA value for current pixel */
)
{
    if (l2_qa_pix & LASRC_FILL_VAL)
        return true;
    else
        return false;
}


/******************************************************************************
MODULE:  lasrc_qa_is_valid_aerosol_retrieval

PURPOSE: Determines if the aerosol retrievel for the LaSRC QA pixel value is
valid

RETURN VALUE:
Type = boolean
Value           Description
-----           -----------
true            Pixel is valid aerosol retrieval
false           Pixel is not valid aerosol retrieval

NOTES:
1. This is an inline function so it should be fast as the function call overhead
   is eliminated by dropping the code inline with the original application.
******************************************************************************/
static inline bool lasrc_qa_is_valid_aerosol_retrieval
(
    uint8_t l2_qa_pix      /* I: Level-2 QA value for current pixel */
)
{
    if (l2_qa_pix & LASRC_VALID_AEROSOL_RET_VAL)
        return true;
    else
        return false;
}


/******************************************************************************
MODULE:  lasrc_qa_is_water

PURPOSE: Determines if the LaSRC aerosol QA pixel is flagged as water, which
changes the way the aerosols are retrieved

RETURN VALUE:
Type = boolean
Value           Description
-----           -----------
true            Pixel is water
false           Pixel is not water

NOTES:
1. This is an inline function so it should be fast as the function call overhead
   is eliminated by dropping the code inline with the original application.
******************************************************************************/
static inline bool lasrc_qa_is_water
(
    uint8_t l2_qa_pix      /* I: Level-2 QA value for current pixel */
)
{
    if (l2_qa_pix & LASRC_WATER_VAL)
        return true;
    else
        return false;
}


/******************************************************************************
MODULE:  lasrc_qa_is_aerosol_interp

PURPOSE: Determines if aerosol value for the LaSRC aerosol QA pixel is
interpolated

RETURN VALUE:
Type = boolean
Value           Description
-----           -----------
true            Pixel is aerosol interpolation
false           Pixel is not aerosol interpolation

NOTES:
1. This is an inline function so it should be fast as the function call overhead
   is eliminated by dropping the code inline with the original application.
******************************************************************************/
static inline bool lasrc_qa_is_aerosol_interp
(
    uint8_t l2_qa_pix      /* I: Level-2 QA value for current pixel */
)
{
    if (l2_qa_pix & LASRC_AEROSOL_INTERP_VAL)
        return true;
    else
        return false;
}


/******************************************************************************
MODULE:  lasrc_qa_aerosol_level

PURPOSE: Returns the aerosol level (0-3) for the current LaSRC aerosol QA
pixel.

RETURN VALUE:
Type = uint8_t
Value           Description
-----           -----------
0               Aerosol level bits are 00 (none)
1               Aerosol level bits are 01 (low)
2               Aerosol level bits are 10 (moderate)
3               Aerosol level bits are 11 (high)

NOTES:
1. This is an inline function so it should be fast as the function call overhead
   is eliminated by dropping the code inline with the original application.
******************************************************************************/
static inline uint8_t lasrc_qa_aerosol_level
(
    uint8_t l2_qa_pix      /* I: Level-2 QA value for current pixel */
)
{
    return ((l2_qa_pix >> LASRC_AEROSOL_LEVEL_BIT) & ESPA_L2_DOUBLE_BIT);
}

#endif
