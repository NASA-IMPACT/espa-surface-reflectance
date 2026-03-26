#ifndef _READ_LEVEL1_QA_H_
#define _READ_LEVEL1_QA_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

/* Defines for the Level-1 pixel QA band */
/* Define the constants used for shifting bits and ANDing with the bits to
   get to the desire quality bits */
#define ESPA_L1_SINGLE_BIT 0x01             /* 00000001 */
#define ESPA_L1_DOUBLE_BIT 0x03             /* 00000011 */
#define ESPA_L1_DESIGNATED_FILL_BIT 0       /* one bit */

/******************************************************************************
MODULE:  level1_qa_is_fill

PURPOSE: Determines if the current Level-1 QA pixel is fill

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
static inline bool level1_qa_is_fill
(
    uint16_t l1_qa_pix      /* I: Level-1 QA value for current pixel */
)
{
    if (((l1_qa_pix >> ESPA_L1_DESIGNATED_FILL_BIT) & ESPA_L1_SINGLE_BIT) == 1)
        return true;
    else
        return false;
}


/******************************************************************************
MODULE:  level1_qa_is_saturated

PURPOSE: Determines if the current Level-1 QA pixel is saturated

RETURN VALUE:
Type = boolean
Value           Description
-----           -----------
true            Pixel is saturated
false           Pixel is not saturated

NOTES:
1. This is an inline function so it should be fast as the function call overhead
   is eliminated by dropping the code inline with the original application.
******************************************************************************/
static inline bool level1_qa_is_saturated
(
    uint16_t l1_qa_pix,     /* I: Level-1 QA value for current pixel */
    uint8_t band            /* I: Band to check for saturation */
)
{
    /* The saturation bit is 0-based for each band */
    if (((l1_qa_pix >> (band-1)) & ESPA_L1_SINGLE_BIT) == 1)
        return true;
    else
        return false;
}

#endif
