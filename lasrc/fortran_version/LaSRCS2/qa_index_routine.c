#include <stdio.h>
#include <string.h>
#include "mfhdf.h"

/* Compile with:

   cc -O0 -c qa_index_routine.c

 */

static char QAMAP_1KM[2600] =
                 {"\n\tBits are listed from the MSB (bit 7) to the LSB (bit 0):\n"
                  "\tBit    Description\n"
                  "\t7      internal test; \n"
                  "\t       1 -- yes\n"
                  "\t       0 -- no\n"
                  "\t6      unused; \n"
                  "\t       1 -- yes\n"
                  "\t       0 -- no\n"
                 "\t4-5     aerosol;\n"
                  "\t       00 -- climatology\n"
                  "\t       01 -- low\n"
		  "\t       10 -- average\n"
		  "\t       11 -- high\n"
                  "\t3      cloud shadow; \n"
                  "\t       1 -- yes\n"
                  "\t       0 -- no\n"
                  "\t2      adjacent to cloud; \n"
                  "\t       1 -- yes\n"
                  "\t       0 -- no\n"
                  "\t1      cloud; \n"
                  "\t       1 -- yes\n"
                  "\t       0 -- no\n"
                  "\t0      cirrus cloud; \n"
                  "\t       1 -- yes\n"
                  "\t       0 -- no\n\0"
		  
};

int set_qamap_(int32 *sds) {
char message[2600] = {"\0"};
strcpy(message, QAMAP_1KM);  
if(SDsetattr(*sds, "QA index\0", 4, strlen(message), (void *)message)){
       printf("SDS %d: Error: Unable to write QA map indices for QA\n", *sds);
       return(1);
     }  
return(0);
}

