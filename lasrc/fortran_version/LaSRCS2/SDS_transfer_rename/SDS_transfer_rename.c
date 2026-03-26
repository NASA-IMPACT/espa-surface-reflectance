/* compile with:
   gcc -O0 -o SDS_transfer_rename SDS_transfer_rename.c -I/usr/local/TOOLKIT/SDPToolkit_gfortran44-5.2.17/hdf/linux64/hdf-4.2.5/include/ \
                           -L/usr/local/TOOLKIT/SDPToolkit_gfortran44-5.2.17/hdf/linux64/hdf-4.2.5/lib/ -lmfhdf -ldf -ljpeg -lz -lsz -lm 
 */
#include <ctype.h>
#include <stdio.h>
#include <strings.h>
#include <math.h>
#include "mfhdf.h"
int main(int argc, char **argv)
{
FILE *fd;
int32 i,j,k;
int32 start[5];
int32 sd_id3, sds_out;
int32 sd_id1, n_sets1, n_gattr1, sds_id1, rank1, dims1[5], number_type1, nattr1, dcount, dnt, dattr;
int32 sd_id2, sds_id2;
char name1[128];
char dimnames[5][128];
int8 transferred[512];
void *image;
int size;
void *get_image(int32 id, int32 *dims, int32 type, int *size);
void transfer_attributes(int32 id1, int32 id_out);

if (argc<5) {
       printf("\nProgram for reading an SDS from one HDF file and copying it to another HDF file -- with renaming\n");
       printf("usage: %s <donor HDF file> <recipient HDF file> <SDSname in donor> <SDSname in recipient>\n", argv[0]);
       printf("program also copies all local (SDS) metadata.\n");
       printf("Jim Ray, SSAI, %s\n\n", __DATE__);
       exit(0);
	}

if (!strcmp(argv[1], argv[2])) {
      printf("error: cannot use same file, '%s', as both input and output.\n", argv[1]);
      exit(-1);
                               }
if ((sd_id1 = SDstart(argv[1], DFACC_RDONLY)) == -1) {
      printf("error: file '%s' can't be opened with SDstart().\n", argv[1]);
      exit(-1);
                                                     }

/* get some OBJECTIVE idea of whether/not file exists... */
fd = fopen(argv[2], "r");
if (fd == NULL) {    /* it doesn't exist, try to create it. */
   if ((sd_id2 = SDstart(argv[2], DFACC_CREATE)) == -1) {
         printf("error: file '%s' can't be opened with SDstart().\n", argv[2]);
         SDend(sd_id1);
         exit(-1);
                                                        }
                }
else {               /* it exists! Open for read/writing... */
   fclose(fd);
   if ((sd_id2 = SDstart(argv[2], DFACC_RDWR)) == -1) {
         printf("error: file '%s' can't be opened with SDstart().\n", argv[2]);
         SDend(sd_id1);
         exit(-1);
                                                      }
     }

SDfileinfo(sd_id1, &n_sets1, &n_gattr1);
if (n_sets1 == 0) {
    printf("HDF file '%s' has no SDSs.  This program, at present, ", argv[1]);
    printf("can only process SDSs; exiting...\n");
    SDend(sd_id1);
    SDend(sd_id2);
    exit(0);
                  }



      for (j=0;j<n_sets1;j++) {
         sds_id1 = SDselect(sd_id1, j);
         for (k=0;k<5;k++) dims1[k] = 0;
         SDgetinfo(sds_id1, name1, &rank1, dims1, &number_type1, &nattr1);

         for (k=0;k<rank1;k++) {
	     SDdiminfo(SDgetdimid(sds_id1, k), dimnames[k], &dcount, &dnt, &dattr);
	    }
	 if (!strcmp(name1, argv[3])) {
               image = get_image(sds_id1,dims1,number_type1,&size);
	        

               if ((sds_out = SDcreate(sd_id2, argv[4], number_type1, rank1, dims1)) == -1) {
                  printf("Error creating SDS in output file '%s', cannot continue\n", argv[2]);
                  exit(-3);
                   }
	    
               for (k=0;k<5;k++) start[k] = 0;
               printf("Transferring SDS '%s' as '%s'...\n", argv[3], argv[4]);
               if ((SDwritedata(sds_out, start, NULL, dims1, (VOID *)image)) == -1) {
                  printf("Error writing SDS in output file '%s', cannot continue\n", argv[2]);
                  exit(-3);
                   }
               for (k=0;k<rank1;k++) {
	          SDsetdimname(SDgetdimid(sds_out, k), dimnames[k]);
	          }
               transfer_attributes(sds_id1, sds_out);
               SDendaccess(sds_out); 
 
	                              }
          SDendaccess(sds_id1); 
	                       }


SDend(sd_id2);
SDend(sd_id1);
}


void *get_image(int32 id, int32 *dims, int32 type, int *size)
{
int32 *i32_array;
uint32 *ui32_array;
float64 *f64_array;
float32 *f32_array;
int8 *i8_array;
uint8 *ui8_array;
int16 *i16_array;
uint16 *ui16_array;
char *char_array;
uchar8 *uchar_array;
int i,count;
int32 start[5];

start[0] = start[1] = start[2] = start[3] = start[4] = *size = 0;
i=0;
count = dims[i++];
while(dims[i] != 0) count *= dims[i++];
    

       switch(type) {
           case 4:
	      *size = count;
              char_array = (char *)malloc(*size);
              SDreaddata(id, start, NULL, dims, char_array);
	      return((void *)char_array);
              break;
           case 3:
	      *size = count;
              uchar_array = (uchar8 *)malloc(*size);
              SDreaddata(id, start, NULL, dims, uchar_array);
	      return((void *)uchar_array);
              break;
           case 5:
	      *size = count*sizeof(float);
              f32_array = (float *)malloc(*size);
              SDreaddata(id, start, NULL, dims, f32_array);
	      return((void *)f32_array);
              break;
           case 6:
	      *size = count*sizeof(double);
              f64_array = (double *)malloc(*size);
              SDreaddata(id, start, NULL, dims, f64_array);
	      return((void *)f64_array);
              break;
           case 20:
	      *size = count*sizeof(int8);
              i8_array = (int8 *)malloc(*size);
              SDreaddata(id, start, NULL, dims, i8_array);
	      return((void *)i8_array);
              break;
           case 21:
	      *size = count*sizeof(uint8);
              ui8_array = (uint8 *)malloc(*size);
              SDreaddata(id, start, NULL, dims, ui8_array);
	      return((void *)ui8_array);
              break;
           case 22:
	      *size = count*sizeof(int16);
              i16_array = (int16 *)malloc(*size);
              SDreaddata(id, start, NULL, dims, i16_array);
	      return((void *)i16_array);
              break;
           case 23:
	      *size = count*sizeof(uint16);
              ui16_array = (uint16 *)malloc(*size);
              SDreaddata(id, start, NULL, dims, ui16_array);
	      return((void *)ui16_array);
              break;
           case 24:
	      *size = count*sizeof(int32);
              i32_array = (int32 *)malloc(*size);
              SDreaddata(id, start, NULL, dims, i32_array);
	      return((void *)i32_array);
              break;
           case 25:
	      *size = count*sizeof(uint32);
              ui32_array = (uint32 *)malloc(*size);
              SDreaddata(id, start, NULL, dims, ui32_array);
	      return((void *)ui32_array);
              break;
                         }
}

void transfer_attributes(int32 id1, int32 id_out)
{
#define MAXLENGTH 512
int32 i, j, jj;
int32 n_attr1, n_sets1, count1, rank1, dims1[5], number_type1;
char name1[MAXLENGTH];
char attrib[MAXLENGTH], attrib1[MAXLENGTH];
char label[MAXLENGTH], tag[MAXLENGTH];
short *contents;
short ii;
char *newstring;
char *charattr;
uchar8 *ucharattr;
int16 *shortattr;
uint16 *ushortattr;
int32 *intattr;
uint32 *uintattr;
float32 *floatattr;
float64 *doubleattr;
char charatt;
uchar8 ucharatt;
int16 shortatt;
uint16 ushortatt;
int32 intatt;
uint32 uintatt;
float32 floatatt;
float64 doubleatt;

SDgetinfo(id1, name1, &rank1, dims1, &number_type1, &n_attr1);	   		   

for (j=0;j<n_attr1;j++) {
     SDattrinfo(id1, j, attrib1, &number_type1, &count1);

     switch(number_type1) {   
       case DFNT_CHAR8:      
       case DFNT_INT8: 
           if (count1 == 1) {
             SDreadattr(id1, j, &charatt);  
	     SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)&charatt);
	                    }
	   else {
	     if ((charattr = (char *)malloc((count1+1)*sizeof(char))) == NULL) {
	         printf("Out of memory, array 'charattr'\n");
		 return;
	                                                                       }
             SDreadattr(id1, j, charattr);  
	     charattr[count1] = '\0';
	     SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)charattr);
 	     free(charattr);
	        }
           break;
       case DFNT_UCHAR8:    /* treat them like integers */
       case DFNT_UINT8:   
           if (count1 == 1) {
             SDreadattr(id1, j, &ucharatt);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)&ucharatt);
 	                    }
	   else {
	     if ((ucharattr = 
	       (unsigned char *)malloc(count1*sizeof(unsigned char))) == NULL) {
	         printf("Out of memory, array 'ucharattr'\n");
		 return;
	                                                                       }
             SDreadattr(id1, j, ucharattr);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)ucharattr);
	     free(ucharattr);
	        }
           break;
       case DFNT_INT16: 
           if (count1 == 1) {
             SDreadattr(id1, j, &shortatt);  
	     SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)&shortatt);
	                    }
	   else {
	     if ((shortattr = (short *)malloc(count1*sizeof(short))) == NULL) {
	         printf("Out of memory, array 'shortattr'\n");
		 return;
	                                                                      }
             SDreadattr(id1, j, shortattr);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)shortattr);
             free(shortattr);
	        }
           break;
       case DFNT_UINT16: 
           if (count1 == 1) {
             SDreadattr(id1, j, &ushortatt);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)&ushortatt);
 	                    }
	   else {
	     if ((ushortattr = (unsigned short *)malloc(count1*sizeof(unsigned short))) == NULL) {
	         printf("Out of memory, array 'ushortattr'\n");
		 return;
	                                                                      }
             SDreadattr(id1, j, ushortattr);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)ushortattr);
	     free(ushortattr);
	        }
           break;
       case DFNT_INT32: 
           if (count1 == 1) {
             SDreadattr(id1, j, &intatt);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)&intatt);
	                   }
	   else {
	     if ((intattr = (int32 *)malloc(count1*sizeof(int32))) == NULL) {
	         printf("Out of memory, array 'intattr'\n");
		 return;
	                                                                    }
             SDreadattr(id1, j, intattr);  
	     SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)intattr);
	     free(intattr);
	        }
           break;
       case DFNT_UINT32: 
           if (count1 == 1) {
             SDreadattr(id1, j, &uintatt);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)&uintatt);
 	                     }
	   else {
	     if ((uintattr = (uint32 *)malloc(count1*sizeof(uint32))) == NULL) {
	         printf("Out of memory, array 'uintattr'\n");
		 return;
	                                                                       }
             SDreadattr(id1, j, uintattr);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)uintattr);
	     free(uintattr);
	        }
           break;
       case DFNT_FLOAT: 
           if (count1 == 1) {
             SDreadattr(id1, j, &floatatt);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)&floatatt);
 	                     }
	   else {
	     if ((floatattr = (float *)malloc(count1*sizeof(float))) == NULL) {
	         printf("Out of memory, array 'floatattr'\n");
		 return;
	                                                                      }
	                                                                      
             SDreadattr(id1, j, floatattr);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)floatattr);
	     free(floatattr);
	        }
           break;
       case DFNT_DOUBLE: 
           if (count1 == 1) {
             SDreadattr(id1, j, &doubleatt);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)&doubleatt);
 	                     }
	   else {
	     if ((doubleattr = (double *)malloc(count1*sizeof(double))) == NULL) {
	         printf("Out of memory, array 'doubleattr'\n");
		 return;
	                                                                      }
             SDreadattr(id1, j, doubleattr);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)doubleattr);
	     free(doubleattr);
	        }
           break;
                   }  /* switch(..) */
                        }  /* for (j=0;j<n_gattr;j++) */


}

