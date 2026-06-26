#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <libgen.h>
#include <time.h>

#define MAXLENGTH 1000 

int ad_hoc_XML_reader(char *acquisition_date)
{

/* "ad-hoc" with respect to reading only PRODUCT_START_TIME from MTD_MSIL1C.xml files -- 
 * to be discarded when (if) the XML-handling in espa_product_formatter starts reading it.
 * 26-FEB-26,  JPR */

FILE *fd;
long k, l;
int i, j, ii, n, n_lines, INGM, start = 0;
void get_a_line(char *text, int lengthoftext, int *start, char line[MAXLENGTH]);
char *text;
char line[MAXLENGTH];
char xxxxxxxxxx[MAXLENGTH*10];  /* need some kind of space between line[] and subline[] on the cluster...  */
char tmpstr[MAXLENGTH], subline[50];
int globstart;

/* Global stuff to look for */
int global_stuff = 1;
char sgreads[1][30] = {
        "<PRODUCT_START_TIME>",
	};
char egreads[1][30] = {
        "</PRODUCT_START_TIME>",
 	};

if ((fd = fopen("MTD_MSIL1C.xml", "r")) == NULL) {
   printf("Error opening MTD_MSIL1C.xml; exiting...\n");
   exit(-1);
   }
   
fseek(fd,0L,SEEK_SET);
fseek(fd,0L,SEEK_END);
k = ftell(fd);
/*printf("File is %ld bytes long\n", k);*/
fseek(fd,0L,SEEK_SET);

text = (char *)malloc(k*sizeof(char));
if ( fread(text, 1, k, fd)  != k) {
    printf("Error reading from file MTD_MSIL1C.xml, cannot continue\n");
    free(text);
    exit(-4);
    }
fclose(fd);
    
/* get number of lines ('\n') in file */
n_lines = 0;
for(l=0L;l<k;l++) if ( text[l] == '\n') n_lines++;
    
INGM = 0;
for (i=0;i<n_lines;i++) {
   for(j=0;j<MAXLENGTH;j++) line[j] = '\0';
   get_a_line(text, (int)k, &start, line);
   n = strlen(line);
   line[n-1] = '\0';  /* get rid of carriage return */
   if ( strstr(line, "<Product_Info>")) INGM = 1;
   else if ( strstr(line, "</Product_Info>")) INGM = 0;
   
   if (INGM) {
      /*printf("%d -- %s -- \n", i, line);*/
      /* Simplest ever parse -- crude as hell! */
      for(j=0;j<global_stuff;j++) {
         if ( (strstr(line, sgreads[j])) && (strstr(line, egreads[j]))) {
	    /* strchr(), strrchr() don't seem to work...? */
	    globstart = 0;
	    for(ii=0;ii<MAXLENGTH;ii++) tmpstr[ii] = '\0'; 
	    
	    for(ii=0;ii<strlen(line);ii++) {
	       if (line[ii] == '>') {
	          globstart = ii+1;
		  break;
		  }
	      }
	    for(ii=globstart;ii<strlen(line);ii++) {
	       tmpstr[ii-globstart] = line[ii];	       
	       if (line[ii] == '<') {
	          tmpstr[ii-globstart] = '\0';
		  break;
		  }
	      }
	    switch (j) {
	       case 0: strcpy(acquisition_date, tmpstr); break;
	      }  
	   }
        }  /* for each global metadatum we want */
      }  /* if in global metadata */
      
   }  /* for (i=0;i<n_lines;i++) */

return(0);
}

void get_a_line(char *text, int lengthoftext, int *start, char line[MAXLENGTH])
{
int i=0;
int where;
int getout = 0;

where = *start;
line[i] = '\0';
if (where >= lengthoftext) return;
if (text[where] == '\0') return;
else {
   while (getout == 0) {  
	 if ((text[where] == '\0')||(text[where] == '\n')||(i>=MAXLENGTH)||(where>=lengthoftext)) getout=1;
         line[i++] = text[where++];
      }
      
    *start = where;
    line[i] = '\0';
    return;
     }
}

