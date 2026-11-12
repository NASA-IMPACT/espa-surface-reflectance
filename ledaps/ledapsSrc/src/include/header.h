#ifndef HEADER_H
#define HEADER_H

#include "edips.h"
#include "jsc.h"
#include "ltwg.h"
#include "tipsp.h"
#include "siat.h"

/*
 *	ANSI X3.27-1978  American National Standard Magnetic
 *		Tape Labels and File Structure for Information Interchange
 *			pp 13-22
 */

/*
 *	Label #1 Information
 */
typedef struct
	{
	char	file_ident[17];		/* file identifier */
	char	file_set_ident[6];	/* file-set identifier */
	char	file_sect_no[4];	/* file section number */
	char	file_seq_no[4];		/* file sequence number */
	char	gener_no[4];		/* generation number */
	char	gener_ver_no[2];	/* generation version number */
	char	create_date[6];		/* creation date " yyddd" */
	char	expire_date[6];		/* expiration date " yyddd" */
	char	accessibility;		/* accessability */
	char	blk_count[6];		/* block count */
	char	system_code[13];	/* system code */
	char	hdr1_spaces[7];		/* spaces */
	} INFO1;			/* each field is dimensioned to be
					   one char longer than is written to
					   to tape, in order to allow them to
					   be null-terminated strings	      */
/*
 *	Label #2 Information
 */
typedef struct
	{
	char	rec_format;		/* "F" - fixed length */
					/* "D" - variable length */
					/* "S" - spanned */
	char	blk_length[5];		/* block length in characters */
	char	rec_length[5];		/* record length in characters */
	char	timestamp[11];		/* file create date YYDDDHHMMSS */
	char	hdr2_spaces1[24];	/* spaces */
	char	buffer_offset_length[2];/* buffer-offset length */
	char	hdr2_spaces2[28];	/* spaces */
	} INFO2;
/*
 *	File-Header #1 Layout
 */
struct HDR1
	{
	char	label_id[3];		/* "HDR" */
	char	label_no;		/* "1" */
	INFO1	info;
	} ;
/*
 *	File-Header #2 Layout
 */
struct HDR2
	{
	char	label_id[3];		/* "HDR" */
	char	label_no;		/* "2" */
	INFO2	info;
	} ;
/*
 *	File-Header #n Layout
 */
struct HDRN
	{
	char	label_id[3];		/* "HDR" */
	char	label_no;		/* "3", "4", ... "9" */
	char	hdrn_data[76];		/* reserved for system use */
	} ;
/*
 *	User Header Labels
 */
struct UHLA
	{
	char	label_id[3];		/* "UHL" */
	char	label_no;		/* "a" character */
	char	uhla_data[76];		/* reserved for user application */
	} ;


int c_edips_anc
(
    char *outh,                 /* image host file name */
    long outsys,                /* output system */
    struct EDPSAC *anc,         /* EDIPS ancillary structure */
    struct EDPSAC_EPH *anceph,  /* EDIPS ancillary structure--Ephemeris Data */
    struct EDPSAC_ATT *ancatt,  /* EDIPS ancillary structure--Attitude Data */
    struct EDPSAC_RS *ancutm,   /* EDIPS ancillary structure--Horizontal and
                                   Vertical resampling for UTM/PS data     */
    struct EDPSAC_RS *ancsom    /* EDIPS ancillary structure--Horizontal and
                                   Vertical resampling for SOM/HOM data    */
);


int c_convertab
(
   char  *header_a,                  /* fast format a header bfer     */
   char  *header_b                   /* fast format b header buffer     */
);

int c_edips_ann
(
    char *outh,
    struct EDIPSA *annot
);

int c_edips_hdr
(
    char *outh,                 /* host image name                      */
    long outsys,                /* output system                        */
    struct EDIPSH *hdr          /* EDIPS header structure               */
);

int c_edips_trl
(
    char *outh,
    struct EDIPST *trl
);

int c_fhdrrd
(
    char *header,               /* the pointer to the fast fmt header       */
    long field,                 /* the field to be read                     */
    char *data                  /* address to place the data read           */
);

int c_fhdrup
(
    char *header,               /* the pointer to the fast fmt header       */
    long field,                 /* the field to be updated                  */
    char *data                  /* pointer to the data to place in the field*/
);

int c_jsc_anc
(
    struct JSCIN *jscin,        /* TIPSP sturcture declaration          */
    char *name,                 /* image name                           */
    long *nl                    /* line number to seek to               */
);

int c_jsc_lhdr
(
    char *outh,         /* image host file name                         */
    struct JSCLH *lhdr  /* JSC Landsat header structure                 */
);

int c_ltwg_eaa
(
    struct LTWGEAA *eaa, /* OUTPUT - LTWG ephomeris volume descriptor */
                         /* attitude ancillary rec structure 	*/
    FILE   *fin,         /* INPUT - File pointer			*/
    int    *size         /* INPUT - record size in bytes		*/
);

int c_ltwg_fpr
(
    struct LTWGFPR *fpr,    /* OUTPUT - LTWG volume descriptor rec structure*/
    FILE *fin               /* INPUT - File pointer			*/
);

int c_ltwg_idr
(
    struct LTWGIDR *idr,    /* OUTPUT - LTWG image data rec structure	*/
    FILE *fin               /* INPUT - File pointer				*/
);

int c_ltwg_ifd
(
    struct LTWGIFD *ifd,    /* OUTPUT - LTWG imagery file descriptor  	*/
    FILE *fin               /* INPUT - File pointer				*/
);

int c_ltwg_lfd
(
    struct LTWGLFD *lfd,    /* OUTPUT - LTWG leader file descriptor  	*/
    FILE *fin               /* INPUT - File pointer				*/
);

int c_ltwg_mpr
(
    struct LTWGMPR *mpr,    /* OUTPUT - LTWG map projection rec structure	*/
    FILE *fin,              /* INPUT - File pointer				*/
    long *size              /* INPUT - number of bytes to read for file	*/
);

int c_ltwg_mtr
(
    struct LTWGMTR *mtr, /* OUTPUT - LTWG mission telemetry anc rec structure */
    FILE *fin,           /* INPUT - File pointer				*/
    long *size           /* INPUT - record size in bytes			*/
);

int c_ltwg_rcr
(
    struct LTWGRCR *rcr, /* OUTPUT - LTWG radiometric calibration rec struct */
    FILE *fin,           /* INPUT - File pointer				*/
    long *size           /* INPUT - Record size in bytes			*/
);

int c_ltwg_rjr
(
    struct LTWGRJR *rjr,    /* OUTPUT - LTWG raw jitter meas rec structure	*/
    FILE *fin,              /* INPUT - File pointer				*/
    long *size              /* INPUT - record size in bytes			*/
);

int c_ltwg_sfd
(
    struct LTWGSFD *sfd,    /* OUTPUT - LTWG supplemental file descriptor    */
    FILE *fin               /* INPUT - File pointer				*/
);

int c_ltwg_shr
(
    struct  LTWGSHR *shr,   /* OUTPUT - LTWG volume descriptor rec structure */
    FILE    *fin,           /* INPUT - File pointer				*/
    long    *size           /* INPUT - number of bytes to read in for file	*/
);

int c_ltwg_sih
(
    struct LTWGSIH *sih,    /* OUTPUT - LTWG interval header rec structure	*/
    FILE   *fin,            /* INPUT - File pointer				*/
    long   *size            /* INPUT - record size in bytes			*/
);

int c_ltwg_str
(
    struct LTWGSTR *str,    /* OUTPUT - LTWG supplemental text rec structure */
    FILE *fin               /* INPUT - File pointer				*/
);

int c_ltwg_tfd
(
    struct LTWGTFD *tfd,    /* OUTPUT - LTWG trailer file des rec structure	*/
    FILE *fin               /* INPUT - File pointer				*/
);

int c_ltwg_tmhk
(
    struct LTWGTMHK *tmhk, /* OUTPUT - LTWG TM house keeping rec structure	*/
    FILE   *fin,           /* INPUT - File pointer				*/
    long   *size           /* INPUT - record size in bytes			*/
);

int c_ltwg_tr
(
    struct LTWGTR *tr,  /* OUTPUT - LTWG trailer rec structure		*/
    FILE *fin,          /* INIPUT - File pointer			*/
    long *size          /* INIPUT - record size in bytes		*/
);

int c_ltwg_triv
(
    struct LTWGTRIV *triv, /* OUTPUT - LTWG volume descriptor rec structure */
    FILE *fin              /* INPUT - File pointer				*/
);

int c_ltwg_vdr
(
    struct LTWGVDR *vdr, /* OUTPUT - LTWG volume descriptor rec structure */
    FILE *fin            /* INPUT - File pointer				*/
);

int c_jsc_radtran
(
    char *outh,          /* image host file name                         */
    struct JSCRT *jscrt4,/* JSC Landsat header Radiometric Transformation
                            band 4                                       */
    struct JSCRT *jscrt5,/* JSC Landsat header Radiometric Transformation
                            band 5                                       */
    struct JSCRT *jscrt6,/* JSC Landsat header Radiometric Transformation
                            band 6                                       */
    struct JSCRT *jscrt7,/* JSC Landsat header Radiometric Transformation
                            band 7                                       */
    struct JSCRT *jscrt8 /* JSC Landsat header Radiometric Transformation
                            band 8                                       */
);

int c_jsc_uhdr
(
    char *outh,         /* image host file name                         */
    struct JSCUH *uhdr  /* JSC Universal header structure               */
);

long get_long_var
(
    unsigned char *val
);

int c_open_header
(
    char *buf,                  /* header buffer                        */
    char *name,                 /* image name                           */
    long *flag                  /* flag to determine read or write acc  */
);

int c_read_annot
(
    struct TIPSP *tipsp,        /* TIPSP sturcture declaration          */
    char *name                  /* image name                           */
);

int c_read_head
(
    struct TIPSP *tipsp,
    char *name                  /* image name                           */
);

int c_read_siat
(
    struct SIAT *siat,      /* SIAT sturcture declaration		*/
    char  *name             /* image name				*/
);

int c_read_trail
(
    struct TIPSP *tipsp,        /* TIPSP sturcture declaration          */
    char *name                  /* image name                           */
);

int c_readfpg
(
    unsigned char *buf,     /* value in FPG format			*/
    double	*tempdbl        /* converted value of double		*/
);

int cdbl2ftn
(
    char *buffer,               /* the buffer to print the string to */
    double value                /* the value to convert */
);

int fstr2cdbl
(
    char *buffer,               /* the buffer to with the string */
    double *value               /* the return address */
);

int convert_dms 
(
    char **strings,    /* pointer to strings     */
    long proj_number,  /* USGS projection number */
    long method        /* converstion method.  0 = old->new,  !0 = new->old */
);

int deg2dms
(
    char *buffer,               /* the buffer to write the string to */
    double value,               /* the value to convert */
    long mode                   /* latitude = 0, longitude = 1 */
);

int dms2deg
(
    char *buffer,               /* the buffer that contains the string */
    double *value,              /* the address to place the result */
    long mode                   /* latitude = 0, longitude = 1 */
);

int c_fhdfld
(
    long field,                 /* field to get information about     */
    long fdata[5]               /* array to fill with the information */
);

int meta_init
(
    char *file
);

int meta_insert
(
    int field,
    char *data
);

char *meta_getfield
(
    int field
);

void meta_setfile
(
    char *file
);

int meta_write(void);

void meta_writearchv
(
    char *tapeloc[],
   char *format,
   int file_number
);

void meta_read(void);

int meta_cleanup(void);

/* NOTE: template_header used to be template, but since template is a 
         reserved word it was changed */
void template_header
(
    char *buf                   /* header buffer                        */
);

#endif
