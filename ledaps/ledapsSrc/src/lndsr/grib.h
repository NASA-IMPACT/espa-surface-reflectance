/* version 1.2.1 of grib headers  w. ebisuzaki */

#ifndef INT2
#define INT2(a,b)   ((1-(int) ((unsigned) (a & 0x80) >> 6)) * (int) (((a & 0x7f) << 8) + b))
#endif

#define BDS_LEN(bds)		((int) ((bds[0]<<16)+(bds[1]<<8)+bds[2]))
#define BDS_Flag(bds)		(bds[3])

#define BDS_Grid(bds)		((bds[3] & 128) == 0)
#define BDS_Harmonic(bds)	(bds[3] & 128)

#define BDS_Packing(bds)	((bds[3] & 64) != 0)
#define BDS_SimplePacking(bds)	((bds[3] & 64) == 0)
#define BDS_ComplexPacking(bds)	((bds[3] & 64) != 0)

#define BDS_OriginalType(bds)	((bds[3] & 32) != 0)
#define BDS_OriginalFloat(bds)	((bds[3] & 32) == 0)
#define BDS_OriginalInt(bds)	((bds[3] & 32) != 0)

#define BDS_MoreFlags(bds)      ((bds[3] & 16) != 0)
#define BDS_UnusedBits(bds)	((int) (bds[3] & 15))

#define BDS_BinScale(bds)	INT2(bds[4],bds[5])

#define BDS_RefValue(bds)	(ibm2flt(bds+6))
#define BDS_NumBits(bds)	((int) bds[10])

#define BDS_DataStart(bds)      ((int) (11 + BDS_MoreFlags(bds)*3))

/* breaks if BDS_NumBits(bds) == 0 */
#define BDS_NValues(bds)        (((BDS_LEN(bds) - BDS_DataStart(bds))*8 - \
				BDS_UnusedBits(bds)) / BDS_NumBits(bds))
/*
#define BDS_NValues(bds)        ((BDS_NumBits(bds) == 0) ? 0 : \
				(((BDS_LEN(bds) - BDS_DataStart(bds))*8 - \
				BDS_UnusedBits(bds)) / BDS_NumBits(bds)))
*/

/* undefined value -- if bitmap */
#define UNDEFINED		9.999e20

/* version 1.2 of grib headers  w. ebisuzaki */

#define BMS_LEN(bms)		((bms) == NULL ? 0 : (bms[0]<<16)+(bms[1]<<8)+bms[2])
#define BMS_UnusedBits(bms)	((bms) == NULL ? 0 : bms[3])
#define BMS_StdMap(bms)		((bms) == NULL ? 0 : ((bms[4]<<8) + bms[5]))
#define BMS_bitmap(bms)		((bms) == NULL ? NULL : (bms)+6)
#define BMS_nxny(bms)		((((bms) == NULL) || BMS_StdMap(bms)) \
	? 0 : (BMS_LEN(bms)*8 - 48 - BMS_UnusedBits(bms)))
/* cnames_file.c */

/* search order for parameter names
 *
 * #define P_TABLE_FIRST
 * look at external parameter table first
 *
 * otherwise use builtin NCEP-2 or ECMWF-160 first
 */
/* #define P_TABLE_FIRST */

/* search order for external parameter table
 * 1) environment variable GRIBTAB
 * 2) environment variable gribtab
 * 3) the file 'gribtab' in current directory
 */


/* cnames.c */
/* then default values */
char *k5toa(unsigned char *pds);
char *k5_comments(unsigned char *pds);
int setup_user_table(int center, int subcenter, int ptable);


struct ParmTable {
	char *name, *comment;
};

/* version 1.4.1 of grib headers  w. ebisuzaki */
/* this version is incomplete */

#ifndef INT3
#define INT3(a,b,c) ((1-(int) ((unsigned) (a & 0x80) >> 6)) * (int) (((a & 127) << 16)+(b<<8)+c))
#endif
#ifndef INT2
#define INT2(a,b)   ((1-(int) ((unsigned) (a & 0x80) >> 6)) * (int) (((a & 127) << 8) + b))
#endif

#ifndef UINT3
#define UINT3(a,b,c) ((int) ((a << 16) + (b << 8) + (c)))
#endif

#ifndef UINT2
#define UINT2(a,b) ((int) ((a << 8) + (b)))
#endif


#define GDS_Len1(gds)		(gds[0])
#define GDS_Len2(gds)		(gds[1])
#define GDS_Len3(gds)		(gds[2])
#define GDS_LEN(gds)		((int) ((gds[0]<<16)+(gds[1]<<8)+gds[2]))

#define GDS_NV(gds)		(gds[3])
#define GDS_DataType(gds)	(gds[5])

#define GDS_LatLon(gds)		(gds[5] == 0)
#define GDS_Mercator(gds)	(gds[5] == 1)
#define GDS_Gnomonic(gds)	(gds[5] == 2)
#define GDS_Lambert(gds)	(gds[5] == 3)
#define GDS_Gaussian(gds)	(gds[5] == 4)
#define GDS_Polar(gds)		(gds[5] == 5)
#define GDS_RotLL(gds)		(gds[5] == 10)
#define GDS_Harmonic(gds)	(gds[5] == 50)
#define GDS_ssEgrid(gds)	(gds[5] == 201)	/* semi-staggered E grid */
#define GDS_fEgrid(gds)		(gds[5] == 202) /* filled E grid */

#define GDS_LatLon_nx(gds)	((int) ((gds[6] << 8) + gds[7]))
#define GDS_LatLon_ny(gds)	((int) ((gds[8] << 8) + gds[9]))
#define GDS_LatLon_La1(gds)	INT3(gds[10],gds[11],gds[12])
#define GDS_LatLon_Lo1(gds)	INT3(gds[13],gds[14],gds[15])
#define GDS_LatLon_mode(gds)	(gds[16])
#define GDS_LatLon_La2(gds)	INT3(gds[17],gds[18],gds[19])
#define GDS_LatLon_Lo2(gds)	INT3(gds[20],gds[21],gds[22])

#define GDS_LatLon_dx(gds)      INT2(gds[23],gds[24])
#define GDS_LatLon_dy(gds)      INT2(gds[25],gds[26])
#define GDS_Gaussian_nlat(gds)  ((gds[25]<<8)+gds[26])

#define GDS_LatLon_scan(gds)	(gds[27])

#define GDS_Polar_nx(gds)	((gds[6] << 8) + gds[7])
#define GDS_Polar_ny(gds)	((gds[8] << 8) + gds[9])
#define GDS_Polar_La1(gds)	(INT3(gds[10],gds[11],gds[12]))
#define GDS_Polar_Lo1(gds)	(INT3(gds[13],gds[14],gds[15]))
#define GDS_Polar_Lov(gds)	(INT3(gds[17],gds[18],gds[19]))
#define GDS_Polar_scan(gds)	(gds[27])
#define GDS_Polar_Dx(gds)	(INT3(gds[20], gds[21], gds[22]))
#define GDS_Polar_Dy(gds)	(INT3(gds[23], gds[24], gds[25]))
#define GDS_Polar_pole(gds)	((gds[26] & 128) == 128)

#define GDS_Lambert_nx(gds)	((gds[6] << 8) + gds[7])
#define GDS_Lambert_ny(gds)	((gds[8] << 8) + gds[9])
#define GDS_Lambert_La1(gds)	INT3(gds[10],gds[11],gds[12])
#define GDS_Lambert_Lo1(gds)	INT3(gds[13],gds[14],gds[15])
#define GDS_Lambert_mode(gds)	(gds[16])
#define GDS_Lambert_Lov(gds)	INT3(gds[17],gds[18],gds[19])
#define GDS_Lambert_dx(gds)	INT3(gds[20],gds[21],gds[22])
#define GDS_Lambert_dy(gds)	INT3(gds[23],gds[24],gds[25])
#define GDS_Lambert_NP(gds)	((gds[26] & 128) == 0)
#define GDS_Lambert_scan(gds)   (gds[27])
#define GDS_Lambert_Latin1(gds)	INT3(gds[28],gds[29],gds[30])
#define GDS_Lambert_Latin2(gds)	INT3(gds[31],gds[32],gds[33])
#define GDS_Lambert_LatSP(gds)	INT3(gds[34],gds[35],gds[36])
#define GDS_Lambert_LonSP(gds)	INT3(gds[37],gds[37],gds[37])

#define GDS_ssEgrid_n(gds)	UINT2(gds[6],gds[7])
#define GDS_ssEgrid_n_dum(gds)  UINT2(gds[8],gds[9])
#define GDS_ssEgrid_La1(gds)	INT3(gds[10],gds[11],gds[12])
#define GDS_ssEgrid_Lo1(gds)	INT3(gds[13],gds[14],gds[15])
#define GDS_ssEgrid_mode(gds)	(gds[16])
#define GDS_ssEgrid_La2(gds)	UINT3(gds[17],gds[18],gds[19])
#define GDS_ssEgrid_Lo2(gds)	UINT3(gds[20],gds[21],gds[22])
#define GDS_ssEgrid_di(gds)	INT2(gds[23],gds[24])
#define GDS_ssEgrid_dj(gds)	INT2(gds[25],gds[26])
#define GDS_ssEgrid_scan(gds)	(gds[27])

#define GDS_fEgrid_n(gds)	UINT2(gds[6],gds[7])
#define GDS_fEgrid_n_dum(gds)   UINT2(gds[8],gds[9])
#define GDS_fEgrid_La1(gds)	INT3(gds[10],gds[11],gds[12])
#define GDS_fEgrid_Lo1(gds)	INT3(gds[13],gds[14],gds[15])
#define GDS_fEgrid_mode(gds)	(gds[16])
#define GDS_fEgrid_La2(gds)	UINT3(gds[17],gds[18],gds[19])
#define GDS_fEgrid_Lo2(gds)	UINT3(gds[20],gds[21],gds[22])
#define GDS_fEgrid_di(gds)	INT2(gds[23],gds[24])
#define GDS_fEgrid_dj(gds)	INT2(gds[25],gds[26])
#define GDS_fEgrid_scan(gds)	(gds[27])


#define GDS_Merc_nx(gds)	UINT2(gds[6],gds[7])
#define GDS_Merc_ny(gds)	UINT2(gds[8],gds[9])
#define GDS_Merc_La1(gds)	INT3(gds[10],gds[11],gds[12])
#define GDS_Merc_Lo1(gds)	INT3(gds[13],gds[14],gds[15])
#define GDS_Merc_mode(gds)	(gds[16])
#define GDS_Merc_La2(gds)	INT3(gds[17],gds[18],gds[19])
#define GDS_Merc_Lo2(gds)	INT3(gds[20],gds[21],gds[22])
#define GDS_Merc_Latin(gds)	INT3(gds[23],gds[24],gds[25])
#define GDS_Merc_scan(gds)	(gds[27])
#define GDS_Merc_dx(gds)        INT3(gds[28],gds[29],gds[30])
#define GDS_Merc_dy(gds)        INT3(gds[31],gds[32],gds[33])

/* rotated Lat-lon grid */

#define GDS_RotLL_nx(gds)	UINT2(gds[6],gds[7])
#define GDS_RotLL_ny(gds)	UINT2(gds[8],gds[9])
#define GDS_RotLL_La1(gds)	INT3(gds[10],gds[11],gds[12])
#define GDS_RotLL_Lo1(gds)	INT3(gds[13],gds[14],gds[15])
#define GDS_RotLL_mode(gds)	(gds[16])
#define GDS_RotLL_La2(gds)	INT3(gds[17],gds[18],gds[19])
#define GDS_RotLL_Lo2(gds)	INT3(gds[20],gds[21],gds[22])
#define GDS_RotLL_dx(gds)       INT2(gds[23],gds[24])
#define GDS_RotLL_dy(gds)       INT2(gds[25],gds[26])
#define GDS_RotLL_scan(gds)	(gds[27])
#define GDS_RotLL_LaSP(gds)	INT3(gds[32],gds[33],gds[34])
#define GDS_RotLL_LoSP(gds)	INT3(gds[35],gds[36],gds[37])
#define GDS_RotLL_RotAng(gds)	ibm2flt(&(gds[38]))

/* index of NV and PV */
#define GDS_PV(gds)		((gds[3] == 0) ? -1 : (int) gds[4] - 1)
#define GDS_PL(gds)		((gds[4] == 255) ? -1 : (int) gds[3] * 4 + (int) gds[4] - 1)

/* version 3.4 of grib headers  w. ebisuzaki */
/* this version is incomplete */

#ifndef INT2
#define INT2(a,b)   ((1-(int) ((unsigned) (a & 0x80) >> 6)) * (int) (((a & 0x7f) << 8) + b))
#endif

#define PDS_Len1(pds)		(pds[0])
#define PDS_Len2(pds)		(pds[1])
#define PDS_Len3(pds)		(pds[2])
#define PDS_LEN(pds)		((int) ((pds[0]<<16)+(pds[1]<<8)+pds[2]))
#define PDS_Vsn(pds)		(pds[3])
#define PDS_Center(pds)		(pds[4])
#define PDS_Model(pds)		(pds[5])
#define PDS_Grid(pds)		(pds[6])
#define PDS_HAS_GDS(pds)	((pds[7] & 128) != 0)
#define PDS_HAS_BMS(pds)	((pds[7] & 64) != 0)
#define PDS_PARAM(pds)		(pds[8])
#define PDS_L_TYPE(pds)		(pds[9])
#define PDS_LEVEL1(pds)		(pds[10])
#define PDS_LEVEL2(pds)		(pds[11])

#define PDS_KPDS5(pds)		(pds[8])
#define PDS_KPDS6(pds)		(pds[9])
#define PDS_KPDS7(pds)		((int) ((pds[10]<<8) + pds[11]))

/* this requires a 32-bit default integer machine */
#define PDS_Field(pds)		((pds[8]<<24)+(pds[9]<<16)+(pds[10]<<8)+pds[11])

#define PDS_Year(pds)		(pds[12])
#define PDS_Month(pds)		(pds[13])
#define PDS_Day(pds)		(pds[14])
#define PDS_Hour(pds)		(pds[15])
#define PDS_Minute(pds)		(pds[16])
#define PDS_ForecastTimeUnit(pds)	(pds[17])
#define PDS_P1(pds)		(pds[18])
#define PDS_P2(pds)		(pds[19])
#define PDS_TimeRange(pds)	(pds[20])
#define PDS_NumAve(pds)		((int) ((pds[21]<<8)+pds[22]))
#define PDS_NumMissing(pds)	(pds[23])
#define PDS_Century(pds)	(pds[24])
#define PDS_Subcenter(pds)	(pds[25])
#define PDS_DecimalScale(pds)	INT2(pds[26],pds[27])
#define PDS_Year4(pds)   (pds[12] + 100*(pds[24] - (pds[12] != 0))) 
/* new but no good: 24-AUG-00 #define PDS_Year4(pds)     (pds[12] + 100*(pds[24] - 1)) */

/* various centers */
#define NMC			7
#define ECMWF			98

/* ECMWF Extensions */

#define PDS_EcLocalId(pds)	(PDS_LEN(pds) >= 41 ? (pds[40]) : 0)
#define PDS_EcClass(pds)	(PDS_LEN(pds) >= 42 ? (pds[41]) : 0)
#define PDS_EcType(pds)		(PDS_LEN(pds) >= 43 ? (pds[42]) : 0)
#define PDS_EcStream(pds)	(PDS_LEN(pds) >= 45 ? (INT2(pds[43], pds[44])) : 0)

#define PDS_EcENS(pds)		(PDS_LEN(pds) >= 52 && pds[40] == 1 && \
				pds[43] * 256 + pds[44] == 1035 && pds[50] != 0)
#define PDS_EcFcstNo(pds)	(pds[50])
#define PDS_EcNoFcst(pds)	(pds[51])


/* NCEP Extensions */

#define PDS_NcepENS(pds)	(PDS_LEN(pds) >= 44 && pds[25] == 2 && pds[40] == 1)
#define PDS_NcepFcstType(pds)	(pds[41])
#define PDS_NcepFcstNo(pds)	(pds[42])
#define PDS_NcepFcstProd(pds)	(pds[43])


#define VERSION "v1.7.0b1 (8-24-98) Wesley Ebisuzaki"

#define CHECK_GRIB

#define MSEEK 1024
#define BUFF_ALLOC0	40000


#ifndef min
#define min(a,b)  ((a) < (b) ? (a) : (b))
#define max(a,b)  ((a) < (b) ? (b) : (a))
#endif

enum Def_NCEP_Table {rean, opn, rean_nowarn, opn_nowarn};

#ifndef DEF_T62_NCEP_TABLE
#define DEF_T62_NCEP_TABLE	rean
#endif
//enum Def_NCEP_Table def_ncep_table = DEF_T62_NCEP_TABLE;


/* Function prototypes */
int read_grib_array(FILE *input, char *what, char *where, int *nx, int *ny,
    float **narray);
int read_grib_date(FILE *input, char *what, char *where, char *date);

