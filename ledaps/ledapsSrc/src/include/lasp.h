
/*******************************************************************************
NAME:     lasp.h 
PURPOSE:  prototypes of LAS support functions.  Must be #include'd after las.h
          The bulk of this file was created by running makep.csh and makep.C
	  makep.C is under $CVSROOT/bin/prototypes/


FILE HISTORY
--------------------------------------------------------------------------------
  Date       Programmer                          Changes 
---------  ---------------  ----------------------------------------------------
 Aug 1994  David Zokaites   original development
 Jan 1995  David Zokaites   revised as html files were updated
 Apr 1995  David Zokaites   prototyped some TAE supports 
 Jul 1995  David Zokaites   added prototype of c_lasmsg()
 Jun 1996  David Zokaites   used const keyword, converted to C++
*******************************************************************************/




extern "C"
{


// prototype the commonly used TAE support functions
   struct VARIABLE *p_find (const struct PARBLK *vptr, const char name[]); 
   int q_intg (struct PARBLK *vptr, char name[], int count, int intg[], 
       int mode); 
   int q_out (struct PARBLK *vptr); 
   int q_real (struct PARBLK *vptr, char name[], int count, double real[], 
       int mode); 
   int q_string (struct PARBLK *vptr, char name[], int count, char *vector[],
       int mode); 




// prototypes of LAS supports by category


// alias
   int c_getalias (const char *instring, char *outstring);

// ddr
   void c_calcor (const double *sl, const double *ss, const double *nl,
                 const double *ns, const struct DDR *ddr, double *upleft,
                 double *upright, double *loleft, double *loright);
   int c_comsys (const char *hname, long *flag);
   int c_getbdf (const struct BDDR *bddr, const long *key,
                unsigned char *value);
   int c_getbdr (const char *hname, struct BDDR *bddr, const long *band);
   int c_getddf (const struct DDR *ddr, const long *key, unsigned char *value);
   int c_getddr (const char *hname, struct DDR *ddr);
   void c_intbdr (struct BDDR *bddr);
   void c_intddr (struct DDR *ddr);
   int c_lused (const char *hname);
   int c_minmax (const char hname[][CMLEN], const long *nimg,
                const long window[][4], const long bands[][MAXBND + 1],
                const long *nbands, double *minval, double *maxval);
   int c_putbdf (struct BDDR *bddr, const long *key,
                const unsigned char *value);
   int c_putbdr (const char *hname, const struct BDDR *bddr);
   int c_putddf (struct DDR *ddr, const long *key, const unsigned char *value);
   int c_putddr (const char *hname, const struct DDR *ddr);
   int c_rotcor (const double proj_coord[8], const long image_coord[8],
                const double window[4], double window_proj_coord[8]);
   int c_rtcoef (const double proj_coord[8], const long image_coord[8],
                double coef[6], const long *flag);
   int c_rtiddr (struct DDR *inddr, const long *rot);
   int c_svmnmx (const char *hname, const long *totbnd, const double minval[],
                const double maxval[], const char *m_flag);
   int c_upddr (const struct PARBLK *parblk, const char host_in[][CMLEN],
               const char *host_out, const long *nmrimg, const long window[][4],
               const long bands[][MAXBND +1], const long nbands[],
               const double *line_inc, const double *samp_inc,
               const long *out_nbands, const long *cflag);

// dltf
   int c_ckdltf (char *host, char *key, long *fd);
   int c_cldltf (const long *fd, const char *host);
   int c_r1dltf (long *fd, const char *host, char *key, char *descr,
                struct ASDLTF *asdltf, long *dtlen, char *dtype);
   int c_r2dltf (const long *fd, const char *key, const long *band,
                long *startpt, long *endpt, const long *dlen,
                unsigned char *lutdata, double *slope, double *offset);
   int c_w1dltf (const long *fd, const char *key, const char *descr,
                const char *type, const struct ASDLTF *asdltf);
   int c_w2dltf (const long *fd, const char *type, const char *key,
                const double *slope, double *offset, const long *startpt,
                const long *endpt, const long *lutnum,
                const unsigned char *lutdata);

// error
   void c_errmsg (const char *message, const char *key, const long *vrity);

// filenm
   int c_cmhost (const char *lname, char *hname);
   int c_cmtae (const char *hname, char *lname);

// georeg
#ifdef GEOMPAK_H  // these prototypes won't compile unless geompak.h is included
   void for_init (const long outsys, const long outzone,
                 const double outparm[15], const long outdatum,
                 const char fn27[], const char fn83[], long *iflg,
                 long (*for_trans[MAXPROJ + 1])());
   void inv_init (const long insys, const long inzone, const double inparm[15],
                 const long indatum, const char fn27[], const char fn83[],
                 long *iflg, long (*inv_trans[MAXPROJ + 1])());
   int c_steplr (const long *npts, const long *indep,
                const double (*data)[MAX_COEF], double *coefs,
                const double *alpha);
#endif 
   int c_ckptid (const char *ptid);
   int c_ckunit (const char *inunit, long *outunit);
   int c_decdeg (double *angle, const char *coform, const char *type);
   int c_degdms (const double *deg, double *dms, const char *code,
                const char *check);
   int c_digpt (double *yx, const long *flag);
   double c_eval (const long *degree, const double *coef, const double *x,
                 const double *y);
   int c_file_corner (const double in_proj_corners[8],
                     const double ls_corners[8], const long window[4],
                     double out_proj_corners[8]);
   void gctp (const double incoor[2], const long *insys, const long *inzone,
             const double inparm[15], const long *inunit, const long *indatum,
             const long *ipr, const char efile[], const long *jpr,
             const char pfile[], double outcoor[2], const long *outsys,
             const long *outzone, const double outparm[15], const long *outunit,
             const long *outdatum, const char fn27[], const char fn83[],
             long *iflg);
   int c_getgrd (const char *hname, struct GEOGRID *grid);
   int c_getrwt (const char *hname, struct RWTAB *table);
   int c_grderr (struct GEOGRID *grid, const double *tol);
   int c_grdred (struct GEOGRID *grid, const double *tol);
   int c_invert (const double *incoef, double *outcoef);
   int c_mapedg (const double *ul_long, const double *lr_long,
                const double *ul_lat, const double *lr_lat, const long *inunits,
                const long *proj_code, const long *zone,
                const double parout[15], const long *units, double *pxmin,
                double *pxmax, double *pymin, double *pymax);
   int c_proj (const long *inproj, const long *inunit, const long *inzone,
              const double inparm[15], long *outproj, long *outunit,
              long *outzone, double outparm[15], const double *inx,
              const double *iny, double *outx, double *outy);
   int c_projon (const long *inproj, const long *inunit, const long *inzone,
                const long *indatum, const double inparm[15], long *outproj,
                long *outunit, long *outzone, long *outdatum,
                double outparm[15], const long prtprm[2], const char *fname);
   void c_prostr (const long *proj_code, char *string);
   int c_ptclse (int *ptdesc);
   int c_ptio (int *ptdesc, long *ptdata);
   int c_ptopen (int *ptdesc, const char *hname, const long *pttype,
                const long *acc, long *pthdr);
   int c_putgrd (const char *hname, const struct GEOGRID *grid);
   int c_putrwt (const char *hname, const struct RWTAB *table);
   int c_regmap (const char *mapfile, char *mapid, char *units);
   int xxxfor (const double lon, const double lat, double *x, double *y);
   int xxxinv (const double x, const double y, double *lon, double *lat);

// gofio
#ifdef GOF_H // these prototypes won't compile unless gof.h is included
   int c_attchk (const struct HEADER *hdr, const char inatt[][ATTIN],
                char outatt[][ATTVAL], int index[MAXATT], const int *search);
#endif 
   int c_gofcls (const int *fd);
   int c_gofopn (int *fd, char *hostname, char *goflasnm, const char *type,
                const int *access, int *cflag);
   int c_gofrd (const int *fd, int *obyte, int *nbyte, const int *type,
               union GOF *gofptr);
   int c_gofwr (const int *fd, int *nbytes, const int *type,
               const union GOF *gofptr, const int *addrec);

// imageio
   unsigned char * c_ealloc (const int *fdarray, long *buffersz);
   int c_eclose (const int *fdesc);
   int c_egclse (const int *gdesc);
   int c_egroup (const int *fdarray, int *gdesc, const unsigned char *buffer,
                const long *buffersz);
   int c_eopenr (int *fdesc, const char *fname, long *nbands, const long *sl,
                const long *ss, long *nl, long *ns, const long *acc,
                long *dtype, const long opt[], const long *nline);
   int c_eopens (int *fdesc, const char *fname, long bands[], const long *sl,
                const long *ss, long *nl, long *ns, const long offset[],
                const long *acc, long *dtype, const long *nbuf,
                const long opt[]);
   int c_eread (const int *fdesc, const long *eband, const long *eline,
               const unsigned char *euserbuf, const long *nl_blk);
   int c_estep (const int *gdesc);
   int c_ewrite (const int *fdesc, const long *eband, const long *eline,
                const unsigned char *euserbuf, const long *nl_blk);
   int get_buf (unsigned char *buffer, const long *nlbuf, const long *nsbuf,
               long *nlav, long *nsav, const long *sli, const long *ssi,
               const long *nli, const long *nsi, const long *nlb,
               const long *nsb, long *slcb, long *sscb, const long *nextflg,
               const long *band, const long *fill, int *fdesc, const char *host,
               long *dtype);
   int put_buf (const unsigned char *buffer, const long *nlbuf,
               const long *nsbuf, const long *nlav, const long *nsav,
               long *slbuf, long *ssbuf, const long *nlti, const long *nsti,
               const long *nextflg, const long *band, const long *bands,
               const int *fdesc, const char *host, const long *dtype,
               const long *wtype);

// lablserv
   int c_lsclos (const long *fd, const char *hostname, const long *action);
   int c_lsdel (const long *fd, const char *key);
   void c_lsmknm (const char *inname, const char *ext, char *outname);
   int c_lsopen (long *fd, const char *hostname, const long *access);
   int c_lsread (const long *fd, char *key, long *clen, long *dlen, char *cbuf,
                unsigned char *dbuf, char *dtype);
   int c_lsrepl (const long *fd, const char *key, const long *clen,
                const long *dlen, const char *cbuf, const unsigned char *dbuf,
                const char *dtype);
   int c_lsrew (const long *fd);
   int c_lsstat (const char *hostname, const long *mode);
   int c_lswrit (const long *fd, const char *key, const long *clen,
                const long *dlen, const char *cbuf, const unsigned char *dbuf,
                const char *dtype);

// labltabl
   void close_tab (struct TAB_DEF *tab);
   struct VECTOR * get_field_ptr (struct TAB_DEF *tab, char *name);
   void get_lt_date (struct VECTOR *vec, struct LT_DATE *date);
   void get_lt_time (struct VECTOR *vec, struct LT_TIME *time);
   void get_matrix (struct VECTOR *vec, long odtype, unsigned char *matrix);
   int get_record (struct TAB_DEF *tab);
   void get_vector (struct TAB_DEF *tab);
   char * make_lt_name (const char *innam, const long num);
   void open_tab (struct TAB_DEF *tab, char *hname, long access);
   void put_lt_date (struct VECTOR *vec, struct LT_DATE *date);
   void put_lt_time (struct VECTOR *vec, struct LT_TIME *time);
   void put_matrix (struct VECTOR *vec, long idtype, struct MATRIX *size,
                   unsigned char *matrix);
   void put_record (struct TAB_DEF *tab);
   void put_vector (struct TAB_DEF *tab, char *name, long dtype, char *desc,
                   long nrow, long ncol, char sep);
   void upcase (char *buf);

// math
   void c_qrdecomp (double *matrix, const long n_row, const long n_col,
                   double *work, const long flag);
   void c_qrsolve (double *matrix, const long n_row, const long n_col,
                  double *work, double *result, long flag);
   double c_randnm (long *seed);

// pixman
   int c_pxadd (const long *ns, const unsigned char *in1,
               const unsigned char *in2, unsigned char *out,
               const float *scale1, const float *scale2, const float *offset,
               const long *dtype);
   int c_pxbtan (const unsigned char *in, const unsigned char *map,
                unsigned char *out, const long *ns, const long *dtype);
   int c_pxbtor (const unsigned char *in, const unsigned char *map,
                unsigned char *out, const long *ns, const long *dtype);
   int c_pxbxor (const unsigned char *in, const unsigned char *map,
                unsigned char *out, const long *ns, const long *dtype);
   int c_pxcomp (const long *in, long *out, const long *ns, const long *dtype);
   int c_pxconv (const long *ftyp, const long *totyp,
                const unsigned char *frombuf, unsigned char *tobuf,
                const long *size);
   void c_pxcopy (const unsigned char *frombuf, unsigned char *tobuf,
                 const long *dtype, const long *ns);
   int c_pxdiv (const long *ns, const unsigned char *in1,
               const unsigned char *in2, unsigned char *out, const float *scale,
               const float *offset, const long *dtype, const long *errval);
   int c_pxe10 (const long *in, double *out, const long *ns, const long *dtype);
   int c_pxe255 (const unsigned char *in, double *out, const long *ns,
                const long *dtype);
   int c_pxe32k (const short *in, double *out, const long *ns,
                const long *dtype);
   void c_pxfill (unsigned char *tobuf, const long *dtype, const long *ns,
                 const unsigned char *fillval);
   int c_pxhiss (const long *ns, const unsigned char *buf, const long *dtype,
                long *hist, long *count, const long *skip);
   int c_pxhist (const long *ns, const unsigned char *buf, const long *dtype,
                long *hist, long *count);
   int c_pxindr (const long *in, const long *val, const long *ns, long *index,
                const long *dtype);
   int c_pxinx (const long *in, const long *val, const long *ns, long *index,
               const long *dtype);
   int c_pxinxd (const long *in, const long *val, const long *ns, long *index,
                const long *dtype);
   int c_pxinxr (const long *in, const long *val, const long *ns, long *index,
                const long *dtype);
   int c_pxl10 (const long *in, double *out, const long *ns, const long *dtype);
   int c_pxl255 (const unsigned char *in, double *out, const long *ns,
                const long *dtype);
   int c_pxland (const long *in, const long *map, long *out, const long *ns,
                const long *dtype);
   int c_pxline (const long *in, long *out, const long *ns, const long *dtype,
                const float *scale, const float *offset);
   int c_pxlkup (const long *in, long *out, const long *ns, const long *dtype,
                const long *lut, const long *luttype);
   int c_pxlor (const long *in, const long *map, long *out, const long *ns,
               const long *dtype);
   int c_pxlxor (const long *in, const long *map, long *out, const long *ns,
                const long *dtype);
   int c_pxmax (const unsigned char *in, const unsigned char *map,
               unsigned char *out, const long *ns, const long *dtype);
   int c_pxmin (const unsigned char *in, const unsigned char *map,
               unsigned char *out, const long *ns, const long *dtype);
   int c_pxmnmx (const unsigned char *in, const long *ns, double *min,
                double *max, const long *dtype);
   int c_pxmul (const long *ns, const unsigned char *in1,
               const unsigned char *in2, unsigned char *out, const float *scale,
               const float *offset, const long *dtype);
   void c_pxsamp (const long *ns, long *dtype, long *inc, unsigned char *buf);
   int c_pxsub (const long *ns, const unsigned char *in1,
               const unsigned char *in2, unsigned char *out,
               const float *scale1, const float *scale2, const float *offset,
               const long *dtype);
   int c_pxswap (unsigned char *buf, const long *nentity, long *size);
   int c_pxsys (const long *insys, const long *outsys, unsigned char *buf,
               long *dtype, const long *size, const unsigned char *under,
               const unsigned char *over, const long *flag);

// statio
   int c_stclse (const long *fd);
   int c_stdldt (const long *fd, const char *data_type, const long *ic,
                const long *is);
   int c_stdlnd (const long *fd, const long *ic, const long *is);
   int c_stfdcl (const long *fd, const char *class_name, long *ic);
   int c_stfdsi (const long *fd, const char *class_name, const char *site_name,
                long *ic, long *is);
   int c_stget (const long *fd, const char *data_type, const long *ic,
               const long *is, long shape[9], unsigned char *buf);
   int c_stnode (const long *fd, const long *ic, const long *is, const long *id,
                char *data_type, long shape[9]);
   int c_stopen (long *fd, const char *hos, const long *access, long *nclass);
   long c_stput (const long *fd, const char *data_type, const long *ic,
                const long *is, const long shape[9], const unsigned char *buf);

// tapeio
   int c_tpbsf (const long *unit, long *count, const long *mode);
   int c_tpbsr (const long *unit, long *count, const long *mode);
   int c_tpclos (const long *unit, const long *mode, const long *rewopt,
                const long *eofopt);
   int c_tpfsf (const long *unit, long *count, const long *mode);
   int c_tpfsr (const long *unit, long *count, const long *mode);
   int c_tplblpr (char *printer);
   int c_tplsrd (const long *unit, char *key, long *clen, long *dlen,
                char *cbuf, unsigned char *dbuf, char *dtype);
   int c_tplswr (const long *unit, const char *key, const long *clen,
                const long *dlen, const char *cbuf, const unsigned char *dbuf,
                const char *dtype);
   int c_tpopen (long *unit, const long *skip, const long *mode,
                const long *opt1, const char *dens, const long *inout,
                char *message, char *tapeid);
   int c_tpread (const long *unit, const unsigned char *buffer,
                const long *bufsz, const long *mode);
   int c_tpreop (long *unit, const long *skip, const long *mode,
                const long *opt1, const char *dens, const long *inout,
                const char *drivename);
   int c_tprew (const long *unit, const long *mode, const long *opt1);
   int c_tpweof (const long *unit, const long *count, const long *mode);
   int c_tpwrit (const long *unit, const unsigned char *buffer,
                const long *bufsz, const long *mode);

// unixdisk
   int c_dkbuf (long *bufsiz);
   int c_dkcopy (const char *from, const char *to,
                const char directories[][CMLEN]);
   long c_dkcre (const char directories[][CMLEN], const char *filename,
                long filesize, int wait);
   long c_dkdel (const char *filename);
   int c_dkfree (const char *path);
   int c_dkmove (const char *from, const char *to,
                const char directories[][CMLEN]);
   int c_dkrspac (const char *path);

// utility
   void c_asc2eb (char *buf, const long *nbytes);
   int c_bitand (const int *in1, const int *in2);
   int c_bitor (const int *in1, const int *in2);
   int c_bitxor (const int *in1, const int *in2);
   void c_index (const char *string1, const char *string2, long *pos);
   void c_clean (const struct PARBLK *vblock, const char *name,
                const char hout[][CMLEN], const long *count);
   void c_complt (const char *name, const char out [][CMLEN],
                 const long *count);
   void c_delfil (const char *hname);
   void c_eb2asc (char *buf, const long *nbytes);
   int c_getcnt (const struct PARBLK *block, const char *param, long *count);
   int c_getin (const struct PARBLK *block, const char *param,
               const long *element, char tnames[][CMLEN], char hnames[][CMLEN],
               long *nmrimg, long (*window)[4], long (*bands)[MAXBND+1],
               long *nbands, long *totbnd, const long *imglimit,
               long *same_window, long (*qflag)[2]);
   int c_getin2 (const char inparm[], char tnames[][CMLEN],
                char hnames[][CMLEN], long *images, long (*window)[4],
                long (*bands)[MAXBND+1], long *nbands, long *totbnd,
                const long *imglimit, long *same_window, long (*qflag)[2]);
   int c_getnam (const struct PARBLK *block, const char *param,
                const long *element, char *hname, const char *ext,
                const long *access);
   int c_getnum (const struct PARBLK *block, const char *param,
                const long *count, const long *dtype, unsigned char *buffer);
   int c_getout (const struct PARBLK *block, const char *param, char *tname,
                char *hname);
   int c_getstr (const struct PARBLK *block, const char *param,
                const long *count, const long *len, char *buffer);
   int c_getsys (char *cursys);
   int c_gettyp (const char hname[][CMLEN], const long *nmrimg, long *dtype);
   int c_hstout (const struct PARBLK *vblock, const char *hosout,
                const char hosin[][CMLEN], const long *numin, const char *text);
   void c_index (const char *string1, const char *string2, long *pos);
   void c_inittm (const char *subject, const char *verb, const char *entity,
                 const char *category, const char *text, const long *max_entity,
                 const long *max_category, const long *msgtim,
                 const long *flag);
   void c_inlas (struct PARBLK *vblock, long *stdlun);
   int c_inmnmx (const long *dtype, const long *totbnd, double *minval,
                double *maxval);
   void c_lasmsg (const long fatalflg, const char *msgkey, const char *format,
                 ...);
   void c_low2up (char *buf, const long *nbytes);
   int c_mkimg (const char *inname, const char *append, char *tname,
               char *hname, const long *flag);
   int c_mknam (const char *filename, char *hname, const char *ext,
               const long *access);
   void c_pterm (const char *buf, const long *flag);
   int c_setprt (const struct PARBLK *vblock, long *term, long *file, long lp,
                char *prtfile, char *filename);
   char squeeze (const char *buf, const long length);
   void c_startm (void);
   void c_stoptm (void);
   int c_sysset (char *csys, long *nsys, const long *sysflg);
   void c_up2low (char *buf, const long *nbytes);



} // end extern "C"



