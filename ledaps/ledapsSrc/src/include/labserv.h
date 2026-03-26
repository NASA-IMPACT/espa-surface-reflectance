#ifndef LABSERV_H
#define LABSERV_H


int c_lsclos
(
    long *fd,           /* File descriptor of opened file (input)     */
    char *hostname,     /* host file name (input)                     */
    long *action        /* action flag:                               */
);

int c_lsdel
(
    long *fd,           /* file descriptor (input)                   */
    char *inkey         /* key of record to be deleted (input)       */
);

void c_lsmknm
(
    char *inname,       /* input host image name                      */
    char *ext,          /* input extension for label services file    */
    char *outname       /* output host label services file name       */
);

int c_lsopen
(
    long *fd,           /* File descriptor of open file (output)      */
    char *hname,        /* TAE host file name   (input)               */
    long *rwmode        /* File access:         (input)               */
);

int c_lsread
(
    long *fd,           /* File descriptor              (input)       */
    char *key,          /* input key of record to be read (in/out)    */
    long *clen,         /* number of char bytes to read (in/out)      */
    long *dlen,         /* number of data bytes to read (in/out)      */
    char *cbuf,         /* buffer to read char data into  (output)    */
    unsigned char *dbuf,/* data buffer to read data into              */
    char *dtype         /* datatype of data record (I2,B,R4,...) (out)*/
);

int c_lsrepl
(
    long *fd,           /* file descriptor                      (input)      */
    char *key,          /* key to be written                    (input)      */
    long *clen,         /* number of char bytes to be written   (input)      */
    long *dlen,         /* number of data bytes to be written   (input)      */
    char *cbuf,         /* char buffer to be written            (input)      */
    char *dbuf,         /* data buffer to be written            (input)      */
    char *dtype         /* data type of data record (I2,B,R4,...) (input)    */
);

int c_lsrew
(
    long *fd            /* file descriptor              (input)       */
);

int c_lsstat
(
    char *hostname,     /* host file name               (input)       */
    long *mode          /* type of access mode to check (input)       */
);

int c_lswrit
(
    long *fd,           /* file descriptor                      (input)       */
    char *key,          /* key of record to be written          (input)       */
    long *clen,         /* number of char bytes to be written   (input)       */
    long *dlen,         /* number of data bytes to be written   (input)       */
    char *cbuf,         /* character buffer to be written       (input)       */
    char *dbuf,         /* data buffer to be written            (input)       */
    char *dtype         /* datatype of data record (I2,B,R4,...)(input)       */
);

#endif
