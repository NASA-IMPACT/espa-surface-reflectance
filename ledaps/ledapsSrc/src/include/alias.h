#ifndef ALIAS_H
#define ALIAS_H

#include <stdio.h>
#include "desc.h"

void bldalias
(
    char *aliasfile          /* name of user's alias file    */
);

int getalias
(
    struct DESC *instring,
    struct DESC *outstring
);

int c_getalias
(
    char *instring,          /* input string                 */
    char *outstring          /* output string                */
);

int delalias
(
    char *instring           /* input string                 */
);

int findalias 
(
    FILE *aliasf,            /* alias file  (input)          */
    char *aliasname,         /* alias name  (input)          */
    char *aliastext          /* alias text  (output)         */
);

int listalias 
(
    FILE *pfile,             /* print file                    */
    char *instring           /* input string                  */
);

int putalias
(
    char *aliasname,         /* name of alias add  (input)   */
    char *aliastext,         /* text of alias add  (input)   */
    char *aliasfile,         /* name of user's alias file (o)*/
    long *newfile            /* new file created? y=1,n=0 (o)*/
);

#endif
