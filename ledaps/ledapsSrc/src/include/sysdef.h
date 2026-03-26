/*
  Defines the valid systems
*/
#define NMR_SYSTEMS      13             /* Number of systems defined         */

#define SYS_GOULD_UTX    0              /* Gould UTX system                  */
#define SYS_SUN_BSD      1              /* Sun BSD system                    */
#define SYS_VAX_ULT      2              /* VAX Ultrix system                 */
#define SYS_VAX_VMS      3              /* VAX/VMS system                    */
#define SYS_IEEE_STD     4              /* IEEE standard (ref 754)           */
#define SYS_HP_MPE       5              /* Hewlert Packard system            */
#define SYS_PRIME_PRIMOS 6              /* PRIME/PRIMOS system               */
#define SYS_IBM_AIX	 7		/* IBM AIX			     */
#define SYS_IBM_MVS	 8		/* IBM MVS			     */
#define SYS_CRAY_UNICOS  9		/* Cray Y-MP Unicos                  */
#define SYS_SGI_IRX	 10		/* SGI IRX			     */
#define SYS_DG_UX 	 11		/* Data General   		     */
#define SYS_OTHER_MSC	 12		/* Misc. Other systems not covered   */

/*
  Strings for the same systems as above
*/
#define UTX     "gould-utx"
#define BSD     "sun-bsd"
#define ULT     "vax-ult"
#define VMS     "vax-vms"
#define IEEE    "ieee-std"
#define HP      "hp-mpe"
#define PRIMOS  "prime-primos"
#define AIX     "ibm-aix"
#define MVS     "ibm-mvs"
#define UNICOS  "cray-unicos"
#define IRX     "sgi-irx"
#define DG_UX 	"dg-ux"
#define MSC     "other-msc"
#define LINUX   "linux"

/*
  Underflow and Overflow flags
*/
#define SYS_ERR         0       /* Error on underflow and overflow           */
#define SYS_SET         1       /* Set to min/max for system                 */
#define SYS_USE         2       /* Use specified underflow & overflow        */

/*
  System code conversions
*/
#define SETN    0       /* Convert system string to system numerical code    */
#define SETS    1       /* Convert system numerical code to system string    */
