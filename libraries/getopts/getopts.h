/* getopts.h */

#ifndef VEGAFEM_GETOPTS_H
#define VEGAFEM_GETOPTS_H

#define OPTINT 1
#define OPTSTR 2
#define OPTBOOL 3
#define OPTLONG 4
#define OPTSUBOPT 5

namespace vegafem
{
typedef struct {const char *sw;
               int opttyp;
               void *var;} opt_t;
int getopts(int argc, char **argv, opt_t opttable[]);
}//namespace vegafem
#endif

