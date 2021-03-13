#include <stdio.h>
#include <stdlib.h>

extern "C"
{
   double randomval_(void) {
     double ret;
     ret = ((double)rand())/((double)(RAND_MAX));
     return(ret);
   }
} 
