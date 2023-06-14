#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>
#include <time.h>



/*--------------------------------------------------------------------------*/
/*                                                                          */
/*            MONOTONE QUINTIC SPLINE INTERPOLATION                         */
/*                                                                          */
/*            Bachelor studies by Moritz van Recum                          */
/*                                                                          */ 
/*            an implementation of                                          */
/*            'AN ALGORITHM FOR CONSTRUCTING MONOTONE                       */
/*            QUINTIC INTERPOLATING SPLINES'                                */
/*            by Thomas C. H. Lux et al. (2019)                             */
/*                                                                          */
/*                                                                          */
/*            Supervisor: Joachim Weickert                                  */  
/*                                                                          */
/*--------------------------------------------------------------------------*/


#define EPS 1e-14          /* small constant */


/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*        functions for piecewise quintic hermite interpolation             */
/*--------------------------------------------------------------------------*/


double eval_hermite_poly

   (double z,             /* evaluation point (in [0,1])*/
    double *coeff)        /* coefficients of size 6 */

/*
   evaluates a hermite polynomial at point z  
   (coeff[i][5] + coeff[i][4] z + coeff[i][3] z² + coeff[i][2] z³ + ...)
*/

{
   
/* evaluation at z using Horner's method */
return ((((coeff[0] 
    * z + coeff[1]) 
    * z + coeff[2]) 
    * z + coeff[3]) 
    * z + coeff[4]) 
    * z + coeff[5];
         
}  /* eval_hermite_poly */




int main() {

    double *coeff = malloc(sizeof(double) * 6);
    for(int i = 0; i < 6; i++) {
        coeff[i] = i;
    }    
    double result = eval_hermite_poly(1.345, coeff);
    printf("result = %lf\n", result);
    return 0;
}