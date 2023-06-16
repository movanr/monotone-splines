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


double eval_quintic_hermite_poly

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
         
}  /* eval_quintic_hermite_poly */

/*--------------------------------------------------------------------------*/

void quintic_hermite_coeffs

   (double  h,          /* interval length */
    double  *f,         /* breakpoint values */
    double  *fx,        /* breakpoint approximated first derivatives */
    double  *fxx,       /* breakpoint approximated second derivatives */
    double  *coeff,     /* coefficient vector of size (n-1)*6 (output) */ 
    long    n)          /* number of breakpoints */

/*
    calculates coefficients of the piecewise quintic hermite for each interval
    and stores them in the coefficient matrix;
    the polynomials for each interval all have the domain [0,1]
    first derivative values d_i are implicitly mapped to h*d_i due to the 
    mapping of the intervals to [0,1] and the chain rule
    second derivative values d_i are implicitly mapped to h^2*d_i due to the 
    mapping of the intervals to [0,1] and the chain rule
    *equidistand version*
*/

{

long  i;    /* loop variable */

for (i=0; i<n-1; i++)
    {   
    coeff[6*i] = 
        -6.0    * f[i] 
        -3.0    * h*fx[i] 
        -0.5    * h*h*fxx[i] 
        +0.5    * h*h*fxx[i+1] 
        -3.0    * h*fx[i+1] 
        +6.0    * f[i+1];
    coeff[6*i+1] = 
        +15.0   * f[i] 
        +8.0    * h*fx[i] 
        +3.0/2  * h*h*fxx[i] 
        -1.0    * h*h*fxx[i+1] 
        +7.0    * h*fx[i+1] 
        -15     * f[i+1];
    coeff[6*i+2] = 
        -10.0   * f[i] 
        -6.0    * h*fx[i] 
        -3.0/2  * h*h*fxx[i] 
        +0.5    * h*h*fxx[i+1] 
        -4.0    * h*fx[i+1] 
        +10.0   * f[i+1];
    coeff[6*i+3] = 0.5 * h*h*fxx[i];
    coeff[6*i+4] = h*fx[i];
    coeff[6*i+5] = f[i];
    }

return;

}  /* quintic_hermite_coeffs */

/*--------------------------------------------------------------------------*/

void tridiagonal_matrix_algorithm
   
   (double  *s,   /* solution vector (output) */
    double  *a,   /* lower diagonal entries */
    double  *b,   /* diagonal entries */
    double  *c,   /* upper diagonal entries */
    double  *d,   /* right side of the tridiagonal system */
    long    n)    /* dimension of the system */

/* 
   Tridiagonal matrix algorithm to solve tridiagonal system 
*/

{
long  i; /* loop variable */

/* forward sweep */
c[0] = c[0]/b[0];
d[0] = d[0]/b[0];

for (i=1; i<n-1; i++) 
   {
   c[i] = c[i]/(b[i]-c[i-1]*a[i]);
   d[i] = (d[i]-d[i-1]*a[i])/(b[i]-c[i-1]*a[i]);
   }
d[n-1] = (d[n-1]-d[n-2]*a[n-1])/(b[n-1]-c[n-2]*a[n-1]);

/* back substitution */
s[n-1] = d[n-1];
for (i=n-2; i>=0; i--)
   s[i] = d[i] - c[i]*s[i+1]; 
   
return;

}  /* tridiagonal_matrix_algorithm */

/*--------------------------------------------------------------------------*/

void cubic_spline_derivs

   (double  h,          /* interval length */ 
    double  *f,         /* breakpoint values */      
    double  *fx,        /* approximated derivatives (output) */
    long    n)          /* number of breakpoints */


/*
   calculates derivatives of the cubic C²-spline
   * equidistant version *
*/

{
long   i;               /* loop variable */
double *d;              /* right side of the tridiagonal system */
double *a, *b, *c;      /* diagonal entries of the tridiagonal system */

/* allocate memory */
a = malloc (sizeof(double) * n);
b = malloc (sizeof(double) * n);
c = malloc (sizeof(double) * n-1);
d = malloc (sizeof(double) * n);

/* calculate entries of the tridiagonal system */
/* ("not-a-knot" boundary conditions) */

b[0] = 2;
c[0] = 4;
d[0] = (4*f[1] - 5*f[0] + f[2])/h;
for (i=1; i<n-1; i++)   
   {        
   a[i] = 1;
   b[i] = 4;
   c[i] = 1;
   d[i] = 3.0*(f[i+1] - f[i-1])/h;
   }
a[n-1] = -4;
b[n-1] = -2;
d[n-1] = (4*f[n-2] + f[n-3] - 5*f[n-1])/h;


/* solve tridiagonal system */
tridiagonal_matrix_algorithm (fx, a, b, c, d, n);

/* free variables */
free(a);
free(b);
free(c);
free(d);

return;

}  /* cubic_spline_derivs */

/*--------------------------------------------------------------------------*/

void cubic_spline_second_derivs

   (double  h,          /* interval length */ 
    double  *f,         /* breakpoint values */      
    double  *fx,        /* cubic spline derivatives (input) */
    double  *fxx,       /* cubic spline second derivatives (output) */
    long    n)          /* number of breakpoints */


/*
   calculates second derivatives of the cubic C²-spline. 
   First calculate the cubic spline derivatives and then call this function 
   to obtain the second derivatives 
   * equidistant version *
*/

{
long   i;               /* loop variable */

/* coefficients a_0[i] and a_1[i] are expressed in terms of f and fx */
/* second derivative is given by 2a_1[i] where a_1[i] is the coefficient
   of x^2 in the standard basis {x^3, x^2, x, 1};*/
for(i = 0; i < n-1; i++) {
    fxx[i] = 2 * (-3*f[i] + 3*f[i+1] - 2*h*fx[i] - h*fx[i+1]);
}
/* the last derivative is given by 6*a_0[n-2] + 2a_1[n-2]; */
fxx[n-1] = 6 * (2*f[n-2] + h*fx[n-2] - 2*f[n-1] + h*fx[n-1])
            + 2 * (-3*f[n-2] + 3*f[n-1] - 2*h*fx[n-2] - h*fx[n-1]);

return;

}  /* cubic_spline_second_derivs */

/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*       Implementation of MQSI by Lux et al. (2019)                        */
/*--------------------------------------------------------------------------*/
/*       ALGORITHM 1 is_monotone                                            */
/*--------------------------------------------------------------------------*/

int is_monotone
    (double x0,     /* first x-value */
     double x1,     /* second x-value */
     double f0,     /* value at x0 */
     double f1,     /* value at x1 */
     double fx0,    /* derivative at x0 */
     double fx1,    /* derivative at x1 */
     double fxx0,   /* second derivative at x0 */
     double fxx1)   /* second derivative at x1 */

/* 
    f is an order six polynomial defined by f0, f1, fx0, fx1, fxx0, fxx1;
    returns 1 (true) if f is monotone on [x0,x1], 0 (false) otherwise
*/
    
{
    /* variables of the algorithm */
    double w;
    double v;
    double A;
    double B;
    double gamma0;
    double gamma1;
    double alpha0;
    double alpha1;
    double beta0;
    double beta1;
    double gamma;
    double alpha;
    double beta;


    if(f0 == f1) { 
        if(!(fabs(fx0) < EPS)){
            //printf("false at condition 0\n");
            return 0; // false
        }
        if(!(fabs(fx1) < EPS)) {
            //printf("false at condition 1\n");
            return 0; // false
        }
        if(!(fabs(fxx0) < EPS)) {
            //printf("false at condition 2\n");
            return 0; // false
        }
        if(!(fabs(fxx1) < EPS)) {
            //printf("false at condition 3\n");
            return 0; // false
        }
    }

    if((fx0 < 0) || (fx1 < 0)) {
        //printf("false at condition 4\n");
        return 0; // false
    }

    if((fabs(fx0) < EPS) || (fabs(fx1) < EPS)) {
        w = x0 - x1;
        v = f0 - f1;
        if(fxx1 > -4*fx1/w) {
            //printf("false at condition 5\n");
            return 0; // false
        }
        if(fxx1 < (3*w*fxx0 - 24*fx0 - 32*fx1 + 60*v/w)/(5*w)) {
            //printf("false at condition 6\n");
            return 0; // false
        }
        if(fxx0 < 3*fx0/w) {
            //printf("false at condition 7\n");
            return 0; // false
        }
        //printf("true at condition 8\n");
        return 1; // true
    }

    // remaining case: fx0 != 0 and fx1 != 0
    
    A = fx0 * (x1-x0)/(f1-f0);
    B = fx1 * (x1-x0)/(f1-f0);

    gamma0 = 4 * fx0/fx1 * pow(B/A, 0.75);
    gamma1 = (x1-x0)/fx1 * pow(B/A, 0.75);
    alpha0 = 4 * pow(B/A, 0.25);
    alpha1 = (x1-x0)/fx1 * pow(B/A, 0.25);
    beta0 = 30 - (12*(fx0+fx1)*(x1-x0))/((f1-f0)*sqrt(A)*sqrt(B));
    beta1 = (-3*pow((x1-x0),2.0))/(2*(f1-f0)*sqrt(A)*sqrt(B));

    gamma = gamma0 + gamma1*fxx0;
    alpha = alpha0 + alpha1*fxx1;
    beta = beta0 + beta1*(fxx0-fxx1);

    if(beta <= 6.0) {
        //printf("condition 12a\n");
        return (alpha > -(beta+2)/2);
    }
    else {
        //printf("condition 12b\n");
        return (gamma > -2 * sqrt(beta-2));
    }
    //printf("condition end\n");
    return 0; // false

}  /* is_monotone */

/*--------------------------------------------------------------------------*/

int main() {

    long N = 100;
    long n = 3;
    long i, j;
    double h = 1.0;
    double H = h/N;
    FILE *file;

    double *f = malloc(sizeof(double) * n);
    double *fx = malloc(sizeof(double) * n);
    double *fxx = malloc(sizeof(double) * n);
    double *coeff = malloc(sizeof(double) * 6*(n-1));

    double *X = malloc(sizeof(double) * (N-1)*n);
    double *Y = malloc(sizeof(double) * (N-1)*n);

/*
    f[0] = 1;
    f[1] = 3;
    f[2] = 5;

    fx[0] = 0;
    fx[1] = 0;
    fx[2] = 0;

    fxx[0] = 1;
    fxx[1] = -1;
    fxx[2] = 1;
*/

/*
    f[0] = 10;
    f[1] = 10;
    f[2] = 10;
    f[3] = 10;
    f[4] = 10;
    f[5] = 10;
    f[6] = 10.5;
    f[7] = 15;
    f[8] = 50;
    f[9] = 60;
    f[10] = 85;
*/

    f[0] = 0;
    f[1] = 1;
    f[2] = 0;
    fx[0] = 0;
    fx[1] = 0;
    fx[2] = 0;
    fxx[0] = 2;
    fxx[1] = -10;
    fxx[2] = 2;

    //cubic_spline_derivs(h, f, fx, n);
    
    //cubic_spline_second_derivs(h, f, fx, fxx, n);

    quintic_hermite_coeffs(h, f, fx, fxx, coeff, n);
  
    // check monotonicity
    /*for(i = 0; i < n-1; i++) {
        printf("%d\n", is_monotone(i*h, (i+1)*h, f[i], f[i+1], fx[i], fx[i+1], fxx[i], fxx[i+1]));
    }*/


    for(i = 0; i < n-1; i++) {
        for(j = 0; j < N; j++) {
            X[i*N+j] = i*h + j*H;
            Y[i*N+j] = eval_quintic_hermite_poly(j*1.0/N, coeff + 6*i);
        }
    }




    file = fopen("interpolant.dat", "w");
    for(i = 0; i < N*(n-1); i++) {
        fprintf(file, "%lf %lf\n", X[i], Y[i]);
    }
    fclose(file);

    file = fopen("points.dat", "w");
    for(i = 0; i < n; i++) {
        fprintf(file, "%lf %lf\n", i*h, f[i]);
    }
    fclose(file);


    free(f);
    free(fx);
    free(fxx);
    free(coeff);
    free(X);
    free(Y);
    return 0;
}