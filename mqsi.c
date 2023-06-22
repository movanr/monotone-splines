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

#define EPS 1e-14 /* small constant */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/* https://www.geeksforgeeks.org/introduction-and-array-implementation-of-queue/ */
// C program for array implementation of queue
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
 
// A structure to represent a queue (with positive integers as values)
struct Queue {
    long front, rear, size;
    unsigned capacity;
    long* array;
};
 
// function to create a queue
// of given capacity.
// It initializes size of queue as 0
struct Queue* createQueue(unsigned capacity)
{
    struct Queue* queue = (struct Queue*)malloc(
        sizeof(struct Queue));
    queue->capacity = capacity;
    queue->front = queue->size = 0;
 
    // This is important, see the enqueue
    queue->rear = capacity - 1;
    queue->array = (long*)malloc(
        queue->capacity * sizeof(long));
    return queue;
}
 
// Queue is full when size becomes
// equal to the capacity
int isFull(struct Queue* queue)
{
    return (queue->size == queue->capacity);
}
 
// Queue is empty when size is 0
int isEmpty(struct Queue* queue)
{
    return (queue->size == 0);
}
 
// Function to add an item to the queue.
// It changes rear and size
void enqueue(struct Queue* queue, long item)
{
    if (isFull(queue))
        return;
    queue->rear = (queue->rear + 1)
                  % queue->capacity;
    queue->array[queue->rear] = item;
    queue->size = queue->size + 1;
    //printf("%d enqueued to queue\n", item);
}
 
// Function to remove an item from queue.
// It changes front and size
long dequeue(struct Queue* queue)
{
    if (isEmpty(queue))
        return -1;
    long item = queue->array[queue->front];
    queue->front = (queue->front + 1)
                   % queue->capacity;
    queue->size = queue->size - 1;
    return item;
}
 
// Function to get front of queue
long front(struct Queue* queue)
{
    if (isEmpty(queue))
        return -1;
    return queue->array[queue->front];
}
 
// Function to get rear of queue
long rear(struct Queue* queue)
{
    if (isEmpty(queue))
        return -1;
    return queue->array[queue->rear];
}

void destroy(struct Queue* queue) {
    free(queue->array);
    free(queue);
}


/*--------------------------------------------------------------------------*/
/*        functions for piecewise quintic hermite interpolation             */
/*--------------------------------------------------------------------------*/

double eval_quintic_hermite_poly

    (double z,      /* evaluation point (in [0,1])*/
     double *coeff) /* coefficients of size 6 */

/*
   evaluates a hermite polynomial at point z
   (coeff[i][5] + coeff[i][4] z + coeff[i][3] z² + coeff[i][2] z³ + ...)
*/

{

    /* evaluation at z using Horner's method */
    return ((((coeff[0] * z + coeff[1]) * z + coeff[2]) * z + coeff[3]) * z + coeff[4]) * z + coeff[5];

} /* eval_quintic_hermite_poly */

/*--------------------------------------------------------------------------*/

void quintic_hermite_coeffs

    (double h,      /* interval length */
     double *f,     /* breakpoint values */
     double *fx,    /* breakpoint approximated first derivatives */
     double *fxx,   /* breakpoint approximated second derivatives */
     double *coeff, /* coefficient vector of size (n-1)*6 (output) */
     long n)        /* number of breakpoints */

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

    long i; /* loop variable */

    for (i = 0; i < n - 1; i++)
    {
        coeff[6 * i] =
            -6.0 * f[i] - 3.0 * h * fx[i] - 0.5 * h * h * fxx[i] + 0.5 * h * h * fxx[i + 1] - 3.0 * h * fx[i + 1] + 6.0 * f[i + 1];
        coeff[6 * i + 1] =
            +15.0 * f[i] + 8.0 * h * fx[i] + 3.0 / 2 * h * h * fxx[i] - 1.0 * h * h * fxx[i + 1] + 7.0 * h * fx[i + 1] - 15 * f[i + 1];
        coeff[6 * i + 2] =
            -10.0 * f[i] - 6.0 * h * fx[i] - 3.0 / 2 * h * h * fxx[i] + 0.5 * h * h * fxx[i + 1] - 4.0 * h * fx[i + 1] + 10.0 * f[i + 1];
        coeff[6 * i + 3] = 0.5 * h * h * fxx[i];
        coeff[6 * i + 4] = h * fx[i];
        coeff[6 * i + 5] = f[i];
    }

    return;

} /* quintic_hermite_coeffs */

/*--------------------------------------------------------------------------*/

void tridiagonal_matrix_algorithm

    (double *s, /* solution vector (output) */
     double *a, /* lower diagonal entries */
     double *b, /* diagonal entries */
     double *c, /* upper diagonal entries */
     double *d, /* right side of the tridiagonal system */
     long n)    /* dimension of the system */

/*
   Tridiagonal matrix algorithm to solve tridiagonal system
*/

{
    long i; /* loop variable */

    /* forward sweep */
    c[0] = c[0] / b[0];
    d[0] = d[0] / b[0];

    for (i = 1; i < n - 1; i++)
    {
        c[i] = c[i] / (b[i] - c[i - 1] * a[i]);
        d[i] = (d[i] - d[i - 1] * a[i]) / (b[i] - c[i - 1] * a[i]);
    }
    d[n - 1] = (d[n - 1] - d[n - 2] * a[n - 1]) / (b[n - 1] - c[n - 2] * a[n - 1]);

    /* back substitution */
    s[n - 1] = d[n - 1];
    for (i = n - 2; i >= 0; i--)
        s[i] = d[i] - c[i] * s[i + 1];

    return;

} /* tridiagonal_matrix_algorithm */

/*--------------------------------------------------------------------------*/

void cubic_spline_derivs

    (double h,   /* interval length */
     double *f,  /* breakpoint values */
     double *fx, /* approximated derivatives (output) */
     long n)     /* number of breakpoints */

/*
   calculates derivatives of the cubic C²-spline
   * equidistant version *
*/

{
    long i;            /* loop variable */
    double *d;         /* right side of the tridiagonal system */
    double *a, *b, *c; /* diagonal entries of the tridiagonal system */

    /* allocate memory */
    a = malloc(sizeof(double) * n);
    b = malloc(sizeof(double) * n);
    c = malloc(sizeof(double) * n - 1);
    d = malloc(sizeof(double) * n);

    /* calculate entries of the tridiagonal system */
    /* ("not-a-knot" boundary conditions) */

    b[0] = 2;
    c[0] = 4;
    d[0] = (4 * f[1] - 5 * f[0] + f[2]) / h;
    for (i = 1; i < n - 1; i++)
    {
        a[i] = 1;
        b[i] = 4;
        c[i] = 1;
        d[i] = 3.0 * (f[i + 1] - f[i - 1]) / h;
    }
    a[n - 1] = -4;
    b[n - 1] = -2;
    d[n - 1] = (4 * f[n - 2] + f[n - 3] - 5 * f[n - 1]) / h;

    /* solve tridiagonal system */
    tridiagonal_matrix_algorithm(fx, a, b, c, d, n);

    /* free variables */
    free(a);
    free(b);
    free(c);
    free(d);

    return;

} /* cubic_spline_derivs */

/*--------------------------------------------------------------------------*/

void cubic_spline_second_derivs

    (double h,    /* interval length */
     double *f,   /* breakpoint values */
     double *fx,  /* cubic spline derivatives (input) */
     double *fxx, /* cubic spline second derivatives (output) */
     long n)      /* number of breakpoints */

/*
   calculates second derivatives of the cubic C²-spline.
   First calculate the cubic spline derivatives and then call this function
   to obtain the second derivatives
   * equidistant version *
*/

{
    long i; /* loop variable */

    /* coefficients a_0[i] and a_1[i] are expressed in terms of f and fx */
    /* second derivative is given by 2a_1[i] where a_1[i] is the coefficient
       of x^2 in the standard basis {x^3, x^2, x, 1};*/
    for (i = 0; i < n - 1; i++)
    {
        fxx[i] = 2 * (-3 * f[i] + 3 * f[i + 1] - 2 * h * fx[i] - h * fx[i + 1]);
    }
    /* the last derivative is given by 6*a_0[n-2] + 2a_1[n-2]; */
    fxx[n - 1] = 6 * (2 * f[n - 2] + h * fx[n - 2] - 2 * f[n - 1] + h * fx[n - 1]) + 2 * (-3 * f[n - 2] + 3 * f[n - 1] - 2 * h * fx[n - 2] - h * fx[n - 1]);

    return;

} /* cubic_spline_second_derivs */


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*       MONOTONE QUINTIC SPLINE INTERPOLATION                              */
/*                                                                          */
/*--------------------------------------------------------------------------*/

int is_monotone_old

    (double x0,   /* first x-value */
     double x1,   /* second x-value */
     double f0,   /* value at x0 */
     double f1,   /* value at x1 */
     double fx0,  /* derivative at x0 */
     double fx1,  /* derivative at x1 */
     double fxx0, /* second derivative at x0 */
     double fxx1) /* second derivative at x1 */

/*
    f is an order six polynomial defined by f0, f1, fx0, fx1, fxx0, fxx1;
    returns 1 (true) if f is monotone on [x0,x1], 0 (false) otherwise
*/

{
    /* variables */
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

    if (f0 == f1)
    {
        if (!(fabs(fx0) < EPS))
        {
            printf("is_monotone: false at condition 0\n");
            return 0; // false
        }
        if (!(fabs(fx1) < EPS))
        {
            printf("is_monotone: false at condition 1\n");
            return 0; // false
        }
        if (!(fabs(fxx0) < EPS))
        {
            printf("is_monotone: false at condition 2\n");
            return 0; // false
        }
        if (!(fabs(fxx1) < EPS))
        {
            printf("is_monotone: false at condition 3\n");
            return 0; // false
        }
    }

    if ((fx0 < 0) || (fx1 < 0))
    {
        printf("is_monotone: false at condition 4\n");
        return 0; // false
    }

    if ((fabs(fx0) < EPS) || (fabs(fx1) < EPS))
    {
        w = x0 - x1;
        v = f0 - f1;
        if (fxx1 > -4 * fx1 / w)
        {
            printf("is_monotone: false at condition 5\n");
            return 0; // false
        }
        if (fxx1 < (3 * w * fxx0 - 24 * fx0 - 32 * fx1 + 60 * v / w) / (5 * w))
        {
            printf("is_monotone: false at condition 6\n");
            return 0; // false
        }
        if (fxx0 < 3 * fx0 / w)
        {
            printf("is_monotone: false at condition 7\n");
            return 0; // false
        }
        printf("is_monotone: true at condition 8\n");
        return 1; // true
    }

    // remaining case: fx0 != 0 and fx1 != 0

    A = fx0 * (x1 - x0) / (f1 - f0);
    B = fx1 * (x1 - x0) / (f1 - f0);

    gamma0 = 4 * fx0 / fx1 * pow(B / A, 0.75);
    gamma1 = (x1 - x0) / fx1 * pow(B / A, 0.75);
    alpha0 = 4 * pow(B / A, 0.25);
    alpha1 = (x1 - x0) / fx1 * pow(B / A, 0.25);
    beta0 = 30 - (12 * (fx0 + fx1) * (x1 - x0)) / ((f1 - f0) * sqrt(A) * sqrt(B));
    beta1 = (-3 * pow((x1 - x0), 2.0)) / (2 * (f1 - f0) * sqrt(A) * sqrt(B));

    gamma = gamma0 + gamma1 * fxx0;
    alpha = alpha0 + alpha1 * fxx1;
    beta = beta0 + beta1 * (fxx0 - fxx1);
    printf("beta = %lf\n", beta);

    if (beta <= 6.0)
    {
        printf("is_monotone: %s at condition 12a\n", alpha > -(beta + 2) / 2 ? "true" : "false");
        return (alpha > -(beta + 2) / 2);
    }
    else
    {
        printf("is_monotone: %s at condition 12b\n", gamma > -2 * sqrt(beta - 2) ? "true" : "false");
        return (gamma > -2 * sqrt(beta - 2));
    }

} /* is_monotone */

/*--------------------------------------------------------------------------*/

int sgn 

   (double x) 
   
/*
   returns the sign of a double format value
*/ 

{
if (x > 0) return 1;
   
if (x < 0) return -1;
   
return 0;
   
}  /* sgn */

/*--------------------------------------------------------------------------*/

int is_monotone

    (double x0,   /* first x-value */
     double x1,   /* second x-value */
     double f0,   /* value at x0 */
     double f1,   /* value at x1 */
     double fx0,  /* derivative at x0 */
     double fx1,  /* derivative at x1 */
     double fxx0, /* second derivative at x0 */
     double fxx1) /* second derivative at x1 */

/*
    f is an order six polynomial defined by f0, f1, fx0, fx1, fxx0, fxx1;
    returns 1 (true) if f is monotone increasing on [x0,x1],
    0 (false) otherwise
*/

{

    /* variables */
    double w, z, alpha, beta, gamma, t, sign;

    if (fabs(f0 - f1) < EPS)
    {
        if (!(fabs(fx0) < EPS))
        {
            printf("is_monotone: false at condition 0\n");
            return 0; // false
        }
        if (!(fabs(fx1) < EPS))
        {
            printf("is_monotone: false at condition 1\n");
            return 0; // false
        }
        if (!(fabs(fxx0) < EPS))
        {
            printf("is_monotone: false at condition 2\n");
            return 0; // false
        }
        if (!(fabs(fxx1) < EPS))
        {
            printf("is_monotone: false at condition 3\n");
            return 0; // false
        }
        printf("true at condition 0\n");
        return 1; // true
    }

    if(f1 > f0) {
        sign = 1.0;
    }
    else {
        sign = -1.0;
    }

    if ((sign*fx0 < 0.0) || (sign*fx1 < 0.0))
    {
        printf("is_monotone: false at condition 4\n");
        return 0; // false
    }

    w = x1 - x0;
    z = f0 - f1;

    if ((fabs(fx0) < EPS) || (fabs(fx1) < EPS))
    {
        if (sign*fxx1 * w > sign*4.0 * fx1)
        {
            printf("is_monotone: false at condition 8\n");
            return 0; // false;
        }
        t = fx0 * (4*fx1 - fxx1 * w);
        if(t > 0.0) {
            t = 2.0 * sqrt(t);
        }
        if (t + sign*(3.0 * fx0 + fxx0 * w) < 0.0)
        {
            printf("is_monotone: false at condition 10\n");
            return 0; // false
        }
        if (-60.0 * z*sign - w * (sign*(24.0 * fx0 + 32.0 * fx1) - 2.0 * t + w*sign * (3.0 * fxx0 - 5.0 * fxx1)) < 0.0)
        {
            printf("is_monotone: false at condition 11\n");
            return 0; // false
        }
        printf("is_monotone: true at condition 1\n");
        return 1; // true
    }

    if (w * (2.0 * sqrt(fx0 * fx1) - sign*3.0 * (fx0 + fx1)) - sign*24.0 * z <= 0.0)
    {
        printf("is_monotone: false at condition 14\n");
        return 0; // false;
    }
    t = pow(fx0 * fx1, 0.75);
    alpha = sign * (4 * fx1 - fxx1 * w) * sqrt(sign*fx0) / t;
    gamma = sign * (4 * fx0 + fxx0 * w) * sqrt(sign*fx1) / t;
    beta = sign * (-60 * z / w + 3 * (w * (fxx1 - fxx0) - 8 * (fx0 + fx1))) / (2 * sqrt(fx0 * fx1));
    //printf("alpha = %lf\nbeta=%lf\ngamma=%lf\n", alpha, beta, gamma);

    if (beta <= 6)
    {
        printf("is_monotone: %s at condition 19\n", fmin(alpha, gamma) > -(beta + 2) / 2 ? "true" : "false");
        return (fmin(alpha, gamma) > -(beta + 2) / 2);
    }
    printf("is_monotone: %s at condition 20\n", fmin(alpha, gamma) > -2 * sqrt(beta - 2) ? "true" : "false");
    if(!(fmin(alpha, gamma) > -2 * sqrt(beta - 2))) {
        //printf("left = %lf, right = %lf\n", fmin(alpha, gamma), -2 * sqrt(beta - 2));
        printf("alpha = %lf, beta = %lf, gamma = %lf\n", alpha, beta, gamma);
    }
    return (fmin(alpha, gamma) > -2 * sqrt(beta - 2));

    return 0;
}

/*--------------------------------------------------------------------------*/

int is_monotone_index

    (long index,
     double h, 
     double *f, 
     double *fx, 
     double *fxx)

{

/*
    helper function to check if a piecewise order six polynomial given by 
    f, fx, fxx is monotone at interval given by index;
*/

return is_monotone(
    h*index, 
    h*(index+1), 
    f[index], 
    f[index+1], 
    fx[index], 
    fx[index+1], 
    fxx[index], 
    fxx[index+1]);
}  /* is_monotone_index */

/*--------------------------------------------------------------------------*/

double median

    (double a,
     double b,
     double c)

/*
    returs the median of a, b, c
*/

{
    /* auxiliary variables */
    double x;
    double y;
    double z;

    x = a - b;
    y = b - c;
    z = a - c;

    if (x * y > 0)
        return b;
    if (x * z > 0)
        return c;
    return a;

} /* median */

/*--------------------------------------------------------------------------*/

void make_monotone

    (double x0,    /* first x-value */
     double x1,    /* second x-value */
     double f0,    /* value at x0 */
     double f1,    /* value at x1 */
     double *fx0,  /* derivative at x0 (changed) */
     double *fx1,  /* derivative at x1 (changed) */
     double *fxx0, /* second derivative at x0 (changed) */
     double *fxx1) /* second derivative at x1 (changed) */

/*
    f is an order six polynomial defined by f0, f1, fx0, fx1, fxx0, fxx1;
    returns f monotone on [x0,x1]
*/

{
    /* variables */
    double v;
    double w;
    double A;
    double B;
    double *eta;
    double *eta0;

    /* allocate memory */
    eta = malloc(sizeof(double) * 2);
    eta0 = malloc(sizeof(double) * 2);

    if ((fabs(f1 - f0) < EPS))
    {
        (*fx0) = 0.0;
        (*fx1) = 0.0;
        (*fxx0) = 0.0;
        (*fxx1) = 0.0;
        printf("make_monotone: case 0\n");
        return;
    }
    (*fx0) = median(0, (*fx0), 14 * (f1 - f0) / (x1 - x0));
    (*fx1) = median(0, (*fx1), 14 * (f1 - f0) / (x1 - x0));

    if ((fabs(*fx0) < EPS) || (fabs(*fx1) < EPS))
    {
        w = x1 - x0;
        v = f0 - f1;
        if (5 * (*fx0) + 4 * (*fx1) > 20 * v / w)
        {
            (*fx0) = (*fx0) * fmax(0, (20 * v) / (w * (5 * (*fx0) + 4 * (*fx1))));
            (*fx1) = (*fx1) * fmax(0, (20 * v) / (w * (5 * (*fx0) + 4 * (*fx1))));
        }

        (*fxx0) = fmin((*fxx0), (4 * (2 * (*fx0) + (*fx1)) + 20 * v / w) / w);
        (*fxx0) = fmax((*fxx0), 3 * (*fx0) / w);
        (*fxx1) = fmin((*fxx1), -4 * (*fx1) / w);
        (*fxx1) = fmax((*fxx1), (3 * w * (*fxx0) - 24 * (*fx0) - 32 * (*fx1) + 60 * v / w) / (5 * w));

        printf("make_monotone: case 1\n");
        return;
    }

    /* remaining case: fx0 != 0 and fx1 != 0 */
    A = (*fx0) * (x1 - x0) / (f1 - f0);
    B = (*fx1) * (x1 - x0) / (f1 - f0);
    if (fmax(A, B) > 6.0)
    {
        (*fx0) = 6 * (*fx0) / (fmax(A, B));
        (*fx1) = 6 * (*fx1) / (fmax(A, B));
    }

    eta[0] = (*fxx0);
    eta[1] = (*fxx1);

    eta0[0] = -(sqrt(A) / 4) * (7 * sqrt(A) + 3 * sqrt(B));
    eta0[1] = (sqrt(B) / 4) * (3 * sqrt(A) + 7 * sqrt(B));

    *fxx0 = eta0[0];
    *fxx1 = eta0[1];

    /* free memory */
    free(eta);
    free(eta0);

    printf("make_monotone: case 2\n");
    return;
}

/*--------------------------------------------------------------------------*/

void mqsi_old

    (double h,    /* constant step size */
     double *f,   /* signal values */
     double *fx,  /* for testing TODO remove */
     double *fxx, /* for testing TODO remove */
     long n)      /* number of values */

/*
    Returns monotone quintic spline interpolant Q(x) such that 
    Q(i*h) = f_i and is monotone increasing (decreasing) on 
    all intervals that f_i is monotone increasing (decreasing)
    * equidistant version*
*/

{
long i;                 /* loop variable */
//double *fx;             /* first derivatives */
//double *fxx;            /* second derivatives */
double *fx_copy;
double *fxx_copy;
long shrink;            /* shrink factor */
struct Queue* queue;    /* queue for indices */

/* allocate memory */
//fx = malloc(sizeof(double) * n); // TODO decomment
//fxx = malloc(sizeof(double) * n); // TODO decomment 
fx_copy =  malloc(sizeof(double) * n);
fxx_copy = malloc(sizeof(double) * n);
/* Approximate first and second derivatives  at X_i = i*h */
//cubic_spline_derivs(h, f, fx, n); // TODO decomment
//cubic_spline_second_derivs(h, f, fx, fxx, n); // TODO decomment

/* initialize queue and shrink factor */
queue = createQueue(n);
shrink = 10;

/* initialize derivative copies */
for(i = 0; i < n; i++) {
    fx_copy[i] = fx[i];
    fxx_copy[i] = fxx[i];
}

for(i = 0; i < n-1; i++) {
    if(!is_monotone_index(i, h, f, fx, fxx)) {
        enqueue(queue, i);
    }
} 

long counter = 0;

while(!isEmpty(queue) && counter < 11) {
    counter++;
    i = dequeue(queue);

    if(is_monotone_index(i, h, f, fx, fxx)) {
        continue;
    }
    fx[i] = fx[i] - fx_copy[i]/shrink;
    fx[i+1] = fx[i+1] - fx_copy[i+1]/shrink;
    fxx[i] = fxx[i] - fxx_copy[i]/shrink;
    fxx[i+1] = fxx[i+1] - fxx_copy[i+1]/shrink;

    // sign change
    if(fx[i]*sgn(fx_copy[i]) < 0) {
        fx[i] = EPS/2;
    }
    if(fx[i+1]*sgn(fx_copy[i+1]) < 0) {
        fx[i+1] = EPS/2;
    }
    if(fxx[i]*sgn(fxx_copy[i]) < 0) {
        fxx[i] = EPS/2;
    }
    if(fxx[i+1]*sgn(fxx_copy[i+1]) < 0) {
        fxx[i+1] = EPS/2;
    }


    if(!is_monotone_index(i, h, f, fx, fxx)) {
        enqueue(queue, i);
    }
    if(!is_monotone_index((long)fmax(i-1, 0), h, f, fx, fxx)) {
        enqueue(queue, (long)fmax(i-1, 0));
    }
    if(!is_monotone_index((long)fmin(i+1, n-2), h, f, fx, fxx)) {
        enqueue(queue, (long)fmin(i+1, n-2));
    }
}




/* free memory */
//free(fx); // TODO decomment
//free(fxx); // TODO decomment
free(fx_copy);
free(fxx_copy);
destroy(queue);

return;
}  /* mqsi_old */

/*--------------------------------------------------------------------------*/

void mqsi

    (double h,
     double *f,
     double *fx,
     double *fxx,
     long   n)

/*
    Computes monotone quintic spline interpolant Q(x) such that 
    Q(i*h) = f_i and is monotone increasing (decreasing) on 
    all intervals that f_i is monotone increasing (decreasing);
    The piecewise quintic spline at interval i is given by the 
    interval values and first and second derivatives
    f[i], f[i+1], fx[i], fx[i+1], fxx[i], fxx[i+1], which are
    being modified by the algorihm to achieve monotonicity
    * equidistant version*
*/

{
int *checking;
int *growing;
int *shrinking;
int *to_check;
int *to_grow;
int *to_shrink;
int *flats;
int *extrema;

double A, B, direction, dx, scale, step_size;
long i, j, nc, ng, ns;
int searching;

/* allocate memory */
checking = malloc(sizeof(int)*n);
growing = malloc(sizeof(int)*n);
shrinking = malloc(sizeof(int)*n);
to_check = malloc(sizeof(int)*n);
to_grow = malloc(sizeof(int)*n);
to_shrink = malloc(sizeof(int)*n);
flats = malloc(sizeof(int)*n);
extrema = malloc(sizeof(int)*n);

/* Identify local extreme points and flat points */
flats[0] = fabs(f[0] - f[1]) < EPS * (1.0 + fabs(f[0] + fabs(f[1])));
flats[n-1] = fabs(f[n-2] - f[n-1]) < EPS * (1.0 + fabs(f[n-2]) + fabs(f[n-1]));
extrema[0] = 0;
extrema[1] = 0;

for(i = 1; i < n-1; i++) {
    if(fabs(f[i-1]-f[i]) < EPS*(1.0+fabs(f[i])+fabs(f[i+1]))
    || fabs(f[i]-f[i+1]) < EPS*(1.0+fabs(f[i])+fabs(f[i+1]))) {
        flats[i] = 1;
        extrema[i] = 0;
    } 
    else {
        flats[i] = 0;
        extrema[i] = (f[i]-f[i-1]) * (f[i+1]-f[i]) < 0.0;
    }
}

/*
    Use local quadratic interpolants to estimate slopes and second
    derivatives at all points. Use zero-sloped quadratic interpolants
    at extrema with left/right neighbors to estimate curvature.
*/
for(i = 0; i < n; i++) {
    /* Initialize the curvature to be maximally large */
    fxx[i] = 1e54;
    /* If this is locally flat, then first and second derivatives are zero 
       valued */
    if(flats[i]) {
        fx[i] = 0.0;
        fxx[i] = 0.0;
    }
    /* If this is an extreme point (local minimum or maximum f),
       construct quadratic interpolants that have zero slope here and
       interpolate left/right neighbors */
    else if(extrema[i]) {
        /* set the first derivative to zero */
        fx[i] = 0.0;
        /* Compute the coefficient A in  Ax^2+Bx+C  that interpolates 
           at h(i-1) */
        fxx[i] = (f[i-1] - f[i]) / (h*h);
        /* Compute the coefficient A in  Ax^2+Bx+C  that 
           interpolates at h(i+1) */
        A = (f[i+1] - f[i]) / (h*h);
        if(fabs(A) < fabs(fxx[i])) {
            fxx[i] = A;
        }
        /* Compute the actual second derivative (instead of coefficient A) */
        fxx[i] *= 2.0;
        printf("fx[%ld] = %lf, fxx[%ld] = %lf\n", i, fx[i], i, fxx[i]);
    }
    else {
        /* Determine the direction of change at the point I */
        if(i==0) {
            if(f[i] < f[i+1]) {
                direction = 1.0;
            }
            else if(f[i] > f[i+1]) {
                direction = -1.0;
            }
        }
        else {
            if(f[i-1] < f[i]) {
                direction = 1.0;
            }
            else if (f[i-1] > f[i]) {
                direction = -1.0;
            }
        }
    
        /* ------------------- */
        /* Quadratic left of i */
        if(i-1 > 0) {
            /* If a zero derivative at left point, use its right interpolant */
            if(flats[i-1] || extrema[i-1]) {
                A = (f[i]-f[i-1])/(h*h);
                B = 0.0;
            }
            /* Otherwise use the standard quadratic on the left */
            else {
                A = (f[i-2] + f[i] - 2*f[i-1])/(2*h*h);
                B = (f[i] - f[i-2])/(2*h);
            }
        }
        dx = 2.0 * A * h + B;
        if(dx*direction >= 0.0) {
            fx[i] = dx;
            fxx[i] = A;
        }
        /* --------------------------- */
        /* Quadratic centered at i (require that it has at least one
        neighbor that is not forced to zero slope) */
        if((i > 0) && (i < n-1)) {
            if(!((extrema[i-1] || flats[i-1]) && (extrema[i+1] || flats[i+1]))) {
                /* Construct quadratic interpolant through this point and 
                neighbors */
                A = (f[i-1] + f[i+1] - 2*f[i])/(2*h*h);
                B = (f[i+1] - f[i-1])/(2*h);
                dx = B;
                if((dx*direction >= 0.0) && (fabs(A) < fabs(fxx[i]))) {
                    fx[i] = dx;
                    fxx[i] = A;
                }
            }
        }
        /* -------------------- */
        /* Quadratic right of i */
        if(i+1 < n-1) {
            /* If a zero derivative at right point, use its left interpolant */
            if(extrema[i+1] || flats[i+1]) {
                A = (f[i]-f[i+1])/(h*h);
                B = 0.0;
            }
            /* Otherwise use the standard quadratic on the right */
            else {
                A = (f[i] + f[i+2] - 2*f[i+1])/(2*h*h);
                B = (f[i+2] - f[i])/(2*h);
            }
            dx = 2.0 * A * (-h) + B;
            /* Keep this new quadratic if it has less curvature */
            if((dx*direction >= 0.0) && (fabs(A) < fabs(fxx[i]))) {
                fx[i] = dx;
                fxx[i] = A;
            }
        }
        /* Set the final quadratic */
        if(fxx[i] == 1e54) {
            fx[i] = 0.0;
            fxx[i] = 0.0;
        }
        /* Compute curvature of quadratic from coefficient of x^2 */
        else {
            fxx[i] = 2.0*fxx[i];
        }
        printf("fx[%ld] = %lf, fxx[%ld] = %lf\n", i, fx[i], i, fxx[i]);
    }
} 

/* free memory */
free(checking);
free(growing);
free(shrinking);
free(to_check);
free(to_grow);
free(to_shrink);
free(flats);
free(extrema);

return;

}  /* mqsi */

/*--------------------------------------------------------------------------*/

int main()
{
    long N = 100;
    long n = 11;
    long i, j;
    double h = 1.0;
    double H = h / N;
    FILE *file;

    double *f = malloc(sizeof(double) * n);
    double *fx = malloc(sizeof(double) * n);
    double *fxx = malloc(sizeof(double) * n);
    double *coeff = malloc(sizeof(double) * 6 * (n - 1));

    double *X = malloc(sizeof(double) * (N - 1) * n);
    double *Y = malloc(sizeof(double) * (N - 1) * n);

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
        f[0] = 10;
        f[1] = 12;
        f[2] = 10;
        f[3] = 15;
        f[4] = 10;
        f[5] = 10;
        f[6] = 10.5;
        f[7] = 15;
        f[8] = 50;
        f[9] = 60;
        f[10] = 85;


/*
    f[0] = 0;
    f[1] = 1;
    fx[0] = -1;
    fx[1] = 1;
    fxx[0] = -1;
    fxx[1] = -5;
*/
    //cubic_spline_derivs(h, f, fx, n);
    //cubic_spline_second_derivs(h, f, fx, fxx, n);

    /* test call mqsi */
    //printf("%lf\n%lf\n%lf\n%lf\n", fx[0], fx[1], fxx[0], fxx[1]);    
    mqsi(h, f, fx, fxx, n);
    //printf("%lf\n%lf\n%lf\n%lf\n", fx[0], fx[1], fxx[0], fxx[1]);
    //is_monotone_index(0, h, f, fx, fxx);

    //printf("%lf, %lf, %lf, %lf, %lf, %lf\n", f[5], f[6], fx[5], fx[6], fxx[5], fxx[6]);


    for (i = 0; i < n - 1; i++)
    {
        if (!is_monotone(i * h, (i + 1) * h, f[i], f[i + 1], fx[i], fx[i + 1], fxx[i], fxx[i + 1]))
        {
            // make_monotone(i*h, (i+1)*h, f[i], f[i+1], fx+i, fx+i+1, fxx+i, fxx+i+1);
        }
    }
    quintic_hermite_coeffs(h, f, fx, fxx, coeff, n);
    //printf("%lf %lf %lf %lf %lf %lf\n", coeff[0], coeff[1], coeff[2], coeff[3], coeff[4], coeff[5]);
    //printf("%lf*x^5+%lf*x^4+%lf*x^3+%lf*x^2+%lf*x+%lf\n", coeff[6*5+0], coeff[6*5+1], coeff[6*5+2], coeff[6*5+3], coeff[6*5+4], coeff[6*5+5]);


    for (i = 0; i < n - 1; i++)
    {
        for (j = 0; j < N; j++)
        {
            X[i * N + j] = i * h + j * H;
            Y[i * N + j] = eval_quintic_hermite_poly(j * 1.0 / N, coeff + 6 * i);
        }
    }


    file = fopen("interpolant.dat", "w");
    for (i = 0; i < N * (n - 1); i++)
    {
        fprintf(file, "%lf %lf\n", X[i], Y[i]);
    }
    fclose(file);

    file = fopen("points.dat", "w");
    for (i = 0; i < n; i++)
    {
        fprintf(file, "%lf %lf\n", i * h, f[i]);
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