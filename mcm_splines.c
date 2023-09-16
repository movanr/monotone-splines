#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>
#include <time.h>
#include <unistd.h>


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*        MEAN CURVATURE MOTION, EXPLICIT SCHEMES WITH (MONOTONE)           */
/*                 CUBIC AND QUINTIC SPLINE INTERPOLATION                   */
/*                                                                          */
/*            Bachelor studies by Moritz van Recum                          */
/*                                                                          */
/*                Supervisor: Joachim Weickert                              */
/*                                                                          */
/*            (Copyright by Joachim Weickert, 11/2021)                      */
/*--------------------------------------------------------------------------*/

#define EPS 1e-14         
#define FALSE 0
#define TRUE 1
#define PI 3.14159
#define SAMPLE_RATE 100

/*--------------------------------------------------------------------------*/

void alloc_double_vector

     (double **vector,   /* vector */
      long   n1)         /* size */

/*
  allocates memory for a double format vector of size n1
*/

{
*vector = (double *) malloc (n1 * sizeof(double));

if (*vector == NULL)
   {
   printf("alloc_double_vector: not enough memory available\n");
   exit(1);
   }

return;

}  /* alloc_double_vector */

/*--------------------------------------------------------------------------*/

void alloc_double_matrix

     (double ***matrix,  /* matrix */
      long   n1,         /* size in direction 1 */
      long   n2)         /* size in direction 2 */

/*
  allocates memory for a double format matrix of size n1 * n2 
*/

{
long i;    /* loop variable */

*matrix = (double **) malloc (n1 * sizeof(double *));

if (*matrix == NULL)
   {
   printf("alloc_double_matrix: not enough memory available\n");
   exit(1);
   }

for (i=0; i<n1; i++)
    {
    (*matrix)[i] = (double *) malloc (n2 * sizeof(double));
    if ((*matrix)[i] == NULL)
       {
       printf("alloc_double_matrix: not enough memory available\n");
       exit(1);
       }
    }

return;

}  /* alloc_double_matrix */

/*--------------------------------------------------------------------------*/

void free_double_vector

     (double  *vector,    /* vector */
      long    n1)         /* size */

/*
  frees memory for a double format vector of size n1
*/

{

free(vector);
return;

}  /* free_double_vector */

/*--------------------------------------------------------------------------*/

void free_double_matrix

     (double  **matrix,   /* matrix */
      long    n1,         /* size in direction 1 */
      long    n2)         /* size in direction 2 */

/*
  frees memory for a double format matrix of size n1 * n2
*/

{
long i;   /* loop variable */

for (i=0; i<n1; i++)
   free(matrix[i]);

free(matrix);

return;

}  /* free_double_matrix */

/*--------------------------------------------------------------------------*/

void read_string

     (char *v)         /* string to be read */

/*
  reads a string v
*/

{
if (fgets (v, 80, stdin) == NULL)
   {
   printf ("read_string: cannot read value, aborting\n");
   exit(1);
   }

if (v[strlen(v)-1] == '\n')
   v[strlen(v)-1] = 0;

return;

}  /* read_string */

/*--------------------------------------------------------------------------*/

void read_long

     (long *v)         /* value to be read */

/*
  reads a long value v
*/

{
char   row[80];    /* string for reading data */

if (fgets (row, 80, stdin) == NULL)
   {
   printf("read_long: cannot read value\n");
   exit(1);
   }
if (row[strlen(row)-1] == '\n')
   row[strlen(row)-1] = 0;
sscanf(row, "%ld", &*v);

return;

}  /* read_long */


/*--------------------------------------------------------------------------*/

void read_double

     (double *v)         /* value to be read */

/*
  reads a double value v
*/

{
char   row[80];    /* string for reading data */

if (fgets (row, 80, stdin) == NULL)
   {
   printf ("read_double: cannot read value\n");
   exit(1);
   }
if (row[strlen(row)-1] == '\n')
   row[strlen(row)-1] = 0;
sscanf(row, "%lf", &*v);

return;

}  /* read_double */

/*--------------------------------------------------------------------------*/

void skip_white_space_and_comments 

     (FILE *inimage)  /* input file */

/*
  skips over white space and comments while reading the file
*/

{

int   ch = 0;   /* holds a character */
char  row[80];  /* for reading data */

/* skip spaces */
while (((ch = fgetc(inimage)) != EOF) && isspace(ch));
  
/* skip comments */
if (ch == '#')
   {
   if (fgets(row, sizeof(row), inimage))
      skip_white_space_and_comments (inimage);
   else
      {
      printf("skip_white_space_and_comments: cannot read file\n");
      exit(1);
      }
   }
else
   fseek (inimage, -1, SEEK_CUR);

return;

} /* skip_white_space_and_comments */

/*--------------------------------------------------------------------------*/

void read_pgm_to_double

     (const char  *file_name,    /* name of pgm file */
      long        *nx,           /* image size in x direction, output */
      long        *ny,           /* image size in y direction, output */
      double      ***u)          /* image, output */

/*
  reads a greyscale image that has been encoded in pgm format P5 to
  an image u in double format;
  allocates memory for the image u;
  adds boundary layers of size 1 such that
  - the relevant image pixels in x direction use the indices 1,...,nx
  - the relevant image pixels in y direction use the indices 1,...,ny
*/

{
char  row[80];      /* for reading data */
long  i, j;         /* image indices */
long  max_value;    /* maximum color value */
FILE  *inimage;     /* input file */

/* open file */
inimage = fopen (file_name, "rb");
if (inimage == NULL)
   {
   printf ("read_pgm_to_double: cannot open file '%s'\n", file_name);
   exit(1);
   }

/* read header */
if (fgets(row, 80, inimage) == NULL)
   {
   printf ("read_pgm_to_double: cannot read file\n");
   exit(1);
   }

/* image type: P5 */
if ((row[0] == 'P') && (row[1] == '5'))
   {
   /* P5: grey scale image */
   }
else
   {
   printf ("read_pgm_to_double: unknown image format\n");
   exit(1);
   }

/* read image size in x direction */
skip_white_space_and_comments (inimage);
if (!fscanf (inimage, "%ld", nx))
   {
   printf ("read_pgm_to_double: cannot read image size nx\n");
   exit(1);
   }

/* read image size in x direction */
skip_white_space_and_comments (inimage);
if (!fscanf (inimage, "%ld", ny))
   {
   printf ("read_pgm_to_double: cannot read image size ny\n");
   exit(1);
   }

/* read maximum grey value */
skip_white_space_and_comments (inimage);
if (!fscanf (inimage, "%ld", &max_value))
   {
   printf ("read_pgm_to_double: cannot read maximal value\n");
   exit(1);
   }
fgetc(inimage);

/* allocate memory */
alloc_double_matrix (u, (*nx)+2, (*ny)+2);

/* read image data row by row */
for (j=1; j<=(*ny); j++)
 for (i=1; i<=(*nx); i++)
     (*u)[i][j] = (double) getc(inimage);

/* close file */
fclose (inimage);

return;

}  /* read_pgm_to_double */

/*--------------------------------------------------------------------------*/

void comment_line

     (char* comment,       /* comment string (output) */
      char* lineformat,    /* format string for comment line */
      ...)                 /* optional arguments */

/* 
  Adds a line to the comment string comment. The string line can contain 
  plain text and format characters that are compatible with sprintf.
  Example call: 
  print_comment_line(comment, "Text %lf %ld", double_var, long_var).
  If no line break is supplied at the end of the input string, it is 
  added automatically.
*/

{
char     line[80];
va_list  arguments;

/* get list of optional function arguments */
va_start (arguments, lineformat);

/* convert format string and arguments to plain text line string */
vsprintf (line, lineformat, arguments);

/* add line to total commentary string */
strncat (comment, line, 80);

/* add line break if input string does not end with one */
if (line[strlen(line)-1] != '\n')
   sprintf (comment, "%s\n", comment);

/* close argument list */
va_end (arguments);

return;

}  /* comment_line */

/*--------------------------------------------------------------------------*/

void write_double_to_pgm

     (double  **u,          /* image, unchanged */
      long    nx,           /* image size in x direction */
      long    ny,           /* image size in y direction */
      char    *file_name,   /* name of pgm file */
      char    *comments)    /* comment string (set 0 for no comments) */

/*
  writes a greyscale image in double format into a pgm P5 file
*/

{
FILE           *outimage;  /* output file */
long           i, j;       /* loop variables */
double         aux;        /* auxiliary variable */
unsigned char  byte;       /* for data conversion */

/* open file */
outimage = fopen (file_name, "wb");
if (NULL == outimage)
   {
   printf("could not open file '%s' for writing, aborting\n", file_name);
   exit(1);
   }

/* write header */
fprintf (outimage, "P5\n");                  /* format */
if (comments != 0)
   fprintf (outimage, "%s", comments);       /* comments */
fprintf (outimage, "%ld %ld\n", nx, ny);     /* image size */
fprintf (outimage, "255\n");                 /* maximal value */

/* write image data */
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
     {
     aux = u[i][j] + 0.499999;    /* for correct rounding */
     if (aux < 0.0)
        byte = (unsigned char)(0.0);
     else if (aux > 255.0)
        byte = (unsigned char)(255.0);
     else
        byte = (unsigned char)(aux);
     fwrite (&byte, sizeof(unsigned char), 1, outimage);
     }

/* close file */
fclose (outimage);

return;

}  /* write_double_to_pgm */

/*--------------------------------------------------------------------------*/

void dummies_double

     (double **u,        /* image */
      long   nx,         /* size in x direction */
      long   ny)         /* size in y direction */

/* 
  creates dummy boundaries for a double format image u by mirroring 
*/

{
long i, j;  /* loop variables */

for (i=1; i<=nx; i++)
    {
    u[i][0]    = u[i][1];
    u[i][ny+1] = u[i][ny];
    }

for (j=0; j<=ny+1; j++)
    {
    u[0][j]    = u[1][j];
    u[nx+1][j] = u[nx][j];
    }

return;

}  /* dummies_double */

/*--------------------------------------------------------------------------*/

void analyse_grey_double

     (double  **u,         /* image, unchanged */
      long    nx,          /* pixel number in x direction */
      long    ny,          /* pixel number in y direction */
      double  *min,        /* minimum, output */
      double  *max,        /* maximum, output */
      double  *mean,       /* mean, output */
      double  *std)        /* standard deviation, output */

/*
  computes minimum, maximum, mean, and standard deviation of a greyscale
  image u in double format
*/

{
long    i, j;       /* loop variables */
double  help1;      /* auxiliary variable */
double  help2;      /* auxiliary variable */

/* compute maximum, minimum, and mean */
*min  = u[1][1];
*max  = u[1][1];
help1 = 0.0;
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     if (u[i][j] < *min) *min = u[i][j];
     if (u[i][j] > *max) *max = u[i][j];
     help1 = help1 + u[i][j];
     }
*mean = help1 / (nx * ny);

/* compute standard deviation */
*std = 0.0;
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     help2  = u[i][j] - *mean;
     *std = *std + help2 * help2;
     }
*std = sqrt (*std / (nx * ny));

return;

}  /* analyse_grey_double */



/*--------------------------------------------------------------------------*/
/*        functions for piecewise cubic hermite interpolation               */
/*--------------------------------------------------------------------------*/

double eval_cubic_hermite_poly

   (double z,             /* evaluation point (in [0,1])*/
    double *coeff)        /* coefficients of size 4 */

/*
   evaluates a hermite polynomial at point z  
   (coeff[i][3] + coeff[i][2] z + coeff[i][1] z^2 + coeff[i][0] z^3)
*/

{
   
/* evaluation at z using Horner's method */
return ((coeff[0] * z + coeff[1]) * z + coeff[2]) * z + coeff[3];
         
}  /* eval_cubic_hermite_poly */

/*--------------------------------------------------------------------------*/

double eval_quintic_hermite_poly

    (double z,      /* evaluation point (in [0,1])*/
     double *coeff) /* coefficients of size 6 */

/*
   evaluates a hermite polynomial at point z
   (coeff[i][5] + coeff[i][4] z + coeff[i][3] z^2 + coeff[i][2] z^3 + ...)
*/

{

    /* evaluation at z using Horner's method */
    return ((((coeff[0] * z + coeff[1]) * z + coeff[2]) * z 
            + coeff[3]) * z + coeff[4]) * z + coeff[5];

} /* eval_quintic_hermite_poly */

/*--------------------------------------------------------------------------*/

double eval_piecewise_poly

   (double z,              /* evaluation point */
    double *coeff,         /* coefficients of size (n-1) * 4 */
    long   p,              /* spline order */ 
    long   n)              /* number of breakpoints */
    
/*
   evaluates a piecewise polynomial at point z in the i-th interval 
   (is given by c[i+0] + c[i+1] (z-i) + c[i+2] (z-i)^2 + c[i+3] (z-i)^3...
   * equidistant version (h=1) *
*/

{
long i; /* indices */             
/* first compute interval index i for z */
/* left boundary */
if (z < 1)
   i = 0;
/* right boundary */
else if (z >= n-2)
   i = n-2;
/* z lies between 1 and n-2 */
i = (long)z;
   
/* evaluation at z-i in the i-th interval */
if(p == 4) {
   return eval_cubic_hermite_poly(z-i, coeff+(i*p));
}
else if(p == 6) {
   return eval_quintic_hermite_poly(z-i, coeff+(i*p));
}

printf("eval_piecewise_poly: invalid spline order: %ld\n", p);
printf("aborting\n");
exit(1);

}  /* eval_piecewise_poly */

/*-------------------------------------------------------------------------*/

void hermite_coeffs

   (double  *Y,         /* breakpoint values */
    double  *d,         /* breakpoint approximated derivatives */
    double  *coeff,     /* coefficient vector of size (n-1)*4 */
    long    n)          /* number of breakpoints */
    
/*
   calculates coefficients of the piecewise cubic hermite for each interval
   and stores them in the coefficient vector;
   uses the standard basis {x^3, x^2, x, 1};
   the polynomials for each interval all have the domain [0,1]
   * equidistant version (h=1) *
*/  

{
long  i;    /* loop variable */

for (i=0; i<n-1; i++)
   {
   /* coefficient for x^3 */
   coeff[4*i] = 2*Y[i] + d[i] - 2*Y[i+1] + d[i+1];   
   /* coefficient for x^2 */
   coeff[4*i+1] = -3*Y[i] + 3*Y[i+1] - 2*d[i] - d[i+1];
   /* coefficient for x^1 */
   coeff[4*i+2] = d[i];
   /* coefficient for x^0 */
   coeff[4*i+3] = Y[i];
   }

return;

}  /* hermite_coeffs */

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

for (i=1; i<n-1; i++) {
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

void cubic_spline_derivs_eq

   (double  *Y,         /* breakpoint values */      
    double  *s,         /* slopes (output) */
    long    n)          /* number of breakpoints */


/*
   calculates derivatives of the cubic C^2-spline
   * equidistant version (h=1) *
*/

{
long   i;               /* loop variable */
double *d;              /* right side of the tridiagonal system */
double *a, *b, *c;      /* diagonal entries of the tridiagonal system */

/* allocate memory */
alloc_double_vector(&a, n);
alloc_double_vector(&b, n);
alloc_double_vector(&c, n-1);
alloc_double_vector(&d, n);

/* calculate entries of the tridiagonal system */
/* ("not-a-knot" boundary conditions) */

b[0] = 2;
c[0] = 4;
d[0] = (4*Y[1] - 5*Y[0] + Y[2]);
for (i=1; i<n-1; i++)   
   {        
   a[i] = 1;
   b[i] = 4;
   c[i] = 1;
   d[i] = 3.0*(Y[i+1] - Y[i-1]);
   }
a[n-1] = -4;
b[n-1] = -2;
d[n-1] = (4*Y[n-2] + Y[n-3] - 5*Y[n-1]);


/* solve tridiagonal system */
tridiagonal_matrix_algorithm (s, a, b, c, d, n);

/* free variables */
free_double_vector(a, n);
free_double_vector(b, n);
free_double_vector(c, n-1);
free_double_vector(d, n);


return;

}  /* cubic_spline_derivs_eq */


/*--------------------------------------------------------------------------*/

void cubic_spline_derivs_zero_bc

    (double *f,  /* breakpoint values */
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
    /* ("zero" boundary conditions) */

    b[0] = 1;
    c[0] = 0;
    d[0] = 0;
    for (i = 1; i < n - 1; i++)
    {
        a[i] = 1;
        b[i] = 4;
        c[i] = 1;
        d[i] = 3.0 * (f[i + 1] - f[i - 1]);
    }
    a[n - 1] = 0;
    b[n - 1] = 1;
    d[n - 1] = 0;

    /* solve tridiagonal system */
    tridiagonal_matrix_algorithm(fx, a, b, c, d, n);

    /* free variables */
    free(a);
    free(b);
    free(c);
    free(d);

    return;

} /* cubic_spline_derivs_zero_bc */

/*--------------------------------------------------------------------------*/

void cubic_spline_second_derivs

    (double *f,   /* breakpoint values */
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
        fxx[i] = 2 * (-3 * f[i] + 3 * f[i + 1] - 2 * fx[i] - fx[i + 1]);
    }
    /* the last derivative is given by 6*a_0[n-2] + 2a_1[n-2]; */
    fxx[n - 1] = 6 * (2 * f[n - 2] + fx[n - 2] - 2 * f[n - 1] 
        + fx[n - 1]) + 2 * (-3 * f[n - 2] + 3 * f[n - 1] 
        - 2 * fx[n - 2] -  fx[n - 1]);

    return;

} /* cubic_spline_second_derivs */

/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*             functions for cubic monotonicity constraints                 */
/*             order 3 accurate monotone cubic interpolation                */
/*--------------------------------------------------------------------------*/

void modify_slopes_basic
   
   (double  *d,            /* slopes, length n */
    double  *delta,        /* secant lines, length n-1 */
    long    n)             /* length */

/*
   Fritsch-Carlson method (1980): 
   modifies the slopes d[i] for monotonicity preservation.  
   Let S2 be the circle with raduis 3 intersected with the first quadrant.
   The values d[i] are modified such that (d[i], d[i+1])
   lies in S2*delta[i] by projecting them in a straight line to the origin
   onto S2*delta[i].
*/

{
long  i;             /* index variable */
double x, y;         /* normalized slopes */
double r;            /* euclidian length of (x,y) */

for (i=0; i<n-1; i++)   
   {
   /* if secant line is zero set derivatives to zero */
   if (fabs(delta[i]) < EPS)
      {
      d[i] = 0.0;
      d[i+1] = 0.0;
      continue;
      }
      
   /* calcuate normalized slope values */
   x = d[i] / delta[i];
   y = d[i+1] / delta[i];
   r = sqrt(x*x + y*y);
   /* modify slopes if necessary */
   if (r > 3.0)
      {
      d[i] *= 3.0/r;
      d[i+1] *= 3.0/r;
      }
   }
}   /* modify_slopes_basic */

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

void enforce_sign

   (double  *d,      /* slopes */
    double  *delta,  /* secant lines */
    long    n)       /* length */
    
/*
   limits derivative values d[i] by setting them to zero if their sign does 
   not match the correct sign.
*/

{
long i;     /* loop variable */

/* ensure derivatives have the same sign as the adj. secant lines or zero */
for (i=1; i<n-1; i++) 
   {
   if ((sgn(d[i]) != sgn(delta[i-1])) | (sgn(d[i]) != sgn(delta[i])))
      d[i] = 0.0;
   }
if(sgn(d[0]) != sgn(delta[0]))
   d[0] = 0;
if(sgn(d[n-1]) != sgn(delta[n-2]))
   d[n-1] = 0;

return;

}  /* enforce_sign */


/*--------------------------------------------------------------------------*/
/*             order 4 accurate monotone cubic interpolation                */
/*                 functions for ext. two-sweep method                      */
/*--------------------------------------------------------------------------*/

double F

   (double x, 
    double y) 
    
/* 
   - used for the definition of the region M (see below)
   - zero set of F is an ellipse 
*/

{

return x*x + y*y + x*y - 6*x - 6*y + 9;

}  /* F */

/*--------------------------------------------------------------------------*/

/*
   Regions M,A,B,C,D,E divide the first quadrant into seperate regions;
   represented as boolean functions;
   precondition for all regions: x>0 and y>0
*/
  
/*--------------------------------------------------------------------------*/

int M

   (double x,        /* >= 0 */
    double y)        /* >= 0 */
    
/* 
   region M 
   consists of the area of an ellipse described by the zero set of 
   F or x + y <= 3;
   returns 1 if (x,y) lies in M and 0 otherwise
*/

{

return (x + y <= 3.0) | (F (x, y) <= 0.0);

}  /* M */

/*--------------------------------------------------------------------------*/

int A

   (double x,        /* >= 0 */
    double y)        /* >= 0 */
    
/*
   region A
   x + y <= 4 and x <= 1 and (x,y) not in M;
   returns 1 if (x,y) lies in A and 0 otherwise
*/

{

return (x + y <= 4) && (x <= 1) && !M (x, y);

}  /* A */

/*--------------------------------------------------------------------------*/

int B

   (double x,        /* >= 0 */
    double y)        /* >= 0 */

/*
   region B 
   x <= 3 and (x,y) not in A and (x,y) not in M;
   returns 1 if (x,y) lies in B and 0 otherwise
*/

{

return (x <= 3) && (x + y >= 4) && (F (x, y) >= 0);

}  /* B */ 

/*--------------------------------------------------------------------------*/

int C

   (double x,        /* >= 0 */
    double y)        /* >= 0 */

/* 
   region C 
   x >= 3 and y >= 3;
   returns 1 if (x,y) lies in C and 0 otherwise
*/

{

return (x >= 3.0) && (y >= 3.0);

}  /* C */

/*--------------------------------------------------------------------------*/

int D

   (double x,        /* >= 0 */
    double y)        /* >= 0 */

/* 
   region D 
   y <= 3 and (x,y) not in E and (x,y) not in M;
   returns 1 if (x,y) lies in D and 0 otherwise
*/

{

return (y <= 3) && (x + y >= 4) && F (x, y) >= 0;
   
}  /* D */

/*--------------------------------------------------------------------------*/


int E

   (double x,        /* >= 0 */
    double y)        /* >= 0 */

/* 
   region E 
   x + y <= 4 and x >= 1 and (x,y) not in M;
   special case for the extended_backward_sweep function; 
   returns 1 if (x,y) lies in E and 0 otherwise
*/

{

return (x + y <= 4) && (x >= 1) && (!M (x, y));

}  /* E */

/*--------------------------------------------------------------------------*/
/*                     extended sweep method                                */
/*                     Eisentstat et al. (1985)                             */
/*--------------------------------------------------------------------------*/

void extended_forward_sweep

   (double  *d,      /* slopes, length n */
    double  *delta,  /* secant lines, length n-1 */
    long    n)       /* length */

/*
   Implementation of the adapted verion of the extended forward sweep function
   by Eisenstat et al (1985): 
   "The Order Of monotone Piecewise Cubic Interpolation"
   SIAM Journal on Numerical Analysis

   extended forward sweep function 
   modifies the slopes (d[i],d[i+1]); projects them onto the union of 
   delta[i]*M and delta[i]*D for i = 0 ,..., n-2;
   only modifies d[i+1] in each iteration except for region A;
   increases or decreases them accordingly onto the boundary of the union 
   of M and D 
*/
    
{
long  i;             /* index variables */
double x,y;          /* normalized slopes */        

for (i = 0; i < n-1; i++) 
   {
   /* if secant line is zero then set derivatives to zero */
   if (fabs(delta[i]) < EPS)
      {
      d[i] = 0.0;
      d[i+1] = 0.0;
      continue;
      }
      
   /* calcuate normalized slope values */
   x = d[i] / delta[i];
   y = d[i+1] / delta[i];
   
   /* check if slopes have to be modified */
   if ((M(x,y)) | (y<=3.0))
      continue;
      
   /* project (x,y) onto the boundary of M or D */
   if (C(x,y))
      {
      y = 3.0;                                              /* decrease y */
      d[i+1] = delta[i] * y;                                /* update d[i+1] */
      }
   else if (B(x,y)) 
      {
      y = 3.0 - 0.5*x + 0.5*sqrt(3.0) * sqrt((4.0-x)*x);    /* decrease y */
      d[i+1] = delta[i] * y;                                /* update d[i+1] */
      }
   else if (A(x,y))
      {
      /* special case for local extrema */
      if( ((i>0 && (fabs(delta[i-1])<EPS) | (delta[i-1]*delta[i]<0))) )
         {
         y = 3.0 - 0.5*x + 0.5*sqrt(3.0) * sqrt((4.0-x)*x); /* decrease y */
         d[i+1] = delta[i] * y;                             /* update d[i+1] */
         continue;
         }
      x = 3.0 - 0.5*y - 0.5*sqrt(3.0) * sqrt((4.0-y)*y);    /* increase x */
      d[i] = delta[i] * x;                                  /* update d[i] */
   
      /* check if (d[i-1], d[i]) still lies in M */
      if (i==0)
         continue;
      x = d[i-1] / delta[i-1];
      y = d[i] / delta[i-1];
      if ( (M(x,y)) | (y<=3.0) )   
         continue;
      if (C(x,y))
         {
         y = 3.0;                                           /* decrease y */
         d[i] = delta[i-1] * y;                             /* update d[i] */
         }
      else if (A(x,y) | B(x,y))
         {
         y = 3.0 - 0.5*x + 0.5*sqrt(3.0) * sqrt((4.0-x)*x); /* decrease y */
         d[i] = delta[i-1] * y;                             /* update d[i] */
         }
         
      /* check if (d[i], d[i+1]) still lies in M */
      x = d[i] / delta[i];
      y  = d[i+1] / delta[i];
      if (!M(x,y))
         {
         y = 3.0 - 0.5*x + 0.5*sqrt(3.0) * sqrt((4.0-x)*x); /* decrease y */
         d[i+1] = delta[i] * y;                             /* update d[i+1] */
         }
      }
   }

return;

}   /* extended_forward_sweep */

/*--------------------------------------------------------------------------*/

void extended_backward_sweep
   
   (double  *d,            /* slopes, length n */
    double  *delta,        /* secant lines, length n-1 */
    long    n)             /* length */
/*
     Implementation of the adapted verion of the extended backward sweep function
     by Eisenstat et al (1985): 
     "The Order Of monotone Piecewise Cubic Interpolation"
     SIAM Journal on Numerical Analysis

     extended backward sweep function 
     modifies the slopes (d[i],d[i+1]); projects them onto delta[i]*M 
     for i = n-2 ,..., 0;
     only modifies d[i] in each iteration except for region E;
     If (x,y) lies in E then the (x,y) is not necessarily projected in a 
     straight (horizontal) line onto the boundary of M;  
     increases or decreases them accordingly onto the boundary of M;
*/

{
long  i;     /* index variables */
double x,y;  /* normalized slopes */      

for (i=n-2; i>=0; i--) 
   {
   /* if secant line is zero then skip */
   if (fabs(delta[i]) < EPS)
      continue;

   /* calcuate normalized slope values */
   x = d[i] / delta[i];
   y = d[i+1] / delta[i];
   /* check if slopes have to be modified */
   if ( (M(x,y)) | (x<=3.0) )
      continue;
   /* project (x,y) onto the boundary of M */
   if (D(x,y) | C(x,y))
      {
      x = 3.0 - 0.5*y + 0.5*sqrt(3.0) * sqrt((4.0-y)*y);    /* decrease x */
      d[i] = delta[i] * x;                                  /* update d[i] */
      }
   else if (E(x,y))
      {
      /* special case for boundaries and zero secant lines */
      if ( ((i<n-2 && (fabs(delta[i+1]) < EPS) | (delta[i]*delta[i+1]<0))) )
         {
         x = 3.0 - 0.5*y + 0.5*sqrt(3.0) * sqrt((4.0-y)*y); /* decrease x */
         d[i] = delta[i] * x;                               /* update d[i] */
         continue;
         }
      y = 3.0 - 0.5*x - 0.5*sqrt(3) * sqrt((4.0-x)*x);      /* increase y */
      d[i+1] = delta[i] * y;                                /* update d[i+1] */
      /* check if (d[i+1], d[i+2]) still lies in M */
      if (i==n-2)
         continue;
      x = d[i+1]/delta[i+1];
      y = d[i+2]/delta[i+1];
      if (M(x,y))
         continue;
      x = 3.0 - 0.5*y + 0.5*sqrt(3.0) * sqrt((4.0-y)*y);    /* decrease x */
      d[i+1] = delta[i+1] *x;                               /* update d[i+1] */
      x = d[i] / delta[i];
      y = d[i+1] / delta[i];
      
      /* check if (d[i], d[i+1]) still lies in M */
      if (!M(x,y))
         {
         x = 3.0 - 0.5*y + 0.5*sqrt(3.0) * sqrt((4.0-y)*y); /* decrease x */
         d[i] = delta[i] * x;                               /* update d[i+1] */
         }  
      }
      
   }
   
return;

}  /* extended_backward_sweep */


/*--------------------------------------------------------------------------*/
/*        functions for piecewise quintic hermite interpolation             */
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/

void quintic_hermite_coeffs

    (double *f,     /* breakpoint values */
     double *fx,    /* breakpoint approximated first derivatives */
     double *fxx,   /* breakpoint approximated second derivatives */
     double *coeff, /* coefficient vector of size (n-1)*6 (output) */
     long n)        /* number of breakpoints */

/*
    calculates coefficients of the piecewise quintic hermite for each interval
    and stores them in the coefficient matrix;
    the polynomials for each interval all have the domain [0,1]
    *equidistand version (h=1) *
*/

{

    long i; /* loop variable */

    for (i = 0; i < n - 1; i++)
    {
        coeff[6 * i] =
            -6.0 * f[i] - 3.0  * fx[i] - 0.5 * fxx[i] 
            + 0.5 * fxx[i + 1] - 3.0 * fx[i + 1] 
            + 6.0 * f[i + 1];
        coeff[6 * i + 1] =
            +15.0 * f[i] + 8.0 * fx[i] + 3.0 / 2 * fxx[i] 
            - 1.0 * fxx[i + 1] + 7.0  * fx[i + 1] - 15 * f[i + 1];
        coeff[6 * i + 2] =
            -10.0 * f[i] - 6.0 * fx[i] - 3.0 / 2  * fxx[i] 
            + 0.5 * fxx[i + 1] - 4.0 * fx[i + 1] + 10.0 * f[i + 1];
        coeff[6 * i + 3] = 0.5 * fxx[i];
        coeff[6 * i + 4] =  fx[i];
        coeff[6 * i + 5] = f[i];
    }

    return;

} /* quintic_hermite_coeffs */


/*--------------------------------------------------------------------------*/
/*                MONOTONE QUINTIC SPLINE INTERPOLATION                     */
/*                           (Lux 2020)                                     */
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
    returns TRUE if f is monotone increasing on [x0,x1],
    FALSE otherwise

    This function is a direct portation of the Fortran 2003 routine by Thomas Lux.
    Most comments have been adopted from original source code documentation.
    Original code:
    Algorithm 1031: MQSI - Monotone Quintic Spline Interpolation
    Thomas C.H. Lux, Layne T. Watson, William Thacker, Tyler H. Chang. 
    ACM Transactions on Mathematical Software. Submitted August, 2020
    https://github.com/tchlux/papers/tree/master/%5B2020-08%5D_ACMTOMS_(MQSI)
*/

{
    /* local variables */
    double w, z, alpha, beta, gamma, t, sign;

    /* When the function values are flat, everything *must* be 0 */
    if (fabs(f0 - f1) < EPS) {
        if (!(fabs(fx0) < EPS))
        {
            return FALSE;
        }
        if (!(fabs(fx1) < EPS))
        {
            return FALSE;
        }
        if (!(fabs(fxx0) < EPS))
        {
            return FALSE;
        }
        if (!(fabs(fxx1) < EPS))
        {
            return FALSE; 
        }
        return TRUE; 
    }
    /* Identify the direction of change of the function 
       (increasing or decreasing) */
    if(f1 > f0) {
        sign = 1.0;
    }
    else {
        /* try something else to enforce symmetry */
        sign = 1.0;
        t = f0;
        f0 = f1;
        f1 = t;
        t = fx0;
        fx0 = fx1;
        fx1 = t;
        t = fxx0;
        fxx0 = fxx1;
        fxx1 = t;
        fx0 = -fx0;
        fx1 = -fx1;
    }
    /* Determine which set of monotonicity conditions to use based on the
       assigned first derivative values at either end of the interval. */
    if ((sign*fx0 < 0.0) || (sign*fx1 < 0.0)) {
        return FALSE;
    }

    w = x1 - x0;
    z = f0 - f1;

    if ((fabs(fx0) < EPS) || (fabs(fx1) < EPS)) {
        /* Simplified monotone case, which reduces to a test of cubic
        positivity studied in
        
        J. W. Schmidt and W. He{\ss}, ``Positivity of cubic polynomials on
        intervals and positive spline interpolation'', {\sl BIT Numerical
        Mathematics}, 28 (1988) 340--352.
        
        Notably, monotonicity results when the following conditions hold:
            alpha >= 0,
            delta >= 0,
            beta  >= alpha - 2 * sqrt{ alpha delta },
            gamma >= delta - 2 * sqrt{ alpha delta },
        
        where alpha, beta, delta, and gamma are defined in the paper. The
        equations that follow are the result of algebraic simplifications
        of the terms as they are defined by Schmidt and He{\ss}.
        
        The condition delta >= 0 was enforced when first estimating
        derivative values (with correct sign). Next check for alpha >= 0 */
        if (sign*fxx1 * w > sign*4.0 * fx1) {
            return FALSE; 
        }
        /* Else compute a simplification of 2 * sqrt{ alpha delta } */
        t = fx0 * (4*fx1 - fxx1 * w);
        if(t > 0.0) {
            t = 2.0 * sqrt(t);
        }
        /* Check for gamma >= delta - 2 * sqrt{ alpha delta } */
        if (t + sign*(3.0 * fx0 + fxx0 * w) < 0.0) {
            return FALSE; 
        }
        /* Check for beta >= alpha - 2 * sqrt{ alpha delta } */
        if (-60.0 * z*sign - w * (sign*(24.0 * fx0 + 32.0 * fx1) 
            - 2.0 * t + w*sign * (3.0 * fxx0 - 5.0 * fxx1)) < 0.0) {
            return FALSE; 
        }
        return TRUE; 
    }
    /* Full quintic monotonicity case related to the theory in 
    G. Ulrich and L. T. Watson, ``Positivity conditions for quartic 
    polynomials'', {\sl SIAM J. Sci. Comput.}, 15 (1994) 528--544.
    
    Monotonicity results when the following conditions hold:
    tau_1 > 0,
    if (beta <= 6) then
        alpha > -(beta + 2) / 2,
        gamma > -(beta + 2) / 2,
    else
        alpha > -2 sqrt(beta - 2),
        gamma > -2 sqrt(beta - 2),
    end if
    
    where alpha, beta, gamma, and tau_1 are defined in the paper. The
    following conditions are the result of algebraic simplifications
    of the terms as defined by Ulrich and Watson. */
        
    /* Check for tau_1 > 0.*/
    if (w * (2.0 * sqrt(fx0 * fx1) - sign*3.0 * (fx0 + fx1)) 
        - sign*24.0 * z <= 0.0) {
        return FALSE; 
    }
    /* Compute alpha, gamma, beta from theorems to determine monotonicity */
    t = pow(fx0 * fx1, 0.75);
    alpha = sign * (4 * fx1 - fxx1 * w) * sqrt(sign*fx0) / t;
    gamma = sign * (4 * fx0 + fxx0 * w) * sqrt(sign*fx1) / t;
    beta = sign * (-60 * z / w + 3 * (w * (fxx1 - fxx0) 
           - 8 * (fx0 + fx1))) / (2 * sqrt(fx0 * fx1));

    if (beta <= 6) {
        return (fmin(alpha, gamma) > -(beta + 2) / 2);
    }
    if(!(fmin(alpha, gamma) > -2 * sqrt(beta - 2))) {
    }
    return (fmin(alpha, gamma) > -2 * sqrt(beta - 2));

    return FALSE;

}  /* is_monotone */


/*--------------------------------------------------------------------------*/

int is_monotone_index

    (long   index,/* index to check for monotonicity*/
     double *f,   /* discrete function values */
     double *fx,  /* discrete derivatives */
     double *fxx) /* discrete second derivatives*/

{

/*
    function to make call of is_monotone easier
*/

return is_monotone(
    index, 
    index+1, 
    f[index], 
    f[index+1], 
    fx[index], 
    fx[index+1], 
    fxx[index], 
    fxx[index+1]);

}  /* is_monotone_index */

/*--------------------------------------------------------------------------*/

void mqsi 
    
    (double *f,             /* discrete function values to interpolate */
     double *coeff,         /* quintic hermite spline coefficients */
     long   n)              /* number of function values */

/*
    This function is a direct portation of the Fortran 2003 routine by Thomas Lux.
    Most comments have been adopted from original source code documentation.
    Original code and comments:
    Algorithm 1031: MQSI - Monotone Quintic Spline Interpolation
    Thomas C.H. Lux, Layne T. Watson, William Thacker, Tyler H. Chang. 
    ACM Transactions on Mathematical Software. 
    https://github.com/tchlux/papers/tree/master/%5B2020-08%5D_ACMTOMS_(MQSI)
*/

{
double accuracy; 
/* Estimated first and second derivatives by quadratic
   facet model at all data points" */
double *fx;         /* estimated first derivatives */
double *fxx;        /* estimated second derivatives */
double *fx_copy;    /* estimated first derivatives copy */
double *fxx_copy;   /* estimated second derivatives copy */
/* Execution queues holding intervals to check for monotonicity
   checking, to_check, derivative values to grow (after shrinking)
   growing, to_grow, and derivative values to shrink (because of
   nonmonotonicity) shrinking, to_shrink." The arrays flats and
   extrema are used to identify locations of local maxima and
   minima of provided Y values early in the routine */
int *checking;      
int *growing;
int *shrinking;
int *flats;
int *extrema;
long *to_check;
long *to_grow;
long *to_shrink;


/* Coefficients for a quadratic interpolant A and B, the direction of
   function change direction, a derivative value dx, and the current
   bisection search (relative) step size step_size */
double A, B, direction, dx, step_size;
/* Iteration indices i and j, number of data points nd, counters for
   execution queues, checking, nc, growing, ng, and shrinking, ns */
long i, j, nc, ng, ns;
/* Boolean indicating whether the bisection search is in progress */
int searching;

accuracy = sqrt(EPS);

// to avoid complier warnings
direction = 0.0;
B = 0.0;
A = 0.0;

/* allocate memory */
fx = malloc(sizeof(double)*n);
fxx = malloc(sizeof(double)*n);
checking = malloc(sizeof(int)*n);
growing = malloc(sizeof(int)*n);
shrinking = malloc(sizeof(int)*n);
flats = malloc(sizeof(int)*n);
extrema = malloc(sizeof(int)*n);
to_check = malloc(sizeof(long)*n);
to_grow = malloc(sizeof(long)*n);
to_shrink = malloc(sizeof(long)*n);
fx_copy = malloc(sizeof(double)*n);
fxx_copy = malloc(sizeof(double)*n);


/*---------------------------------------------------------------------------*/
/*          Algorithm 1: Estimate initial derivatives by                     */
/*         using a minimum curvature quadratic facet model                   */
/*---------------------------------------------------------------------------*/

/* Identify local extreme points and flat points */
flats[0] = fabs(f[0] - f[1]) < EPS * (1.0 + fabs(f[0] + fabs(f[1])));
flats[n-1] = fabs(f[n-2] - f[n-1]) < EPS * (1.0 + fabs(f[n-2]) + fabs(f[n-1]));
extrema[0] = 0;
extrema[n-1] = 0;

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
if(TRUE) {
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
        fxx[i] = (f[i-1] - f[i]);
        /* Compute the coefficient A in  Ax^2+Bx+C  that 
           interpolates at h(i+1) */
        A = (f[i+1] - f[i]);
        if(fabs(A) < fabs(fxx[i])) {
            fxx[i] = A;
        }
        /* Compute the actual second derivative (instead of coefficient A) */
        fxx[i] *= 2.0;
    }
    else {
      /* Determine the direction of change at the point i */
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
    
        /*-------------------------------------------------------------------*/
        /* Quadratic left of i */
        if((i-1 > 0)) {
            /* If a zero derivative at left point, use its right interpolant */
            if(flats[i-1] || extrema[i-1]) {
                A = (f[i]-f[i-1]);
                B = 0.0;
            }
            /* Otherwise use the standard quadratic on the left */
            else {
                A = (f[i-2] + f[i] - 2*f[i-1])/2;
                B = (f[i] - f[i-2])/2;
            }
        
            dx = 2.0 * A  + B;
            if(dx*direction >= 0.0) {
                fx[i] = dx;
                fxx[i] = A;
            }
        }
        
        /*-------------------------------------------------------------------*/
        /* Quadratic centered at i (require that it has at least one
        neighbor that is not forced to zero slope) */
        if((i > 0) && (i < n-1)) {
            if(!((extrema[i-1] || flats[i-1]) 
            && (extrema[i+1] || flats[i+1]))) {
                /* Construct quadratic interpolant through this point and 
                neighbors */
                A = (f[i-1] + f[i+1] - 2*f[i])/2;
                B = (f[i+1] - f[i-1])/2;
                dx = B;
                /* Keep this new quadratic if it has less curvature */
                /* add a case to original code for symmetry: 
                   if i is on the left of the center of the signal f
                   then check for "<", if i is on the right of the center
                   then check for "<=". This elminates the problem of asymmetric
                   derivative choices at the "=" case */
                if(i <= n/2-1) { // left side of the signal f
                    if((dx*direction >= 0.0) && (fabs(A) < fabs(fxx[i]))) {
                        fx[i] = dx;
                        fxx[i] = A;
                    }
                }
                else { // right side of the signal f
                    if((dx*direction >= 0.0) && (fabs(A) <= fabs(fxx[i]))) {
                        fx[i] = dx;
                        fxx[i] = A;
                    }
                }
            }
        }
        /* ----------------------------------------------------------------- */
        /* Quadratic right of i */
        if((i+1 < n-1)) { 
            /* If a zero derivative at right point, use its left interpolant */
            if(extrema[i+1] || flats[i+1]) {
                A = (f[i]-f[i+1]);
                B = 0.0;
            }
            /* Otherwise use the standard quadratic on the right */
            else {
                A = (f[i] + f[i+2] - 2*f[i+1])/2;
                B = (f[i+2] - f[i])/2;
            }
            dx = 2.0 * A * (-1) + B;
            /* Keep this new quadratic if it has less curvature */
            /* add a case to original code for symmetry: 
            if i is on the left of the center of the signal f
            then check for "<", if i is on the right of the center
            then check for "<=". This elminates the problem of asymmetric
            derivative choices at the "=" case */
            if(i <= n/2-1) { // left side of the signal f
                if((dx*direction >= 0.0) && (fabs(A) < fabs(fxx[i]))) {
                    fx[i] = dx;
                    fxx[i] = A;
                }
            }
            else { // right side of the signal f
                if((dx*direction >= 0.0) && (fabs(A) <= fabs(fxx[i]))) {
                    fx[i] = dx;
                    fxx[i] = A;
                }
            }
        }
        /* "Set the final quadratic" */
        if(fxx[i] == 1e54) {
            fx[i] = 0.0;
            fxx[i] = 0.0;
        }
        /* Compute curvature of quadratic from coefficient of x^2 */
        else {
            fxx[i] = 2.0*fxx[i];
        }
    }
}
}

/*---------------------------------------------------------------------------*/
/*          TEST: Estimate initial derivatives by                            */
/*         using a piecewise cubic splines in monotone intervals             */
/*          with zero-derivative boundary conditions                         */
/*---------------------------------------------------------------------------*/
if(FALSE) {
/* cubic spline derivs */
/* start with first interval */
i = 0;
/* initialize curvature to be maximally large */
A = 1e54;
/* find each monotone section of the signal */
for(j = 1; j < n; j++) {
    /* loop until an extremum or flat point */
    if(!((flats[j]) || (extrema[j]))) {
        continue;
    }
    /* check if interval [i,j] is monotone and not flat */
    if((extrema[i]) || (extrema[j]) || (j == n-1)) {
        /* calculate cubic spline derivatives for the monotone interval */
        /* use zero-derivative boundary conditions */
        cubic_spline_derivs_zero_bc(f+i, fx+i, j-i+1);
        cubic_spline_second_derivs(f+i, fx+i, fxx+i, j-i+1);
        /* keep the smaller curvature */
        if((fabs(A) < fabs(fxx[i]))) {
            fxx[i] = A;
        }
        /* keep track of the current curvature to compare next step */
        A = fxx[j];
    } 
    /* start at the end of the current monotone section */
    i = j;
}

/* set flats to zero */
for(i = 0; i < n; i++) {
    if(flats[i]) {
        fx[i] = 0.0;
        fxx[i] = 0.0;
    }
}}


/*---------------------------------------------------------------------------*/
/*         Algorithm 3: Identify viable piecewise monotone                   */
/*        derivative values by doing a quasi-bisection search                */
/*---------------------------------------------------------------------------*/

if(TRUE) { // hard coded switch for testing initial derivatives

/* Store the initially estimated first and second derivative values */
for(i = 0; i < n; i++) {
    fx_copy[i] = fx[i];
    fxx_copy[i] = fxx[i];
}

/* Identify which spline segments are not monotone and need to be fixed */
for(i = 0; i < n; i++) {
    checking[i] = 0;
    shrinking[i] = 0;
}
nc = 0;
ng = 0;
ns = 0;
for(i = 0; i < n-1; i++) {
    /* Check for monotonicity on all segments that are not flat */
    if((!(flats[i] && flats[i+1])) 
    && (!(is_monotone_index(i, f, fx, fxx)))) {
        /* Store points bounding nonmonotone segments in the to_shrink queue */
        if(!shrinking[i]) {
            shrinking[i] = TRUE; 
            to_shrink[ns] = i;
            ns++;
        }
        if(!shrinking[i+1]) {
            shrinking[i+1] = TRUE;
            to_shrink[ns] = i+1;
            ns++;
        }
    }
}

/* Initialize step size to 1.0 (will be halved at beginning of loop) */
step_size = 1.0;
searching = TRUE;
for(i = 0; i < n; i++) {
    growing[i] = FALSE;
}

/* Loop until the accuracy is achieved and *all* intervals are monotone */
while((searching || (ns > 0))) {
    /* Compute the step size for this iteration */
    if(searching) {
        step_size = step_size/2.0;
        if(step_size < accuracy) {
            searching = FALSE;
            step_size = accuracy;
            ng = 0;
        }
    }
    /* Grow the step size (at a slower rate than the step size reduction
        rate) if there are still intervals to fix */
    else {
        step_size = step_size * 1.5;
    }

    /* Grow all those first and second derivatives that were previously 
        shrunk, and correspond to currently monotone spline pieces */

    /* grow_values */
    for(j = 0; j < ng; j++) {
        i = to_grow[j];
        /* Do not grow values that are actively related to a nonmonotone
            spline segment */
        if(shrinking[i]) {
            continue;
        }
        /* Otherwise, grow those values that have been modified 
            previously */
        fx[i] += step_size * fx_copy[i];
        fxx[i] += step_size * fxx_copy[i];
        /* Make sure the first derivative does not exceed its original 
            value */
        if((fx_copy[i] < 0.0) && (fx[i] < fx_copy[i])) {
            fx[i] = fx_copy[i];
        }
        else if((fx_copy[i] > 0.0) && (fx[i] > fx_copy[i])) {
            fx[i] = fx_copy[i];
        }
        /* Make sure the second derivative does not exceed its original 
            value */
        if((fxx_copy[i] < 0.0) && (fxx[i] < fxx_copy[i])) {
            fxx[i] = fxx_copy[i];
        }
        else if((fxx_copy[i] > 0.0) && (fxx[i] > fxx_copy[i])) {
            fxx[i] = fxx_copy[i];
        }
        /* Set this point and its neighboring intervals to be checked for
            monotonicity. Use sequential if statements for 
            short-circuiting */
        if(i > 0) if(!(checking[i-1])) {
            checking[i-1] = TRUE;
            to_check[nc] = i-1;
            nc++;
        }
        if(i < n-1) if(!checking[i]) {
            checking[i] = TRUE;
            to_check[nc] = i;
            nc++;
        }
    }  /* grow_values */
    
    /* Shrink the first and second derivatives that cause 
        nonmonotonicity */
    /* shrink_values */
    for(j = 0; j < ns; j++) {
        i = to_shrink[j];
        shrinking[i] = FALSE;
        if(searching && (!growing[i])) {
            growing[i] = TRUE;
            to_grow[ng] = i;
            ng++;
        }
        /* Shrink the values that are causing nonmonotonicity */
        fx[i] -= step_size * fx_copy[i];
        fxx[i] -= step_size * fxx_copy[i];
        /* Make sure the first derivative does not pass zero */
        if((fx_copy[i] < 0.0) && (fx[i] > 0.0)) {
            fx[i] = 0.0;
        }
        else if((fx_copy[i] > 0.0) && (fx[i] < 0.0)) {
            fx[i] = 0.0;
        }
        /* Make sure the second derivative does not pass zero */
        if((fxx_copy[i] < 0.0) && (fxx[i] > 0.0)) {
            fxx[i] = 0.0;
        }
        else if((fxx_copy[i] > 0.0) && (fxx[i] < 0.0)) {
            fxx[i] = 0.0;
        }
        /* Set this point and its neighboring intervals to be checked for
            monotonicity */
        if(i > 0) {
            if(!checking[i-1]) {
                checking[i-1] = TRUE;
                to_check[nc] = i-1;
                nc++;
            }
        }
        if(i < n-1) {
            if(!checking[i]) {
                checking[i] = TRUE;
                to_check[nc] = i;
                nc++;
            }
        }
    }  /* shrink_values*/
    /* Identify which spline segments are nonmonotone after the updates */
    ns = 0;
    /* check_monotonicity */
    for(j = 0; j < nc; j++) {
        i = to_check[j];
        checking[i] = FALSE;
        if(!is_monotone_index(i, f, fx, fxx)) {
            if(!shrinking[i]) {
                shrinking[i] = TRUE;
                to_shrink[ns] = i;
                ns++;
            }
            if(!shrinking[i+1]) {
                shrinking[i+1] = TRUE;
                to_shrink[ns] = i+1;
                ns++;
            }
        }
    }  /* check_monotonicity */
    nc = 0;
}
}
/* calculate quintic spline representation in terms of hermite coefficients */
quintic_hermite_coeffs(f, fx, fxx, coeff, n);

/* free memory */
free(fx);
free(fxx);
free(checking);
free(growing);
free(shrinking);
free(to_check);
free(to_grow);
free(to_shrink);
free(flats);
free(extrema);
free(fx_copy);
free(fxx_copy);

return;

}  /* mqsi */


/*--------------------------------------------------------------------------*/
/*                      interpolation methods                               */        
/*--------------------------------------------------------------------------*/

void linear_interpolation

   (double  *Y,         /* y-values */ 
    double  *coeff,     /* hermite spline coefficients (output) */
    long    n)          /* number of (x,y) values */

/*
   sets the hermite spline coefficients for linear interpolation
*/

{
long  i;    /* loop variale */
for (i=0; i<n-1; i++)
   {
   coeff[4*i] = 0;
   coeff[4*i+1] = 0;
   coeff[4*i+2] = Y[i+1]-Y[i];     /* set slope */
   coeff[4*i+3] = Y[i];            /* set shift */
   }

return;

}  /* linear_interpolation */

/*--------------------------------------------------------------------------*/

void cubic_spline_interpolation
   
   (double  *Y,         /* breakpoint values */ 
    double  *coeff,     /* hermite spline coefficients (output) */
    long    n)          /* number of (x,y) values */

/*
   interpolation with cubic spline derivatives;
   calculates the hermite coefficients; 
   nonmonotone, order 4 accurate
   * equidistant version (h=1) *
*/

{
long     i;      /* index variable */
double   *d;     /* approximate derivatives */
double   *delta; /* secant lines */

/* allocate memory */
alloc_double_vector(&d,n);
alloc_double_vector(&delta,n-1);

/* initialize delta */
for (i=0; i<n-1; i++) 
   delta[i] = (Y[i+1]-Y[i]);

/* calculate approximate derivatives */
cubic_spline_derivs_eq (Y, d, n); 

/* calculate cubic hermite coefficients */
hermite_coeffs(Y, d, coeff, n);

/* free variables */
free_double_vector(d,n);
free_double_vector(delta,n-1);

}  /* cubic_spline */         

/*--------------------------------------------------------------------------*/

void mpci_basic
   
   (double  *Y,         /* breakpoint values */ 
    double  *coeff,     /* hermite spline coefficients (output) */
    long    n)          /* number of (x,y) values */
/*
   monotone piecewise cubic interpolation with Fritsch-Carlson method;
   calculates the hermite coefficients; 
   monotone, order 3 accurate
   * equidistant version (h=1) *
*/

{
long     i;      /* index variable */
double   *d;     /* approximate derivatives */
double   *delta; /* secant lines */

/* allocate memory */
alloc_double_vector(&d,n);
alloc_double_vector(&delta,n-1);

/* initialize delta */
for (i=0; i<n-1; i++) 
   delta[i] = (Y[i+1]-Y[i]);

/* STEP 1: calculate approximate derivatives */
cubic_spline_derivs_eq (Y, d, n); 

/* STEP 2: ensure that d[i] have the right sign */
enforce_sign(d, delta, n);

/* STEP 3: modify derivatives to preserve monotonicity */
modify_slopes_basic(d, delta, n);

/* calculate cubic hermite coefficients */
hermite_coeffs(Y, d, coeff, n);

/* free variables */
free_double_vector(d,n);
free_double_vector(delta,n-1);

}  /* mpci_basic */

/*--------------------------------------------------------------------------*/

void mpci_extended_sweep
   
   (double  *Y,         /* breakpoint values */ 
    double  *coeff,     /* hermite spline coefficients (output) */
    long    n)          /* number of (x,y) values */

/*
   monotone piecewise cubic interpolation with extended sweep method;
   calculates the hermite coefficients;
   monotone, order 4 accurate
   * equidistant version (h=1) *
*/

{
long     i;      /* index variable */
double   *d;     /* approximate derivatives */
double   *delta; /* secant lines */

/* allocate memory */
alloc_double_vector(&d,n);
alloc_double_vector(&delta,n-1);

/* initialize delta */
for (i=0; i<n-1; i++) 
   delta[i] = (Y[i+1]-Y[i]);

/* STEP 1: calculate approximate derivatives */
cubic_spline_derivs_eq (Y, d, n); 
   
/* STEP 2: ensure that d[i] have the right sign */
enforce_sign(d, delta, n);

/* STEP 3: modify derivatives to preserve monotonicity */
extended_forward_sweep (d, delta, n);
extended_backward_sweep (d, delta, n);

/* calculate cubic hermite coefficients */
hermite_coeffs(Y, d, coeff, n);

/* free variables */
free_double_vector(d,n);
free_double_vector(delta,n-1);

}  /* mpci_extended_sweep */


/*--------------------------------------------------------------------------*/

void interpolate_piecewise_Hermite

   (double  *f,      /* input samples */
    long    n,       /* number of sample points */
    long    method,  /* interpolation method */
    double  *coeff) /* coefficients (output) */

/*
   interpolates f using the specified piecwise Hermite interpolation method
*/

{
void     (*interpolation)(double*, double*, long); /* interpolation method */

/* initialize interpolation method */
switch(method) {
case 1:
   interpolation = &linear_interpolation;
   break;
case 2:
   interpolation = &cubic_spline_interpolation;
   break;
case 3:
   interpolation = &mpci_basic;
   break;
case 4: 
   interpolation = &mpci_extended_sweep;
   break;
case 5: 
   interpolation = &mqsi;
   break;
default:
   printf("interpolate_image_row_to_file: invalid method: %ld", method);
   exit(1);
}

/* interpolate data with given method */
interpolation(f, coeff, n);

return; 

}  /* interpolate_piecewise_Hermite */


/*--------------------------------------------------------------------------*/
/*                      image interpolation                                 */
/*--------------------------------------------------------------------------*/

void interpolate_image_horizontal
     
     (long     nx,              /* image dimension in x direction */
      long     ny,              /* image dimension in y direction */
      double   **u,             /* input image */
      double   **coeff,         /* coefficients in x direction, output */
      long     method)

/*
   interpolates image in x and y direction along pixels;
   returns coefficients of piecewise cubic polynomial interpolant
   * equidistant version (h=1) *
*/

{
long     i, j;      /* loop variables */
double   **Y;       /* values in x-direction (size ny * nx) */

void (*interpolation)(double*, double*, long);

switch(method) {
case 1:
   interpolation = &linear_interpolation;
   break;
case 2:
   interpolation = &cubic_spline_interpolation;
   break;
case 3:
   interpolation = &mpci_basic;
   break;
case 4: 
   interpolation = &mpci_extended_sweep;
   break;
case 5: 
   interpolation = &mqsi;
   break;
default:
   printf("interpolate_image_horizontal: invalid method: %ld", method);
   exit(1);
}

alloc_double_matrix (&Y, ny+2, nx+2);

/* initialize breakpoints */
for (j=0; j<=ny+1; j++)
   for (i=0; i<=nx+1; i++)
      Y[j][i] = u[i][j];      /* Y-values of j-th image row */

      
/* ---- loop for x-direction ---- */

for (j=0; j<ny+2; j++) {
   interpolation(Y[j], (coeff)[j], nx+2);
}

 
/* free memory */
free_double_matrix (Y, ny+2, nx+2);

return;

}  /* interpolate_image_horizontal */

/*--------------------------------------------------------------------------*/

void interpolate_image_vertical

     (long     nx,              /* image dimension in x direction */
      long     ny,              /* image dimension in y direction */
      double   **u,             /* input image */
      double   **coeff,         /* coefficients in y direction, output */
      long     method)
/*
   interpolates image in y direction along pixels;
   returns coefficients of piecewise cubic polynomial interpolant
   * equidistant version *
*/

{
long     i;         /* loop variable */

void (*interpolation)(double*, double*, long);

switch(method) {
case 1:
   interpolation = &linear_interpolation;
   break;
case 2:
   interpolation = &cubic_spline_interpolation;
   break;
case 3:
   interpolation = &mpci_basic;
   break;
case 4: 
   interpolation = &mpci_extended_sweep;
   break;
case 5: 
   interpolation = &mqsi;
   break;
default:
   printf("interpolate_image_vertical: invalid method: %ld", method);
   exit(1);
}

/* ---- loop for y-direction ---- */
for (i=0; i<nx+2; i++) {
   interpolation(u[i], (coeff)[i], ny+2);
}


return;

}  /* interpolate_image_vertical */


/*--------------------------------------------------------------------------*/

void interpolate_image

     (long     nx,              /* image dimension in x direction */
      long     ny,              /* image dimension in y direction */
      double   **u,             /* input image */
      double   **coeff_x,       /* coefficients in x direction, output */
      double   **coeff_y,       /* coefficients in y direction, output */
      long     method)          /* interpolation method */

/*
   interpolates image in x and y direction along pixels;
   returns coefficients of piecewise cubic polynomial interpolant
*/

{

interpolate_image_horizontal (nx, ny, u, coeff_x, method);
interpolate_image_vertical (nx, ny, u, coeff_y, method);

return;

}  /* interpolate_image */
      
/*--------------------------------------------------------------------------*/
/*      get interpolated adjacent image values in a given direcion          */
/*--------------------------------------------------------------------------*/

void approx_neighbours

    (double    *neighbour1,    /* value of neighbour 1, output */
     double    *neighbour2,    /* value of neighbour 2, output */
     double    *h_sqr,         /* squared distance to a neighbour, output*/
     long      i,              /* row of the pixel in the image */
     long      j,              /* column of the pixel in the image */
     double    dir_x,          /* x component of the normed direction */
     double    dir_y,          /* y component of the normed direction */
     double    h,              /* pixel size in x and y direction */
     long      nx,             /* image dimension in x direction */
     long      ny,             /* image dimension in y direction */
     double    **u,            /* image */
     double    **coeff_x,      /* quintic coefficients in x direction */
     double    **coeff_y,      /* quintic coefficients in x direction */
     long      method)         /* interpolation method */

/*
   approximates the neighbouring image values of u at pixel (i,j) in the 
   given direction;
   use interpolated image to obtain values between pixels
*/

{
double   offs;                             /* offset variable, time saver */
long     p;                                /* 1 + degree of the used spline */
double   (*eval_hermite)(double, double*); /* hermite evaluation function */

/* initialize the spline order depending on the used method */  
p = method == 5 ? 6 : 4;

/* initialize hermite interpolation method*/
eval_hermite = p == 6 ? &eval_quintic_hermite_poly : &eval_cubic_hermite_poly;

/* let alpha be the angle of the direction vector on the unit circle */
/* case 1: |cos(alpha)| <= sqrt(2)/2 (intersection with horizontal line)*/
if (fabs(dir_x) <= fabs(dir_y)) 
  {
  /* calculate offset from starting position (i+1)*h */
  offs = 1.0 + dir_x / dir_y; // (0 < offs < 2 since |dir_x| <= |dir_y|)
  
  /* calculate h_sqr = h^2 * || (xi_x/xi_y, 1)^T ||^2 */
  *h_sqr = h * h * (1.0 + dir_x * dir_x / (dir_y * dir_y));
  
  /* evaluate interpolated image at intersection points */
  /* evaluate corresponding hermite polynomial */ 
  if (offs < 1.0) 
     {
     *neighbour1 = eval_hermite(1-offs, coeff_x[j-1] + p * i); 
     *neighbour2 = eval_hermite(offs, coeff_x[j+1] + p * (i-1)); 
     }
  else 
     {
     *neighbour1 = eval_hermite(2-offs, coeff_x[j-1] + p *(i-1)); 
     *neighbour2 = eval_hermite(offs-1, coeff_x[j+1] + p *i); 
     }
  }

  
/* case 2: |cos(alpha)| > sqrt(2)/2 (intersection with vertical line) */  
else 
  {
  /* calculate offset from starting position (j+1)*h */
  offs = 1.0 + dir_y / dir_x; // (0 < offs < 2 since |dir_y| < |dir_x|)
  
  /* calculate h_sqr = h^2 * || (1, dir_y/dir_x)^T ||^2 */
  *h_sqr = h * h * (1.0 + dir_y * dir_y / (dir_x * dir_x));
  
  /* evaluate interpolated image at intersection points */
  /* evaluate corresponding hermite polynomial */      
 if (offs < 1.0) 
     {
     *neighbour1 = eval_hermite (1-offs, coeff_y[i-1] + p * j); 
     *neighbour2 = eval_hermite (offs, coeff_y[i+1] + p * (j-1)); 
     }
  else
     {
     *neighbour1 = eval_hermite (2-offs, coeff_y[i-1] + p * (j-1)); 
     *neighbour2 = eval_hermite (offs-1, coeff_y[i+1] + p * j); 
     }  
  } 

return;

}  /* approx_neighbours */



/*--------------------------------------------------------------------------*/
/*              MEAN CURVATURE MOTION PDE EXPLICIT TIME SCHEME              */
/*--------------------------------------------------------------------------*/

void mcm_delta 

     (double   tau,       /* time step size */
      double   alpha,     /* dissipativity parameter in delta stencil */
      double   gamma,     /* nonnegativity parameter in delta stencil */
      long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      double   h,         /* pixel size in x and y direction */
      double   **u)       /* image, changed */

/*
   original source code by Joachim Weickert (2021)
   for the explicit delta scheme
   explicit scheme for MCM with delta stencil
*/

{
long     i, j;                  /* loop variables */
double   **f;                   /* image u at old time level */
double   fx, fy;                /* Sobel derivatives */
double   f00, f11, f22, f33;    /* 2nd derivatives in principal directions */
double   inv_4h, inv_8h;        /* time savers */
double   inv_hh, inv_2hh;       /* time savers */
double   inv_grad_sqr, aux;     /* time savers */
double   a, b, c;               /* weights for u_{xx}, u_{xy}, u_{yy} */
double   delta;                 /* free parameter in the delta stencil */


/* ---- allocate memory ---- */

alloc_double_matrix (&f, nx+2, ny+2);


/* ---- copy u to f ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     f[i][j] = u[i][j];


/* ---- reflecting dummy boundaries ---- */

dummies_double (f, nx, ny);


/* ---- compute time savers ---- */

inv_4h  = 1.0 / (4.0 * h);
inv_8h  = 1.0 / (8.0 * h);
inv_hh  = 1.0 / (h * h);
inv_2hh = 1.0 / (2.0 * h * h);
aux     = gamma * (1 - 2.0 * alpha);


/* ---- loop ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     /* compute Sobel derivatives */
     fx  =   (f[i+1][j+1] - f[i-1][j+1]) * inv_8h
           + (f[i+1][j]   - f[i-1][j]  ) * inv_4h
           + (f[i+1][j-1] - f[i-1][j-1]) * inv_8h;
     fy  =   (f[i+1][j+1] - f[i+1][j-1]) * inv_8h
           + (f[i]  [j+1] - f[i]  [j-1]) * inv_4h
           + (f[i-1][j+1] - f[i-1][j-1]) * inv_8h;

     /* compute (slightly regularised) inverse squared gradient */
     inv_grad_sqr = 1.0 / (fx * fx + fy * fy + 1.0e-10);
    
     /* compute important intermediate expressions */ 
     a =   fy * fy * inv_grad_sqr;
     b = - fx * fy * inv_grad_sqr;
     c =   fx * fx * inv_grad_sqr; 
     delta = alpha * (a + c) + aux * fabs(b);

     /* compute second order derivatives in the 4 principal directions */
     f00 = (f[i+1][j]   - 2.0 * f[i][j] + f[i-1][j])   * inv_hh; 
     f22 = (f[i][j+1]   - 2.0 * f[i][j] + f[i][j-1])   * inv_hh;
     f11 = (f[i+1][j+1] - 2.0 * f[i][j] + f[i-1][j-1]) * inv_2hh;   
     f33 = (f[i-1][j+1] - 2.0 * f[i][j] + f[i+1][j-1]) * inv_2hh;   
     
     /* explicit update step for MCM */
     u[i][j]  = f[i][j] + tau * ( (a - delta) * f00 + (delta + b) * f11
                                + (c - delta) * f22 + (delta - b) * f33 );
     }


/* ---- free memory ---- */
 
free_double_matrix (f, nx+2, ny+2);

return;

}  /* mcm_delta */
/*--------------------------------------------------------------------------*/

void mcm_step

     (double   tau,           /* time step size */
      long     nx,            /* image dimension in x direction */
      long     ny,            /* image dimension in y direction */
      double   h,             /* pixel size in x and y direction */
      double   **u,           /* image, changed */
      double   **coeff_x,     /* hermite coefficients in x direction */
      double   **coeff_y,     /* hermite coefficients in y direction */
      long     method)        /* interpolation method */

/*
   Adapted version of original source code by Joachim Weickert (2021)
   for the explicit delta scheme. 
   explicit scheme for MCM;
   approximates u_{\xi \xi} direcltly using interpolated neighbouring
   image values in \xi-direction (direction of the isophote)
*/

{
long     i, j;                /* loop variables */
double   fx, fy;              /* Sobel derivatives */
double   inv_4h, inv_8h;      /* time savers */
double   grad_sqr;            /* time saver */
double   inv_grad;            /* time saver */
double   **f;                 /* image u at old time level */
double   xi_x, xi_y;          /* xi vector */
double   n1;                  /* value of neighbour 1 */
double   n2;                  /* value of neighbour 2 */
double   h_sqr;               /* squared distance to a neighbour */

/* ---- alloc memory ---- */
alloc_double_matrix (&f, nx+2, ny+2);

/* ---- copy u to f ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     f[i][j] = u[i][j];
    
/* ---- reflecting dummy boundaries ---- */

dummies_double (f, nx, ny);

/* ---- interpolate image ---- */
     
interpolate_image (nx, ny, f, coeff_x, coeff_y, method); 

/* ---- compute time savers ---- */

inv_4h  = 1.0 / (4.0 * h);
inv_8h  = 1.0 / (8.0 * h);

/* ---- loop ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     /* compute Sobel derivatives */
     fx  =   (f[i+1][j+1] - f[i-1][j+1]) * inv_8h
           + (f[i+1][j]   - f[i-1][j]  ) * inv_4h
           + (f[i+1][j-1] - f[i-1][j-1]) * inv_8h;
     fy  =   (f[i+1][j+1] - f[i+1][j-1]) * inv_8h
           + (f[i]  [j+1] - f[i]  [j-1]) * inv_4h
           + (f[i-1][j+1] - f[i-1][j-1]) * inv_8h;
     
     /* compute squared gradient */
     grad_sqr = fx * fx + fy * fy;
     
     
     /* no changes if gradient magnitude is zero */
     if (grad_sqr < EPS)
        continue;
     
        
     /* compute inverse gradient magnitude */
     inv_grad = 1.0 / sqrt (grad_sqr);

     
     /* compute normed xi vector (direction of the isophote) */
     xi_x = -fy * inv_grad;
     xi_y = fx * inv_grad;
  
     
     /* get interpolated adjacent image values in xi-direction */
     
     approx_neighbours 
        (&n1, &n2, &h_sqr, i, j, xi_x, xi_y, h, nx, ny, u, 
            coeff_x, coeff_y, method);
     
   
     /* explicit update step for MCM */
     u[i][j] = f[i][j] + tau * (n1 - 2*u[i][j] + n2) / h_sqr;
     
     }
     
free_double_matrix (f, nx+2, ny+2);

return;

}  /* mcm_step */

/*--------------------------------------------------------------------------*/

void mcm

     (double   tau,       /* time step size */
      long     kmax,      /* number of iterations */
      long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      double   h,         /* pixel size in x and y direction */
      double   **u,       /* image, changed */
      long     method)    /* interpolation method */
/*
   iterative explicit scheme for MCM with interpolation;
*/
   
{
long     k;               /* loop variable */
long     p;               /* 1 + degree of the used spline */
double   **coeff_x;       /* hermite coefficients in x direction */
double   **coeff_y;       /* hermite coefficients in y direction */

/* initialize the spline order depending on the used method */  
p = method == 5 ? 6 : 4;

/* ---- allocate memory ---- */
alloc_double_matrix (&coeff_y, nx+2, (ny+2-1) * p);
alloc_double_matrix (&coeff_x, ny+2, (nx+2-1) * p);

/* ---- loop ---- */
for (k=1; k<=kmax; k++) {
   printf ("iteration:   %ld / %ld \n", k, kmax);
   mcm_step (tau, nx, ny, h, u, coeff_x, coeff_y, method);
    
}  
   
/* ---- free memory ---- */
free_double_matrix (coeff_x, ny+2, (nx+2-1) * p);
free_double_matrix (coeff_y, nx+2, (ny+2-1) * p);

return;

}  /* mcm */

/*--------------------------------------------------------------------------*/
/*                 FUNCTIONS FOR EXPERIMENTS / EVALUATION                   */
/*--------------------------------------------------------------------------*/

long* create_random_zeros 

   (long  m,     /* number of zeros */
     long min,   /* minimum zero */
     long max)   /* maxumum zero (min < max) */

/*
   returns a random zero set between min and max 
*/

{

   long i;            // loop variable 
   long n;            // width of the zero set domain  
   long *z;           // zeros of the polynomial
   
   n = max-min;

   /* allocate memory */
   z = malloc(sizeof(long)*(m+1));
   

   /* generate random zeros and store them in z */
   for(i = 1; i <= m; i++) {
      // pseudo-random integer between min and max
      z[i] = (rand() % (n+1)) + min; 
   }

   return z;

}  /* create_random_zeros */

/*--------------------------------------------------------------------------*/

void expand_zero
   
   (long *c,   /* polynomial in expanded form of degree < deg */
    long z,    /* linear factor (x-z) */
    long deg)  /* max degree */

/* 
   subroutine for expand_zeros
   expands the polynomial given by p = (sum_c_i)*(x-z) 
   modifies the coefficients c_i
*/

{
long  i;    /* loop variable */
long  *c_copy; 

/* allocate memory */
c_copy = malloc(sizeof(long)*(deg+1));

/* initialize copy */
for(i = 0; i <= deg; i++) {
   c_copy[i] = c[i];
}

/* p = p-z */
for(i = 0; i < deg; i++) {
   c[i] *= -z;
}

/* p += p*x */
for(i = 1; i <= deg; i++) {
   c[i] += c_copy[i-1];
}

/* free memory */
free(c_copy);

return;

}  /* expand_zero */

/*--------------------------------------------------------------------------*/

long* expand_zeros

   (long *z,    /* zero set of the polynomial */
    long deg)   /* degree of the polynomial */

/*
   returns integer coefficients of a  
   polynomial with zero set domain between min and max 
*/

{

   long i;            // loop variable
   long *c;           // coefficients of the polynomial

   /* allocate memory */
   c = malloc(sizeof(long)*(deg+1));

   /* the polylomial p is given by prod_i (x-z[i]) */
   /* goal: coefficients s.t. p = sum_i c[i]x^i */
   c[0] = 1;
   for(i = 1; i <= deg; i++) {
      expand_zero(c, z[i], deg);
   }

   return c;

}  /* expand_zeros */

/*--------------------------------------------------------------------------*/

double* integrate_poly
   
   (long *c,    /* coefficients */
    long deg)   /* degree */

/*
   integrates the polynomial sum_i c_i x^i and returns the new coefficients
*/

{
   long    i;       /* loop variable */
   double  *c_integ; 

   /* allocate memory */
   c_integ = malloc(sizeof(double)*(deg+2));


   c_integ[0] = 0; /* integration constant */

   
   for(i = 1; i <= deg+1; i++) {
      c_integ[i] = (double)(1.0/i * c[i-1]);
   } 

   return c_integ;

}  /* integrate_poly */

/*--------------------------------------------------------------------------*/

double eval_poly

   (double  x,       /* evaluation point */
    double  *c,      /* coefficients */
    long    deg)     /* degree */
 
 /*
   returns sum_i c_i x^i
 */

{
   long   i;       /* loop variable */ 
   double  res;   /* result */

   res = 0.0;
   for(i = 0; i <= deg; i++) {
      res += c[i] * pow(x,i);
   }

   return res;

}  /* eval_poly */

/*--------------------------------------------------------------------------*/

void print_poly

   (double *c, /* coefficients */
    long   deg)   /* degree */

/*
   prints polynomial sum_i c[i]c^i to console
*/

{
long i;  /* loop variable */

printf("%+lf", c[0]);
for(i = 1; i < deg; i++) {
   if(c[i] != 0) {
      printf("%+lfx^%ld", c[i], i);
   }
}
if(deg > 0) {
   if(c[deg] != 0) {
      printf("%+lfx^%ld", c[deg], deg);
   }
}
printf("\n");

return;

}  /* print_poly */

/*--------------------------------------------------------------------------*/

void sample_poly

   (double  *c,
    double  deg,
    long    start,
    long    end,
    double  *Y)
/*
   samples the polynomial sum_i c_i x^i stores sample breakpoints in Y
   samples at sample distance 0.5 over the interval [start, end]
*/

{

long i, j, k;  /* loop variables */
long   n;      /* sample size for each interval */
double x;      /* sample point x-value*/
n = 2;
k = 0;
for(i = start; i < end; i++) {
   for(j = 0; j < n; j++) {
      x = i + 1.0*j/n;
      Y[k] = eval_poly(x, c, deg);
      k++;
   }
}
Y[k] = eval_poly(end, c, deg);

return;

}  /* sample_poly */

/*--------------------------------------------------------------------------*/

void sample_poly_to_file

   (double *c,     /* coefficients */
    long    deg,   /* degree */
    long    start, /* sample start point */
    long    end)   /* sample end point */

/*
   samples the polynomial sum_i c_i x^i and writes samples to file
   samples at sample distance 0.5 over the interval [start, end]
*/

{

   long   i, j;   /* loop variables */
   long   n;      /* sample size for each interval */
   long   m;      /* sample size for exact data */
   double x;      /* sample point x-value */
   FILE  *f;      /* output file */
   FILE  *g;      /* output file 2 */

   n = 2;
   m = 100;
   f = fopen("sample.dat", "w");
   if (NULL == f) {
      printf("could not open file 'sample.dat' for writing, aborting\n");
      exit(1);
   }
   for(i = start; i < end; i++) {
      for(j = 0; j < n; j++) {
         x = i + 1.0*j/n;
         fprintf(f, "%lf %lf\n", x, eval_poly(x, c, deg));
      }
   }
   fprintf(f, "%lf %lf\n", (double)end, eval_poly((double)end, c, deg));

   fclose(f);

   
   g = fopen("sample_exact.dat", "w");
   if (NULL == g) {
      printf("could not open file 'sample_exact.dat' for writing, aborting\n");
      exit(1);
   }
   for(i = start; i < end; i++) {
      for(j = 0; j < m; j++) {
         x = i + 1.0*j/m;
         fprintf(g, "%lf %lf\n", x, eval_poly(x, c, deg));
      }
   }
   fprintf(g, "%lf %lf\n", (double)end, eval_poly((double)end, c, deg));
   fclose(g);

   printf(
      "output files 'sample.dat' and 'sample_exact.dat' successfully written\n"
   );

   return;

}  /* sample_poly_to_file */

/*--------------------------------------------------------------------------*/

void sample_piecewise_poly_to_file

   (char   *fname,   /* output file name */
    double h,        /* interval length */
    double *Y,       /* breakpoint values */
    double *coeff,   /* Hermite coefficients */
    long   n,        /* number of breakpoitns */
    long   p)        /* spline order */

/*
   samples piecewise polynomial and writes samples Y and interpolant to file
   * equidistant version *
*/
   
{
long i, j;  /* loop variable */
long N;     /* sample rate */
double x;   /* sample point x-value*/
FILE *f;    /* output file */

N = SAMPLE_RATE;

f = fopen(fname, "w");
if (NULL == f) {
   printf("could not open file '%s' for writing, aborting\n", fname);
   exit(1);
}

for(i = 0; i < n-1; i++) {
   for(j = 0; j < N; j++) {
      x = i + j*1.0/N;
      fprintf(f, "%lf %lf\n", x*h, eval_piecewise_poly(x, coeff, p, n));
   }
}

fclose(f);

printf("output file %s successfully written\n", fname);

return; 

}  /* sample_piecewise_poly_to_file */

/*--------------------------------------------------------------------------*/

void specific_poly_to_file

   (long    deg,     /* degree */
    double  *c)      /* coefficients of size degree+1 */

/*
   samples given polynomial at 10 equidistant
   positions in [0,5], interpolates it with the different Hermite
   interpolation methods and writes the results to files for plotting
*/

{
long i;        /* loop variable */
double  *Y;          /* sampled values */
double  *coeff;      /* coefficients of Hermite inerpolant */
char    fname[80];   /* plot file name */

/* allocate memory */
alloc_double_vector(&coeff, 10*6);
alloc_double_vector(&Y, 11);

/* sample polynomial for plotting  */
sample_poly(c, deg, 0, 5, Y);

/* sample polynomial for plotting */
sample_poly_to_file(c, deg, 0, 5);

for(i = 1; i <= 5; i++) {
   /* interpolate */
   interpolate_piecewise_Hermite(Y, 11, i, coeff);

   /* write file name */
   snprintf(fname, sizeof(char) * 32, "interpolant_%ld.dat", i);

   /* sample interpolant for plotting */
   sample_piecewise_poly_to_file(fname, 0.5, Y, coeff, 11, i < 5 ? 4 : 6);
}

/* free memory */
free(Y);
free(coeff);

return;

}  /* specific_poly_to_file */

/*--------------------------------------------------------------------------*/

void random_poly_to_file

   (long deg)  /* degree */

/*
   creates a random polynomial of degree deg and samples it at 10 equidistant
   positions in [0,5], interpolates it with the different Hermite
   interpolation methods and writes the results to files for plotting
*/

{
long i;        /* loop variable */
long    *z;          /* zero set */
long    *z_expand;   /* expanded coefficients */
double  *c;          /* integrated polynomial coefficients */
double  *Y;          /* sampled values */
double  *coeff;      /* coefficients of Hermite inerpolant */
char    fname[80];   /* plot file name */

/* allocate memory */
z = malloc(sizeof(long)*(deg-1));
alloc_double_vector(&coeff, 10*6);
alloc_double_vector(&Y, 11);

/* create integer zero set (defines derivative, thus deg-1) */
z = create_random_zeros(deg-1, 0, 5);

printf("zeros of derivative:\n");
for(i = 1; i <= deg-1; i++) {
   printf("%ld\n", z[i]);
}

/* expand in order to be able to integrate analytic */
z_expand = expand_zeros(z, deg-1);

/* integrate */
c = integrate_poly(z_expand, deg);

/* sample polynomial for plotting  */
sample_poly(c, deg, 0, 5, Y);

/* sample polynomial for plotting */
sample_poly_to_file(c, deg, 0, 5);

for(i = 1; i <= 5; i++) {
   /* interpolate */
   interpolate_piecewise_Hermite(Y, 11, i, coeff);

   /* write file name */
   snprintf(fname, sizeof(char) * 32, "interpolant_%ld.dat", i);

   /* sample interpolant for plotting */
   sample_piecewise_poly_to_file(fname, 0.5, Y, coeff, 11, i < 5 ? 4 : 6);
}

/* free memory */
free(z);
free(Y);
free(coeff);

return;

}  /* random_poly_to_file */

/*--------------------------------------------------------------------------*/

void test_reproduction_order
   
   (long    deg,     /* degree */
    long    method)  /* method */

/*
   test if interpolant of given method accurately reproduces
   polynomials p of degree deg on [0,5] that have extrema in {0,...,5} and 
   have a representation p = int(prod_i=1^{deg-1} (x-z[i]))
*/

{
   long    i, j, k;     /* loop variables */
   long    *z;          /* zero set */
   long    *z_expand;   /* expanded coefficients */
   double  *c;          /* integrated polynomial coefficients */
   long    N;           /* number of tests */
   long    p;           /* spline order of the Hermite interpolant method */
   double  *Y;          /* sampled values */
   double  x, y1, y2;   /* evaluation points */
   int     exact;       /* boolean variable */
   double  *coeff;      /* coefficients of Hermite inerpolant */


   /* initialize the spline order depending on the used method */  
   p = method == 5 ? 6 : 4;

   /* 10 * number of possible zero sets */
   N = 10* (long)pow(5, deg-1); 

   exact = TRUE;

   /* allocate memory */
   z = malloc(sizeof(long)*(deg-1));
   alloc_double_vector(&coeff, 10*p);
   alloc_double_vector(&Y, 11);
   
   for(k = 0; k < N; k++) {
      /* create integer zero set (defines derivative, thus deg-1) */
      z = create_random_zeros(deg-1, 0, 5);

      /* expand in order to be able to integrate analytic */
      z_expand = expand_zeros(z, deg-1);

      /* integrate */
      c = integrate_poly(z_expand, deg);

      /* sample polynomial for plotting  */
      sample_poly(c, deg, 0, 5, Y);

      /* interpolate */
      interpolate_piecewise_Hermite(Y, 11, method, coeff);

      /* evaluate */
      for(i = 0; i < 5; i++) {
         /* check 10 samples for two sample interavals [i,i+1/2], [i+1/2,i+1]*/
         /* results in 6 checks per sample interval */
         /* sufficient to check for exact approx. in cubic and quintic case */
         for(j = 0; j < 10; j++) {
            x = i + j * 1.0/10;
            y1 = eval_poly(x, c, deg);
            y2 = eval_piecewise_poly(2*x, coeff, p, 11);
            if(fabs(y1-y2) > EPS) {
               exact = FALSE;
            }
         }
      }
   }
   if(exact) {
      printf("interpolants are exact.\n");
   }
   else {
      printf("interpolants are not exact.\n");
   }

   /* free memory */
   free(z);
   free(Y);
   free(coeff);

   return;
}  /* test_reproduction_order */

/*--------------------------------------------------------------------------*/

void interpolate_image_row_to_file 

   (char    in[80],  /* image file name (for comments) */
    double  **u,     /* image */
    long    nx,      /* image size in x-direction */
    long    method,  /* method (method_parameter) */
    long    row)     /* image row to interpolate */
    
/*
   interpolates the specified row of the image with given method 
   and writes the sampled interpolant to a file for plotting purposes
*/

{
long     i, j;                /* loop variables */
double   *Y;                  /* sampled image values */
double   *coeff;              /* coefficients of the piecwise Hermite interpolant*/
long     N;                   /* number of samples per interval */
long     p;                   /* spline order of the Hermite interpolant method */
char     comments[1600];      /* string for comments */
char     file_name[32];       /* output file name */
FILE     *f;                  /* file 1 */
FILE     *g;                  /* file 2 */

/* method descriptions */
char description[5][80] = {
   "Linear Spline", 
   "Cubic Spline", 
   "Fritsch-Carlson (1980)", 
   "Extended Two-Sweep (Eisenstat 1985)", 
   "MQSI (Lux 2020)"
};

N =  SAMPLE_RATE;       

/* initialize the spline order depending on the used method */  
p = method == 5 ? 6 : 4;

/* allocate memory */
alloc_double_vector(&Y, nx+2);
alloc_double_vector(&coeff, (nx+1)*p);

/* write file name which depends on the method */
snprintf(file_name, sizeof(char) * 32, "interpolant_%ld.dat", method);

/* initialize row to interpolate */
for (i=0; i<=nx+1; i++) {
   Y[i] = u[row][i];
}

/* interpolate data with given method */
interpolate_piecewise_Hermite(Y, nx+2, method, coeff);

/* write comments */
comments[0] = 0;
comment_line (comments, "# data for plotting purposes\n");
comment_line (comments, "# interpolant of single row of initial pgm image\n");
comment_line (comments, 
"#----------------------------------------------------------------------#\n");
comment_line (comments, "# initial image:  %s\n", in);
comment_line (comments, "# method:         %s\n", description[method-1]);
comment_line (comments, "# row:                          %ld\n", row);
comment_line (comments, 
"#----------------------------------------------------------------------#\n");

/* write data to file */
f = fopen("breakpoints.dat", "w");
if (NULL == f) {
   printf("could not open file 'breakpoints.dat' for writing, aborting\n");
   exit(1);
}
/* write comments */
fprintf (f, "%s", comments);  
for (i = 1; i <= nx; i++) {
   fprintf(f, "%ld %lf\n", i, Y[i]);
}
fclose(f);

g = fopen(file_name, "w");
if (NULL == f) {
   printf("could not open file '%s' for writing, aborting\n", file_name);
   exit(1);
}
fprintf (g, "%s", comments); 
for(i = 1; i < nx; i++) {
   for(j = 0; j < N; j++) {
      fprintf(
         g, 
         "%lf %lf\n", 
         i + j*1.0/N, 
         eval_piecewise_poly(i + j*1.0/N, coeff, p, nx+2)
      );
   }
}
fclose(g);

printf (
   "breakpoints.dat and interpolant_%ld.dat successfully written\n\n", 
   method
);

/* free memory */
free_double_vector(Y, nx+2);
free_double_vector(coeff, (nx+1)*p);

return;

}  /* interpolate_image_row_to_file */

/*--------------------------------------------------------------------------*/


void resize_image
     
     (long     method,          /* interpolation method parameter */
      long     nx,              /* image dimension in x direction */
      long     ny,              /* image dimension in y direction */
      double   **u,             /* input image */
      long     nx_new,           /* new x-dimension */
      long     ny_new,           /* new y-dimension */
      double   ***w)             /* resized image, output */

/*
   resizes the image with the specified interpolation 
   method to match the new dimensions 
*/

{
long     i, j;                       /* loop variables */
double   **coeff_x, **coeff_y;       /* coefficients of interpolants */
double   **v;                        /* horizontally resized image */
double   aux;                        /* auxiliary variable */
long     p;
   
/* initialize the spline order depending on the used method */  
p = method == 5 ? 6 : 4;

/* ---- allocate memory ---- */
alloc_double_matrix (&coeff_x, ny+2, (nx+2-1) * p);
alloc_double_matrix (&coeff_y, nx_new+2, (ny+2-1) * p);
alloc_double_matrix (&v, nx_new+2, ny+2);
alloc_double_matrix (w, nx_new+2, ny_new+2);

interpolate_image_horizontal (nx, ny, u, coeff_x, method);

for (j=0; j<=ny+1; j++) {
   for (i=0; i<=nx_new+1; i++) {
      aux = (double)i/(nx_new+1) * (nx+1);
      v[i][j] = eval_piecewise_poly (aux, coeff_x[j], p, nx+2);
   }
}

interpolate_image_vertical (nx_new, ny, v, coeff_y, method);

for (i=0; i<=nx_new+1; i++) {
   for (j=0; j<=ny_new+1; j++) {
      aux = (double)j/(ny_new+1) * (ny+1);
      (*w)[i][j] = eval_piecewise_poly (aux, coeff_y[i], p, ny+2);
   }
}

/* free memory */
free_double_matrix (coeff_x, ny+2, (nx+2-1) * p);
free_double_matrix (coeff_y, nx_new+2, (ny+2-1) * p);
free_double_matrix (v, nx_new+2, ny+2);

return;
}

/*--------------------------------------------------------------------------*/

void generate_disk_image

   (long    nx,      /* image dimension in x direction */
    long    ny,      /* image dimension in y direction */
    double  sigma,   /* radius of the disk */
    double  **u)     /* image, ouput */

/*
   generate disk shaped image of size nx*ny
*/ 
 
{
long i, j;       /* loop variables */
double mx, my;   /* middle point coordinates */

/* initialize middle point coordinates */
mx = (double)(nx+1) / 2.0;
my = (double)(ny+1) / 2.0;

/* set image values */
for(i = 0; i < nx+2; i++) {
   for(j = 0; j < ny+2; j++) {
      u[i][j] = 0.0;
      if((mx-i)*(mx-i) + (my-j)*(my-j) <= sigma * sigma) {
         u[i][j] = 255.0;
      }
   }
}  

return;  

}  /* generate_disk_image */

/*--------------------------------------------------------------------------*/

void generate_square_image

   (long    nx,      /* image dimension in x direction */
    long    ny,      /* image dimension in y direction */
    double  len,     /* side length of the square */
    double  **u)     /* image, ouput */

/*
   generate disk shaped image of size nx*ny
*/ 
 
{
long i, j;       /* loop variables */
double mx, my;   /* middle point coordinates */

/* initialize middle point coordinates */
mx = (double)(nx+1) / 2.0;
my = (double)(ny+1) / 2.0;

/* set image values */
for(i = 0; i < nx+2; i++) {
   for(j = 0; j < ny+2; j++) {
      u[i][j] = 30.0;
      if(fmax(fabs(mx-i),fabs(my-j)) <= len/2) {
         u[i][j] = 225.0;
      }
   }
}  

return;  

}  /* generate_disk_image */

/*--------------------------------------------------------------------------*/

void diff_image

   (long    nx,   /* image dimension in x direction */
    long    ny,   /* image dimension in y direction */
    double  **u,  /* first image */
    double  **v,  /* second image */
    double  **w)  /* absolute difference of the images, output */

/*
   absolute difference of u and v, stored in w;
*/

{
   long i, j;     /* loop variables */
   for(i = 1; i <= nx; i++) {
      for(j = 1; j <= ny; j++) {
         w[i][j] =  fabs(u[i][j] - v[i][j]);
      }
   }
   
   return; 
   
}  /* diff_image */

/*--------------------------------------------------------------------------*/

void threshold_image
   
     (long   nx,     /* image dimension in x direction */
      long   ny,     /* image dimension in y direction */
      double **u,    /* image, edited */
      double thresh) /* threshold */

/*
   sets image values below given threshold to 0 and 
   sets image values above given threshold to 255
*/

{
long i, j;     /* loop variables */

for(i = 0; i < nx+2; i++) {
   for(j = 0; j < nx+2; j++) {
      if(u[i][j] < thresh) {
         u[i][j] = 0.0;
      }
      else {
         u[i][j] = 255;
      }
   }
}

return;

}  /* threshold_image */

/*--------------------------------------------------------------------------*/

double disk_image_radius

   (double  h,       /* pixel size */
    long    nx,      /* image dimension in x direction */
    long    ny,      /* image dimension in y direction */
    double  **u)     /* disk shaped image */
    
/*
   returns the radius of a centered disk shaped image
*/

{
double radius; /* result */
long i;        /* loop variable */

radius = 0.0;
/* distance between the first pixel in the middle row which is 
   255 and the center of the image */
for(i = 1; i <= nx; i++) {
    if((long)(u[i][(ny+1)/2]) == 255) {
      radius = fabs(h * ((double)(nx+1) / 2.0 - (i + 0.499)));
    }
}

return radius;

}

/*--------------------------------------------------------------------------*/

double mse
   (double  **u,  /* error image */
    long    nx,   /* image size in x direction */
    long    ny)   /* image size in y direction */

/*
   returns the mean squared error of the image u in comparison to zero
*/

{
   long     i, j;     /* loop variables */
   double   mse;      /* mean squared error */
   mse = 0.0;
   for(i = 1; i <= nx; i++) {
      for(j = 1; j <= ny; j++) {
         mse += u[i][j] * u[i][j];
      }
   }
   return mse / (nx*ny);
}  /* mse */

/*--------------------------------------------------------------------------*/

double partial_disk_area
   
   (double sigma,    /* radius of the disk centred at (0,0) */
    double posX,     /* x position of the pixel's centre */
    double posY,     /* y position of the pixel's centre */
    double N)        /* number of samples in one direction */

/*
   returns the approximate area of the disk of radius sigma that covers the 
   pixel of size 1 and distance r from the centre of the disk.
*/

{
   long   i, j;      /* loop variables */
   double startX;    /* start x-value of the pixel as a square */
   double startY;    /* start y-value of the pixel as a square */ 
   double step;      /* step size of the numerical integration  */
   double sampleX;   /* sample x-position of the numerical integration */
   double sampleY;   /* sample y-position of the numerical integration */
   double sigma_sqr; /* squared radius for calculations */

   N = 100;
   step = 1.0/N;
   
   sigma_sqr = sigma*sigma;

   startX = posX - 0.5;
   startY = posY - 0.5;
   
   sampleX = startX;
   sampleY = startY;

   double area = 0.0;

   for(i = 0; i <= N; i++) {
      sampleY = startY;
      sampleX += step;
      for(j = 0; j <= N; j++) {
         sampleY += step;
         if(sampleX*sampleX + sampleY*sampleY <= sigma_sqr) {
            area++;
         }
      }
   }
   area = area/(N*N);

   return area;

}  /* partial_area */

/*--------------------------------------------------------------------------*/

void generate_disk_image_partial_area

   (long    nx,      /* image dimension in x direction */
    long    ny,      /* image dimension in y direction */
    double  sigma,   /* radius of the disk */
    double  **u)     /* image, ouput */

/*
   generate disk shaped image of size nx*ny with partial area effect
*/ 
 
{
long    i, j;              /* loop variables */
double  mx, my;            /* middle point coordinates */
double  sigma_sqr_lower;   /* (sigma-1)^2 */
double  sigma_sqr_upper;   /* (sigma+1)^2 */

/* initialize middle point coordinates */
mx = (double)(nx+1) / 2.0;
my = (double)(ny+1) / 2.0;

sigma_sqr_lower = (sigma-1)*(sigma-1);
sigma_sqr_upper = (sigma+1)*(sigma+1);

/* set image values */
for(i = 0; i < nx+2; i++) {
   for(j = 0; j < ny+2; j++) {
      u[i][j] = 0.0;
      if((mx-i)*(mx-i) + (my-j)*(my-j) <= sigma_sqr_lower) {
         u[i][j] = 255.0;
      }
      /* partial area effect kicks in */
      else if((mx-i)*(mx-i) + (my-j)*(my-j) <= sigma_sqr_upper) {
         u[i][j] = 255.0 * partial_disk_area(sigma, mx-i, my-j, 1000);
      }
      /* outside of the range of the partial area effect */
      else  {
         u[i][j] = 0.0;
      }
   }
}  

return;  

}  /* generate_disk_image_partial_area */

/*--------------------------------------------------------------------------*/

void disk_shrinkage

     (double   tau,       /* time step size */
      long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      double   h,         /* pixel size in x and y direction */
      double   sigma)     /* radius of the disk */
/*
   iterative explicit scheme for MCM with interpolation on a disk shaped image;
   evaluation of the scheme with different interpolation methods
   analytic extiction time of the disk is given by 1/2*sigma^2 where sigma is
   the radius of the disk
*/
   
{
long     i, j;            /* loop variables */
double   **u;             /* image */
double   **v;             /* copy of the image u (remains constant) */;
double   **disk;          /* shrinking disk image solution */
double   **diff;          /* difference to numerical solution */
double   **coeff_x;       /* hermite coefficients in x direction */
double   **coeff_y;       /* hermite coefficients in y direction */
double   max, min;        /* largest, smallest grey value */
double   mean;            /* average grey value */
double   std;             /* standard deviation */
double   cutoff;          /* zero threshold value */
long     method;          /* interpolation method parameter */
long     k[6];            /* iteration indices for different methods */
double   t;               /* auxiliary variable for current time step */
double   t_ext;           /* auxiliary variable, extinction time */
double   t_ext_appr;      /* auxiliary variable, approximated extinction time */
long     sample_dist;     /* distance between two sampled time steps */
long     sample_rate;     /* sample rate for plots */
long     sample_count;    /* number of samples positions */
double   area;
double   slope;           /* estimated slope of the line A(t) = disk area */
double   slope_error;     /* area between estimated and analytic solution */
double   sum_t_sqr;       /* sum_{i=1}^N t_i^2 */
double   error;           /* mean squared error of approx. solution */

/* ---- allocate memory ---- */

alloc_double_matrix (&v, nx+2, ny+2);
alloc_double_matrix (&u, nx+2, ny+2);
alloc_double_matrix (&disk, nx+2, ny+2);
alloc_double_matrix (&diff, nx+2, ny+2);
alloc_double_matrix (&coeff_y, nx+2, (ny+2-1) * 6);
alloc_double_matrix (&coeff_x, ny+2, (nx+2-1) * 6);

/* initialize disk image with partial area effect */
generate_disk_image_partial_area(nx, ny, sigma, u);

/* extinction time 1/2 sigma^2 h^2 */
t_ext = 0.5 * sigma * sigma * h * h;

/* initialize sample rate */
sample_rate = sigma*sigma/64;  

/* copy the image u to v */
for(i = 0; i < nx+2; i++) {
   for(j = 0; j < nx+2; j++) {
      v[i][j] = u[i][j];
   }
}

/* initialize cutoff */
cutoff = 0.499;

printf("simga = %3.1lf\n", sigma);

/* ---- loop through interpolation methods ---- */
for(method = 0; method <= 5; method++) { 
   if(method == 1)
      continue;
   
   /* fresh counter for snapshots */
   //sc = 0;
   /* fresh copy of the original image */
   for(i = 0; i < nx+2; i++) {
      for(j = 0; j < nx+2; j++) {
         u[i][j] = v[i][j];
      }
   }

   max = 255.0;      /* maximal grey value */
   k[method] = 0;    /* start with method 0 */
   slope = 0;        /* set slope of A(t) to zero */
   sample_count = 0; /* set sample count to zero */
   error = 0.0;      /* set error to zero */
   slope_error = 0.0;/* set slope error to zero */
   t = 0.0;          /* initialize time */
   sum_t_sqr = 0.0;  /* initialize squared time steps sum */
   /* set sample distance (1/tau for k-units) */
   sample_dist = sample_rate/tau;
   
   
   /* run evolution with current method until image vanished 
   (or ext. time is reached) */
   while((max > cutoff) && (t <= t_ext)){

      /* increment step counter for current method */
      k[method]++;

      /* update current time value */
      t = k[method]*tau;
      
      /* --- write plot sample data: area of shrunk disk --- */
      if((k[method] % sample_dist) == 0) { 
         /* get copy of the original image for editing */
         for(i = 0; i < nx+2; i++) {
            for(j = 0; j < nx+2; j++) {
               disk[i][j] = u[i][j];
            }
         }
         /* get area as sum of grey values divided by max value */
         area = 0.0;
         for(i = 1; i <= nx; i++) {
            for(j = 1; j <= ny; j++) {
               area += u[i][j];
            }
         }
         area = area / 255.0;
         /* update slope */
         slope += t * (area - PI * sigma * sigma);
         sum_t_sqr += t*t;
         sample_count++;
         /*analytic solution with partial area effect*/
         generate_disk_image_partial_area(nx, ny, sqrt(sigma*sigma-2*t), disk);
         diff_image(nx, ny, u, disk, diff);  
         /* update mse error */
         error += mse(diff, nx, ny);
      }      
 
      /* hard break condition */
      if(t > 10 * t_ext) {
         /* factor 10 above theoretical extiction time, abort */
         printf(
            "theoretical extiction time: %5.3lf\n",
            t_ext
         );
         printf("exceeded extinction time by factor 10, aborting\n");
         break;
      }

      /* --- apply mcm method step --- */
      if(method == 0) {
         /* delta stencil */
         mcm_delta (tau, 0.0, 0.0, nx, ny, h, u);
      }
      else {
         /* interpolation scheme */
         mcm_step (tau, nx, ny, h, u, coeff_x, coeff_y, method);
      }

      /* obtain max image value for stop condition */
      analyse_grey_double (u, nx, ny, &min, &max, &mean, &std); 
      

   }

   /* finalize slope estimation */
   slope /= sum_t_sqr;
   //slope_error = 0.5*fabs(slope+2*PI) * t_ext * t_ext;#
   t_ext_appr = -PI*sigma*sigma/slope;
   slope_error = 100*fabs(t_ext-t_ext_appr)/t_ext;
   /* finalize mse error estimation */
   error /= sample_count;
   /* console output */
   printf(
      "%ld %10.3lf %10.3lf %10.3lf\n", 
      method,
      slope,
      slope_error, 
      //fabs(slope + 2*PI)/(2*PI)*100.0,
      error
   );

}

/* ---- free memory ---- */
free_double_matrix (u, nx+2, ny+2);
free_double_matrix (v, nx+2, ny+2);
free_double_matrix (diff, nx+2, ny+2);
free_double_matrix (disk, nx+2, ny+2);

free_double_matrix (coeff_x, ny+2, (nx+2-1) * 6);
free_double_matrix (coeff_y, nx+2, (ny+2-1) * 6);

return;

}  /* disk_shrinkage */

/*--------------------------------------------------------------------------*/
/*                           MAIN MENU FUNCTIONS                            */
/*--------------------------------------------------------------------------*/

void print_interpolation_methods() 

/*
   print interpolation methods table 
*/ 

{
int i;              /* loop variable */
char divider[81];   /* string for list divider */

/* initialize list divider */
for(i = 0; i < 80; i++) {
   divider[i] = '-';
}
divider[80] = '\0';

/* print table */
printf("%s\n", divider);
printf("interpolation method                    ");
printf("monotone  order of   spline   smoothness\n");
printf("                                                  "); 
printf("accuracy   degree   degree\n");
printf("%s\n", divider);
printf("1: Linear Spline                        ");         
printf("yes       2          1        0\n");
printf("%s\n", divider);
printf("2: Cubic Spline                         ");
printf("no        4          3        2\n");
printf("%s\n", divider);
printf("3: Fritsch-Carlson (1980)               ");
printf("yes       3          3        1\n");
printf("%s\n", divider);
printf("4: Ext. Two-Sweep (Eisenstat 1985)      ");
printf("yes       4          3        1\n");
printf("%s\n", divider);
printf("5: MQSI (Lux 2020)                      ");
printf("yes       3          5        2\n");
printf("%s\n", divider);

return;

}  /* print_interpolation_methods */

/*--------------------------------------------------------------------------*/

int main (int argc, char* argv[])

{
double  **u;                  /* image */
long    nx, ny;               /* image size in x, y direction */ 
double  h;                    /* pixel size in x and y direction */
double  max, min;             /* largest, smallest grey value */
double  mean;                 /* average grey value */
double  std;                  /* standard deviation */
char    comments[1600];       /* string for comments */


/* main menu parameters */
char    in[80];               /* for reading data */
char    out[128];             /* for writing data */
double  tau;                  /* time step size */
long    kmax;                 /* number of iterations */
long    method;               /* inteprolation method parameter */

/* testing menu parameters (testing mode) */
long    i, j, k;              /* loop variables */
long    test_method;          /* test method parameter */
double  sigma;                /* raduis of disk */
double  len;                  /* side length of square */
long    nx_new, ny_new;       /* new image dimensions */
double  **w;                  /* image copy */
long    row;                  /* image row */
long    deg;                  /* degree of the polynomials */
double  *Y;                   /* sample values */
double  *coeff;               /* Hermite coefficients */
long    t[8];                 /* snapshot time steps */
double  coeffs[9];            /* coefficients */

srand(time(NULL));            /* rng initialization */

/* method descriptions */
char description[5][80] = {
   "Linear Spline", 
   "Cubic Spline", 
   "Fritsch-Carlson (1980)", 
   "Extended Two-Sweep (Eisenstat 1985)", 
   "MQSI (Lux 2020)"
};

/* parsing arguments */
int c;            /* variable for character arguments */
int tflag = 0;    /* flag for argment 't' for testing */
//int index;        /* index variable for looping through arguments */

opterr = 0;

while((c = getopt (argc, argv, "t")) != -1) {
   switch(c) {
   case 't':
      tflag = 1;
      break;
   case '?':
      if (isprint (optopt)) {
         fprintf (stderr, "Unknown option `-%c'.\n", optopt);
      }
      else {
         fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
      }
      return 1;
   default:
      abort ();
   }
}

printf("\nMEAN CURVATURE MOTION FOR GREYSCALE (PGM) IMAGES\n");
printf("EXPLICIT SCHEMES WITH (MONOTONE) SPLINE INTERPOLATION\n\n");
printf ("**************************************************\n\n");
printf ("    Bachelor studies by Moritz van Recum (2023)   \n\n");
printf ("    Supervisor: Joachim Weickert                  \n\n");
printf ("    Copyright 2021 by Joachim Weickert            \n");
printf ("    Dept. of Mathematics and Computer Science     \n");
printf ("    Saarland University, Saarbruecken, Germany    \n\n");
printf ("    All rights reserved. Unauthorised usage,      \n");
printf ("    copying, hiring, and selling prohibited.      \n\n");
printf ("**************************************************\n\n");

/* initialize pixel size h */
h = 1.0;

/* testing/evaluation mode */
if(tflag) {
   printf("EXPERIMENTS/EVALUATION\n\n");
   printf("1: disk shrinkage\n");
   printf("2: generate disk image\n");
   printf("3: generate square image\n");
   printf("4: resize image\n");
   printf("5: interpolate image row to file\n");
   printf("6: test accurate polynomial reproduction order\n");
   printf("7: plot interpolated polynomials for all methods\n");
   printf("8: mcm image sequence at 7 time steps (t=900 to t=3000)\n");
   printf("9: polynomial to file for plotting");
   printf("10: generate disk image with partial area effect\n");
   printf("\n");
   printf("test method number:          ");
   
   /* read test method parameter */
   read_long(&test_method);
   printf("\n");
  
   switch(test_method) {
   case 1:
      /* --------------------------------------------------*/
      /*                 disk shrinkage                    */
      
      /* -- read input image (pgm format P5) -- */
      

      /* read time step size */
      printf ("time step size tau (<%5.3lf):               ", h*h/2.0);
      read_double (&tau);

      printf("0: Delta scheme\n");
      for(method = 2; method <= 5; method++) {
         printf("%ld: %s\n", method, description[method-1]);
      }
      printf("method, slope, area between lines, avg mse of samples\n");
      for(sigma = 16; sigma <= 128; sigma += 16) {
         nx = 4*sigma;
         ny = 4*sigma;
         disk_shrinkage (tau, nx, nx, h, sigma);
      }      
            
      break; 
      /* --------------------------------------------------*/
      
   case 2:
      /* --------------------------------------------------*/
      /*          generate disk shaped image               */
      
      /* read disk radius */
      printf("disk radius:                               ");
      read_double (&sigma);
      
      /* read image dimensions */
      printf("image dimensions (nx=ny)                   ");
      read_long(&nx);
      ny = nx;
      
      /* read output file name */
      printf ("output image (pgm):                        ");
      read_string (out);
      
      /* allocate memory */
      alloc_double_matrix (&u, nx+2, ny+2);
      
      /* generate image of size nx*ny with disk of radius sigma */
      generate_disk_image (nx, ny, sigma, u);
      
      /* ---- write output image (pgm format) ---- */

      /* write parameter values in comment string */
      comments[0] = '\0';
      comment_line (comments, "# disk\n");
      comment_line (comments, "# radius: %5.3lf\n", sigma); 
      
      

      /* write image */
      write_double_to_pgm (u, nx, ny, out, comments);
      printf ("output image %s successfully written\n\n", out);
     
      /* free memory */
      free_double_matrix (u, nx+2, ny+2);
      
      break;
      /* --------------------------------------------------*/

   case 3:
      /* --------------------------------------------------*/
      /*          generate square shaped image               */
      
      /* read disk radius */
      printf("square side length:                        ");
      read_double (&len);
      
      /* read image dimensions */
      printf("image dimensions (nx=ny)                   ");
      read_long(&nx);
      ny = nx;
      
      /* read output file name */
      printf ("output image (pgm):                        ");
      read_string (out);
      
      /* allocate memory */
      alloc_double_matrix (&u, nx+2, ny+2);
      
      /* generate image of size nx*ny with disk of radius sigma */
      generate_square_image (nx, ny, len, u);
      
      /* ---- write output image (pgm format) ---- */

      /* write parameter values in comment string */
      comments[0] = '\0';
      comment_line (comments, "# square\n");
      comment_line (comments, "# side length: %5.3lf\n", len); 
      
   
      /* write image */
      write_double_to_pgm (u, nx, ny, out, comments);
      printf ("output image %s successfully written\n\n", out);
     
      /* free memory */
      free_double_matrix (u, nx+2, ny+2);
      
      break;
      /* --------------------------------------------------*/
      

   case 4:
      /* --------------------------------------------------*/
      /*                 resize image                      */

      printf ("input image (pgm):           ");
      read_string (in);

      printf ("output image (pgm):          ");
      read_string (out);

      /* print interpolation methods table */
      print_interpolation_methods();

      printf ("method (1-5):                ");
      read_long (&method);

      printf ("new x-dimension:             ");
      read_long (&nx_new);

      printf ("new y-dimension:             ");
      read_long (&ny_new);

      read_pgm_to_double (in, &nx, &ny, &u);

      /* reflecting dummy boundaries */
      dummies_double (u, nx, ny);

      resize_image (method, nx, ny, u, nx_new, ny_new, &w);

      /* write image */
      comments[0] = '\0';
      comment_line(comments, "#resize image\n");
      comment_line(comments, "#original image: %s\n", in);
      write_double_to_pgm (w, nx_new, ny_new, out, comments);
      printf ("output image %s successfully written\n\n", out);

      /* free memory */
      free_double_matrix (u, nx+2, ny+2);
      free_double_matrix (w, nx_new+2, ny_new+2);

      break;

      /* --------------------------------------------------*/

   case 5:
      /* --------------------------------------------------*/
      /*       interpolate single row from a pgm image     */
   

      /* read input image (pgm format P5) */
      printf ("input image (pgm):\n");
      read_string (in);
      printf ("\n");
      read_pgm_to_double (in, &nx, &ny, &u);
         
      /*  read row parameter */
      printf("row to interpolate\n");
      printf("(integer between 1 and %ld)\n", ny);
      read_long (&row);
      if( (row < 1) || (row > ny))
         {
         printf("inavid parameter, aborting\n");
         exit(1);
         }
      printf("\n");
      
      /* print interpolation methods table */
      print_interpolation_methods();

      printf ("method (1-5):                ");
      read_long (&method);
      
      /* reflecting dummy boundaries */
      dummies_double (u, nx, ny);

      /* interpolate and write breakpoints and sample values to file */
      interpolate_image_row_to_file (in, u, nx, method, row);


      /* free memory */
      free_double_matrix (u, nx+2, ny+2);
      
      break;
      /* --------------------------------------------------*/

   case 6:
      /* --------------------------------------------------*/
      /*       polynomial accurate reproduction test       */

      /* print interpolation methods table */
      print_interpolation_methods();

      printf ("method (1-5):                ");
      read_long (&method);

      printf("polynomial degree (>=2):     ");
      read_long (&deg);

      /* allocate memory */
      alloc_double_vector(&Y, 11);
      alloc_double_vector(&coeff, 10*6);
   
      test_reproduction_order(deg, method);

      /* free memory */
      free(Y);
      free(coeff);

      break;

      /* --------------------------------------------------*/

   case 7: 
      /* --------------------------------------------------*/  
      /*   plot interpolated polynomials for all methods   */

      printf("polynomial degree (>=2):     ");
      read_long (&deg);

      random_poly_to_file(deg);


      break;
      /* --------------------------------------------------*/

   case 8:
      /* mcm image sequence at 7 time steps (t=900 to t=3000) */

      /* initialize time steps */
      t[0] = 0;
      t[1] = 900;
      t[2] = 3100;
      t[3] = 5300;
      t[4] = 7800;
      t[5] = 13000;
      t[6] = 22000;
      t[7] = 30000;

      /* -- read input image (pgm format P5) -- */
      printf ("input image (pgm):                         ");
      read_string (in);

      read_pgm_to_double (in, &nx, &ny, &u);

      /* alloc memory */
      alloc_double_matrix(&w, nx+2, ny+2);

      printf ("\n");

      /* ---- read parameters ---- */

      /* read time step size */
      printf ("time step size tau (<%5.3lf):               ", h*h/2.0);
      read_double (&tau);

      /* copy image */
      for (i=1; i<=nx; i++) {
          for (j=1; j<=ny; j++) {
              w[i][j] = u[i][j];
          }
      }

      for(method = 1; method <= 5; method++) {
         /* reset image to original image */
         for (i=1; i<=nx; i++) {
             for (j=1; j<=ny; j++) {
                 u[i][j] = w[i][j];
             }
         }
         /* reflecting dummy boundaries */
         dummies_double (u, nx, ny);

         for(k = 1; k <= 7; k++) {
            /* iterations until the next snapshot time step */
            kmax = (t[k]-t[k-1])/tau;

            /* print output file name */
            snprintf(out, sizeof(char) * 128, "m%ld-t%ld-%s",method, t[k], in);

            /* ---- process image ---- */

            /* process image */
            mcm(tau, kmax, nx, ny, h, u, method);

            /* analyse processed image */   
            analyse_grey_double (u, nx, ny, &min, &max, &mean, &std);

            /* ---- write output image (pgm format) ---- */

            /* write parameter values in comment string */
            comments[0] = '\0';
            comment_line (comments, "# MCM, explicit scheme with spline interpolation\n");
            comment_line (comments, "# initial image:        %s\n", in);
            comment_line (comments, "# interpolation method: %s\n", description[method-1]);
            comment_line (comments, "# tau:       %8.2lf\n", tau);
            comment_line (comments, "# kmax:      %8ld\n",   (long)(t[k]/tau));
            comment_line (comments, "# min:       %8.2lf\n", min);
            comment_line (comments, "# max:       %8.2lf\n", max);
            comment_line (comments, "# mean:      %8.2lf\n", mean);
            comment_line (comments, "# st. dev.:  %8.2lf\n", std);
          
            /* write image */
            write_double_to_pgm (u, nx, ny, out, comments);
            printf ("output image %s successfully written\n\n", out);
         }
      }

      /* ---- free memory ---- */
      free_double_matrix (u, nx+2, ny+2);
      free_double_matrix (w, nx+2, ny+2);

      break;
      /* --------------------------------------------------*/

   case 9: 
      /* print specific sampled polynomial to file for plotting */
      
      /* deg 2: int (x-2) */
      coeffs[0] = 0.0;
      coeffs[1] = -2.0;
      coeffs[2] = 0.5;
      /* deg 3: int (x-1)(x-2) */
      coeffs[0] = 0.0;
      coeffs[1] = 2.0;
      coeffs[2] = -3.0/2.0;
      coeffs[3] = 1.0/3.0;
      /* deg 4: int (x-1)(x-2)(x-3) */
      coeffs[0] = 0.0;
      coeffs[1] = -6.0;
      coeffs[2] = 11.0/2;
      coeffs[3] = -2.0;
      coeffs[4] = 1.0/4.0;
      /* deg 7: int (x-1)(x-2)(x-3)(x-4)(x-5)(x) */
      coeffs[0] = 0.0;
      coeffs[1] = 0.0;
      coeffs[2] = -60.0;
      coeffs[3] = 274.0/3;
      coeffs[4] = -225.0/4.0;
      coeffs[5] = 17.0;
      coeffs[6] = -5.0/2;
      coeffs[7] = 1.0/7;

      specific_poly_to_file(7,coeffs);

      break;
   /* --------------------------------------------------*/

   case 10:
      /* partial area test */

      nx = 256;
      ny = 256;

      printf("output image:                            ");
      read_string(out);

      alloc_double_matrix(&u, nx+2, ny+2);
      generate_disk_image_partial_area(nx, ny, 64, u);

      comments[0] = '\0';
      write_double_to_pgm (u, nx, ny, out, comments);
      printf ("output image %s successfully written\n\n", out);
      
      free_double_matrix(u, nx, ny);
      

      break;

   default: 
      printf("invalid test method number, aborting\n"); 
      exit(1);
      /* --------------------------------------------------*/
   }
   return 0;
}


/* --------------------------------------------------*/
/*                 main program                      */

/* -- read input image (pgm format P5) -- */
printf ("input image (pgm):                         ");
read_string (in);

read_pgm_to_double (in, &nx, &ny, &u);
printf ("\n");

/* print interpolation methods table */
print_interpolation_methods();

/* ---- read parameters ---- */
/* read interpolation method */
printf ("\n");
printf ("interpolation method (1-5):                ");
read_long (&method);

/* read time step size */
printf ("time step size tau (<%5.3lf):               ", h*h/2.0);
read_double (&tau);

/* read number of iterations */
printf ("number of iterations (>0):                 ");
read_long (&kmax);

/* read output file name */
printf ("output image (pgm):                        ");
read_string (out);

printf ("\n");


/* ---- analyse initial image ---- */

analyse_grey_double (u, nx, ny, &min, &max, &mean, &std);


/* ---- process image ---- */

/* reflecting dummy boundaries */
dummies_double (u, nx, ny);

/* process image */
mcm(tau, kmax, nx, ny, h, u, method);

printf ("\ninitial image\n");
printf ("minimum:          %8.2lf \n", min);
printf ("maximum:          %8.2lf \n", max);
printf ("mean:             %8.2lf \n", mean);
printf ("standard dev.:    %8.2lf \n\n", std);

/* analyse processed image */   
analyse_grey_double (u, nx, ny, &min, &max, &mean, &std);
printf ("processed image:\n");
printf ("minimum:          %8.2lf \n", min);
printf ("maximum:          %8.2lf \n", max);
printf ("mean:             %8.2lf \n", mean);
printf ("standard dev.:    %8.2lf \n\n", std);


/* ---- write output image (pgm format) ---- */

/* write parameter values in comment string */
comments[0] = '\0';
comment_line (comments, "# MCM, explicit scheme with spline interpolation\n");
comment_line (comments, "# initial image:        %s\n", in);
comment_line (comments, "# interpolation method: %s\n", description[method-1]);
comment_line (comments, "# tau:       %8.2lf\n", tau);
comment_line (comments, "# kmax:      %8ld\n", kmax);
comment_line (comments, "# min:       %8.2lf\n", min);
comment_line (comments, "# max:       %8.2lf\n", max);
comment_line (comments, "# mean:      %8.2lf\n", mean);
comment_line (comments, "# st. dev.:  %8.2lf\n", std);


/* write image */
write_double_to_pgm (u, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);

/* ---- free memory ---- */
free_double_matrix (u, nx+2, ny+2);

return(0);

}  /* main */
