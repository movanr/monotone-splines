#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>
#include <time.h>


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


#define EPS 1e-14          /* small constant */
#define FALSE 0
#define TRUE 1


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
   (coeff[i][3] + coeff[i][2] z + coeff[i][1] z² + coeff[i][0] z³)
*/

{
   
/* evaluation at z using Horner's method */
return ((coeff[0] * z + coeff[1]) * z + coeff[2]) * z + coeff[3];
         
}  /* eval_cubic_hermite_poly */

/*-------------------------------------------------------------------------*/

void hermite_coeffs

   (double  h,          /* interval length */
    double  *Y,         /* breakpoint values */
    double  *d,         /* breakpoint approximated derivatives */
    double  *coeff,     /* coefficient vector of size (n-1)*4 */
    long    n)          /* number of breakpoints */
    
/*
   calculates coefficients of the piecewise cubic hermite for each interval
   and stores them in the coefficient matrix;
   uses the standard basis {x^3, x^2, x, 1};
   the polynomials for each interval all have the domain [0,1]
   derivative values d_i are implicitly mapped to h*d_i due to the 
   mapping of the intervals to [0,1] and the chain rule
   * equidistant version *
*/  

{
long  i;    /* loop variable */

for (i=0; i<n-1; i++)
   {
   coeff[4*i] = 2*Y[i] + h*d[i] - 2*Y[i+1] + h*d[i+1];
   coeff[4*i+1] = -3*Y[i] + 3*Y[i+1] - 2*h*d[i] - h*d[i+1];
   coeff[4*i+2] = h*d[i];
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

void cubic_spline_derivs_eq

   (double  h,          /* interval length */ 
    double  *Y,         /* breakpoint values */      
    double  *s,         /* slopes (output) */
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
alloc_double_vector(&a, n);
alloc_double_vector(&b, n);
alloc_double_vector(&c, n-1);
alloc_double_vector(&d, n);

/* calculate entries of the tridiagonal system */
/* ("not-a-knot" boundary conditions) */

b[0] = 2;
c[0] = 4;
d[0] = (4*Y[1] - 5*Y[0] + Y[2])/h;
for (i=1; i<n-1; i++)   
   {        
   a[i] = 1;
   b[i] = 4;
   c[i] = 1;
   d[i] = 3.0*(Y[i+1] - Y[i-1])/h;
   }
a[n-1] = -4;
b[n-1] = -2;
d[n-1] = (4*Y[n-2] + Y[n-3] - 5*Y[n-1])/h;


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

/*--------------------------------------------------------------------------*/
/*             functions for cubic monotonicity constraints                 */
/*                                                                          */
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
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
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
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

double eval_quintic_hermite_poly

    (double z,      /* evaluation point (in [0,1])*/
     double *coeff) /* coefficients of size 6 */

/*
   evaluates a hermite polynomial at point z
   (coeff[i][5] + coeff[i][4] z + coeff[i][3] z² + coeff[i][2] z³ + ...)
*/

{

    /* evaluation at z using Horner's method */
    return ((((coeff[0] * z + coeff[1]) * z + coeff[2]) * z 
            + coeff[3]) * z + coeff[4]) * z + coeff[5];

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
            -6.0 * f[i] - 3.0 * h * fx[i] - 0.5 * h * h * fxx[i] 
            + 0.5 * h * h * fxx[i + 1] - 3.0 * h * fx[i + 1] 
            + 6.0 * f[i + 1];
        coeff[6 * i + 1] =
            +15.0 * f[i] + 8.0 * h * fx[i] + 3.0 / 2 * h * h * fxx[i] 
            - 1.0 * h * h * fxx[i + 1] + 7.0 * h * fx[i + 1] - 15 * f[i + 1];
        coeff[6 * i + 2] =
            -10.0 * f[i] - 6.0 * h * fx[i] - 3.0 / 2 * h * h * fxx[i] 
            + 0.5 * h * h * fxx[i + 1] - 4.0 * h * fx[i + 1] + 10.0 * f[i + 1];
        coeff[6 * i + 3] = 0.5 * h * h * fxx[i];
        coeff[6 * i + 4] = h * fx[i];
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

    Original code:
    Algorithm XXXX: MQSI – Monotone Quintic Spline Interpolation
    Thomas C.H. Lux, Layne T. Watson, William Thacker, Tyler H. Chang. 
    ACM Transactions on Mathematical Software. Submitted August, 2020
    https://tchlux.github.io/
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

    (long index,  /* index to check for monotonicity*/
     double h,    /* grid size */
     double *f,   /* discrete function values */
     double *fx,  /* discrete derivatives */
     double *fxx) /* discrete second derivatives*/

{

/*
    function to make call of is_monotone easier
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

void mqsi 
    
    (double h,              /* grid size */
     double *f,             /* discrete function values to interpolate */
     double *coeff,         /* quintic hermite spline coefficients */
     long   n)              /* number of function values */

/*
    Original code and comments:
    Algorithm XXXX: MQSI – Monotone Quintic Spline Interpolation
    Thomas C.H. Lux, Layne T. Watson, William Thacker, Tyler H. Chang. 
    ACM Transactions on Mathematical Software. Submitted August, 2020
    https://tchlux.github.io/
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
                A = (f[i]-f[i-1])/(h*h);
                B = 0.0;
            }
            /* Otherwise use the standard quadratic on the left */
            else {
                A = (f[i-2] + f[i] - 2*f[i-1])/(2*h*h);
                B = (f[i] - f[i-2])/(2*h);
            }
        
            dx = 2.0 * A * h + B;
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
                A = (f[i-1] + f[i+1] - 2*f[i])/(2*h*h);
                B = (f[i+1] - f[i-1])/(2*h);
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

/*---------------------------------------------------------------------------*/
/*         Algorithm 3: Identify viable piecewise monotone                   */
/*        derivative values by doing a quasi-bisection search                */
/*---------------------------------------------------------------------------*/

exit(1);
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
    && (!(is_monotone_index(i, h, f, fx, fxx)))) {
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
        if(!is_monotone_index(i, h, f, fx, fxx)) {
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

/* calculate quintic spline representation in terms of hermite coefficients */
quintic_hermite_coeffs(h, f, fx, fxx, coeff, n);

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

   (double  h,
    double  *Y,         /* y-values */ 
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
   
   (double  h,          /* interval length */ 
    double  *Y,         /* breakpoint values */ 
    double  *coeff,     /* hermite spline coefficients (output) */
    long    n)          /* number of (x,y) values */

/*
   interpolation with cubic spline derivatives;
   calculates the hermite coefficients; 
   nonmonotone, order 4 accurate
   * equidistant version *
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
   delta[i] = (Y[i+1]-Y[i]) / h;

/* calculate approximate derivatives */
cubic_spline_derivs_eq (h, Y, d, n); 

/* calculate cubic hermite coefficients */
hermite_coeffs(h, Y, d, coeff, n);

/* free variables */
free_double_vector(d,n);
free_double_vector(delta,n-1);

}  /* cubic_spline */         

/*--------------------------------------------------------------------------*/

void mpci_basic
   
   (double  h,          /* interval length */ 
    double  *Y,         /* breakpoint values */ 
    double  *coeff,     /* hermite spline coefficients (output) */
    long    n)          /* number of (x,y) values */
/*
   monotone piecewise cubic interpolation with Fritsch-Carlson method;
   calculates the hermite coefficients; 
   monotone, order 3 accurate
   * equidistant version *
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
   delta[i] = (Y[i+1]-Y[i]) / h;

/* STEP 1: calculate approximate derivatives */
cubic_spline_derivs_eq (h, Y, d, n); 

/* STEP 2: ensure that d[i] have the right sign */
enforce_sign(d, delta, n);

/* STEP 3: modify derivatives to preserve monotonicity */
modify_slopes_basic(d, delta, n);

/* calculate cubic hermite coefficients */
hermite_coeffs(h, Y, d, coeff, n);

/* free variables */
free_double_vector(d,n);
free_double_vector(delta,n-1);

}  /* mpci_basic */

/*--------------------------------------------------------------------------*/

void mpci_extended_sweep
   
   (double  h,          /* interval length */ 
    double  *Y,         /* breakpoint values */ 
    double  *coeff,     /* hermite spline coefficients (output) */
    long    n)          /* number of (x,y) values */

/*
   monotone piecewise cubic interpolation with extended sweep method;
   calculates the hermite coefficients;
   monotone, order 4 accurate
   * equidistant version *
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
   delta[i] = (Y[i+1]-Y[i]) / h;

/* STEP 1: calculate approximate derivatives */
cubic_spline_derivs_eq (h, Y, d, n); 
   
/* STEP 2: ensure that d[i] have the right sign */
enforce_sign(d, delta, n);

/* STEP 3: modify derivatives to preserve monotonicity */
extended_forward_sweep (d, delta, n);
extended_backward_sweep (d, delta, n);

/* calculate cubic hermite coefficients */
hermite_coeffs(h, Y, d, coeff, n);

/* free variables */
free_double_vector(d,n);
free_double_vector(delta,n-1);

}  /* mpci_extended_sweep */


/*--------------------------------------------------------------------------*/
/*                      image interpolation                                 */
/*--------------------------------------------------------------------------*/

void interpolate_image_horizontal
     
     (long     nx,              /* image dimension in x direction */
      long     ny,              /* image dimension in y direction */
      double   h,               /* pixel size in x and y direction */
      double   **u,             /* input image */
      double   **coeff,         /* coefficients in x direction, output */
      long     method)

/*
   interpolates image in x and y direction along pixels;
   returns coefficients of piecewise cubic polynomial interpolant
   * equidistant version *
*/

{
long     i, j;      /* loop variables */
double   **Y;       /* values in x-direction (size ny * nx) */

void (*interpolation)(double, double*, double*, long);

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
   interpolation(h, Y[j], (coeff)[j], nx+2);
}

 
/* free memory */
free_double_matrix (Y, ny+2, nx+2);

return;

}  /* interpolate_image_horizontal */

/*--------------------------------------------------------------------------*/

void interpolate_image_vertical

     (long     nx,              /* image dimension in x direction */
      long     ny,              /* image dimension in y direction */
      double   h,               /* pixel size in x and y direction */
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

void (*interpolation)(double, double*, double*, long);

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
   interpolation(h, u[i], (coeff)[i], ny+2);
}


return;

}  /* interpolate_image_vertical */


/*--------------------------------------------------------------------------*/

void interpolate_image

     (long     nx,              /* image dimension in x direction */
      long     ny,              /* image dimension in y direction */
      double   h,               /* pixel size in x and y direction */
      double   **u,             /* input image */
      double   **coeff_x,       /* coefficients in x direction, output */
      double   **coeff_y,       /* coefficients in y direction, output */
      long     method)          /* interpolation method */

/*
   interpolates image in x and y direction along pixels;
   returns coefficients of piecewise cubic polynomial interpolant
*/

{

interpolate_image_horizontal (nx, ny, h, u, coeff_x, method);
interpolate_image_vertical (nx, ny, h, u, coeff_y, method);

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
  
  /* calculate h_sqr = || (xi_x/xi_y, 1)^T ||² */
  *h_sqr = h * (1.0 + dir_x * dir_x / (dir_y * dir_y));
 
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
  
  /* calculate h_sqr = || (1, dir_y/dir_x)^T ||² */
  *h_sqr = h * (1.0 + dir_y * dir_y / (dir_x * dir_x));
  
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

interpolate_image (nx, ny, h, f, coeff_x, coeff_y, method); 

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
   printf ("iteration:   %ld / %ld \r", k, kmax);
   mcm_step (tau, nx, ny, h, u, coeff_x, coeff_y, method);
   }  
   
/* ---- free memory ---- */
free_double_matrix (coeff_x, ny+2, (nx+2-1) * p);
free_double_matrix (coeff_y, nx+2, (ny+2-1) * p);

return;

}  /* mcm */

/*--------------------------------------------------------------------------*/
/*                           MAIN MENU                                      */
/*--------------------------------------------------------------------------*/

int main ()

{
long    i;                    /* loop variable */
char    in[80];               /* for reading data */
char    out[80];              /* for reading data */
double  **u;                  /* image */
long    nx, ny;               /* image size in x, y direction */ 
double  tau;                  /* time step size */
long    kmax;                 /* number of iterations */
long    method;               /* method parameter */
double  h;                    /* pixel size in x and y direction */
double  max, min;             /* largest, smallest grey value */
double  mean;                 /* average grey value */
double  std;                  /* standard deviation */
char    comments[1600];       /* string for comments */
char    divider[81];          /* string for list divider */

/* method descriptions */
char description[5][80] = {
   "Linear Spline", 
   "Cubic Spline", 
   "Fritsch-Carlson (1980)", 
   "Extended Two-Sweep (Eisenstat 1985)", 
   "MQSI (Lux 2020)"
};

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


/* ---- read input image (pgm format P5) ---- */

printf ("input image (pgm):                         ");
read_string (in);

read_pgm_to_double (in, &nx, &ny, &u);
printf ("\n");

/* ---- read parameters ---- */

/* initialize list divider */
for(i = 0; i < 80; i++) {
   divider[i] = '-';
}
divider[80] = '\0';

/* print interpolation methods table */
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

/* read interpolation method */
printf ("\n");
printf ("interpolation method (1-5):                ");
read_long (&method);

/* read time step size */
printf ("time step size tau (<0.500):               ");
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
printf ("initial image\n");
printf ("minimum:          %8.2lf \n", min);
printf ("maximum:          %8.2lf \n", max);
printf ("mean:             %8.2lf \n", mean);
printf ("standard dev.:    %8.2lf \n\n", std);


/* ---- process image ---- */

/* initialize h */
h = 1.0;

/* reflecting dummy boundaries */

dummies_double (u, nx, ny);

/* process image */
mcm(tau, kmax, nx, ny, h, u, method);

/* analyse processed image */   
analyse_grey_double (u, nx, ny, &min, &max, &mean, &std);
printf ("\n\nprocessed image:\n");
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
