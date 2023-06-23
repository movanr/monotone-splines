#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>
#include <time.h>


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*        MEAN CURVATURE MOTION, EXPLICIT SCHEME WITH MONOTONE              */
/*                 PIECEWISE CUBIC INTERPOLATION                            */
/*                                                                          */
/*            Bachelor studies by Moritz van Recum                          */
/*                                                                          */
/*                Supervisor: Joachim Weickert                              */
/*                                                                          */
/*            (Copyright by Joachim Weickert, 11/2021)                      */
/*--------------------------------------------------------------------------*/


#define EPS 1e-14          /* small constant */


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
/*             functions for monotonicity constraints                       */
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
long i;
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
/*                   functions for ext. two-sweep method                    */
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
/*                     extended sweep methods                               */
/*--------------------------------------------------------------------------*/

void extended_forward_sweep

   (double  *d,      /* slopes, length n */
    double  *delta,  /* secant lines, length n-1 */
    long    n)       /* length */

/*
   extended forward sweep function 
   - modifies the slopes (d[i],d[i+1]); projects them onto the union of 
     delta[i]*M and delta[i]*D for i = 0 ,..., n-2 
   - only modifies d[i+1] in each iteration except for region A
   - increases or decreases them accordingly onto the boundary of the union of
     M and D 
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
   - modifies the slopes (d[i],d[i+1]); projects them onto delta[i]*M 
     for i = n-2 ,..., 0
   - only modifies d[i] in each iteration except for region E
   - If (x,y) lies in E then the (x,y) is not necessarily projected in a 
     straight (horizontal) line onto the boundary of M.  
   - increases or decreases them accordingly onto the boundary of M
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
/*                 monotone piecewise cubic interpolation                   */                   
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
      double   **coeff)        /* coefficients in x direction, output */

/*
   interpolates image in x and y direction along pixels;
   returns coefficients of piecewise cubic polynomial interpolant
   * equidistant version *
*/

{
long     i, j;      /* loop variables */
double   **Y;       /* values in x-direction (size ny * nx) */

alloc_double_matrix (&Y, ny+2, nx+2);

/* initialize breakpoints */
for (j=0; j<=ny+1; j++)
   for (i=0; i<=nx+1; i++)
      Y[j][i] = u[i][j];      /* Y-values of j-th image row */

      
/* ---- loop for x-direction ---- */

/* extended two-sweep method */
for (j=0; j<ny+2; j++)
   mpci_extended_sweep (h, Y[j], (coeff)[j], nx+2);
 
 
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
      double   **coeff)        /* coefficients in y direction, output */
/*
   interpolates image in y direction along pixels;
   returns coefficients of piecewise cubic polynomial interpolant
   * equidistant version *
*/

{
long     i;         /* loop variable */

/* ---- loop for y-direction ---- */

/* extended two-sweep method */
for (i=0; i<nx+2; i++)
   mpci_extended_sweep (h, u[i], (coeff)[i], ny+2);

return;

}  /* interpolate_image_vertical */


/*--------------------------------------------------------------------------*/

void interpolate_image

     (long     nx,              /* image dimension in x direction */
      long     ny,              /* image dimension in y direction */
      double   h,               /* pixel size in x and y direction */
      double   **u,             /* input image */
      double   **coeff_x,      /* coefficients in x direction, output */
      double   **coeff_y)      /* coefficients in y direction, output */

/*
   interpolates image in x and y direction along pixels;
   returns coefficients of piecewise cubic polynomial interpolant
*/

{

interpolate_image_horizontal (nx, ny, h, u, coeff_x);
interpolate_image_vertical (nx, ny, h, u, coeff_y);

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
     double    h, 
     long      nx,
     long      ny,
     double    **u,
     double    **coeff_x,
     double    **coeff_y)
     
/*
   approximates the neighbouring image values of u at pixel (i,j) in the 
   given direction;
   use interpolated image to obtain values between pixels
*/

{
double   offs;                /* offset variable, time saver */

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
     *neighbour1 = eval_cubic_hermite_poly (1-offs, coeff_x[j-1]+4*i); 
     *neighbour2 = eval_cubic_hermite_poly (offs, coeff_x[j+1]+4*(i-1)); 
     }
  else 
     {
     *neighbour1 = eval_cubic_hermite_poly (2-offs, coeff_x[j-1] + 4 * (i-1) ); 
     *neighbour2 = eval_cubic_hermite_poly (offs-1, coeff_x[j+1] + 4 * i); 
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
     *neighbour1 = eval_cubic_hermite_poly (1-offs, coeff_y[i-1] + 4 * j); 
     *neighbour2 = eval_cubic_hermite_poly (offs, coeff_y[i+1] + 4 * (j-1)); 
     }
  else
     {
     *neighbour1 = eval_cubic_hermite_poly (2-offs, coeff_y[i-1] + 4 * (j-1)); 
     *neighbour2 = eval_cubic_hermite_poly (offs-1, coeff_y[i+1] + 4 * j); 
     }  
  }

return;

}  /* approx_neighbours */


/*--------------------------------------------------------------------------*/
/*                              MCM                                         */
/*--------------------------------------------------------------------------*/

void mcm_step

     (double   tau,           /* time step size */
      long     kmax,          /* number of iterations */
      long     nx,            /* image dimension in x direction */
      long     ny,            /* image dimension in y direction */
      double   h,             /* pixel size in x and y direction */
      double   **u,           /* image, changed */
      double   **coeff_x,     /* hermite coefficients in x direction */
      double   **coeff_y)     /* hermite coefficients in y direction */

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

interpolate_image (nx, ny, h, f, coeff_x, coeff_y); 

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
        (&n1, &n2, &h_sqr, i, j, xi_x, xi_y, h, nx, ny, u, coeff_x, coeff_y);

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
      double   **u)       /* image, changed */
/*
   iterative explicit scheme for MCM with interpolation;
*/
   
{
long     k;               /* loop variable */
double   **coeff_x;       /* hermite coefficients in x direction */
double   **coeff_y;       /* hermite coefficients in y direction */

/* ---- allocate memory ---- */
/* (allocate memory for coefficients here to save computation time) */
alloc_double_matrix (&coeff_y, nx+2, (ny+2-1) * 4);
alloc_double_matrix (&coeff_x, ny+2, (nx+2-1) * 4);

/* ---- loop ---- */
for (k=1; k<=kmax; k++)
   {
   printf ("iteration:   %ld / %ld \r", k, kmax);
   
   mcm_step (tau, kmax, nx, ny, h, u, coeff_x, coeff_y);

   }  
   
/* ---- free memory ---- */
free_double_matrix (coeff_x, ny+2, (nx+2-1) * 4);
free_double_matrix (coeff_y, nx+2, (ny+2-1) * 4);

return;



}  /* mcm */

/*--------------------------------------------------------------------------*/

int main ()

{
char    in[80];               /* for reading data */
char    out[80];              /* for reading data */
double  **u;                  /* image */
long    nx, ny;               /* image size in x, y direction */ 
double  tau;                  /* time step size */
long    kmax;                 /* number of iterations */
double  h;                    /* pixel size in x and y direction */
double  max, min;             /* largest, smallest grey value */
double  mean;                 /* average grey value */
double  std;                  /* standard deviation */
char    comments[1600];       /* string for comments */

printf("\nMEAN CURVATURE MOTION FOR GREYSCALE (PGM) IMAGES\n");
printf("EXPLICIT SCHEME WITH MONOTOINE PIECEWISE CUBIC INTERPOLATION\n\n");
printf ("**************************************************\n\n");
printf ("    Bachelor studies by Moritz van Recum          \n\n");
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


/* ---- read parameters ---- */

printf ("time step size tau (<0.500):               ");
read_double (&tau);

printf ("number of iterations (>0):                 ");
read_long (&kmax);

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
mcm (tau, kmax, nx, ny, h, u);

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
comment_line (comments, "# MCM, explicit scheme with interpolation\n");
comment_line (comments, "# initial image:        %s\n", in);
comment_line (comments, "# interpolation method: %s\n", "monotone, order 4");
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
