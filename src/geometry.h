#include "math.h"

double getDistance(   int, int, double **, double * );
double getDistanceSq( int, int, double **, double * );
double getDistance(   int, int, double *,  double * );
double getDistanceSq( int, int, double *,  double * );
void minimumImageConvention( double &, double &, double &, double * );
void plane2plane( double*,  double*,  double*,  double*,  double*, double*,
                  double*&, double*&, double&, double*&, double& );
void     rotate( double *, double, double *, double *& );
double   getAngle( double *, double *, double * );
double   getDotProduct( double *, double * );
double * getCrossProduct( double *, double *, double *& );
double * getUnitVector( double *, double *& );
double   getMagnitude( double * );


/*-----------------------------------------
    Vector addition template
-----------------------------------------*/
template <typename TYPE>
void vector_add( TYPE *a, TYPE *b, TYPE *&c )
{
  c[0] = a[0] + b[0];
  c[1] = a[1] + b[1];
  c[2] = a[2] + b[2];
  return; 
}
