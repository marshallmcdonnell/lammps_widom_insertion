/*_._._._._._._._._._._._._._._._._._._._._
 
* File Name : distance.cpp

* Purpose :

* Creation Date : 27-07-2015

* Last Modified : Wed 14 Oct 2015 09:18:10 AM EDT

* Created by :

_._._._._._._._._._._._._._._._._._._._._._ */

#include "geometry.h"
#include "stdio.h"
#include "stdlib.h"
#include <iostream>

/*--------------------------------------------------------
    Get distance between iatom and jatom, return r_ij
    uses r as 2-D array
--------------------------------------------------------*/

double getDistance( int iatom, int jatom, double **r, double *bbox ) {
  double dx, dy, dz;
  double dis2, r_ij;
  dx = r[iatom][0] - r[jatom][0];
  dy = r[iatom][1] - r[jatom][1];
  dz = r[iatom][2] - r[jatom][2];

  // Apply minimum image convention
  minimumImageConvention( dx, dy, dz, bbox );

  // Get delta magnitude
  dis2  = dx*dx + dy*dy + dz*dz;
  r_ij = sqrt( dis2 );

  return r_ij;
}

/*--------------------------------------------------------
    Get squared distance between iatom and jatom, return r_ij
    uses r as 2-D array
--------------------------------------------------------*/

double getDistanceSq( int iatom, int jatom, double **r, double *bbox ) {
  double dx, dy, dz;
  double dis2;
  dx = r[iatom][0] - r[jatom][0];
  dy = r[iatom][1] - r[jatom][1];
  dz = r[iatom][2] - r[jatom][2];

  // Apply minimum image convention
  minimumImageConvention( dx, dy, dz, bbox );

  // Get delta magnitude
  dis2  = dx*dx + dy*dy + dz*dz;

  return dis2;
}

/*--------------------------------------------------------
    Get distance between iatom and jatom, return r_ij
    uses r as 1-D vector
--------------------------------------------------------*/

double getDistance( int iatom, int jatom, double *r, double *bbox ) {
  double dx, dy, dz;
  double dis2, r_ij;
  dx = r[3 * iatom + 0] - r[3 * jatom + 0];
  dy = r[3 * iatom + 1] - r[3 * jatom + 1];
  dz = r[3 * iatom + 2] - r[3 * jatom + 2];

  // Apply minimum image convention
  minimumImageConvention( dx, dy, dz, bbox );

  // Get delta magnitude
  dis2  = dx*dx + dy*dy + dz*dz;
  r_ij = sqrt( dis2 );

  return r_ij;
}

/*--------------------------------------------------------
    Get squared distance between iatom and jatom, return r_ij
    uses r as 1-D vector
--------------------------------------------------------*/

double getDistanceSq( int iatom, int jatom, double *r, double *bbox ) {
  double dx, dy, dz;
  double dis2;
  dx = r[3 * iatom + 0] - r[3 * jatom + 0];
  dy = r[3 * iatom + 1] - r[3 * jatom + 1];
  dz = r[3 * iatom + 2] - r[3 * jatom + 2];

  // Apply minimum image convention
  minimumImageConvention( dx, dy, dz, bbox );

  // Get delta magnitude
  dis2  = dx*dx + dy*dy + dz*dz;

  return dis2;
}
/*--------------------------------------------------------
    Apply minimum image convention
--------------------------------------------------------*/
void minimumImageConvention( double & dx, double & dy, double & dz, double * bbox ) {
  if      ( dx >      bbox[0] / 2.0 ) dx -= bbox[0];
  else if ( dx < -1.0*bbox[0] / 2.0 ) dx += bbox[0];

  if      ( dy >      bbox[1] / 2.0 ) dy -= bbox[1];
  else if ( dy < -1.0*bbox[1] / 2.0 ) dy += bbox[1];

  if      ( dz >      bbox[2] / 2.0 ) dz -= bbox[2];
  else if ( dz < -1.0*bbox[2] / 2.0 ) dz += bbox[2];

  return;
}

/*--------------------------------------------------------
    Get shift vector and two angles and axes of rotation
    from one plane to another plane defined by 
    2 sets of 3 points
--------------------------------------------------------*/
void plane2plane( double *r1, double *r2, double *r3, 
                  double *q1, double *q2, double *q3,
                  double *&shift, 
                  double *&axis1, double &theta,
                  double *&axis2, double &phi   )
{
  shift[0] = -1.0*r2[0] + q2[0];
  shift[1] = -1.0*r2[1] + q2[1];
  shift[2] = -1.0*r2[2] + q2[2];

  double * r1_shifted = new double[3];
  double * r2_shifted = new double[3];
  double * r3_shifted = new double[3];
  vector_add( r1, shift, r1_shifted );
  vector_add( r2, shift, r2_shifted );
  vector_add( r3, shift, r3_shifted );

  double origin[3];
  origin[0] =  r2_shifted[0];
  origin[1] =  r2_shifted[1];
  origin[2] =  r2_shifted[2];

  theta = getAngle( r1_shifted, origin, q1 );

  origin[0] *= -1.0;
  origin[1] *= -1.0;
  origin[2] *= -1.0;

  double * v1 = new double[3];
  double * v2 = new double[3];
  vector_add( r1_shifted, origin, v1 );
  vector_add( q1,         origin, v2 );

  getCrossProduct( v1, v2, axis1 );
  getUnitVector( axis1, axis1 );

  double * r1_rotated = new double[3];
  double * r2_rotated = new double[3];
  double * r3_rotated = new double[3];
  rotate( r1_shifted, -1.0*theta, axis1, r1_rotated );
  rotate( r2_shifted, -1.0*theta, axis1, r2_rotated );
  rotate( r3_shifted, -1.0*theta, axis1, r3_rotated );
 

  double * r_zaxis = new double[3];
  double * q_zaxis = new double[3];
  origin[0] = r2_rotated[0];
  origin[1] = r2_rotated[1];
  origin[2] = r2_rotated[2];

  origin[0] *= -1.0;
  origin[1] *= -1.0;
  origin[2] *= -1.0;

  vector_add( r1_rotated, origin, v1 ); 
  vector_add( r3_rotated, origin, v2 ); 

  getCrossProduct( v1, v2, r_zaxis );

  origin[0] = q2[0];
  origin[1] = q2[1];
  origin[2] = q2[2];

  origin[0] *= -1.0;
  origin[1] *= -1.0;
  origin[2] *= -1.0;

  vector_add( q1, origin, v1 );
  vector_add( q3, origin, v2 );

  getCrossProduct( v1, v2, q_zaxis );

  phi = getAngle( r_zaxis, origin, q_zaxis );

  getCrossProduct( r_zaxis, q_zaxis, axis2 );

  origin[0] = r2_rotated[0];
  origin[1] = r2_rotated[1];
  origin[2] = r2_rotated[2];
  
  origin[0] *= -1.0;
  origin[1] *= -1.0;
  origin[2] *= -1.0;

  vector_add( r1_rotated, origin, axis2 );
  getUnitVector( axis2, axis2 );

  delete [] r1_rotated;
  delete [] r2_rotated;
  delete [] r3_rotated;
  delete [] r1_shifted;
  delete [] r2_shifted;
  delete [] r3_shifted;
  delete [] r_zaxis;
  delete [] q_zaxis;
  delete [] v1;
  delete [] v2;

/*
   Now you can call  " ___(r_rotated)___ = rotate( ___(r)___, -1.d0*theta, axis1 ) "
   and               " ___(r_rotated)___ = rotate( ___(r)___, -1.d0*phi,   axis2 ) "
   and         call rotate( ____, -1.d0*phi,   axis2, _____ ) 
   to rotate from coordinate system A (  ___(r_rotated)___ on the left  )
   to coordinate system B (  __(r)__ on the right  )

   NOTE: the negative symbol is due to the fact that we 
         measured the angle in a positive direction but
         we rotate it in the opposite direction
*/

  return;
}

/*--------------------------------------------------------
    Perform a finite rotation of a vector (r) around a
    given axis by a given angle. Return rotated vector.

      This rotation is performed using the rotation formula
      found in "Classical Mechanics" by Goldstein, 2nd Ed.
      pg. 164-165
--------------------------------------------------------*/
void rotate( double *r, double angle, double *axis, double *& rotated_r )
{
  // Get dot product and cross product of vector w/
  // axis of rotation
  double dotprod;
  dotprod = getDotProduct( axis, r );

  double * cross = new double[3];
  getCrossProduct( r, axis, cross );

  // Rotate vector around the axis of rotation
  rotated_r[0] =                  cos( angle )   *     r[0] 
                 + dotprod*(1.0 - cos( angle ) ) *  axis[0]
                 +                sin( angle )   * cross[0];

  rotated_r[1] =                  cos( angle )   *     r[1] 
                 + dotprod*(1.0 - cos( angle ) ) *  axis[1]
                 +                sin( angle )   * cross[1];

  rotated_r[2] =                  cos( angle )   *     r[2] 
                 + dotprod*(1.0 - cos( angle ) ) *  axis[2]
                 +                sin( angle )   * cross[2];

  delete [] cross;
  return;
}

/*--------------------------------------------------------
    Get angle from an input of 3 points 
    w/ the 2nd point being the vertex
--------------------------------------------------------*/
double getAngle( double *r1, double *vertex, double *r2 ) 
{
  // Get vector between r1 / r2 and vertex 
  double dxyz1[3];
  double dxyz2[3];
  dxyz1[0] = r1[0] - vertex[0];
  dxyz1[1] = r1[1] - vertex[1];
  dxyz1[2] = r1[2] - vertex[2];

  dxyz2[0] = r2[0] - vertex[0];
  dxyz2[1] = r2[1] - vertex[1];
  dxyz2[2] = r2[2] - vertex[2];

  // Get magnitude of these two vectors
  double dis1, dis2;
  dis1 = sqrt(   dxyz1[0]*dxyz1[0] 
               + dxyz1[1]*dxyz1[1] 
               + dxyz1[2]*dxyz1[2] );

  dis2 = sqrt(   dxyz2[0]*dxyz2[0] 
               + dxyz2[1]*dxyz2[1] 
               + dxyz2[2]*dxyz2[2] );

  // Get dot product between these two vectors
  double dotprod;
  dotprod = getDotProduct( dxyz1, dxyz2 );

  // Get denominator of equation ( product of magnitudes )
  double denominator;
  denominator = dis1*dis2;
  

  double angle;
  if ( denominator != 0.00 ) 
    angle = acos( dotprod / denominator );
  else 
    angle = acos( dotprod );

  return angle;
}

/*--------------------------------------------------------
    Get dot product of 2 vectors 
--------------------------------------------------------*/
double getDotProduct( double *r1, double *r2 ) 
{
  double dotprod;
  dotprod = r1[0]*r2[0] + r1[1]*r2[1] + r1[2]*r2[2];
  return dotprod;
}

/*--------------------------------------------------------
    Get cross product of 2 vectors
--------------------------------------------------------*/
double * getCrossProduct( double *r1, double *r2, double *&crossprod ) 
{
  crossprod[0] = r1[1]*r2[2] - r1[2]*r2[1];
  crossprod[1] = r1[0]*r2[2] - r1[2]*r2[0];
  crossprod[1] *= -1.0;
  crossprod[2] = r1[0]*r2[1] - r1[1]*r2[0];
}

/*--------------------------------------------------------
    Get unit vector of a given vector
--------------------------------------------------------*/
double * getUnitVector( double *r, double*&unit_r ) 
{
  double magnitude;
  magnitude = getMagnitude( r );

  if ( magnitude != 0.0 ) {
    unit_r[0] = r[0] / magnitude;
    unit_r[1] = r[1] / magnitude;
    unit_r[2] = r[2] / magnitude;
  }
}

/*--------------------------------------------------------
    Get magnitude of vector
--------------------------------------------------------*/
double getMagnitude( double *r ) 
{
  double magnitude = getDotProduct( r, r );
  magnitude = sqrt( magnitude );
  return magnitude;
}
