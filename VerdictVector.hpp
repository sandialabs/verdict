/*=========================================================================

  Module:    VerdictVector.hpp

  Copyright 2003,2006,2019 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
  Under the terms of Contract DE-NA0003525 with NTESS,
  the U.S. Government retains certain rights in this software.

  See LICENSE for details.

=========================================================================*/

/*
 *
 * VerdictVector.hpp contains declarations of vector operations
 *
 * This file is part of VERDICT
 *
 */

// .SECTION Thanks
// Prior to its inclusion within VTK, this code was developed by the CUBIT
// project at Sandia National Laboratories.

#ifndef VERDICTVECTOR_HPP
#define VERDICTVECTOR_HPP

#include "verdict.h"

#include <cassert>
#include <math.h>

namespace VERDICT_NAMESPACE
{
class VerdictVector
{
public:
  //- Heading: Constructors and Destructor
  VERDICT_HOST_DEVICE VerdictVector(); //- Default constructor.

  VERDICT_HOST_DEVICE VerdictVector(const double x, const double y, const double z);
  //- Constructor: create vector from three components

  VERDICT_HOST_DEVICE VerdictVector(const double xyz[3]);
  //- Constructor: create vector from tuple

  VERDICT_HOST_DEVICE VerdictVector(const VerdictVector& tail, const VerdictVector& head);
  VERDICT_HOST_DEVICE VerdictVector(const double *tail, const double *head, int dimension);
  VERDICT_HOST_DEVICE VerdictVector(const double *tail, const double *head);
  //- Constructor for a VerdictVector starting at tail and pointing
  //- to head.

  template <typename ARG1, typename ARG2, typename ARG3> VERDICT_HOST_DEVICE VerdictVector(ARG1, ARG2, ARG3) = delete;
  //- define this template to avoid ambiguity between the (double, double, double) and (double *, double *, int) constructors

  VERDICT_HOST_DEVICE VerdictVector(const VerdictVector& copy_from); //- Copy Constructor

  //- Heading: Set and Inquire Functions
  VERDICT_HOST_DEVICE void set(const double xv, const double yv, const double zv);
  //- Change vector components to {x}, {y}, and {z}

  VERDICT_HOST_DEVICE void set(const double xyz[3]);
  //- Change vector components to xyz[0], xyz[1], xyz[2]

  VERDICT_HOST_DEVICE void set(const VerdictVector& tail, const VerdictVector& head);
  //- Change vector to go from tail to head.

  VERDICT_HOST_DEVICE void set(const VerdictVector& to_copy);
  //- Same as operator=(const VerdictVector&)

  VERDICT_HOST_DEVICE double x() const; //- Return x component of vector

  VERDICT_HOST_DEVICE double y() const; //- Return y component of vector

  VERDICT_HOST_DEVICE double z() const; //- Return z component of vector

  VERDICT_HOST_DEVICE void get_xyz(double& x, double& y, double& z); //- Get x, y, z components
  VERDICT_HOST_DEVICE void get_xyz(double xyz[3]);                   //- Get xyz tuple

  VERDICT_HOST_DEVICE double& r(); //- Return r component of vector, if (r,theta) format

  VERDICT_HOST_DEVICE double& theta(); //- Return theta component of vector, if (r,theta) format

  VERDICT_HOST_DEVICE void x(const double xv); //- Set x component of vector

  VERDICT_HOST_DEVICE void y(const double yv); //- Set y component of vector

  VERDICT_HOST_DEVICE void z(const double zv); //- Set z component of vector

  VERDICT_HOST_DEVICE void r(const double xv); //- Set r component of vector, if (r,theta) format

  VERDICT_HOST_DEVICE void theta(const double yv); //- Set theta component of vector, if (r,theta) format

  VERDICT_HOST_DEVICE double normalize();
  //- Normalize (set magnitude equal to 1) vector - return the magnitude

  VERDICT_HOST_DEVICE VerdictVector& length(const double new_length);
  //- Change length of vector to {new_length}. Can be used to move a
  //- location a specified distance from the origin in the current
  //- orientation.

  VERDICT_HOST_DEVICE double length() const;
  //- Calculate the length of the vector.
  //- Use {length_squared()} if only comparing lengths, not adding.

  VERDICT_HOST_DEVICE double length_squared() const;
  //- Calculate the squared length of the vector.
  //- Faster than {length()} since it eliminates the square root if
  //- only comparing other lengths.

  VERDICT_HOST_DEVICE double interior_angle(const VerdictVector& otherVector);
  //- Calculate the interior angle: acos((a%b)/(|a||b|))
  //- Returns angle in degrees.

  VERDICT_HOST_DEVICE void perpendicular_z();
  //- Transform this vector to a perpendicular one, leaving
  //- z-component alone. Rotates clockwise about the z-axis by pi/2.

  //- Heading: Operator Overloads  *****************************
  VERDICT_HOST_DEVICE VerdictVector& operator+=(const VerdictVector& vec);
  //- Compound Assignment: addition: {this = this + vec}

  VERDICT_HOST_DEVICE VerdictVector& operator-=(const VerdictVector& vec);
  //- Compound Assignment: subtraction: {this = this - vec}

  VERDICT_HOST_DEVICE VerdictVector& operator*=(const VerdictVector& vec);
  //- Compound Assignment: cross product: {this = this * vec},
  //- non-commutative

  VERDICT_HOST_DEVICE VerdictVector& operator*=(const double scalar);
  //- Compound Assignment: multiplication: {this = this * scalar}

  VERDICT_HOST_DEVICE VerdictVector& operator/=(const double scalar);
  //- Compound Assignment: division: {this = this / scalar}

  VERDICT_HOST_DEVICE VerdictVector operator-() const;
  //- unary negation.

  VERDICT_HOST_DEVICE friend VerdictVector operator~(const VerdictVector& vec);
  //- normalize. Returns a new vector which is a copy of {vec},
  //- scaled such that {|vec|=1}. Uses overloaded bitwise NOT operator.

  VERDICT_HOST_DEVICE friend VerdictVector operator+(const VerdictVector& v1, const VerdictVector& v2);
  //- vector addition

  VERDICT_HOST_DEVICE friend VerdictVector operator-(const VerdictVector& v1, const VerdictVector& v2);
  //- vector subtraction

  VERDICT_HOST_DEVICE friend VerdictVector operator*(const VerdictVector& v1, const VerdictVector& v2);
  //- vector cross product, non-commutative

  VERDICT_HOST_DEVICE friend VerdictVector operator*(const VerdictVector& v1, const double sclr);
  //- vector * scalar

  VERDICT_HOST_DEVICE friend VerdictVector operator*(const double sclr, const VerdictVector& v1);
  //- scalar * vector

  VERDICT_HOST_DEVICE friend double operator%(const VerdictVector& v1, const VerdictVector& v2);
  //- dot product

  static VERDICT_HOST_DEVICE double Dot(const VerdictVector& v1, const VerdictVector& v2);
  //- dot product

  VERDICT_HOST_DEVICE friend VerdictVector operator/(const VerdictVector& v1, const double sclr);
  //- vector / scalar

  VERDICT_HOST_DEVICE friend int operator==(const VerdictVector& v1, const VerdictVector& v2);
  //- Equality operator

  VERDICT_HOST_DEVICE friend int operator!=(const VerdictVector& v1, const VerdictVector& v2);
  //- Inequality operator

  VERDICT_HOST_DEVICE VerdictVector& operator=(const VerdictVector& from);

private:
  double xVal; //- x component of vector.
  double yVal; //- y component of vector.
  double zVal; //- z component of vector.
};

VERDICT_HOST_DEVICE inline double VerdictVector::x() const
{
  return xVal;
}
VERDICT_HOST_DEVICE inline double VerdictVector::y() const
{
  return yVal;
}
VERDICT_HOST_DEVICE inline double VerdictVector::z() const
{
  return zVal;
}
VERDICT_HOST_DEVICE inline void VerdictVector::get_xyz(double xyz[3])
{
  xyz[0] = xVal;
  xyz[1] = yVal;
  xyz[2] = zVal;
}
VERDICT_HOST_DEVICE inline void VerdictVector::get_xyz(double& xv, double& yv, double& zv)
{
  xv = xVal;
  yv = yVal;
  zv = zVal;
}
VERDICT_HOST_DEVICE inline double& VerdictVector::r()
{
  return xVal;
}
VERDICT_HOST_DEVICE inline double& VerdictVector::theta()
{
  return yVal;
}
VERDICT_HOST_DEVICE inline void VerdictVector::x(const double xv)
{
  xVal = xv;
}
VERDICT_HOST_DEVICE inline void VerdictVector::y(const double yv)
{
  yVal = yv;
}
VERDICT_HOST_DEVICE inline void VerdictVector::z(const double zv)
{
  zVal = zv;
}
VERDICT_HOST_DEVICE inline void VerdictVector::r(const double xv)
{
  xVal = xv;
}
VERDICT_HOST_DEVICE inline void VerdictVector::theta(const double yv)
{
  yVal = yv;
}
VERDICT_HOST_DEVICE inline VerdictVector& VerdictVector::operator+=(const VerdictVector& vector)
{
  xVal += vector.x();
  yVal += vector.y();
  zVal += vector.z();
  return *this;
}

VERDICT_HOST_DEVICE inline VerdictVector& VerdictVector::operator-=(const VerdictVector& vector)
{
  xVal -= vector.x();
  yVal -= vector.y();
  zVal -= vector.z();
  return *this;
}

VERDICT_HOST_DEVICE inline VerdictVector& VerdictVector::operator*=(const VerdictVector& vector)
{
  double xcross, ycross, zcross;
  xcross = yVal * vector.z() - zVal * vector.y();
  ycross = zVal * vector.x() - xVal * vector.z();
  zcross = xVal * vector.y() - yVal * vector.x();
  xVal = xcross;
  yVal = ycross;
  zVal = zcross;
  return *this;
}

VERDICT_HOST_DEVICE inline VerdictVector::VerdictVector(const VerdictVector& copy_from)
  : xVal(copy_from.xVal)
  , yVal(copy_from.yVal)
  , zVal(copy_from.zVal)
{
}

VERDICT_HOST_DEVICE inline VerdictVector::VerdictVector()
  : xVal(0)
  , yVal(0)
  , zVal(0)
{
}

VERDICT_HOST_DEVICE inline VerdictVector::VerdictVector(const double *tail, const double *head, int dimension)
  : xVal{head[0] - tail[0]}
  , yVal{head[1] - tail[1]}
  , zVal{dimension == 2 ? 0.0 : head[2] - tail[2]}
{
}

VERDICT_HOST_DEVICE inline VerdictVector::VerdictVector(const double *tail, const double *head)
  : xVal{head[0] - tail[0]}
  , yVal{head[1] - tail[1]}
  , zVal{head[2] - tail[2]}
{
}

VERDICT_HOST_DEVICE inline VerdictVector::VerdictVector(const VerdictVector& tail, const VerdictVector& head)
  : xVal(head.xVal - tail.xVal)
  , yVal(head.yVal - tail.yVal)
  , zVal(head.zVal - tail.zVal)
{
}

VERDICT_HOST_DEVICE inline VerdictVector::VerdictVector(const double xIn, const double yIn, const double zIn)
  : xVal(xIn)
  , yVal(yIn)
  , zVal(zIn)
{
}

// This sets the vector to be perpendicular to it's current direction.
// NOTE:
//      This is a 2D function.  It only works in the XY Plane.
VERDICT_HOST_DEVICE inline void VerdictVector::perpendicular_z()
{
  double temp = x();
  x(y());
  y(-temp);
}

VERDICT_HOST_DEVICE inline void VerdictVector::set(const double xv, const double yv, const double zv)
{
  xVal = xv;
  yVal = yv;
  zVal = zv;
}

VERDICT_HOST_DEVICE inline void VerdictVector::set(const double xyz[3])
{
  xVal = xyz[0];
  yVal = xyz[1];
  zVal = xyz[2];
}

VERDICT_HOST_DEVICE inline void VerdictVector::set(const VerdictVector& tail, const VerdictVector& head)
{
  xVal = head.xVal - tail.xVal;
  yVal = head.yVal - tail.yVal;
  zVal = head.zVal - tail.zVal;
}

VERDICT_HOST_DEVICE inline VerdictVector& VerdictVector::operator=(const VerdictVector& from)
{
  xVal = from.xVal;
  yVal = from.yVal;
  zVal = from.zVal;
  return *this;
}

VERDICT_HOST_DEVICE inline void VerdictVector::set(const VerdictVector& to_copy)
{
  *this = to_copy;
}

// Scale all values by scalar.
VERDICT_HOST_DEVICE inline VerdictVector& VerdictVector::operator*=(const double scalar)
{
  xVal *= scalar;
  yVal *= scalar;
  zVal *= scalar;
  return *this;
}

// Scales all values by 1/scalar
VERDICT_HOST_DEVICE inline VerdictVector& VerdictVector::operator/=(const double scalar)
{
  assert(scalar != 0);
  xVal /= scalar;
  yVal /= scalar;
  zVal /= scalar;
  return *this;
}

// Returns the normalized 'this'.
VERDICT_HOST_DEVICE inline VerdictVector operator~(const VerdictVector& vec)
{
  double mag = sqrt(vec.xVal * vec.xVal + vec.yVal * vec.yVal + vec.zVal * vec.zVal);

  VerdictVector temp = vec;
  if (mag != 0.0)
  {
    temp /= mag;
  }
  return temp;
}

// Unary minus.  Negates all values in vector.
VERDICT_HOST_DEVICE inline VerdictVector VerdictVector::operator-() const
{
  return VerdictVector(-xVal, -yVal, -zVal);
}

VERDICT_HOST_DEVICE inline VerdictVector operator+(const VerdictVector& vector1, const VerdictVector& vector2)
{
  double xv = vector1.x() + vector2.x();
  double yv = vector1.y() + vector2.y();
  double zv = vector1.z() + vector2.z();
  return VerdictVector(xv, yv, zv);
  //  return VerdictVector(vector1) += vector2;
}

VERDICT_HOST_DEVICE inline VerdictVector operator-(const VerdictVector& vector1, const VerdictVector& vector2)
{
  double xv = vector1.x() - vector2.x();
  double yv = vector1.y() - vector2.y();
  double zv = vector1.z() - vector2.z();
  return VerdictVector(xv, yv, zv);
  //  return VerdictVector(vector1) -= vector2;
}

// Cross products.
// vector1 cross vector2
VERDICT_HOST_DEVICE inline VerdictVector operator*(const VerdictVector& vector1, const VerdictVector& vector2)
{
  return VerdictVector(vector1) *= vector2;
}

// Returns a scaled vector.
VERDICT_HOST_DEVICE inline VerdictVector operator*(const VerdictVector& vector1, const double scalar)
{
  return VerdictVector(vector1) *= scalar;
}

// Returns a scaled vector
VERDICT_HOST_DEVICE inline VerdictVector operator*(const double scalar, const VerdictVector& vector1)
{
  return VerdictVector(vector1) *= scalar;
}

// Returns a vector scaled by 1/scalar
VERDICT_HOST_DEVICE inline VerdictVector operator/(const VerdictVector& vector1, const double scalar)
{
  return VerdictVector(vector1) /= scalar;
}

VERDICT_HOST_DEVICE inline int operator==(const VerdictVector& v1, const VerdictVector& v2)
{
  return (v1.xVal == v2.xVal && v1.yVal == v2.yVal && v1.zVal == v2.zVal);
}

VERDICT_HOST_DEVICE inline int operator!=(const VerdictVector& v1, const VerdictVector& v2)
{
  return (v1.xVal != v2.xVal || v1.yVal != v2.yVal || v1.zVal != v2.zVal);
}

VERDICT_HOST_DEVICE inline double VerdictVector::length_squared() const
{
  return (xVal * xVal + yVal * yVal + zVal * zVal);
}

VERDICT_HOST_DEVICE inline double VerdictVector::length() const
{
  return sqrt(length_squared())
}

VERDICT_HOST_DEVICE inline double VerdictVector::normalize()
{
  double mag = length();
  if (mag != 0)
  {
    xVal = xVal / mag;
    yVal = yVal / mag;
    zVal = zVal / mag;
  }
  return mag;
}

// Dot Product.
VERDICT_HOST_DEVICE inline double operator%(const VerdictVector& vector1, const VerdictVector& vector2)
{
  return VerdictVector::Dot(vector1, vector2);
}
VERDICT_HOST_DEVICE inline double VerdictVector::Dot(const VerdictVector& vector1, const VerdictVector& vector2)
{
  return (vector1.xVal * vector2.xVal + vector1.yVal * vector2.yVal + vector1.zVal * vector2.zVal);
}
} // namespace verdict

#endif
