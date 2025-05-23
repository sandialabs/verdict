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
#include <utility>
#include <math.h>

namespace VERDICT_NAMESPACE
{
class VerdictVector
{
public:
  //- Heading: Constructors and Destructor
  VERDICT_HOST_DEVICE constexpr VerdictVector(); //- Default constructor.

  VERDICT_HOST_DEVICE constexpr VerdictVector(const double x, const double y, const double z);
  //- Constructor: create vector from three components

  VERDICT_HOST_DEVICE constexpr VerdictVector(const double xyz[3]);
  //- Constructor: create vector from tuple

  VERDICT_HOST_DEVICE constexpr VerdictVector(const VerdictVector& tail, const VerdictVector& head);
  VERDICT_HOST_DEVICE constexpr VerdictVector(const double *tail, const double *head, int dimension);
  VERDICT_HOST_DEVICE constexpr VerdictVector(const double *tail, const double *head);
  //- Constructor for a VerdictVector starting at tail and pointing
  //- to head.

  template <typename ARG1, typename ARG2, typename ARG3> constexpr VerdictVector(ARG1, ARG2, ARG3) = delete;
  //- define this template to avoid ambiguity between the (double, double, double) and (double *, double *, int) constructors

  VERDICT_HOST_DEVICE constexpr VerdictVector(const VerdictVector& copy_from); //- Copy Constructor

  //- Heading: Set and Inquire Functions
  VERDICT_HOST_DEVICE constexpr void set(const double xv, const double yv, const double zv);
  //- Change vector components to {x}, {y}, and {z}

  VERDICT_HOST_DEVICE constexpr void set(const double xyz[3]);
  //- Change vector components to xyz[0], xyz[1], xyz[2]

  VERDICT_HOST_DEVICE constexpr void set(const VerdictVector& tail, const VerdictVector& head);
  //- Change vector to go from tail to head.

  VERDICT_HOST_DEVICE constexpr void set(const VerdictVector& to_copy);
  //- Same as operator=(const VerdictVector&)

  VERDICT_HOST_DEVICE constexpr double x() const; //- Return x component of vector

  VERDICT_HOST_DEVICE constexpr double y() const; //- Return y component of vector

  VERDICT_HOST_DEVICE constexpr double z() const; //- Return z component of vector

  VERDICT_HOST_DEVICE constexpr void get_xyz(double& x, double& y, double& z); //- Get x, y, z components
  VERDICT_HOST_DEVICE constexpr void get_xyz(double xyz[3]);                   //- Get xyz tuple

  VERDICT_HOST_DEVICE constexpr double& r(); //- Return r component of vector, if (r,theta) format

  VERDICT_HOST_DEVICE constexpr double& theta(); //- Return theta component of vector, if (r,theta) format

  VERDICT_HOST_DEVICE constexpr void x(const double xv); //- Set x component of vector

  VERDICT_HOST_DEVICE constexpr void y(const double yv); //- Set y component of vector

  VERDICT_HOST_DEVICE constexpr void z(const double zv); //- Set z component of vector

  VERDICT_HOST_DEVICE constexpr void r(const double xv); //- Set r component of vector, if (r,theta) format

  VERDICT_HOST_DEVICE constexpr void theta(const double yv); //- Set theta component of vector, if (r,theta) format

  VERDICT_HOST_DEVICE double normalize();
  //- Normalize (set magnitude equal to 1) vector - return the magnitude

  VERDICT_HOST_DEVICE VerdictVector& length(const double new_length);
  //- Change length of vector to {new_length}. Can be used to move a
  //- location a specified distance from the origin in the current
  //- orientation.

  VERDICT_HOST_DEVICE double length() const;
  //- Calculate the length of the vector.
  //- Use {length_squared()} if only comparing lengths, not adding.

  VERDICT_HOST_DEVICE constexpr double length_squared() const;
  //- Calculate the squared length of the vector.
  //- Faster than {length()} since it eliminates the square root if
  //- only comparing other lengths.

  VERDICT_HOST_DEVICE double interior_angle(const VerdictVector& otherVector);
  //- Calculate the interior angle: acos((a%b)/(|a||b|))
  //- Returns angle in degrees.

  VERDICT_HOST_DEVICE constexpr void perpendicular_z();
  //- Transform this vector to a perpendicular one, leaving
  //- z-component alone. Rotates clockwise about the z-axis by pi/2.

  //- Heading: Operator Overloads  *****************************
  VERDICT_HOST_DEVICE constexpr VerdictVector& operator+=(const VerdictVector& vec);
  //- Compound Assignment: addition: {this = this + vec}

  VERDICT_HOST_DEVICE constexpr VerdictVector& operator-=(const VerdictVector& vec);
  //- Compound Assignment: subtraction: {this = this - vec}

  VERDICT_HOST_DEVICE constexpr VerdictVector& operator*=(const VerdictVector& vec);
  //- Compound Assignment: cross product: {this = this * vec},
  //- non-commutative

  VERDICT_HOST_DEVICE constexpr VerdictVector& operator*=(const double scalar);
  //- Compound Assignment: multiplication: {this = this * scalar}

  VERDICT_HOST_DEVICE constexpr VerdictVector& operator/=(const double scalar);
  //- Compound Assignment: division: {this = this / scalar}

  VERDICT_HOST_DEVICE constexpr VerdictVector operator-() const;
  //- unary negation.

  VERDICT_HOST_DEVICE friend VerdictVector operator~(const VerdictVector& vec);
  //- normalize. Returns a new vector which is a copy of {vec},
  //- scaled such that {|vec|=1}. Uses overloaded bitwise NOT operator.

  VERDICT_HOST_DEVICE friend constexpr VerdictVector operator+(const VerdictVector& v1, const VerdictVector& v2);
  //- vector addition

  VERDICT_HOST_DEVICE friend constexpr VerdictVector operator-(const VerdictVector& v1, const VerdictVector& v2);
  //- vector subtraction

  VERDICT_HOST_DEVICE friend constexpr VerdictVector operator*(const VerdictVector& v1, const VerdictVector& v2);
  //- vector cross product, non-commutative

  VERDICT_HOST_DEVICE friend constexpr VerdictVector operator*(const VerdictVector& v1, const double sclr);
  //- vector * scalar

  VERDICT_HOST_DEVICE friend constexpr VerdictVector operator*(const double sclr, const VerdictVector& v1);
  //- scalar * vector

  VERDICT_HOST_DEVICE friend constexpr double operator%(const VerdictVector& v1, const VerdictVector& v2);
  //- dot product

  VERDICT_HOST_DEVICE static constexpr double Dot(const VerdictVector& v1, const VerdictVector& v2);
  //- dot product

  VERDICT_HOST_DEVICE friend constexpr VerdictVector operator/(const VerdictVector& v1, const double sclr);
  //- vector / scalar

  VERDICT_HOST_DEVICE friend constexpr int operator==(const VerdictVector& v1, const VerdictVector& v2);
  //- Equality operator

  VERDICT_HOST_DEVICE friend constexpr int operator!=(const VerdictVector& v1, const VerdictVector& v2);
  //- Inequality operator

  VERDICT_HOST_DEVICE constexpr VerdictVector& operator=(const VerdictVector& from);

private:
  double xVal; //- x component of vector.
  double yVal; //- y component of vector.
  double zVal; //- z component of vector.
};

constexpr double VerdictVector::x() const
{
  return xVal;
}
constexpr double VerdictVector::y() const
{
  return yVal;
}
constexpr double VerdictVector::z() const
{
  return zVal;
}
constexpr void VerdictVector::get_xyz(double xyz[3])
{
  xyz[0] = xVal;
  xyz[1] = yVal;
  xyz[2] = zVal;
}
constexpr void VerdictVector::get_xyz(double& xv, double& yv, double& zv)
{
  xv = xVal;
  yv = yVal;
  zv = zVal;
}
constexpr double& VerdictVector::r()
{
  return xVal;
}
constexpr double& VerdictVector::theta()
{
  return yVal;
}
constexpr void VerdictVector::x(const double xv)
{
  xVal = xv;
}
constexpr void VerdictVector::y(const double yv)
{
  yVal = yv;
}
constexpr void VerdictVector::z(const double zv)
{
  zVal = zv;
}
constexpr void VerdictVector::r(const double xv)
{
  xVal = xv;
}
constexpr void VerdictVector::theta(const double yv)
{
  yVal = yv;
}
constexpr VerdictVector& VerdictVector::operator+=(const VerdictVector& vector)
{
  xVal += vector.x();
  yVal += vector.y();
  zVal += vector.z();
  return *this;
}

constexpr VerdictVector& VerdictVector::operator-=(const VerdictVector& vector)
{
  xVal -= vector.x();
  yVal -= vector.y();
  zVal -= vector.z();
  return *this;
}

constexpr VerdictVector& VerdictVector::operator*=(const VerdictVector& vector)
{
  const double xcross = yVal * vector.z() - zVal * vector.y();
  const double ycross = zVal * vector.x() - xVal * vector.z();
  const double zcross = xVal * vector.y() - yVal * vector.x();
  xVal = xcross;
  yVal = ycross;
  zVal = zcross;
  return *this;
}

constexpr VerdictVector::VerdictVector(const VerdictVector& copy_from)
  : xVal(copy_from.xVal)
  , yVal(copy_from.yVal)
  , zVal(copy_from.zVal)
{
}

constexpr VerdictVector::VerdictVector()
  : xVal(0)
  , yVal(0)
  , zVal(0)
{
}

constexpr VerdictVector::VerdictVector(const double *tail, const double *head, int dimension)
  : xVal{head[0] - tail[0]}
  , yVal{head[1] - tail[1]}
  , zVal{dimension == 2 ? 0.0 : head[2] - tail[2]}
{
}

constexpr VerdictVector::VerdictVector(const double *tail, const double *head)
  : xVal{head[0] - tail[0]}
  , yVal{head[1] - tail[1]}
  , zVal{head[2] - tail[2]}
{
}

constexpr VerdictVector::VerdictVector(const VerdictVector& tail, const VerdictVector& head)
  : xVal(head.xVal - tail.xVal)
  , yVal(head.yVal - tail.yVal)
  , zVal(head.zVal - tail.zVal)
{
}

constexpr VerdictVector::VerdictVector(const double xIn, const double yIn, const double zIn)
  : xVal(xIn)
  , yVal(yIn)
  , zVal(zIn)
{
}

constexpr VerdictVector::VerdictVector(const double xyz[3])
  : xVal(xyz[0])
  , yVal(xyz[1])
  , zVal(xyz[2])
{
}

// This sets the vector to be perpendicular to it's current direction.
// NOTE:
//      This is a 2D function.  It only works in the XY Plane.
constexpr void VerdictVector::perpendicular_z()
{
  double temp = x();
  x(y());
  y(-temp);
}

constexpr void VerdictVector::set(const double xv, const double yv, const double zv)
{
  xVal = xv;
  yVal = yv;
  zVal = zv;
}

constexpr void VerdictVector::set(const double xyz[3])
{
  xVal = xyz[0];
  yVal = xyz[1];
  zVal = xyz[2];
}

constexpr void VerdictVector::set(const VerdictVector& tail, const VerdictVector& head)
{
  xVal = head.xVal - tail.xVal;
  yVal = head.yVal - tail.yVal;
  zVal = head.zVal - tail.zVal;
}

constexpr VerdictVector& VerdictVector::operator=(const VerdictVector& from)
{
  xVal = from.xVal;
  yVal = from.yVal;
  zVal = from.zVal;
  return *this;
}

constexpr void VerdictVector::set(const VerdictVector& to_copy)
{
  *this = to_copy;
}

// Scale all values by scalar.
constexpr VerdictVector& VerdictVector::operator*=(const double scalar)
{
  xVal *= scalar;
  yVal *= scalar;
  zVal *= scalar;
  return *this;
}

// Scales all values by 1/scalar
constexpr VerdictVector& VerdictVector::operator/=(const double scalar)
{
#ifndef __HIP_DEVICE_COMPILE__
  assert(scalar != 0);
#endif
  xVal /= scalar;
  yVal /= scalar;
  zVal /= scalar;
  return *this;
}

// Returns the normalized 'this'.
VERDICT_HOST_DEVICE inline VerdictVector operator~(const VerdictVector& vec)
{
  VerdictVector temp = vec;
  temp.normalize();
  return temp;
}

// Unary minus.  Negates all values in vector.
constexpr VerdictVector VerdictVector::operator-() const
{
  return VerdictVector(-xVal, -yVal, -zVal);
}

constexpr VerdictVector operator+(const VerdictVector& vector1, const VerdictVector& vector2)
{
  double xv = vector1.x() + vector2.x();
  double yv = vector1.y() + vector2.y();
  double zv = vector1.z() + vector2.z();
  return VerdictVector(xv, yv, zv);
  //  return VerdictVector(vector1) += vector2;
}

constexpr VerdictVector operator-(const VerdictVector& vector1, const VerdictVector& vector2)
{
  double xv = vector1.x() - vector2.x();
  double yv = vector1.y() - vector2.y();
  double zv = vector1.z() - vector2.z();
  return VerdictVector(xv, yv, zv);
  //  return VerdictVector(vector1) -= vector2;
}

// Cross products.
// vector1 cross vector2
constexpr VerdictVector operator*(const VerdictVector& vector1, const VerdictVector& vector2)
{
  return VerdictVector(vector1) *= vector2;
}

// Returns a scaled vector.
constexpr VerdictVector operator*(const VerdictVector& vector1, const double scalar)
{
  return VerdictVector(vector1) *= scalar;
}

// Returns a scaled vector
constexpr VerdictVector operator*(const double scalar, const VerdictVector& vector1)
{
  return VerdictVector(vector1) *= scalar;
}

// Returns a vector scaled by 1/scalar
constexpr VerdictVector operator/(const VerdictVector& vector1, const double scalar)
{
  return VerdictVector(vector1) /= scalar;
}

constexpr int operator==(const VerdictVector& v1, const VerdictVector& v2)
{
  return (v1.xVal == v2.xVal && v1.yVal == v2.yVal && v1.zVal == v2.zVal);
}

constexpr int operator!=(const VerdictVector& v1, const VerdictVector& v2)
{
  return (v1.xVal != v2.xVal || v1.yVal != v2.yVal || v1.zVal != v2.zVal);
}

constexpr double VerdictVector::length_squared() const
{
  return (xVal * xVal + yVal * yVal + zVal * zVal);
}

VERDICT_HOST_DEVICE inline double VerdictVector::length() const
{
  return (sqrt(xVal * xVal + yVal * yVal + zVal * zVal));
}

// Dot Product.
constexpr double operator%(const VerdictVector& vector1, const VerdictVector& vector2)
{
  return VerdictVector::Dot(vector1, vector2);
}
constexpr double VerdictVector::Dot(const VerdictVector& vector1, const VerdictVector& vector2)
{
  return (vector1.xVal * vector2.xVal + vector1.yVal * vector2.yVal + vector1.zVal * vector2.zVal);
}

struct ElemScale
{
  VerdictVector center;
  double scale;
};

template <typename CoordsContainerType>
VERDICT_HOST_DEVICE constexpr ElemScale elem_scaling(int num_coords, const CoordsContainerType coordinates, int dimension = 3)
{
  VerdictVector min(VERDICT_DBL_MAX, VERDICT_DBL_MAX, dimension == 3 ? VERDICT_DBL_MAX : 0);
  VerdictVector max(-VERDICT_DBL_MAX, -VERDICT_DBL_MAX, dimension == 3 ? -VERDICT_DBL_MAX : 0);
  VerdictVector center(0.0, 0.0, 0.0);

  for (int i = 0; i < num_coords; i++)
  {
    if (coordinates[i][0] < min.x())
      min.x(coordinates[i][0]);
    if (coordinates[i][1] < min.y())
      min.y(coordinates[i][1]);
    if (coordinates[i][0] > max.x())
      max.x(coordinates[i][0]);
    if (coordinates[i][1] > max.y())
      max.y(coordinates[i][1]);
    if (dimension == 3)
    {
      if (coordinates[i][2] < min.z())
        min.z(coordinates[i][2]);
      if (coordinates[i][2] > max.z())
        max.z(coordinates[i][2]);
    }
    center += VerdictVector(coordinates[i][0], coordinates[i][1], dimension == 3 ? coordinates[i][2] : 0.0);
  }
  center /= (double)num_coords;

  double len = (max - min).length();
  if (len < VERDICT_DBL_MIN)
  {
    center = VerdictVector(0.0,0.0,0.0);
    len = 1.0;
  }
  return ElemScale{center, len};
}

template <typename CoordsContainerType>
VERDICT_HOST_DEVICE constexpr double apply_elem_scaling_on_points(int num_coords, const CoordsContainerType coordinates, int num_vec, VerdictVector* v, int dimension = 3)
{
  auto char_size = elem_scaling(num_coords, coordinates, dimension);
  for (int i = 0; i < num_vec; i++)
  {
    v[i] -= char_size.center;
    v[i] /= char_size.scale;
  }
  return char_size.scale;
}

template <typename CoordsContainerType>
VERDICT_HOST_DEVICE constexpr double apply_elem_scaling_on_edges(int num_coords, const CoordsContainerType coordinates, int num_vec, VerdictVector* v, int dimension = 3)
{
  auto char_size = elem_scaling(num_coords, coordinates, dimension);
  for (int i = 0; i < num_vec; i++)
  {
    v[i] /= char_size.scale;
  }
  return char_size.scale;
}


} // namespace verdict

#endif
