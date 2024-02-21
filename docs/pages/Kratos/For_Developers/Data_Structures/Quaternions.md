---
title: Quaternions in Kratos
keywords: 
tags: [C++ Quaternions Tutorial]
sidebar: kratos_for_developers
summary: 
---

## Description
The Quaternion is a Kratos utility class used to represent 3D orientations and rotations, as an efficient and robust alternative to other representations such as Euler angles or rotation matrices. A Quaternion is represented in the form:

$$ w + x\boldsymbol{i} + y\boldsymbol{j} + z\boldsymbol{k} $$  

Where \\(w\\), \\(x\\), \\(y\\) and \\(z\\) are scalar coefficients (the only four members of the Quaternion class), and \\(\boldsymbol{i}\\), \\(\boldsymbol{j}\\), and \\(\boldsymbol{k}\\) are the fundamental quaternion units.

You can find an interactive and funny way to learn about quaternions [here](https://eater.net/quaternions).

## How to construct a Quaternion
To instantiate a quaternion you can either use one of the constructors, or some ad-hoc static functions.  
You can use the default constructor to instantiate a zero quaternion:

```cpp
Quaternion<double> q;
```

Or you can create a quaternion directly from its four coefficients like this:

```cpp
double w,x,y,z;
// ... code to set up the coefficients ...
Quaternion<double> q(w, x,y,z);
```

The identity quaternion (i.e. the quaternion that represents a zero rotation) can be either created using the full constructor or the dedicated static function:
```cpp
Quaternion<double> q(1.0, 0.0,0.0,0.0);
// or
Quaternion<double> q = Quaternion<double>::Identity();
```
You can construct a quaternion that represents a rotation of a given angle about a given axis like this:
```cpp
double x,y,z; // the components of the rotation axis
double angle; // the angle in radians
// ... code to set up the axis and angle ...
Quaternion<double> q = Quaternion<double>::FromAxisAngle(x,y,z, angle);
```

## How to access the coefficients
You can access the four coefficients of the quaternion like this:
```cpp
Quaternion<double> q;
double w = q.W();
double x = q.X();
double y = q.Y();
double z = q.Z();
```

## Conversion from/to other representations:
Quaternion from/to rotation vector (for example the 3 rotation DOFs of a node of beams or shell elements):
```cpp
double rx,ry,rz; // the components of the rotation vector
array_1d<double, 3> rxyz; // or the rotation vector itself
// quaternion from rotation vector
Quaternion<double> q = Quaternion<double>::FromRotationVector(rx, ry, rz); 
// or...
Quaternion<double> q = Quaternion<double>::FromRotationVector(rxyz);
// quaternion to rotation vector
q.ToRotationVector(rx, ry, rz); // rx,ry and rz are now the coefficient of the rotation vector obtained from q
// or...
q.ToRotationVector(rxyz); // rxyz is now the rotation vector obtained from q
```
Quaternion from/to Euler Angles (in Z(-X)Z sequence as in GiD):
```cpp
array_1d<double, 3> ea; // euler angles
// quaternion from euler angles
Quaternion<double> q = Quaternion<double>::FromEulerAngles(ea); 
// quaternion to euler anfles
q.ToEulerAngles(ea); // ea is now a vector containing the 3 euler angles
```
Quaternion from/to rotation matrices:
```cpp
Matrix R; // a 3x3 rotation matrix
// quaternion from rotation matrix
Quaternion<double> q = Quaternion<double>::FromRotationMatrix(R); 
// quaternion to rotation matrix
q.ToRotationMatrix(R); // R is now the rotation matrix obtained from q
```

## Unit Quaternion
Before performing operations that interpret the quaternion as a rotation, you must always normalize the quaternion!. The Quaternion class has functions to get the norm of the quaternion or to normalize the quaternion itself:
```cpp
Quaternion<double> q; // some quaternion...
// get the norm of q
double qnorm = q.norm(); // the norm of q, i.e. ||q||
double qnorm2 = q.squaredNorm(); // the squared norm of q, i.e. ||q||^2
// if you want to normalize q:
q.normalize(); // now q is a unit quaternion
```

## Quaternion as a rotation
A quaternion represents a rotation. A compound rotation can be obtained using quaternion multiplication (as for rotation matrix multiplication):
```cpp
Quaternion<double> q1; // a quaternion representing the first rotation R1
Quaternion<double> q2; // a quaternion representing the second rotation R2
Quaternion<double> q3 = q2*q1; // the quaternion representing the compound rotation R2*R1
```
The inverse rotation can be obtained conjugating the quaternion:
```cpp
Quaternion<double> q; // a quaternion representing a rotation R
Quaternion<double> iq = q.conjugate(); // a quaternion representing the inverse rotation = inverse(R) = transpose(R)
```
You can use the quaternion to rotate a 3D vector like this:
```cpp
array_1d<double, 3> v; // a 3D vector
Quaternion<double> q; // a quaternion representing a rotation R
array_1d<double, 3> Rv; // the rotated vector
q.RotateVector3(v, Rv); // sets Rv as the rotated vector of v, i.e. R*v
// or you can rotate the input vector itself
q.RotateVector3(v); // now v is R*v
```