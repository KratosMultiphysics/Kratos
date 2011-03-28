/*
==============================================================================
KratosOpenCLApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Farshid Mossaiby
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
mossaiby@yahoo.com
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: mossaiby $
//   Date:                $Date: 2011-03-20 00:42:01 $
//   Revision:            $Revision: 1.00 $
//
//

//
// opencl_common.cl
//
// OpenCL common kernels and functions


// Include guard, we do not need this more than once

#ifndef KRATOS_OPENCL_COMMON_CL_INCLUDED

#define KRATOS_OPENCL_COMMON_CL_INCLUDED


#include "opencl_enable_fp64.cl"


//
// Used types

typedef unsigned int IndexType;

typedef double3 VectorType;

typedef double ValueType;

typedef double16 EdgeType;


//
// OpenCL 1.0 adjustment

#if KRATOS_OCL_VERSION < 110

typedef double4 double3;

#endif


//
// Fast math macros

// Currently these are not supported on GPUs

#ifdef __CPU__

	#define KRATOS_OCL_DIVIDE(x, y)			native_divide(x, y)
	#define KRATOS_OCL_RECIP(x)				native_recip(x)
	#define KRATOS_OCL_SQRT(x)				native_sqrt(x)

#else

	#define KRATOS_OCL_DIVIDE(x, y)			((x) / (y))
	#define KRATOS_OCL_RECIP(x)				(1.00 / (x))
	#define KRATOS_OCL_SQRT(x)				sqrt(x)

#endif

//
// OpenCL defines length() as length of the vector, so if we use double4 instead of double3, we have to take care of this

inline double length3(double4 x)
{
	double4 t = x;
	t.s3 = 0.00;

	return KRATOS_OCL_SQRT(dot(t, t));
}

#if KRATOS_OCL_VERSION < 110

	#define KRATOS_OCL_LENGTH3(x)			length3(x)
	#define KRATOS_OCL_VECTOR3(x, y, z)		((VectorType)(x, y, z, 0.00))

#else

	#define KRATOS_OCL_LENGTH3(x)			length(x)
	#define KRATOS_OCL_VECTOR3(x, y, z)		((VectorType)(x, y, z))

#endif


#endif  // KRATOS_OPENCL_COMMON_CL_INCLUDED
