/*
==============================================================================
KratosMultiScaleApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-11-04 12:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#if !defined(RVE_UTILITIES_H_INCLUDED)
#define RVE_UTILITIES_H_INCLUDED

#include <limits>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <ctime>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/define.h"
#include "includes/serializer.h"


namespace Kratos
{

namespace RveUtilities
{

// misc stuff...

inline const double Precision() {
	return 1.0E-10;
}


class RveTimer
{
public:
	RveTimer():t0(0.0),t1(0.0){}
	inline double get_time() {
#ifndef _OPENMP
	return std::clock()/static_cast<double>(CLOCKS_PER_SEC);
#else
	return omp_get_wtime();
#endif
	}
	inline void start() { t0 = get_time(); }
	inline void stop() { t1 = get_time(); }
	inline double value() { return t1-t0; }
private:
	double t0;
	double t1;
};

// model part utilties

template<class TModelPart>
inline size_t CalculateTotalNumberOfDofs(TModelPart& mp)
{
	size_t n(0);
	for (typename TModelPart::NodeIterator node_iter = mp.NodesBegin(); node_iter != mp.NodesEnd(); ++node_iter)
	{
		typename TModelPart::NodeType& iNode = *node_iter;
		n += iNode.GetDofs().size();
	}
	return n;
}

template<class TModelPart>
inline void SaveSolutionVector(TModelPart& mp, Vector& U)
{
	size_t counter = 0;
	for (typename TModelPart::NodeIterator node_iter = mp.NodesBegin(); 
		 node_iter != mp.NodesEnd(); 
		 ++node_iter)
	{
		typename TModelPart::NodeType& iNode = *node_iter;
		for(typename TModelPart::NodeType::DofsContainerType::iterator dof_iter = iNode.GetDofs().begin(); 
			dof_iter != iNode.GetDofs().end(); 
			++dof_iter)
		{
			typename TModelPart::DofType& iDof = *dof_iter;
			U(counter++) = iDof.GetSolutionStepValue();
		}
	}
}

template<class TModelPart>
inline void RestoreSolutionVector(TModelPart& mp, const Vector& U)
{
	size_t counter = 0;
	for (typename TModelPart::NodeIterator node_iter = mp.NodesBegin(); 
		 node_iter != mp.NodesEnd(); 
		 ++node_iter)
	{
		typename TModelPart::NodeType& iNode = *node_iter;
		for(typename TModelPart::NodeType::DofsContainerType::iterator dof_iter = iNode.GetDofs().begin(); 
			dof_iter != iNode.GetDofs().end(); 
			++dof_iter)
		{
			typename TModelPart::DofType& iDof = *dof_iter;
			iDof.GetSolutionStepValue() = U(counter++);
		}
	}
}


// utilities for geometry descriptor object


struct RveBoundarySortXFunctor_2DGeneric
{
	RveBoundarySortXFunctor_2DGeneric(const Matrix& A)
		: mA(A) {}

	bool operator()(Node<3>::Pointer a, Node<3>::Pointer b) {
		double ax = mA(0,0)*a->X0() + mA(0,1)*a->Y0();
		double bx = mA(0,0)*b->X0() + mA(0,1)*b->Y0();
        return ax < bx;
    }

private:
	const Matrix& mA;
};

struct RveBoundarySortYFunctor_2DGeneric
{
	RveBoundarySortYFunctor_2DGeneric(const Matrix& A)
		: mA(A) {}

	bool operator()(Node<3>::Pointer a, Node<3>::Pointer b) {
		double ay = mA(1,0)*a->X0() + mA(1,1)*a->Y0();
		double by = mA(1,0)*b->X0() + mA(1,1)*b->Y0();
        return ay < by;
    }

private:
	const Matrix& mA;
};

struct RveBoundarySortXFunctor_3DGeneric
{
	RveBoundarySortXFunctor_3DGeneric(const Matrix& A)
		: mA(A) {}

	bool operator()(Node<3>::Pointer a, Node<3>::Pointer b) {
		double ax = mA(0,0)*a->X0() + mA(0,1)*a->Y0() + mA(0,2)*a->Z0();
		double bx = mA(0,0)*b->X0() + mA(0,1)*b->Y0() + mA(0,2)*b->Z0();
        return ax < bx;
    }

private:
	const Matrix& mA;
};

struct RveBoundarySortYFunctor_3DGeneric
{
	RveBoundarySortYFunctor_3DGeneric(const Matrix& A)
		: mA(A) {}

	bool operator()(Node<3>::Pointer a, Node<3>::Pointer b) {
		double ay = mA(1,0)*a->X0() + mA(1,1)*a->Y0() + mA(1,2)*a->Z0();
		double by = mA(1,0)*b->X0() + mA(1,1)*b->Y0() + mA(1,2)*b->Z0();
        return ay < by;
    }

private:
	const Matrix& mA;
};

struct RveBoundarySortZFunctor_3DGeneric
{
	RveBoundarySortZFunctor_3DGeneric(const Matrix& A)
		: mA(A) {}

	bool operator()(Node<3>::Pointer a, Node<3>::Pointer b) {
		double az = mA(2,0)*a->X0() + mA(2,1)*a->Y0() + mA(2,2)*a->Z0();
		double bz = mA(2,0)*b->X0() + mA(2,1)*b->Y0() + mA(2,2)*b->Z0();
        return az < bz;
    }

private:
	const Matrix& mA;
};

struct RveBoundarySortXYFunctor_3DGeneric
{
	RveBoundarySortXYFunctor_3DGeneric(const Matrix& A)
		: mA(A) {}

	bool operator()(Node<3>::Pointer a, Node<3>::Pointer b) {
		double ax = mA(0,0)*a->X0() + mA(0,1)*a->Y0() + mA(0,2)*a->Z0();
		double bx = mA(0,0)*b->X0() + mA(0,1)*b->Y0() + mA(0,2)*b->Z0();
        double ay = mA(1,0)*a->X0() + mA(1,1)*a->Y0() + mA(1,2)*a->Z0();
		double by = mA(1,0)*b->X0() + mA(1,1)*b->Y0() + mA(1,2)*b->Z0();
		if(ay < by)
			return true;
		else if(ay > by)
			return false;
		else
			return ax < bx;
    }

private:
	const Matrix& mA;
};

struct RveBoundarySortXZFunctor_3DGeneric
{
	RveBoundarySortXZFunctor_3DGeneric(const Matrix& A)
		: mA(A) {}

	bool operator()(Node<3>::Pointer a, Node<3>::Pointer b) {
		double ax = mA(0,0)*a->X0() + mA(0,1)*a->Y0() + mA(0,2)*a->Z0();
		double bx = mA(0,0)*b->X0() + mA(0,1)*b->Y0() + mA(0,2)*b->Z0();
        double az = mA(2,0)*a->X0() + mA(2,1)*a->Y0() + mA(2,2)*a->Z0();
		double bz = mA(2,0)*b->X0() + mA(2,1)*b->Y0() + mA(2,2)*b->Z0();
		if(az < bz)
			return true;
		else if(az > bz)
			return false;
		else
			return ax < bx;
    }

private:
	const Matrix& mA;
};

struct RveBoundarySortYZFunctor_3DGeneric
{
	RveBoundarySortYZFunctor_3DGeneric(const Matrix& A)
		: mA(A) {}

	bool operator()(Node<3>::Pointer a, Node<3>::Pointer b) {
		double ay = mA(1,0)*a->X0() + mA(1,1)*a->Y0() + mA(1,2)*a->Z0();
		double by = mA(1,0)*b->X0() + mA(1,1)*b->Y0() + mA(1,2)*b->Z0();
        double az = mA(2,0)*a->X0() + mA(2,1)*a->Y0() + mA(2,2)*a->Z0();
		double bz = mA(2,0)*b->X0() + mA(2,1)*b->Y0() + mA(2,2)*b->Z0();
		if(az < bz)
			return true;
		else if(az > bz)
			return false;
		else
			return ay < by;
    }

private:
	const Matrix& mA;
};

struct RveBoundarySortXFunctor
{
    bool operator()(Node<3>::Pointer a, Node<3>::Pointer b) {
        return a->X0() < b->X0();
    }
};

struct RveBoundarySortYFunctor
{
    bool operator()(Node<3>::Pointer a, Node<3>::Pointer b) {
        return a->Y0() < b->Y0();
    }
};

struct RveBoundarySortZFunctor
{
    bool operator()(Node<3>::Pointer a, Node<3>::Pointer b) {
        return a->Z0() < b->Z0();
    }
};

struct OneToOneMasterSlaveMap
{
	OneToOneMasterSlaveMap()
        : Master(), Slave()
	{
	}
	OneToOneMasterSlaveMap(Node<3>::Pointer m, Node<3>::Pointer s)
		: Master(m), Slave(s)
	{
	}
	Node<3>::Pointer Master;
	Node<3>::Pointer Slave;
};

struct TwoToOneMasterSlaveMap
{
	typedef double RealType;
	TwoToOneMasterSlaveMap()
        : Master1(), Master2(), Slave()
	{
	}
	TwoToOneMasterSlaveMap(Node<3>::Pointer m1, Node<3>::Pointer m2, Node<3>::Pointer s)
		: Master1(m1), Master2(m2), Slave(s)
		, C1(0.5), C2(0.5)
	{
	}
	TwoToOneMasterSlaveMap(Node<3>::Pointer m1, Node<3>::Pointer m2, Node<3>::Pointer s,
		                   RealType c1, RealType c2)
		: Master1(m1), Master2(m2), Slave(s)
		, C1(c1), C2(c2)
	{
	}
	Node<3>::Pointer Master1;
	Node<3>::Pointer Master2;
	Node<3>::Pointer Slave;
	RealType C1;
	RealType C2;
};

struct ThreeToOneMasterSlaveMap
{
	typedef double RealType;
	ThreeToOneMasterSlaveMap()
        : Master1(), Master2(), Master3(), Slave()
	{
	}
	ThreeToOneMasterSlaveMap(Node<3>::Pointer m1, Node<3>::Pointer m2, Node<3>::Pointer m3, Node<3>::Pointer s)
		: Master1(m1), Master2(m2), Master3(m3), Slave(s)
		, C1(1.0/3.0), C2(1.0/3.0), C3(1.0/3.0)
	{
	}
	ThreeToOneMasterSlaveMap(Node<3>::Pointer m1, Node<3>::Pointer m2, Node<3>::Pointer m3, Node<3>::Pointer s,
		                     RealType c1, RealType c2, RealType c3)
		: Master1(m1), Master2(m2), Master3(m3), Slave(s)
		, C1(c1), C2(c2), C3(c3)
	{
	}
	Node<3>::Pointer Master1;
	Node<3>::Pointer Master2;
	Node<3>::Pointer Master3;
	Node<3>::Pointer Slave;
	RealType C1;
	RealType C2;
	RealType C3;
};

}

} // namespace Kratos

#endif // RVE_UTILITIES_H_INCLUDED
