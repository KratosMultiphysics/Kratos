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

#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/define.h"
#include "includes/serializer.h"

namespace Kratos
{

namespace RveUtilities
{

inline const double Precision() {
	return 1.0E-8;
}

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
