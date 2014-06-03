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
//   Date:                $Date: 2013-11-05 19:00:00 $
//   Revision:            $Revision: 1.00 $
//
//


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "rve_condition_base.h"
#include "utilities/math_utils.h"
#include "geometries/point_3d.h"
#include "multiscale_application.h"

namespace Kratos
{
	
RveConditionBase::RveConditionBase(IndexType NewId, GeometryType::Pointer pGeometry)
	: Condition(NewId, pGeometry)
	, mpMacroscaleStatus(RveMacroscaleStatus::Pointer())
{
}

RveConditionBase::RveConditionBase(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
	: Condition(NewId, pGeometry, pProperties)
	, mpMacroscaleStatus(RveMacroscaleStatus::Pointer())
{
}
	
RveConditionBase::RveConditionBase(const RveConditionBase& rOther)
	: Condition(rOther)
	, mpMacroscaleStatus(rOther.mpMacroscaleStatus)
{
}

RveConditionBase::~RveConditionBase()
{
}

int RveConditionBase::Check(const ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY
	if(mpMacroscaleStatus == NULL)
		KRATOS_ERROR(std::logic_error, "RveConditionBase - Missing the mpMacroscaleStatus", "");
	return 0;
	KRATOS_CATCH("")
}

bool RveConditionBase::GetMacroStrainVector(VectorType& strainVector)
{
	if(mpMacroscaleStatus == NULL) return false;
	
	const VectorType& macroStrainVector = mpMacroscaleStatus->GetStrainVector();
	if(strainVector.size() != macroStrainVector.size())
		strainVector.resize(macroStrainVector.size(), false);
	noalias(strainVector) = macroStrainVector;
	
	return true;
}

}


