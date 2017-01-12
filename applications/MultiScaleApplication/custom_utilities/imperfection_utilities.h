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
//   Date:                $Date: 2015-14-08 12:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#if !defined(IMPERFECTION_UTILITIES_H_INCLUDED)
#define IMPERFECTION_UTILITIES_H_INCLUDED

#include "includes/model_part.h"

namespace Kratos
{

namespace ImperfectionUtilties
{

	static double CalculateRandomImperfectionScaleFactor(
		const Element::GeometryType& rElementGeometry, 
		const Vector& rShapeFunctionsValues)
	{
		double impf = 1.0;
		bool do_imperf = true;
		for(unsigned int i=0; i<rElementGeometry.PointsNumber(); i++) {
			if(!rElementGeometry[i].SolutionStepsDataHas(RANDOM_IMPERFECTION_FACTOR)) {
				do_imperf = false;
				break;
			}
		}
		if(do_imperf) {
			double noiseval=0.0;
			for(unsigned int i=0; i<rElementGeometry.PointsNumber(); i++) {
				noiseval += rShapeFunctionsValues[i]*(rElementGeometry[i].FastGetSolutionStepValue(RANDOM_IMPERFECTION_FACTOR));
				/*double inoise = rElementGeometry[i].FastGetSolutionStepValue(RANDOM_IMPERFECTION_FACTOR);
				if(std::abs(inoise) > std::abs(noiseval))
					noiseval = inoise;*/
				/*if(std::abs(inoise) < 1.0e-9)
				{
					noiseval = 0.0;
					break;
				}*/
			}
			impf = 1.0-noiseval;
		}
		return impf;
	}

}

} // namespace Kratos

#endif // IMPERFECTION_UTILITIES_H_INCLUDED
