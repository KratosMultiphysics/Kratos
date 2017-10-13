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
//   Last Modified by:    $Author: Stefano Zaghi $
//   Date:                $Date: 2013-11-04 12:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#if !defined(RVE_MACROSCALE_TEMPERATURE_DATA_H_INCLUDED)
#define RVE_MACROSCALE_TEMPERATURE_DATA_H_INCLUDED

#include <limits>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <algorithm>

#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/define.h"
#include "includes/serializer.h"
#include "rve_geometry_descriptor.h"
#include "rve_macroscale_data.h"

namespace Kratos
{

    class RveMacroscaleTemperatureData : public RveMacroscaleData
    {

    public:

        KRATOS_CLASS_POINTER_DEFINITION( RveMacroscaleTemperatureData );

    public:
        
        RveMacroscaleTemperatureData()
			: RveMacroscaleData()
        {
        }
        
        virtual ~RveMacroscaleTemperatureData()
        {
        }

	public: 

		virtual void SetData(ConstitutiveLaw::Parameters& param_macro,
					 ModelPart& modp_micro,
					 const RveGeometryDescriptor& geomdes)
		{
			mStrainVector = param_macro.GetStrainVector();

			// Reed T_Data
			double alpha(0.0);
			double totalVolume(0.0);
			//double T0 = param_macro.GetMaterialProperties().Has(AMBIENT_TEMPERATURE) ? param_macro.GetMaterialProperties()[AMBIENT_TEMPERATURE] : 0.0;
			//const Vector& N = param_macro.GetShapeFunctionsValues();


			//ProcessInfo& procInfo = modp_micro.GetProcessInfo();
			//size_t strain_size = geomdes.Dimension();
			for (ModelPart::ElementIterator it = modp_micro.ElementsBegin(); it != modp_micro.ElementsEnd(); ++it)
			{
				Element& ielem = *it;
				const Properties& props = ielem.GetProperties();
				double ialpha = props[COEFFICIENT_THERMAL_EXPANSION]; // alpha = CTE
				// Homogenize alpha = CTE
				double ivolume = ielem.GetGeometry().DomainSize();
				alpha += ialpha * ivolume;
				totalVolume += ivolume;
			}

			if (totalVolume > 0.0) 
				alpha /= totalVolume;
			else
				alpha = 0.0;

			if (param_macro.GetProcessInfo()[TIME_STEPS] == 1)
			{
				std::stringstream HCTE;
				std::ofstream mHCTE;
				mHCTE.open("Homogenized_Thermal_Parameters.txt");
				HCTE << " CTE Homogenized = " << std::endl;
				HCTE << "=============================================================================" << std::endl;
				HCTE << alpha << std::endl;
				HCTE << "=============================================================================" << std::endl;
				mHCTE << HCTE.str();
				mHCTE.close();
			}
			
			// TODO: CALL HOMOGENIZE - K = [alpha     0  ;
			//								  0   ; alpha]
		}

    protected:

		friend class Serializer;

		virtual void save(Serializer& rSerializer) const
		{
			KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, RveMacroscaleData);
		}

		virtual void load(Serializer& rSerializer)
		{
			KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, RveMacroscaleData);
		}

    };

} // namespace Kratos

#endif // RVE_MACROSCALE_TEMPERATURE_DATA_H_INCLUDED
