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

#if !defined(RVE_MACROSCALE_DATA_H_INCLUDED)
#define RVE_MACROSCALE_DATA_H_INCLUDED

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
#include "multiscale_application_variables.h"

namespace Kratos
{

    class RveMacroscaleData
    {

    public:

        KRATOS_CLASS_POINTER_DEFINITION( RveMacroscaleData );
        typedef double RealType;
        typedef Vector VectorType;
        typedef Matrix MatrixType;
        typedef size_t IndexType;

    public:

        RveMacroscaleData()
			: mStrainVector(0)
			, mMean_Temp(0.0)
			, mhomogen_alpha(0.0)
			, mHasTemperature(false)
			, mHasTemperatureChecked(false)
        {
        }

        virtual ~RveMacroscaleData()
        {
        }

    public:

		virtual void SetData(ConstitutiveLaw::Parameters& param_macro,
					 ModelPart& modp_micro,
					 const RveGeometryDescriptor& geomdes)
		{
			mStrainVector = param_macro.GetStrainVector();

			/**
			we can get also other stuff like:
			TEMPERATURE on nodes (geom[i].FastGet,....) and interpolate them with (params.GetShapeFunc...)
			*/

			// BEGIN MOD STEFANO
			// Reed T_Data
			if(!mHasTemperatureChecked)
			{
				mHasTemperature = true;
				const Element::GeometryType& geom_macro = param_macro.GetElementGeometry();
				for (size_t i = 0; i < geom_macro.size(); i++)
				{
					if(!geom_macro[i].HasDofFor(TEMPERATURE))
					{
						mHasTemperature = false;
						break;
					}
				}
				mHasTemperatureChecked = true;
			}
			if(mHasTemperature)
			{
				double alpha(0.0);
				double totalVolume(0.0);
				double iterpolated_temp = 0.0;
				double T0 = param_macro.GetMaterialProperties().Has(AMBIENT_TEMPERATURE) ? param_macro.GetMaterialProperties()[AMBIENT_TEMPERATURE] : 0.0;
				const Vector& N = param_macro.GetShapeFunctionsValues();

				const Element::GeometryType& geom_macro = param_macro.GetElementGeometry();
				for (size_t i = 0; i < geom_macro.size(); i++)
				{
					iterpolated_temp += geom_macro[i].FastGetSolutionStepValue(TEMPERATURE) * N[i];
				}

				mMean_Temp = iterpolated_temp - T0;

				// Calcutate alpha only for homogenize this properties!! Change it if alpha is function of temperature
				for (ModelPart::ElementIterator it = modp_micro.ElementsBegin(); it != modp_micro.ElementsEnd(); ++it)
				{
					Element& ielem = *it;
					const Properties& props = ielem.GetProperties();
					double ialpha = props.Has(COEFFICIENT_THERMAL_EXPANSION) ? props[COEFFICIENT_THERMAL_EXPANSION] : 0.0; // alpha = CTE
					// Homogenize alpha = CTE
					double ivolume = ielem.GetGeometry().DomainSize();
					mhomogen_alpha += ialpha * ivolume;
					totalVolume += ivolume;
				}

				if (totalVolume > 0.0)
					mhomogen_alpha /= totalVolume;
				else
					mhomogen_alpha = 0.0;
			}
			// END MOD STEFANO

		}

        inline const VectorType& StrainVector()const{return mStrainVector;}
		inline VectorType&       StrainVector()     {return mStrainVector;}

		inline const VectorType& StrainVectorOld()const{return mStrainVectorOld;}
		inline VectorType&       StrainVectorOld()     {return mStrainVectorOld;}

		// BEGIN MOD STEFANO
		inline double Mean_Temp()    const{ return mMean_Temp; }
		inline double homogen_alpha()const{ return mhomogen_alpha; }
		// END MOD STEFANO

        std::string GetInfo()const
        {
            std::stringstream ss;
            ss << "RVE Macroscale Data:" << std::endl;

            ss << "StrainVector: [" << mStrainVector.size() << "] (";
            for(size_t i = 0; i < mStrainVector.size(); i++) {
                ss << mStrainVector(i);
                if(i != mStrainVector.size() - 1)
                    ss << ",";
            }
            ss << ")" << std::endl;

            ss << std::endl;
            return ss.str();
        }

    protected:

        VectorType mStrainVector;
		VectorType mStrainVectorOld;

		// BEGIN MOD STEFANO
		double mMean_Temp;
		double mhomogen_alpha;
		// END MOD STEFANO

		bool mHasTemperature;
		bool mHasTemperatureChecked;

		friend class Serializer;

		virtual void save(Serializer& rSerializer) const
		{
			rSerializer.save("strain", mStrainVector);
		}

		virtual void load(Serializer& rSerializer)
		{
			rSerializer.load("strain", mStrainVector);
		}

    };

    inline std::ostream & operator << (std::ostream& rOStream, const RveMacroscaleData& rThis)
    {
        return rOStream << rThis.GetInfo();
    }

} // namespace Kratos

#endif // RVE_MACROSCALE_DATA_H_INCLUDED
