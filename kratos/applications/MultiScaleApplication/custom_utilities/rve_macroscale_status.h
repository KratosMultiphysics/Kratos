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

#if !defined(RVE_MACROSCALE_STATUS_H_INCLUDED)
#define RVE_MACROSCALE_STATUS_H_INCLUDED

#include <limits>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <algorithm>

#include "includes/define.h"
#include "includes/serializer.h"

namespace Kratos
{

    class RveMacroscaleStatus
    {

    public:

		///@name Type Definitions
		///@{

        KRATOS_CLASS_POINTER_DEFINITION( RveMacroscaleStatus );

        typedef double RealType;

        typedef Vector VectorType;
        
        typedef Matrix MatrixType;
        
        typedef size_t IndexType;

		///@}

    public:
        
		///@name Life Cycle
		///@{

        RveMacroscaleStatus()
			: mStrainVector(0)
        {
        }
        
        virtual ~RveMacroscaleStatus()
        {
        }

		///@}

    public:

		///@name Operations
		///@{

        inline const VectorType& GetStrainVector()const
        {
            return mStrainVector;
        }
        
        inline void SetStrainVector(const VectorType& strainVector)
        {
            if(mStrainVector.size() != strainVector.size())
                mStrainVector.resize(strainVector.size(), false);
            noalias(mStrainVector) = strainVector;
        }

        std::string GetInfo()const
        {
            std::stringstream ss;
            ss << "RVE Macroscale Status:" << std::endl;
            
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

		///@}

    private:
    
		///@name Member Variables
		///@{

        VectorType mStrainVector;
        
		///@}

		///@name Private Operations
		///@{

		friend class Serializer;

		virtual void save(Serializer& rSerializer) const
		{
			rSerializer.save("strain", mStrainVector);
		}

		virtual void load(Serializer& rSerializer)
		{
			rSerializer.load("strain", mStrainVector);
		}

		///@}

    };

	///@name Input and output
	///@{

    inline std::ostream & operator << (std::ostream& rOStream, const RveMacroscaleStatus& rThis)
    {
        return rOStream << rThis.GetInfo();
    }

	///@}

} // namespace Kratos

#endif // RVE_MACROSCALE_STATUS_H_INCLUDED
