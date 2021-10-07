//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Salva Latorre$
//   Date:                $Date: 2015-10-26 09:56:42 $
//   Revision:            $Revision: 1.5 $
//

// Project includes

// System includes
#include <limits>
#include <iostream>
#include <iomanip>

// External includes
#ifdef _OPENMP
#include <omp.h>
#endif

// Project includes
#include "dem_fem_utilities.h"
#include "custom_conditions/RigidFace.h"

namespace Kratos
{

    DEMFEMUtilities::~DEMFEMUtilities() {}

    void DEMFEMUtilities::CreateRigidFacesFromAllElements(ModelPart& r_model_part, PropertiesType::Pointer pProps) {

        ElementsArrayType& all_elements = r_model_part.Elements();

        for (unsigned int i = 0; i < all_elements.size(); i++) {

            ConditionType::Pointer pCondition;
            ElementsArrayType::iterator pElement = all_elements.ptr_begin() + i;
            pCondition = ConditionType::Pointer(new RigidFace3D( pElement->Id(), pElement->pGetGeometry(), pProps));
            r_model_part.Conditions().push_back(pCondition);
        }
    }

    /// Turn back information as a string.
    std::string DEMFEMUtilities::Info() const {
            return "";
    }

    /// Print information about this object.
    void DEMFEMUtilities::PrintInfo(std::ostream& rOStream) const {}

    /// Print object's data.
    void DEMFEMUtilities::PrintData(std::ostream& rOStream) const {}

/// output stream function
// 	template<std::size_t TDim>
// 	inline std::ostream& operator << (std::ostream& rOStream)
// 	{
// 		rThis.PrintInfo(rOStream);
// 		rOStream << std::endl;
// 		rThis.PrintData(rOStream);
//
// 		return rOStream;
// 	}
///@}

} // namespace Kratos
