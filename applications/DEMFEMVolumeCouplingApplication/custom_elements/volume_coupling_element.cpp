// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//                   Alejandro Cornejo Velazquez
//

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "utilities/atomic_utilities.h"

// Application includes
#include "custom_elements/base_solid_element.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{


double VolumeCouplingElement::GetIntegrationWeight(
    const GeometryType::IntegrationPointsArrayType& rThisIntegrationPoints,
    const IndexType PointNumber,
    const double detJ
    ) const 
{

    auto N = row(this->GetGeometry().ShapeFunctionValues(this->GetIntegrationMethod()), PointNumber);
    double interpolated_coupling_weight_at_int_point=0; 
    for (int i=0, i < this->GetGeometry().GetNodes(), i++)
    {
       interpolated_coupling_weight_at_int_point += N[i]*this->GetNodes()[i].GetSolutionStepValue(KM.NODAL_COUPLING_WEIGHT)
    }
    return (1-interpolated_coupling_weight_at_int_point) * TotalLagrangian::GetIntegrationWeight(rThisIntegrationPoints,PointNumber,detJ);
}


} // Namespace Kratos
