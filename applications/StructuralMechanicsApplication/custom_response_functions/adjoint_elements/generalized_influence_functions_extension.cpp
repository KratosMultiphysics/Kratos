//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:   Martin Fusseder, https://github.com/MFusseder
//
//


// Project includes
#include "generalized_influence_functions_extension.h"
#include "custom_response_functions/response_utilities/element_finite_difference_utility.h"
#include "structural_mechanics_application_variables.h"


namespace Kratos
{
    GeneralizedInfluenceFunctionsExtension::GeneralizedInfluenceFunctionsExtension()
    {
    }

    GeneralizedInfluenceFunctionsExtension::~GeneralizedInfluenceFunctionsExtension()
    {
    }

    void GeneralizedInfluenceFunctionsExtension::CalculatePseudoQuantityOnIntegrationPoints(Element& rElement, const Variable<array_1d<double, 3>>& rPseudoQuantityVariable,
                                            std::vector< array_1d<double, 3> >& rOutput, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const double delta = 1e-6;

        if (rPseudoQuantityVariable == PSEUDO_MOMENT)
        {
            ElementFiniteDifferenceUtility::CalculateIntegrationPointsResultsDerivative(rElement, MOMENT,
                I22, delta, rOutput, rCurrentProcessInfo);
        }
        else if (rPseudoQuantityVariable == PSEUDO_FORCE)
        {
            ElementFiniteDifferenceUtility::CalculateIntegrationPointsResultsDerivative(rElement, FORCE,
                I22, delta, rOutput, rCurrentProcessInfo);
        }
        else
            KRATOS_ERROR << "It is possible to provide a pseudo quantity for: " << rPseudoQuantityVariable.Name() << "!" << std::endl;

        KRATOS_CATCH("");
    }

};  // namespace Kratos.

