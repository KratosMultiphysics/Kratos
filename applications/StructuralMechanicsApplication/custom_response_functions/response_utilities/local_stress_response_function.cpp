// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//

// System includes

// External includes

// Project includes
#include "local_stress_response_function.h"

namespace Kratos {

LocalStressResponseFunction::LocalStressResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
: mrModelPart(rModelPart), mResponseSettings(ResponseSettings)
{
    // Get traced element
    const int id_of_traced_element = ResponseSettings["traced_element_id"].GetInt();
    mpTracedElement = rModelPart.pGetElement(id_of_traced_element);

    // Tell traced element the stress type
    mTracedStressType = StressResponseDefinitions::ConvertStringToTracedStressType(ResponseSettings["stress_type"].GetString());
    mpTracedElement->SetValue(TRACED_STRESS_TYPE, static_cast<int>(mTracedStressType) );

    // Get info how and where to treat the stress
    mStressTreatment = StressResponseDefinitions::ConvertStringToStressTreatment( ResponseSettings["stress_treatment"].GetString() );

    if(mStressTreatment == StressTreatment::GaussPoint || mStressTreatment == StressTreatment::Node)
    {
        mIdOfLocation = ResponseSettings["stress_location"].GetInt();
        KRATOS_ERROR_IF(mIdOfLocation < 1) << "Choose a 'stress_location' > 0. Specified 'stress_location': " << mIdOfLocation << std::endl;
    }

    // Check if there are at primal elements, because the primal state is required
    ProcessInfo &r_current_process_info = mrModelPart.GetProcessInfo();
    KRATOS_ERROR_IF( r_current_process_info.Has(IS_ADJOINT) && r_current_process_info[IS_ADJOINT] )
            << "LocalStressResponseFunction: Can not use adjoint model part!" << std::endl;
}

double LocalStressResponseFunction::CalculateValue()
{
    KRATOS_TRY;

    double stress_value = 0.0;

    if(mStressTreatment == StressTreatment::Mean)
        stress_value = CalculateMeanElementStress();
    else if (mStressTreatment == StressTreatment::GaussPoint)
        stress_value = CalculateGaussPointStress();
    else if (mStressTreatment == StressTreatment::Node)
        stress_value = CalculateNodeStress();
    return stress_value;

    KRATOS_CATCH("");
}

double LocalStressResponseFunction::CalculateMeanElementStress()
{
    double stress_value = 0.0;

    Vector element_stress;

    StressCalculation::CalculateStressOnGP(*mpTracedElement, mTracedStressType, element_stress, mrModelPart.GetProcessInfo());

    const SizeType stress_vec_size = element_stress.size();

    for(IndexType i = 0; i < stress_vec_size; ++i)
        stress_value += element_stress[i];

    stress_value /= stress_vec_size;

    return stress_value;
}

double LocalStressResponseFunction::CalculateGaussPointStress()
{
    Vector element_stress;

    StressCalculation::CalculateStressOnGP(*mpTracedElement, mTracedStressType, element_stress, mrModelPart.GetProcessInfo());

    const SizeType stress_vec_size = element_stress.size();

    if(stress_vec_size >= mIdOfLocation)
        return element_stress[mIdOfLocation - 1];

    KRATOS_ERROR << "Chosen Gauss-Point is not available. Chose 'stress_location' between 1 and " <<
                    stress_vec_size  << "!"<< std::endl;
}

double LocalStressResponseFunction::CalculateNodeStress()
{
    Vector element_stress;

    StressCalculation::CalculateStressOnGP(*mpTracedElement, mTracedStressType, element_stress, mrModelPart.GetProcessInfo());

    const SizeType num_ele_nodes = mpTracedElement->GetGeometry().PointsNumber();

    if(num_ele_nodes >= mIdOfLocation)
        return element_stress[mIdOfLocation - 1];

    KRATOS_ERROR << "Chosen Node is not available. The element has only " <<
                    num_ele_nodes  << " nodes."<< std::endl;

}

} // namespace Kratos.

