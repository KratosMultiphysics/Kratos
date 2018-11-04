//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//


// System includes

// External includes

// Project includes
#include "stress_response_definitions.h"
#include "utilities/compare_elements_and_conditions_utility.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

namespace StressResponseDefinitions
{
    TracedStressType ConvertStringToTracedStressType(const std::string& Str)
    {
        const std::map<std::string, TracedStressType> traced_stress_type_map {
            { "FX",  TracedStressType::FX },
            { "FY",  TracedStressType::FY },
            { "FZ",  TracedStressType::FZ },
            { "MX",  TracedStressType::MX },
            { "MY",  TracedStressType::MY },
            { "MZ",  TracedStressType::MZ },
            { "FXX", TracedStressType::FXX },
            { "FXY", TracedStressType::FXY },
            { "FXZ", TracedStressType::FXZ },
            { "FYX", TracedStressType::FYX },
            { "FYY", TracedStressType::FYY },
            { "FYZ", TracedStressType::FYZ },
            { "FZX", TracedStressType::FZX },
            { "FZY", TracedStressType::FZY },
            { "FZZ", TracedStressType::FZZ },
            { "MXX", TracedStressType::MXX },
            { "MXY", TracedStressType::MXY },
            { "MXZ", TracedStressType::MXZ },
            { "MYX", TracedStressType::MYX },
            { "MYY", TracedStressType::MYY },
            { "MYZ", TracedStressType::MYZ },
            { "MZX", TracedStressType::MZX },
            { "MZY", TracedStressType::MZY },
            { "MZZ", TracedStressType::MZZ },
            { "PK2", TracedStressType::PK2 },
        };

        auto stress_type_it = traced_stress_type_map.find(Str);
        if (stress_type_it == traced_stress_type_map.end())
        {
            std::stringstream available_types;
            for(auto& it : traced_stress_type_map)
                available_types << it.first << ",\n";

            KRATOS_ERROR << "Chosen stress type '" <<Str<<"' is not available!" <<
                " Available types are: \n" << available_types.str() << std::endl;
        }

        return stress_type_it->second;
    }

    StressTreatment ConvertStringToStressTreatment(const std::string& Str)
    {
        const std::map<std::string, StressTreatment> stress_treatment_map {
            { "mean", StressTreatment::Mean },
            { "node", StressTreatment::Node },
            { "GP", StressTreatment::GaussPoint }
        };

        auto stress_treatment_it = stress_treatment_map.find(Str);
        if (stress_treatment_it == stress_treatment_map.end())
        {
            std::stringstream available_types;
            for(auto& it : stress_treatment_map)
                available_types << it.first << ",\n";

            KRATOS_ERROR << "Chosen stress treatment '" <<Str<<"' is not available!" <<
                " Available types are: \n" << available_types.str() << std::endl;
        }

        return stress_treatment_it->second;
    }

}  // namespace StressResponseDefinitions.

void StressCalculation::CalculateStressOnNode(Element& rElement,
                        const TracedStressType rTracedStressType,
                        Vector& rOutput,
                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    std::string name_current_element;
    CompareElementsAndConditionsUtility::GetRegisteredName(rElement, name_current_element);

    if(name_current_element == "CrLinearBeamElement3D2N")
        StressCalculation::CalculateStressOnNodeBeam(rElement, rTracedStressType, rOutput, rCurrentProcessInfo);
    else if(name_current_element == "ShellThinElement3D3N")
        KRATOS_ERROR << "Stress calculation on node not yet implemented for " << name_current_element << std::endl;
    else if(name_current_element == "TrussElement3D2N")
        KRATOS_ERROR << "Stress calculation on node not yet implemented for " << name_current_element << std::endl;
    else if(name_current_element == "TrussLinearElement3D2N")
        KRATOS_ERROR << "Stress calculation on node not yet implemented for " << name_current_element << std::endl;
    else
        KRATOS_ERROR << "Stress calculation on node not yet implemented for " << name_current_element << std::endl;

    KRATOS_CATCH("");
}

void StressCalculation::CalculateStressOnGP(Element& rElement,
                        const TracedStressType rTracedStressType,
                        Vector& rOutput,
                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    std::string name_current_element;
    CompareElementsAndConditionsUtility::GetRegisteredName(rElement, name_current_element);

    if(name_current_element == "CrLinearBeamElement3D2N")
        StressCalculation::CalculateStressOnGPBeam(rElement, rTracedStressType, rOutput, rCurrentProcessInfo);
    else if(name_current_element == "ShellThinElement3D3N")
        StressCalculation::CalculateStressOnGPShell(rElement, rTracedStressType, rOutput, rCurrentProcessInfo);
    else if(name_current_element == "TrussElement3D2N")
        StressCalculation::CalculateStressOnGPTruss(rElement, rTracedStressType, rOutput, rCurrentProcessInfo);
    else if(name_current_element == "TrussLinearElement3D2N")
        StressCalculation::CalculateStressOnGPLinearTruss(rElement, rTracedStressType, rOutput, rCurrentProcessInfo);
    else
        KRATOS_ERROR << "Stress calculation on GP not yet implemented for " << name_current_element << std::endl;

    KRATOS_CATCH("");
}

void StressCalculation::CalculateStressOnGPLinearTruss(Element& rElement,
                        const TracedStressType rTracedStressType,
                        Vector& rOutput,
                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const SizeType  GP_num = rElement.GetGeometry().IntegrationPoints().size();
    if (rOutput.size() != GP_num)
        rOutput.resize(GP_num, false);

    switch (rTracedStressType)
    {
        case TracedStressType::FX:
        {
            std::vector< array_1d<double, 3 > > force_vector;
            rElement.CalculateOnIntegrationPoints(FORCE, force_vector, rCurrentProcessInfo);
            for(IndexType i = 0; i < GP_num ; ++i)
                rOutput(i) = force_vector[i][0];
            break;
        }
        default:
            KRATOS_ERROR << "Invalid stress type! Stress type not supported for this element!" << std::endl;
    }

    KRATOS_CATCH("");
}

void StressCalculation::CalculateStressOnGPTruss(Element& rElement,
                        const TracedStressType rTracedStressType,
                        Vector& rOutput,
                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const SizeType  GP_num = (rElement.GetGeometry().IntegrationPoints()).size();
    if (rOutput.size() != GP_num)
        rOutput.resize(GP_num, false);

    switch (rTracedStressType)
    {
        case TracedStressType::FX:
        {
            std::vector< array_1d<double, 3 > > force_vector;
            rElement.CalculateOnIntegrationPoints(FORCE, force_vector, rCurrentProcessInfo);
            for(IndexType i = 0; i < GP_num ; ++i)
                rOutput(i) = force_vector[i][0];
            break;
        }
        case TracedStressType::PK2:
        {
            std::vector<Vector> stress_vector;
            rElement.CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, stress_vector, rCurrentProcessInfo);
            for(IndexType i = 0; i < GP_num ; ++i)
                rOutput(i) = stress_vector[i][0];
            break;
        }
        default:
            KRATOS_ERROR << "Invalid stress type! Stress type not supported for this element!" << std::endl;
    }

    KRATOS_CATCH("");

}

void StressCalculation::CalculateStressOnGPShell(Element& rElement,
                        const TracedStressType rTracedStressType,
                        Vector& rOutput,
                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const SizeType num_gps = rElement.GetGeometry().IntegrationPointsNumber(rElement.GetIntegrationMethod());

    int direction_1 = 0;
    int direction_2 = 0;
    std::vector<Matrix> stress_vector;
    bool stress_is_moment = true;

    switch (rTracedStressType)
    {
        case TracedStressType::MXX:
        {
            direction_1 = 0;
            direction_2 = 0;
            break;
        }
        case TracedStressType::MXY:
        {
            direction_1 = 0;
            direction_2 = 1;
            break;
        }
        case TracedStressType::MXZ:
        {
            direction_1 = 0;
            direction_2 = 2;
            break;
        }
        case TracedStressType::MYX:
        {
            direction_1 = 1;
            direction_2 = 0;
            break;
        }
        case TracedStressType::MYY :
        {
            direction_1 = 1;
            direction_2 = 1;
            break;
        }
        case TracedStressType::MYZ:
        {
            direction_1 = 1;
            direction_2 = 2;
            break;
        }
        case TracedStressType::MZX:
        {
            direction_1 = 2;
            direction_2 = 0;
            break;
        }
        case TracedStressType::MZY:
        {
            direction_1 = 2;
            direction_2 = 1;
            break;
        }
        case TracedStressType::MZZ :
        {
            direction_1 = 2;
            direction_2 = 2;
            break;
        }
        case TracedStressType::FXX :
        {
            direction_1 = 0;
            direction_2 = 0;
            stress_is_moment = false;
            break;
        }
        case TracedStressType::FXY:
        {
            direction_1 = 0;
            direction_2 = 1;
            stress_is_moment = false;
            break;
        }
        case TracedStressType::FXZ:
        {
            direction_1 = 0;
            direction_2 = 2;
            stress_is_moment = false;
            break;
        }
        case TracedStressType::FYX:
        {
            direction_1 = 1;
            direction_2 = 0;
            stress_is_moment = false;
            break;
        }
        case TracedStressType::FYY:
        {
            direction_1 = 1;
            direction_2 = 1;
            stress_is_moment = false;
            break;
        }
        case TracedStressType::FYZ:
        {
            direction_1 = 1;
            direction_2 = 2;
            stress_is_moment = false;
            break;
        }
        case TracedStressType::FZX:
        {
            direction_1 = 2;
            direction_2 = 0;
            stress_is_moment = false;
            break;
        }
        case TracedStressType::FZY:
        {
            direction_1 = 2;
            direction_2 = 1;
            stress_is_moment = false;
            break;
        }
        case TracedStressType::FZZ:
        {
            direction_1 = 2;
            direction_2 = 2;
            stress_is_moment = false;
            break;
        }
        default:
            KRATOS_ERROR << "Invalid stress type! Stress type not supported for this element!" << std::endl;
    }

    if(stress_is_moment)
        rElement.CalculateOnIntegrationPoints(SHELL_MOMENT_GLOBAL, stress_vector, rCurrentProcessInfo);
    else
        rElement.CalculateOnIntegrationPoints(SHELL_FORCE_GLOBAL, stress_vector, rCurrentProcessInfo);

    rOutput.resize(num_gps);
    for(IndexType i = 0; i < num_gps; i++)
    {
        rOutput(i) = stress_vector[i](direction_1, direction_2);
    }

    KRATOS_CATCH("");
}

void StressCalculation::StressCalculation::CalculateStressBeam(Element& rElement,
                        const TracedStressType rTracedStressType,
                        std::vector< array_1d<double, 3 > >& rStressVector,
                        const ProcessInfo& rCurrentProcessInfo,
                        int& rDirection)
{

    rDirection = 0;
    bool stress_is_moment = true;

    switch (rTracedStressType)
    {
        case TracedStressType::MX:
        {
            rDirection = 0;
            break;
        }
        case TracedStressType::MY:
        {
            rDirection = 1;
            break;
        }
        case TracedStressType::MZ:
        {
            rDirection = 2;
            break;
        }
        case TracedStressType::FX:
        {
            rDirection = 0;
            stress_is_moment = false;
            break;
        }
        case TracedStressType::FY:
        {
            rDirection = 1;
            stress_is_moment = false;
            break;
        }
        case TracedStressType::FZ:
        {
            rDirection = 2;
            stress_is_moment = false;
            break;
        }
        default:
            KRATOS_ERROR << "Invalid stress type! Stress type not supported for this element!" << std::endl;
    }

    if(stress_is_moment)
        rElement.CalculateOnIntegrationPoints(MOMENT, rStressVector, rCurrentProcessInfo);
    else
        rElement.CalculateOnIntegrationPoints(FORCE, rStressVector, rCurrentProcessInfo);
}

void StressCalculation::CalculateStressOnGPBeam(Element& rElement,
                        const TracedStressType rTracedStressType,
                        Vector& rOutput,
                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    int direction_1;
    std::vector< array_1d<double, 3 > > stress_vector;
    StressCalculation::CalculateStressBeam(rElement, rTracedStressType,
                                            stress_vector, rCurrentProcessInfo,
                                            direction_1);

    const SizeType GP_num = rElement.GetGeometry().IntegrationPointsNumber(Kratos::GeometryData::GI_GAUSS_3);

    rOutput.resize(GP_num);
    for(IndexType i = 0; i < GP_num ; i++)
    {
        rOutput(i) = stress_vector[i][direction_1];
    }

    KRATOS_CATCH("")
}

void StressCalculation::CalculateStressOnNodeBeam(Element& rElement,
                        const TracedStressType rTracedStressType,
                        Vector& rOutput,
                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    int direction_1;
    std::vector< array_1d<double, 3 > > stress_vector;
    StressCalculation::CalculateStressBeam(rElement, rTracedStressType,
                                            stress_vector, rCurrentProcessInfo,
                                            direction_1);

    rOutput.resize(2);
    rOutput(0) = 2 * stress_vector[0][direction_1] - stress_vector[1][direction_1];
    rOutput(1) = 2 * stress_vector[2][direction_1] - stress_vector[1][direction_1];

    KRATOS_CATCH("")
}

}  // namespace Kratos.

