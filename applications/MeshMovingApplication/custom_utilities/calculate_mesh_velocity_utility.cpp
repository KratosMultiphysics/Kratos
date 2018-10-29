//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//                   Andreas Winterstein
//


// System includes


// External includes


// Project includes
#include "includes/mesh_moving_variables.h"
#include "calculate_mesh_velocity_utility.h"


namespace Kratos
{

CalculateMeshVelocityUtility::CalculateMeshVelocityUtility(ModelPart& rModelPart,
                                                           Parameters Settings)
    : mrModelPart(rModelPart)
{
    KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(MESH_DISPLACEMENT))
        << "The ModelPart does not have the variable \"MESH_DISPLACEMENT\" as "
        << "nodal-solutionstepvariable!" << std::endl;
    KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(MESH_VELOCITY))
        << "The ModelPart does not have the variable \"MESH_VELOCITY\" as "
        << "nodal-solutionstepvariable!" << std::endl;
    KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(MESH_ACCELERATION))
        << "The ModelPart does not have the variable \"MESH_ACCELERATION\" as "
        << "nodal-solutionstepvariable!" << std::endl;

    const std::string integration_method = Settings["integration_method"].GetString();

    const std::map<std::string, CalculateMeshVelocityUtility::IntegrationMethod>
    available_integration_methods {
        {"bdf1",              bdf1},
        {"bdf2",              bdf2},
        {"bdf3",              bdf3},
        {"bdf4",              bdf4},
        {"bdf5",              bdf5},
        {"bdf6",              bdf6},
        {"generalized_alpha", generalized_alpha},
        {"bossak",            generalized_alpha},
        {"newmark",           generalized_alpha}
    };

    if (available_integration_methods.find(integration_method) != available_integration_methods.end()) {
        mIntegrationMethod = available_integration_methods.at(integration_method);
    }
    else {
        std::stringstream err_msg;

        err_msg << "The requested integration-method \"" << integration_method
                << "\" is not available!"
                << "The following methods are available:" << std::endl;

        for (const auto& avail_method : available_integration_methods) {
            err_msg << "\t" << avail_method.first << "\n";
        }

        KRATOS_ERROR << err_msg.str() << std::endl;
    }

    const SizeType min_buffer_size = GetMinimumBufferSize(integration_method);
    const SizeType current_buffer_size = rModelPart.GetBufferSize();

    KRATOS_ERROR_IF(current_buffer_size < min_buffer_size) << "Insufficient buffer size: "
        << current_buffer_size << " < " << min_buffer_size << "!" << std::endl;

    // setting mAlphaF and mAlphaM
    if (integration_method == "generalized_alpha") {

    } else if (integration_method == "bossak") {

    } else if (integration_method == "newmark") {

    }
}

CalculateMeshVelocityUtility::SizeType CalculateMeshVelocityUtility::GetMinimumBufferSize(const std::string& rIntegrationMethod)
{
    const std::map<std::string, SizeType> required_buffer_sizes {
        {"bdf1",              2},
        {"bdf2",              3},
        {"bdf3",              4},
        {"bdf4",              5},
        {"bdf5",              6},
        {"bdf6",              7},
        {"generalized_alpha", 2},
        {"bossak",            2},
        {"newmark",           2}
    };
    return required_buffer_sizes.at(rIntegrationMethod);
}

void CalculateMeshVelocityUtility::CalculateMeshVelocities()
{
    const double delta_time = mrModelPart.GetProcessInfo()[DELTA_TIME];
    KRATOS_ERROR_IF(delta_time < 1.0e-24) << "Detected delta_time = 0! "
        << "check if the time step is created correctly "
        << "for the current time step" << std::endl;

    if ( mIntegrationMethod == bdf1 ||
         mIntegrationMethod == bdf2 ||
         mIntegrationMethod == bdf3 ||
         mIntegrationMethod == bdf4 ||
         mIntegrationMethod == bdf5 ||
         mIntegrationMethod == bdf6 ) {
        CalculateMeshVelocitiesBDF(delta_time);
    }
    else {
        CalculateMeshVelocitiesGeneralizedAlpha(delta_time);
    }
}

void CalculateMeshVelocityUtility::CalculateMeshVelocitiesBDF(const double DeltaTime)
{
    switch(mIntegrationMethod) {
        case bdf1 : {

            break;
        }
        case bdf2 : {

            break;
        }
        case bdf3 : {

            break;
        }
        case bdf4 : {

            break;
        }
        case bdf5 : {

            break;
        }
        case bdf6 : {

            break;
        }
        default : KRATOS_ERROR << "unknown integration-order!" << std::endl;
    }
}

void CalculateMeshVelocityUtility::CalculateMeshVelocitiesGeneralizedAlpha(const double DeltaTime)
{

}


}  // namespace Kratos.


