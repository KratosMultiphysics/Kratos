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
#include<map>


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

    const std::string time_scheme = Settings["time_scheme"].GetString();

    mIntegrationMethod = std::get<0>(GetMethodIterator(time_scheme)->second);

    const SizeType min_buffer_size = GetMinimumBufferSize(time_scheme);
    const SizeType current_buffer_size = rModelPart.GetBufferSize();

    KRATOS_ERROR_IF(current_buffer_size < min_buffer_size) << "Insufficient buffer size: "
        << current_buffer_size << " < " << min_buffer_size << "!" << std::endl;

    // set Generalized-Alpha specific Settings
    if (mIntegrationMethod == generalized_alpha) {
        KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(MESH_ACCELERATION))
            << "The ModelPart does not have the variable \"MESH_ACCELERATION\" as "
            << "nodal-solutionstepvariable!" << std::endl;

        // default (Newmark)
        double alpha_m = 0.0;
        double alpha_f = 0.0;

        // setting mAlphaF and mAlphaM
        if (time_scheme == "generalized_alpha") {
            alpha_m = Settings["alpha_m"].GetDouble();
            alpha_f = Settings["alpha_f"].GetDouble();
        } else if (time_scheme == "bossak") {
            if (Settings.Has("alpha_m")) {
                alpha_m = Settings["alpha_m"].GetDouble();
            }
            else {
                alpha_m = -0.3; // default in Kratos
            }
        }

        mBossakBeta = std::pow((1.0 + alpha_f - alpha_m), 2) * 0.25;
        mBossakGamma = 0.5 + alpha_f - alpha_m;
    }
}

CalculateMeshVelocityUtility::SizeType CalculateMeshVelocityUtility::GetMinimumBufferSize(
    const std::string& rIntegrationMethod)
{
    return std::get<1>(GetMethodIterator(rIntegrationMethod)->second);
}

void CalculateMeshVelocityUtility::CalculateMeshVelocities()
{
    const double delta_time = mrModelPart.GetProcessInfo()[DELTA_TIME];
    KRATOS_ERROR_IF(delta_time < 1.0e-24) << "Detected delta_time = 0! "
        << "check if the time step is created correctly "
        << "for the current time step" << std::endl;

    if ( mIntegrationMethod == bdf1 ||
         mIntegrationMethod == bdf2 ) {
        CalculateMeshVelocitiesBDF(delta_time);
    }
    else {
        CalculateMeshVelocitiesGeneralizedAlpha(delta_time);
    }

    mrModelPart.GetCommunicator().SynchronizeVariable(MESH_VELOCITY);
}

void CalculateMeshVelocityUtility::CalculateMeshVelocitiesBDF(const double DeltaTime)
{
    const int num_local_nodes = mrModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    const auto nodes_begin = mrModelPart.GetCommunicator().LocalMesh().NodesBegin();

    const double coeff = 1 / DeltaTime;

    switch(mIntegrationMethod) {
        case bdf1 : {
            #pragma omp parallel for
            for (int i=0; i<num_local_nodes; i++) {
                const auto it_node  = nodes_begin + i;
                auto& r_mesh_v0       = it_node->FastGetSolutionStepValue(MESH_VELOCITY);
                const auto& r_mesh_u0 = it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT);
                const auto& r_mesh_u1 = it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
                noalias(r_mesh_v0) = r_mesh_u0 - r_mesh_u1;
                r_mesh_v0 *= coeff;
            }
            break;
        }
        case bdf2 : {
            const double c1 = 1.50 * coeff;
            const double c2 = -2.0 * coeff;
            const double c3 = 0.50 * coeff;

            #pragma omp parallel for
            for (int i=0; i<num_local_nodes; i++) {
                const auto it_node  = nodes_begin + i;
                auto& r_mesh_v0 = it_node->FastGetSolutionStepValue(MESH_VELOCITY);
                noalias(r_mesh_v0)  = c1 * it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT);
                noalias(r_mesh_v0) += c2 * it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
                noalias(r_mesh_v0) += c3 * it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT, 2);
            }
            break;
        }
        default : KRATOS_ERROR << "unknown bdf integration-order!" << std::endl;
    }
}

void CalculateMeshVelocityUtility::CalculateMeshVelocitiesGeneralizedAlpha(const double DeltaTime)
{
    const int num_local_nodes = mrModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    const auto nodes_begin = mrModelPart.GetCommunicator().LocalMesh().NodesBegin();
    const double const_u = mBossakGamma / (DeltaTime * mBossakBeta);
    const double const_v = 1.0 - mBossakGamma / mBossakBeta;
    const double const_a = DeltaTime * (1.0 - mBossakGamma / (2.0 * mBossakBeta));

    #pragma omp parallel for
    for (int i=0; i<num_local_nodes; i++) {
        const auto it_node  = nodes_begin + i;
        const auto& r_mesh_u0 = it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT);
        auto&       r_mesh_v0 = it_node->FastGetSolutionStepValue(MESH_VELOCITY);
        auto&       r_mesh_a0 = it_node->FastGetSolutionStepValue(MESH_ACCELERATION);

        const auto& r_mesh_u1 = it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
        const auto& r_mesh_v1 = it_node->FastGetSolutionStepValue(MESH_VELOCITY, 1);
        const auto& r_mesh_a1 = it_node->FastGetSolutionStepValue(MESH_ACCELERATION, 1);

        r_mesh_v0 = const_u * (r_mesh_u0 - r_mesh_u1) + const_v * r_mesh_v1 + const_a * r_mesh_a1;
        r_mesh_a0 = (1.0 / (DeltaTime * mBossakGamma)) * (r_mesh_v0 - r_mesh_v1) - ((1 - mBossakGamma) / mBossakGamma) * r_mesh_a1;
    }

    mrModelPart.GetCommunicator().SynchronizeVariable(MESH_ACCELERATION);
}

const CalculateMeshVelocityUtility::MethodsMapType::const_iterator CalculateMeshVelocityUtility::GetMethodIterator(
    const std::string& rIntegrationMethod)
{
    const auto it_method = msAvailableMethods.find(rIntegrationMethod);

    if (it_method != msAvailableMethods.end()) {
        return it_method;
    }
    else {
        std::stringstream err_msg;

        err_msg << "The requested \"time_scheme\" \"" << rIntegrationMethod
                << "\" is not available!\n"
                << "The following options are available:" << std::endl;

        for (const auto& avail_method : msAvailableMethods) {
            err_msg << "\t" << avail_method.first << "\n";
        }

        KRATOS_ERROR << err_msg.str() << std::endl;
    }
}

// register the available integration-methods
CalculateMeshVelocityUtility::MethodsMapType CalculateMeshVelocityUtility::msAvailableMethods = {
    {"bdf1",              TupleType{bdf1,              2}},
    {"bdf2",              TupleType{bdf2,              3}},
    {"generalized_alpha", TupleType{generalized_alpha, 2}},
    {"bossak",            TupleType{generalized_alpha, 2}},
    {"newmark",           TupleType{generalized_alpha, 2}}
};


}  // namespace Kratos.


