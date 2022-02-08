//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Pablo Becker
//

// System includes

// External includes

// Project includes
#include "includes/convection_diffusion_settings.h"
#include "includes/kratos_components.h"

namespace Kratos {

void ConvectionDiffusionSettings::save(Serializer& rSerializer) const
{
    KRATOS_TRY

    // Save the is defined bool flag for each variable
    rSerializer.save("mis_defined_DensityVar",mis_defined_DensityVar);
    rSerializer.save("mis_defined_DiffusionVar",mis_defined_DiffusionVar);
    rSerializer.save("mis_defined_UnknownVar",mis_defined_UnknownVar);
    rSerializer.save("mis_defined_VolumeSourceVar",mis_defined_VolumeSourceVar);
    rSerializer.save("mis_defined_SurfaceSourceVar",mis_defined_SurfaceSourceVar);
    rSerializer.save("mis_defined_ProjectionVar",mis_defined_ProjectionVar);
    rSerializer.save("mis_defined_ConvectionVar",mis_defined_ConvectionVar);
    rSerializer.save("mis_defined_GradientVar",mis_defined_GradientVar);
    rSerializer.save("mis_defined_MeshVelocityVar",mis_defined_MeshVelocityVar);
    rSerializer.save("mis_defined_TransferCoefficientVar",mis_defined_TransferCoefficientVar);
    rSerializer.save("mis_defined_VelocityVar",mis_defined_VelocityVar);
    rSerializer.save("mis_defined_SpecificHeatVar",mis_defined_SpecificHeatVar);
    rSerializer.save("mis_defined_ReactionVar",mis_defined_ReactionVar);
    rSerializer.save("mIsDefinedReactionGradientVar", mIsDefinedReactionGradientVar);

    // Save the variable names
    // Note that the variable class save method only saves the name of the variables
    if (mpDensityVar != nullptr && mis_defined_DensityVar) {
        rSerializer.save("DensityVarName",mpDensityVar);
    }
    if (mpDiffusionVar != nullptr && mis_defined_DiffusionVar) {
        rSerializer.save("DiffusionVarName",mpDiffusionVar);
    }
    if (mpUnknownVar != nullptr && mis_defined_UnknownVar) {
        rSerializer.save("UnknownVarName",mpUnknownVar);
    }
    if (mpVolumeSourceVar != nullptr && mis_defined_VolumeSourceVar) {
        rSerializer.save("VolumeSourceVarName",mpVolumeSourceVar);
    }
    if (mpSurfaceSourceVar != nullptr && mis_defined_SurfaceSourceVar) {
        rSerializer.save("SurfaceSourceVarName",mpSurfaceSourceVar);
    }
    if (mpProjectionVar != nullptr && mis_defined_ProjectionVar) {
        rSerializer.save("ProjectionVarName",mpProjectionVar);
    }
    if (mpConvectionVar != nullptr && mis_defined_ConvectionVar) {
        rSerializer.save("ConvectionVarName",mpConvectionVar);
    }
    if (mpGradientVar != nullptr && mis_defined_GradientVar) {
        rSerializer.save("GradientVarName",mpGradientVar);
    }
    if (mpMeshVelocityVar != nullptr && mis_defined_MeshVelocityVar) {
        rSerializer.save("MeshVelocityVarName",mpMeshVelocityVar);
    }
    if (mpTransferCoefficientVar != nullptr && mis_defined_TransferCoefficientVar) {
        rSerializer.save("TransferCoefficientVarName",mpTransferCoefficientVar);
    }
    if (mpVelocityVar != nullptr && mis_defined_VelocityVar) {
        rSerializer.save("VelocityVarName",mpVelocityVar);
    }
    if (mpSpecificHeatVar != nullptr && mis_defined_SpecificHeatVar) {
        rSerializer.save("SpecificHeatVarName",mpSpecificHeatVar);
    }
    if (mpReactionVar != nullptr && mis_defined_ReactionVar) {
        rSerializer.save("ReactionVarName",mpReactionVar);
    }
    if (mpReactionGradientVar != nullptr && mIsDefinedReactionGradientVar) {
        rSerializer.save("ReactionGradientVarName",mpReactionGradientVar);
    }

    KRATOS_CATCH("")
}

void ConvectionDiffusionSettings::load(Serializer& rSerializer)
{
    KRATOS_TRY

    // Load the is defined bool flags for each variable
    rSerializer.load("mis_defined_DensityVar",mis_defined_DensityVar);
    rSerializer.load("mis_defined_DiffusionVar",mis_defined_DiffusionVar);
    rSerializer.load("mis_defined_UnknownVar",mis_defined_UnknownVar);
    rSerializer.load("mis_defined_VolumeSourceVar",mis_defined_VolumeSourceVar);
    rSerializer.load("mis_defined_SurfaceSourceVar",mis_defined_SurfaceSourceVar);
    rSerializer.load("mis_defined_ProjectionVar",mis_defined_ProjectionVar);
    rSerializer.load("mis_defined_ConvectionVar",mis_defined_ConvectionVar);
    rSerializer.load("mis_defined_GradientVar",mis_defined_GradientVar);
    rSerializer.load("mis_defined_MeshVelocityVar",mis_defined_MeshVelocityVar);
    rSerializer.load("mis_defined_TransferCoefficientVar",mis_defined_TransferCoefficientVar);
    rSerializer.load("mis_defined_VelocityVar",mis_defined_VelocityVar);
    rSerializer.load("mis_defined_SpecificHeatVar",mis_defined_SpecificHeatVar);
    rSerializer.load("mis_defined_ReactionVar",mis_defined_ReactionVar);
    rSerializer.load("mIsDefinedReactionGradientVar", mIsDefinedReactionGradientVar);

    // If the variables are defined, load their name
    // Note that only the name has been saved to retrieve the already existent variable from the KratosComponents
    if(mis_defined_DensityVar) {
        std::string density_var_name;
        rSerializer.load("DensityVarName", density_var_name);
        mpDensityVar = &(KratosComponents<Variable<double>>::Get(density_var_name));
    }
    if(mis_defined_DiffusionVar) {
        std::string diffusion_var_name;
        rSerializer.load("DiffusionVarName", diffusion_var_name);
        mpDiffusionVar = &(KratosComponents<Variable<double>>::Get(diffusion_var_name));
    }
    if(mis_defined_UnknownVar) {
        std::string unknown_var_name;
        rSerializer.load("UnknownVarName", unknown_var_name);
        mpUnknownVar = &(KratosComponents<Variable<double>>::Get(unknown_var_name));
    }
    if(mis_defined_VolumeSourceVar) {
        std::string volume_source_var_name;
        rSerializer.load("VolumeSourceVarName", volume_source_var_name);
        mpVolumeSourceVar = &(KratosComponents<Variable<double>>::Get(volume_source_var_name));
    }
    if(mis_defined_SurfaceSourceVar) {
        std::string surface_source_var_name;
        rSerializer.load("SurfaceSourceVarName", surface_source_var_name);
        mpSurfaceSourceVar = &(KratosComponents<Variable<double>>::Get(surface_source_var_name));
    }
    if(mis_defined_ProjectionVar) {
        std::string projection_var_name;
        rSerializer.load("ProjectionVarName", projection_var_name);
        mpProjectionVar = &(KratosComponents<Variable<double>>::Get(projection_var_name));
    }
    if(mis_defined_ConvectionVar) {
        std::string convection_var_name;
        rSerializer.load("ConvectionVarName", convection_var_name);
        mpConvectionVar = &(KratosComponents<Variable<array_1d<double,3>>>::Get(convection_var_name));
    }
    if(mis_defined_GradientVar) {
        std::string gradient_var_name;
        rSerializer.load("GradientVarName", gradient_var_name);
        mpGradientVar = &(KratosComponents<Variable<array_1d<double,3>>>::Get(gradient_var_name));
    }
    if(mis_defined_MeshVelocityVar) {
        std::string mesh_velocity_var;
        rSerializer.load("MeshVelocityVarName", mesh_velocity_var);
        mpMeshVelocityVar = &(KratosComponents<Variable<array_1d<double,3>>>::Get(mesh_velocity_var));
    }
    if(mis_defined_TransferCoefficientVar) {
        std::string transfer_coefficient_var_name;
        rSerializer.load("TransferCoefficientVarName", transfer_coefficient_var_name);
        mpTransferCoefficientVar = &(KratosComponents<Variable<double>>::Get(transfer_coefficient_var_name));
    }
    if(mis_defined_VelocityVar) {
        std::string velocity_var_name;
        rSerializer.load("VelocityVarName", velocity_var_name);
        mpVelocityVar = &(KratosComponents<Variable<array_1d<double,3>>>::Get(velocity_var_name));
    }
    if(mis_defined_SpecificHeatVar) {
        std::string specific_heat_var_name;
        rSerializer.load("SpecificHeatVarName", specific_heat_var_name);
        mpSpecificHeatVar = &(KratosComponents<Variable<double>>::Get(specific_heat_var_name));
    }
    if(mis_defined_ReactionVar) {
        std::string reaction_var_name;
        rSerializer.load("ReactionVarName", reaction_var_name);
        mpReactionVar = &(KratosComponents<Variable<double>>::Get(reaction_var_name));
    }
    if(mIsDefinedReactionGradientVar) {
        std::string reaction_gradient_var_name;
        rSerializer.load("ReactionGradientVarName", reaction_gradient_var_name);
        mpReactionGradientVar = &(KratosComponents<Variable<array_1d<double,3>>>::Get(reaction_gradient_var_name));
    }

    KRATOS_CATCH("")
}

}  // namespace Kratos.
