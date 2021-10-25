//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Manuel Messmer

// Project includes
#include "geometries/geometry.h"
#include "containers/pointer_vector.h"
#include "geometries/nurbs_volume_geometry.h"
#include "projection_nurbs_volume_to_embedded_geometry_process.h"

namespace Kratos
{
    typedef Node<3>                                         NodeType;
    typedef NurbsVolumeGeometry<PointerVector<NodeType>>    NurbsVolumeGeometryType;
    typedef NurbsVolumeGeometryType::Pointer                NurbsVolumeGeometryPointerType;

    ProjectionNurbsVolumeToEmbeddedGeometryProcess::ProjectionNurbsVolumeToEmbeddedGeometryProcess(
        Model& rModel, Parameters ThisParameters) : mrModel(rModel), mThisParameters(ThisParameters)
    {
        mThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
        KRATOS_ERROR_IF_NOT( rModel.HasModelPart( mThisParameters["main_model_part_name"].GetString()) )
            << "ProjectionNurbsVolumeToEmbeddedGeometryProcess: Model Part '" <<  mThisParameters["main_model_part_name"].GetString() << "' does not exist." << std::endl;

        KRATOS_ERROR_IF_NOT( rModel.HasModelPart( mThisParameters["embedded_model_part_name"].GetString()) )
            << "ProjectionNurbsVolumeToEmbeddedGeometryProcess: Model Part '" <<  mThisParameters["embedded_model_part_name"].GetString() << "' does not exist." << std::endl;

        ModelPart& main_model_part = mrModel.GetModelPart(mThisParameters["main_model_part_name"].GetString());
        KRATOS_ERROR_IF_NOT( main_model_part.HasGeometry(mThisParameters["nurbs_volume_name"].GetString()) )
            << "ProjectionNurbsVolumeToEmbeddedGeometryProcess: Model Part '" <<  mThisParameters["embedded_model_part_name"].GetString() << "' does not have Geometry: '"
                << mThisParameters["nurbs_volume_name"].GetString() << "'. " << std::endl;
    }

    void ProjectionNurbsVolumeToEmbeddedGeometryProcess::ExecuteBeforeOutputStep()
    {
        ModelPart& main_model_part = mrModel.GetModelPart(mThisParameters["main_model_part_name"].GetString());
        ModelPart& embedded_model_part = mrModel.GetModelPart(mThisParameters["embedded_model_part"].GetString());

        ModelPart::GeometryType::Pointer p_geometry = main_model_part.pGetGeometry(mThisParameters["nurbs_volume_name"].GetString());
        // Todo: Check geometry type!
        NurbsVolumeGeometryPointerType p_nurbs_volume = std::dynamic_pointer_cast<NurbsVolumeGeometryType>(p_geometry);


        for( auto variable_name : mThisParameters["nodal_results"].GetVector() ){
            std::cout << "name" << std::endl;
        }

    }

    // void ProjectNodalValue( Variable<double> const& rVariable,
    //     GidIOType::NodesContainerType& rNodes, double SolutionTag,
    //     std::size_t SolutionStepNumber);

    // void CreateVariablesListScalar(std::string& variable_name)
    // {

    // }

    // void CreateVariablesList(Parameters ThisParameters)
    // {
    //     const std::size_t n_variables = ThisParameters["solution_variables"].size();

    //     // The current dimension
    //     mDomainSize = ThisParameters["domain_size"].GetInt();

    //     const auto variable_names = ThisParameters["solution_variables"].GetStringArray();

    //     for (std::size_t p_var = 0; p_var < n_variables; ++p_var){
    //         const std::string& variable_name = variable_names[p_var];

    //         if(KratosComponents<Variable<double>>::Has(variable_name)){
    //             const auto& r_var = KratosComponents<Variable<double>>::Get(variable_name);
    //             mDoubleVariable.push_back(&r_var);
    //         } else if (KratosComponents< Variable< array_1d< double, 3> > >::Has(variable_name)) {
    //             // Components
    //             const auto& r_var_x = KratosComponents<Variable<double>>::Get(variable_name+"_X");
    //             const auto& r_var_y = KratosComponents<Variable<double>>::Get(variable_name+"_Y");
    //             mDoubleVariable.push_back(&r_var_x);
    //             mDoubleVariable.push_back(&r_var_y);

    //             if (mDomainSize == 3) {
    //                 const auto& r_var_z = KratosComponents<Variable<double>>::Get(variable_name+"_Z");
    //                 mDoubleVariable.push_back(&r_var_z);
    //             }
    //         } else {
    //             KRATOS_ERROR << "Only double and vector variables are allowed in the variables list." ;
    //         }
    //     }
    // }

} // End namespace Kratos
