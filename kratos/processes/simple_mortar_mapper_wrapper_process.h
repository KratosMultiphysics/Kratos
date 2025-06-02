//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "processes/simple_mortar_mapper_process.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// The definition of the size type
    using SizeType = std::size_t;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @ingroup KratosCore
 * @class SimpleMortarMapperProcessWrapper
 * @brief This class wraps automatically the different types mof mappers
 * @author Vicente Mataix Ferrandiz
 */
class SimpleMortarMapperProcessWrapper
        : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SimpleMortarMapperProcessWrapper
    KRATOS_CLASS_POINTER_DEFINITION(SimpleMortarMapperProcessWrapper);

    /// Type definition for sparse space type
    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;

    /// Type definition for local space type
    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

    /// Type definition for matrix
    using MatrixType = typename SparseSpaceType::MatrixType;

    /// Type definition for vector
    using VectorType = typename SparseSpaceType::VectorType;

    /// Type definition for linear solver
    using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;

    /// Index type definition
    using IndexType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @param rOriginModelPart The origin model part to compute
     * @param rDestinationModelPart The destination model part to compute
     * @param rThisVariable The variable to transfer and be transferred
     * @param ThisParameters The configuration parameters
     * @param pThisLinearSolver The pointer to the linear to be used (in case of implicit resolution)
     */
    SimpleMortarMapperProcessWrapper(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        Parameters ThisParameters = Parameters(R"({})" ),
        LinearSolverType::Pointer pThisLinearSolver = nullptr
        )
    {
        // The default parameters
        const Parameters default_parameters = GetDefaultParameters();
        const bool origin_are_conditions_is_defined = ThisParameters.Has("origin_are_conditions");
        const bool destination_are_conditions_is_defined = ThisParameters.Has("destination_are_conditions");
        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        // Automatic detect the entities
        if (!origin_are_conditions_is_defined) {
            if (rOriginModelPart.NumberOfElements() > 0) {
                ThisParameters["origin_are_conditions"].SetBool(false);
                KRATOS_WARNING("SimpleMortarMapperProcessWrapper") << "\'origin_are_conditions\' changed to \'False\'. Mapping from elements." << std::endl;
            } else if (rOriginModelPart.NumberOfConditions() == 0) {
                KRATOS_ERROR << "No conditions defined on origin model part" << std::endl;
            }
        }
        if (!destination_are_conditions_is_defined) {
            if (rDestinationModelPart.NumberOfElements() > 0) {
                ThisParameters["destination_are_conditions"].SetBool(false);
                KRATOS_WARNING("SimpleMortarMapperProcessWrapper") << "\'destination_are_conditions\' changed to \'False\'. Mapping from elements." << std::endl;
            } else if (rDestinationModelPart.NumberOfConditions() == 0) {
                KRATOS_ERROR << "No conditions defined on destination model part" << std::endl;
            }
        }

        // The condition iterators
        auto& r_geometry_origin = ThisParameters["origin_are_conditions"].GetBool() ? rOriginModelPart.Conditions().begin()->GetGeometry() : rOriginModelPart.Elements().begin()->GetGeometry();
        auto& r_geometry_destination = ThisParameters["destination_are_conditions"].GetBool() ? rDestinationModelPart.Conditions().begin()->GetGeometry() : rDestinationModelPart.Elements().begin()->GetGeometry();

        // The dimensions
        const SizeType dimension = r_geometry_origin.WorkingSpaceDimension();
        const SizeType size_1 = r_geometry_origin.size();
        const SizeType size_2 = r_geometry_destination.size();

        // The variable names
        const std::string& r_origin_variable_name = ThisParameters["origin_variable"].GetString();
        const std::string& r_destination_variable_name = ThisParameters["destination_variable"].GetString();

        bool double_variable = true;
        if(KratosComponents<Variable<double>>::Has(r_origin_variable_name)) {
            if (r_destination_variable_name != "" && !(KratosComponents<Variable<double>>::Has(r_destination_variable_name))) {
                KRATOS_ERROR << "The destination variable is not the same type (double) as the origin" << std::endl;
            }
        } else if (KratosComponents< Variable< array_1d< double, 3> > >::Has(r_origin_variable_name)) {
            double_variable = false;
            if (r_destination_variable_name != "" && !(KratosComponents<Variable<array_1d< double, 3>>>::Has(r_destination_variable_name)))
                KRATOS_ERROR << "The destination variable is not the same type (array_1d< double, 3>) as the origin" << std::endl;
        } else {
            KRATOS_ERROR << "The types of the variables are not supported array_1d< double, 3> or double" << std::endl;
        }

        // Creating the mapper
        if (dimension == 2) {
            // 2D
            if (double_variable) {
                mpMapperProcess = Kratos::make_shared<SimpleMortarMapperProcess<2, 2, Variable<double>>>(rOriginModelPart, rDestinationModelPart, ThisParameters, pThisLinearSolver);
            } else {
                mpMapperProcess = Kratos::make_shared<SimpleMortarMapperProcess<2, 2, Variable<array_1d< double, 3>>>>(rOriginModelPart, rDestinationModelPart, ThisParameters, pThisLinearSolver);
            }
        } else {
            // 3D
            if (size_1 == 3 && size_2 == 3) {
                if (double_variable) {
                    mpMapperProcess = Kratos::make_shared<SimpleMortarMapperProcess<3, 3, Variable<double>>>(rOriginModelPart, rDestinationModelPart, ThisParameters, pThisLinearSolver);
                } else {
                    mpMapperProcess = Kratos::make_shared<SimpleMortarMapperProcess<3, 3, Variable<array_1d< double, 3>>>>(rOriginModelPart, rDestinationModelPart, ThisParameters, pThisLinearSolver);
                }
            } else if (size_1 == 4 && size_2 == 4) {
                if (double_variable) {
                    mpMapperProcess = Kratos::make_shared<SimpleMortarMapperProcess<3, 4, Variable<double>>>(rOriginModelPart, rDestinationModelPart, ThisParameters, pThisLinearSolver);
                } else {
                    mpMapperProcess = Kratos::make_shared<SimpleMortarMapperProcess<3, 4, Variable<array_1d< double, 3>>>>(rOriginModelPart, rDestinationModelPart, ThisParameters, pThisLinearSolver);
                }
            } else if (size_1 == 4 && size_2 == 3) {
                if (double_variable) {
                    mpMapperProcess = Kratos::make_shared<SimpleMortarMapperProcess<3, 3, Variable<double>, 4>>(rOriginModelPart, rDestinationModelPart, ThisParameters, pThisLinearSolver);
                } else {
                    mpMapperProcess = Kratos::make_shared<SimpleMortarMapperProcess<3, 3, Variable<array_1d< double, 3>>, 4>>(rOriginModelPart, rDestinationModelPart, ThisParameters, pThisLinearSolver);
                }
            } else if (size_1 == 3 && size_2 == 4) {
                if (double_variable) {
                    mpMapperProcess = Kratos::make_shared<SimpleMortarMapperProcess<3, 4, Variable<double>, 3>>(rOriginModelPart, rDestinationModelPart, ThisParameters, pThisLinearSolver);
                } else {
                    mpMapperProcess = Kratos::make_shared<SimpleMortarMapperProcess<3, 4, Variable<array_1d< double, 3>>, 3>>(rOriginModelPart, rDestinationModelPart, ThisParameters, pThisLinearSolver);
                }
            } else {
                KRATOS_ERROR << "Combination of dimensions and sizes not compatible. Dimension: " << dimension << " Size Origin: " << size_1 << " Size Destination: " << size_2 << std::endl;
            }
        }
    }

    /// Destructor.
    ~SimpleMortarMapperProcessWrapper() override = default;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Execute method is used to execute the Process algorithms.
     */
    void Execute() override
    {
        mpMapperProcess->Execute();
    }

    /**
     * @details This function will be executed at every time step BEFORE performing the solve phase
     */
    void ExecuteInitializeSolutionStep() override
    {
        mpMapperProcess->ExecuteInitializeSolutionStep();
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters = Parameters(R"(
        {
            "echo_level"                       : 0,
            "consider_tessellation"            : false,
            "using_average_nodal_normal"       : true,
            "discontinuous_interface"          : false,
            "discontinous_interface_factor"    : 1.0e-4,
            "absolute_convergence_tolerance"   : 1.0e-9,
            "relative_convergence_tolerance"   : 1.0e-4,
            "max_number_iterations"            : 10,
            "integration_order"                : 2,
            "distance_threshold"               : 1.0e24,
            "zero_tolerance_factor"            : 1.0e0,
            "remove_isolated_conditions"       : false,
            "mapping_coefficient"              : 1.0e0,
            "origin_variable"                  : "TEMPERATURE",
            "destination_variable"             : "",
            "origin_variable_historical"       : true,
            "origin_are_conditions"            : true,
            "destination_variable_historical"  : true,
            "destination_are_conditions"       : true,
            "update_interface"                 : true,
            "search_parameters"                : {
                "allocation_size"                  : 1000,
                "bucket_size"                      : 4,
                "search_factor"                    : 3.5
            }
        })" );

        return default_parameters;
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "SimpleMortarMapperProcessWrapper";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SimpleMortarMapperProcessWrapper";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

    Process::Pointer mpMapperProcess = nullptr; /// The real mapper process

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{
        
    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{

    ///@}
};// class SimpleMortarMapperProcessWrapper
///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.
