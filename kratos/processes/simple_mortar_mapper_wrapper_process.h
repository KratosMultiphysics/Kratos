//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_SIMPLE_MORTAR_MAPPER_WRAPPER_PROCESS)
#define KRATOS_SIMPLE_MORTAR_MAPPER_WRAPPER_PROCESS

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
    typedef std::size_t SizeType;

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

    /// Linear solver
    typedef UblasSpace<double, CompressedMatrix, Vector>    SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector>               LocalSpaceType;
    typedef typename SparseSpaceType::MatrixType                 MatrixType;
    typedef typename SparseSpaceType::VectorType                 VectorType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    /// Index type definition
    typedef std::size_t                                          IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @param rOriginModelPart The origin model part to compute
     * @param rDestinationModelPart The destination model part to compute
     * @param rThisVariable The variable to transfer and be transfered
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
        Parameters default_parameters = GetDefaultParameters();
        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        // The condition iterators
        auto it_cond_origin_begin = rOriginModelPart.Conditions().begin();
        auto it_cond_destination_begin = rDestinationModelPart.Conditions().begin();

        // The dimensions
        const SizeType dimension = it_cond_origin_begin->GetGeometry().WorkingSpaceDimension();
        const SizeType size_1 = it_cond_origin_begin->GetGeometry().size();
        const SizeType size_2 = it_cond_destination_begin->GetGeometry().size();

        // The variable names
        const std::string& origin_variable_name = ThisParameters["origin_variable"].GetString();
        const std::string& destination_variable_name = ThisParameters["destination_variable"].GetString();

        bool double_variable = true;
        if(KratosComponents<Variable<double>>::Has(origin_variable_name)) {
            if (destination_variable_name != "" && !(KratosComponents<Variable<double>>::Has(destination_variable_name))) {
                KRATOS_ERROR << "The destination variable is not the same type (double) as the origin" << std::endl;
            }
        } else if (KratosComponents< Variable< array_1d< double, 3> > >::Has(origin_variable_name)) {
            double_variable = false;
            if (destination_variable_name != "" && !(KratosComponents<Variable<array_1d< double, 3>>>::Has(destination_variable_name)))
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

    void Execute() override
    {
        mpMapperProcess->Execute();
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
protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{
    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    Parameters GetDefaultParameters()
    {
        Parameters default_parameters = Parameters(R"(
        {
            "echo_level"                       : 0,
            "absolute_convergence_tolerance"   : 1.0e-9,
            "relative_convergence_tolerance"   : 1.0e-4,
            "max_number_iterations"            : 10,
            "integration_order"                : 2,
            "distance_threshold"               : 1.0e24,
            "mapping_coefficient"              : 1.0e0,
            "origin_variable"                  : "TEMPERATURE",
            "destination_variable"             : "",
            "origin_variable_historical"       : true,
            "destination_variable_historical"  : true,
            "search_parameters"                : {
                "allocation_size"                  : 1000,
                "bucket_size"                      : 4,
                "search_factor"                    : 3.5
            }
        })" );

        return default_parameters;
    }

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
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
#endif /* KRATOS_SIMPLE_MORTAR_MAPPER_WRAPPER_PROCESS defined */
