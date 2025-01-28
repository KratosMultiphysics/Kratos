//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Juan I. Camarotti

#if !defined(KRATOS_CALCULATE_L2_NORM_ERROR_PROCESS_H_INCLUDED )
#define  KRATOS_CALCULATE_L2_NORM_ERROR_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "processes/process.h"
#include "utilities/function_parser_utility.h"


// Application includes
#include "iga_application_variables.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/* @class CalculateL2NormErrorProcess
 * @ingroup IgaApplication
 **/
class KRATOS_API(IGA_APPLICATION) CalculateL2NormErrorProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;
    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointerType;
    typedef array_1d<double, 3> CoordinatesArrayType;
    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    /// Pointer definition of CalculateL2NormErrorProcess
    KRATOS_CLASS_POINTER_DEFINITION(CalculateL2NormErrorProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    CalculateL2NormErrorProcess(
        Model& rModel,
        Parameters ThisParameters);

    /// Destructor.
    ~CalculateL2NormErrorProcess() = default;

    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

    void ExecuteFinalize() override {
        Execute();
    };


    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters = Parameters(R"(
        {
            "compute_error_model_part_name"                 : "please_specify_model_part_name",
            "analytical_solution"              : "please_specify_analitical_solution",
            "unknown_variable_name"            : "please_specify_unknown_variable_name"
        })" );

        return default_parameters;
    }

    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "CalculateL2NormErrorProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CalculateL2NormErrorProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}

private:
    ///@name Iga functionalities
    ///@{

    Model* mpModel = nullptr;
    Parameters mParameters;
    SizeType mEchoLevel;
    ModelPart* mpComputeErrorModelPart = nullptr;
    std::string mUnknownVariable;
    Kratos::unique_ptr<GenericFunctionUtility> mpAnalyticalSolution;

    ///@}
    ///@name Iga functionalities
    ///@{

    ///@}
    ///@}

    ///@}
    ///@name Utility
    ///@{


    ///@}
    ///@name Input and output
    ///@{


   

    ///@}

}; // Class CalculateL2NormErrorProcess

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  CalculateL2NormErrorProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const CalculateL2NormErrorProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_CALCULATE_L2_NORM_ERROR_PROCESS_H_INCLUDED 
