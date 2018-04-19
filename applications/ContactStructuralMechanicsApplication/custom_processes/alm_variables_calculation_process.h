// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_ALM_VARIABLES_CALCULATION_PROCESS)
#define KRATOS_ALM_VARIABLES_CALCULATION_PROCESS

// System includes

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "geometries/point.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

/**
 * @class ALMVariablesCalculationProcess
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This process calculates the variables related with the ALM
 * @details It is assumed that the YOUNG_MODULUS is provided
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) ALMVariablesCalculationProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ALMVariablesCalculationProcess
    KRATOS_CLASS_POINTER_DEFINITION(ALMVariablesCalculationProcess);

    // General type definitions
    typedef Node<3>                                          NodeType;
    typedef Point                                        PointType;
    typedef Geometry<NodeType>                           GeometryType;
    typedef Geometry<PointType>                     GeometryPointType;
    typedef ModelPart::NodesContainerType              NodesArrayType;
    typedef ModelPart::ConditionsContainerType    ConditionsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ALMVariablesCalculationProcess(
        ModelPart& rThisModelPart,
        Variable<double>& rNodalLengthVariable = NODAL_H,
        Parameters ThisParameters =  Parameters(R"({})")
        ):mrThisModelPart(rThisModelPart),
          mrNodalLengthVariable(rNodalLengthVariable)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "stiffness_factor"                     : 10.0,
            "penalty_scale_factor"                 : 1.0
        })" );

        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        mFactorStiffness = ThisParameters["stiffness_factor"].GetDouble();
        mPenaltyScale = ThisParameters["penalty_scale_factor"].GetDouble();

        if (rThisModelPart.GetNodalSolutionStepVariablesList().Has( rNodalLengthVariable ) == false ) // TODO: Ask Riccardo if it is necessary to use GetNodalSolutionStepVariablesList (REQUIRES BUFFER!!!!)
        {
            KRATOS_ERROR << "Missing variable " << rNodalLengthVariable;
        }

        KRATOS_CATCH("")
    }

    /// Destructor.
    ~ALMVariablesCalculationProcess() override
    = default;

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

    void Execute() override;

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
        return "ALMVariablesCalculationProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ALMVariablesCalculationProcess";
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

    ModelPart& mrThisModelPart;              // The main model part
    Variable<double>& mrNodalLengthVariable; // The variable used to messure the lenght of the element
    double mFactorStiffness;                 // The proportion between stiffness and penalty/scale factor
    double mPenaltyScale;                    // The penalty/scale factor proportion

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
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ALMVariablesCalculationProcess& operator=(ALMVariablesCalculationProcess const& rOther) = delete;

    /// Copy constructor.
    //ALMVariablesCalculationProcess(ALMVariablesCalculationProcess const& rOther);


    ///@}

}; // Class ALMVariablesCalculationProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   ALMVariablesCalculationProcess& rThis);
//
// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const ALMVariablesCalculationProcess& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);
//
//     return rOStream;
// }

}
#endif /* KRATOS_ALM_VARIABLES_CALCULATION_PROCESS defined */
