//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_APPLY_PERTURBATION_FUNCTION_PROCESS_H_INCLUDED
#define KRATOS_APPLY_PERTURBATION_FUNCTION_PROCESS_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/model_part.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"


namespace Kratos
{
///@addtogroup ShallowWaterApplication
///@{

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

///@}
///@name Kratos Classes
///@{

/** 
 * @ingroup ShallowWaterApplication
 * @class ApplyPerturbationFunctionProcess
 * @brief This process assigns a default value or a perturbation if the node is close to an influence area
 */
template<class TVarType>
class KRATOS_API(SHALLOW_WATER_APPLICATION) ApplyPerturbationFunctionProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t                     IndexType;
    typedef Node                         NodeType;
    typedef Node::Pointer                NodePointerType;
    typedef ModelPart::NodesContainerType   NodesArrayType;
    typedef NodesArrayType::iterator        NodeIteratorType;

    /// Pointer definition of ApplyPerturbationFunctionProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyPerturbationFunctionProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with a node
    ApplyPerturbationFunctionProcess(
        ModelPart& rThisModelPart,
        NodePointerType pNode,
        TVarType& rThisVariable,
        Parameters& rThisParameters);

    /// Constructor with an array of nodes
    ApplyPerturbationFunctionProcess(
        ModelPart& rThisModelPart,
        NodesArrayType& rSourcePoints,
        TVarType& rThisVariable,
        Parameters& rThisParameters);

    /// Destructor.
    ~ApplyPerturbationFunctionProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    /*
     * @brief This operator is provided to call the process as a function and simply calls the Execute method.
     */
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
    void Execute() override;

    /**
     * @brief Perform a check with the parameters.
     */
    int Check() override;

    /**
     * @brief this function is designed for being execute once before the solution loop but
     * after all of the solvers where built
     */
    void ExecuteBeforeSolutionLoop() override;

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
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "ApplyPerturbationFunctionProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "ApplyPerturbationFunctionProcess";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {}


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

    ModelPart& mrModelPart;
    NodesArrayType mSourcePoints;
    TVarType& mrVariable;
    double mDefaultValue;
    double mInfluenceDistance;
    double mPerturbation;
    double mHalfWaveNumber;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void ValidateParameters(Parameters& rParameters);

    double ComputeDistance(NodeType& rNode);

    double PointPointSquareDistance(array_1d<double, 3>& rCoordA, array_1d<double, 3>& rCoordB);

    double ComputeInitialValue(double& rDistance);

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
    // ApplyPerturbationFunctionProcess& operator=(ApplyPerturbationFunctionProcess const& rOther);

    /// Copy constructor.
    // ApplyPerturbationFunctionProcess(ApplyPerturbationFunctionProcess const& rOther);


    ///@}

}; // Class ApplyPerturbationFunctionProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

// /// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                 ApplyPerturbationFunctionProcess& rThis);

// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                 const ApplyPerturbationFunctionProcess& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_APPLY_PERTURBATION_FUNCTION_PROCESS_H_INCLUDED  defined
