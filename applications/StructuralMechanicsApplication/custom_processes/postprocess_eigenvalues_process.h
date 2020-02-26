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
//

#if !defined(KRATOS_POSTPROCESS_EIGENVALUES_H_INCLUDED )
#define  KRATOS_POSTPROCESS_EIGENVALUES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"

namespace Kratos {

///@name Kratos Globals
///@{

///@name Kratos Classes
///@{

/// Process to create the animated Eigenvectors
/** This process takes the results of an Eigenvalue Analysis and creates the
 * animated Eigenvectors (Eigenmodes) for GiD using the GidEigenIO, which
 * is customized for this case
 * The Input should be the ComputingModelPart! (Otherwise nodal results migth be messed up)
 * It is of particular importance that all Nodes have the same Dofs!
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) PostprocessEigenvaluesProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of PostprocessEigenvaluesProcess
    KRATOS_CLASS_POINTER_DEFINITION(PostprocessEigenvaluesProcess);

    typedef std::size_t SizeType;

    typedef ModelPart::NodeType::DofsContainerType DofsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    PostprocessEigenvaluesProcess(ModelPart& rModelPart,
                                  Parameters OutputParameters);

    ///@}
    ///@name Operations
    ///@{

    void ExecuteFinalizeSolutionStep() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override {
        return "PostprocessEigenvaluesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {
        rOStream << "PostprocessEigenvaluesProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    Parameters mOutputParameters;

    ///@}
    ///@name Private Operations
    ///@{

    std::string GetLabel(const int NumberOfEigenValue,
                         const double EigenValueSolution) const;

    void GetVariables(std::vector<Variable<double>>& rRequestedDoubleResults,
                      std::vector<Variable<array_1d<double,3>>>& rRequestedVectorResults) const;

    ///@}

}; // Class PostprocessEigenvaluesProcess

///@}

}  // namespace Kratos.

#endif // KRATOS_POSTPROCESS_EIGENVALUES_H_INCLUDED  defined
