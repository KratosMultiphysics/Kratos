// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#pragma once

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

    PostprocessEigenvaluesProcess(
        Model& rModel,
        Parameters OutputParameters);

    ///@}
    ///@name Operations
    ///@{

    void ExecuteFinalizeSolutionStep() override;

    ///@}
    ///@name Inquiry
    ///@{

    const Parameters GetDefaultParameters() const override;

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

    ModelPart* mpModelPart;
    Parameters mOutputParameters;

    ///@}
    ///@name Private Operations
    ///@{

    std::string GetLabel(const int NumberOfEigenValue,
                         const int NumberOfEigenvalues,
                         const double EigenValueSolution) const;

    void GetVariables(std::vector<const Variable<double>*>& rRequestedDoubleResults,
                      std::vector<const Variable<array_1d<double,3>>*>& rRequestedVectorResults) const;

    ///@}

}; // Class PostprocessEigenvaluesProcess

///@}

}  // namespace Kratos.
