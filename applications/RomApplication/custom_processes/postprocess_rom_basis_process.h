//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Altug Emiroglu, http://github.com/emiroglu
//

#if !defined(KRATOS_POSTPROCESS_ROM_BASIS_H_INCLUDED )
#define  KRATOS_POSTPROCESS_ROM_BASIS_H_INCLUDED

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

/// Process to create the animated ROM_BASIS
/** This process takes the results of an Eigenvalue Analysis and creates the
 * animated Eigenvectors (Eigenmodes) for GiD using the GidEigenIO, which
 * is customized for this case
 * The Input should be the ComputingModelPart! (Otherwise nodal results migth be messed up)
 * It is of particular importance that all Nodes have the same Dofs!
 */
class KRATOS_API(ROM_APPLICATION) PostprocessRomBasisProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of PostprocessRomBasisProcess
    KRATOS_CLASS_POINTER_DEFINITION(PostprocessRomBasisProcess);

    typedef std::size_t SizeType;

    typedef ModelPart::NodeType::DofsContainerType DofsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    PostprocessRomBasisProcess(ModelPart& rModelPart,
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
        return "PostprocessRomBasisProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {
        rOStream << "PostprocessRomBasisProcess";
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

    std::string GetLabel(const int NumberOfRomBasis) const;

    void GetVariables(std::vector<Variable<double>>& rRequestedDoubleResults,
                      std::vector<Variable<array_1d<double,3>>>& rRequestedVectorResults) const;

    ///@}

}; // Class PostprocessEigenvaluesProcess

///@}

}  // namespace Kratos.

#endif // KRATOS_POSTPROCESS_ROM_BASIS_H_INCLUDED  defined
