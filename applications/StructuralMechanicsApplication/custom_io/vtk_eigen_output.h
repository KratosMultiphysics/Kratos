// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Philipp Bucher
//

#pragma once

// System includes

// External includes

// Project includes
#include "input_output/vtk_output.h"

namespace Kratos
{
/** \brief VtkEigenOutput
* A simple class to write Eigenresults in Vtk format
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) VtkEigenOutput : public VtkOutput
{
public:

    /// Pointer definition of VtkEigenOutput
    KRATOS_CLASS_POINTER_DEFINITION(VtkEigenOutput);

    ///@name Life Cycle
    ///@{

    explicit VtkEigenOutput(
        ModelPart& rModelPart,
        Parameters EigenOutputParameters,
        Parameters VtkParameters)
            : VtkOutput(rModelPart, VtkParameters),
            mEigenOutputSettings(EigenOutputParameters) {};

    /// Destructor.
    virtual ~VtkEigenOutput() = default;

    ///@}
    ///@name Operations
    ///@{

    void PrintEigenOutput(
        const std::string& rLabel,
        const int AnimationStep,
        const std::vector<const Variable<double>*>& rRequestedDoubleResults,
        const std::vector<const Variable<array_1d<double,3>>*>& rRequestedVectorResults);

    ///@}

    /// Turn back information as a string.
    std::string Info() const override
    {
        return " VtkEigenOutput object ";
    }

    /**
     * @brief Prints information about the class
     * @param rOStream ostream object where output is printed
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << " VtkEigenOutput object " << std::endl;
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

private:
    ///@name Member Variables
    ///@{

    int mLastWrittenAnimationStepIndex = -1;
    Parameters mEigenOutputSettings;

    ///@}
    ///@name Operations
    ///@{

    /// opens the output file with the given openmode
    void OpenOutputFile(
        const std::string& rFileName,
        const std::ios::openmode OpenModeFlags,
        std::ofstream& rOutputFile) const;

    /// returns the name of the output file, depending on the current animation step
    std::string GetEigenOutputFileName(const int AnimationStep) const;

    /// writes a scalar variable of the nodes to the stream
    void WriteScalarEigenVariable(
        const ModelPart::NodesContainerType& rNodes,
        const Variable<double>& rVariable,
        const std::string& rLabel,
        std::ofstream& rFileStream) const;

    /// writes a vector variable of the nodes to the stream
    void WriteVectorEigenVariable(
        const ModelPart::NodesContainerType& rNodes,
        const Variable<array_1d<double, 3>>& rVariable,
        const std::string& rLabel,
        std::ofstream& rFileStream) const;

    ///@}
};

} // namespace Kratos
