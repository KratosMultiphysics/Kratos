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
//

#if !defined( KRATOS_VTK_ROM_BASIS_OUTPUT_H_INCLUDED )
#define KRATOS_VTK_ROM_BASIS_OUTPUT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "input_output/vtk_output.h"

namespace Kratos
{
/** \brief VtkRomBasisOutput
* A simple class to write RomBasis in Vtk format
* This class is adapted from VtkEigenOutput in StructuralMechanicsApplication
*/
class KRATOS_API(ROM_APPLICATION) VtkRomBasisOutput : public VtkOutput
{
public:

    /// Pointer definition of VtkRomBasisOutput
    KRATOS_CLASS_POINTER_DEFINITION(VtkRomBasisOutput);

    ///@name Life Cycle
    ///@{

    explicit VtkRomBasisOutput(
        ModelPart& rModelPart,
        Parameters RomBasisOutputParameters,
        Parameters VtkParameters)
            : VtkOutput(rModelPart, VtkParameters),
            mRomBasisOutputSettings(RomBasisOutputParameters) {};

    /// Destructor.
    virtual ~VtkRomBasisOutput() = default;

    ///@}
    ///@name Operations
    ///@{

    void PrintRomBasisOutput(
        const std::string& rLabel,
        const int AnimationStep,
        const std::vector<Variable<double>>& rRequestedDoubleResults,
        const std::vector<Variable<array_1d<double,3>>>& rRequestedVectorResults);

    ///@}

    /// Turn back information as a string.
    std::string Info() const override
    {
        return " VtkRomBasisOutput object ";
    }

    /**
     * @brief Prints information about the class
     * @param rOStream ostream object where output is printed
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << " VtkRomBasisOutput object " << std::endl;
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

private:
    ///@name Member Variables
    ///@{

    int mLastWrittenAnimationStepIndex = -1;
    Parameters mRomBasisOutputSettings;

    ///@}
    ///@name Operations
    ///@{

    /// opens the output file with the given openmode
    void OpenOutputFile(
        const std::string& rFileName,
        const std::ios::openmode OpenModeFlags,
        std::ofstream& rOutputFile) const;

    /// returns the name of the output file, depending on the current animation step
    std::string GetRomBasisOutputFileName(const int AnimationStep) const;

    /// writes a scalar variable of the nodes to the stream
    void WriteScalarRomBasisVariable(
        const ModelPart::NodesContainerType& rNodes,
        const Variable<double>& rVariable,
        const std::string& rLabel,
        std::ofstream& rFileStream) const;

    /// writes a vector variable of the nodes to the stream
    void WriteVectorRomBasisVariable(
        const ModelPart::NodesContainerType& rNodes,
        const Variable<array_1d<double, 3>>& rVariable,
        const std::string& rLabel,
        std::ofstream& rFileStream) const;

    ///@}
};

} // namespace Kratos

#endif // KRATOS_VTK_ROM_BASIS_OUTPUT_H_INCLUDED
