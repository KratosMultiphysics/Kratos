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
//

#if !defined( KRATOS_VTK_EIGEN_OUTPUT_H_INCLUDED )
#define KRATOS_VTK_EIGEN_OUTPUT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "input_output/vtk_output.h"

namespace Kratos
{
/** \brief VtkEigenOutput
* A simple class that has functionality to write vtk output
* @see : https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
*/
class KRATOS_API(KRATOS_CORE) VtkEigenOutput : public VtkOutput
{
public:

    /// Pointer definition of VtkEigenOutput
    KRATOS_CLASS_POINTER_DEFINITION(VtkEigenOutput);

    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor by passing a ModelPart and Kratos-Parameters
     * @param rModelPart The modelpart which is used for output
     * @param Parameters Parameters including settings for the output
     */
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

    /**
     * @brief Prints mrModelPart in VTK format together with the results
     */
    void PrintEigenOutput(
        const std::string& rLabel,
        const int AnimationStep,
        const std::vector<Variable<double>>& rRequestedDoubleResults,
        const std::vector<Variable<array_1d<double,3>>>& rRequestedVectorResults);

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

    void WriteScalarEigenVariable(
        const ModelPart::NodesContainerType& rNodes,
        const Variable<double>& rVariable,
        const std::string& rLabel,
        std::ofstream& rFileStream) const;

    void WriteVectorEigenVariable(
        const ModelPart::NodesContainerType& rNodes,
        const Variable<array_1d<double, 3>>& rVariable,
        const std::string& rLabel,
        std::ofstream& rFileStream) const;

    void OpenOutputFile(
        const std::string& rFileName,
        const std::ios::openmode OpenModeFlags,
        std::ofstream& rOutputFile) const;

    std::string GetEigenOutputFileName(const int AnimationStep) const;

    ///@}
};

} // namespace Kratos

#endif // KRATOS_VTK_EIGEN_OUTPUT_H_INCLUDED
