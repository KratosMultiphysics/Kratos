//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (based on vtk_output.h)
//
//

#if !defined( KRATOS_VTU_OUTPUT_H_INCLUDED )
#define KRATOS_VTU_OUTPUT_H_INCLUDED

// Project includes
#include "includes/kratos_parameters.h"
#include "includes/io.h"

namespace Kratos
{
/** \brief VtuOutput
* A simple class that has functionality to write vtu output
* @see : https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
*/
class KRATOS_API(KRATOS_CORE) VtuOutput : public IO
{
public:

    /// Definition of the size type
    typedef std::size_t SizeType;

    /// Definition of the index type
    typedef std::size_t IndexType;

    /// Pointer definition of VtuOutput
    KRATOS_CLASS_POINTER_DEFINITION(VtuOutput);

    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor by passing a ModelPart and Kratos-Parameters
     * @param rModelPart The modelpart which is used for output
     * @param Parameters Parameters including settings for the output
     */
    explicit VtuOutput(ModelPart& rModelPart, Parameters ThisParameters);

    /// Destructor.
    virtual ~VtuOutput() = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    static Parameters GetDefaultParameters();

    /**
     * @brief Prints mrModelPart in VTU format together with the results
     */
    void PrintOutput();

    ///@}

    /// Turn back information as a string.
    std::string Info() const override
    {
        return " VtuOutput object ";
    }

    /**
     * @brief Prints information about the class
     * @param rOStream ostream object where output is printed
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << " VtuOutput object " << std::endl;
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    enum class FileFormat {
        VTU_ASCII,
        VTU_BINARY_RAW,
        VTU_BINARY_RAW_COMPRESSED,
        VTU_BINARY_BASE64,
        VTU_BINARY_BASE64_APPENDED
    };

protected:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;                        /// The main model part to post process
    VtuOutput::FileFormat mFileFormat;             /// The file format considered

    Parameters mOutputSettings;                    /// The configuration parameters

    ///@}
    ///@name Operations
    ///@{

    std::string GetOutputFileName(const ModelPart& rModelPart) const;

    ///@}
};

} // namespace Kratos

#endif // KRATOS_VTU_OUTPUT_H_INCLUDED
