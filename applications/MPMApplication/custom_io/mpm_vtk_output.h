//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Crescenzio
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "input_output/vtk_output.h"

namespace Kratos
{
/** \brief MPMVtkOutput
* A simple class that has functionality to write vtk output
* @see : https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
*/
class KRATOS_API(MPM_APPLICATION) MPMVtkOutput : public VtkOutput
{
public:

    /// Definition of the size type
    using SizeType = std::size_t;

    /// Definition of the index type
    using IndexType = std::size_t;

    /// Pointer definition of MPMVtkOutput
    KRATOS_CLASS_POINTER_DEFINITION(MPMVtkOutput);

    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor by passing a ModelPart and Kratos-Parameters
     * @param rModelPart The modelpart which is used for output
     * @param Parameters Parameters including settings for the output
     */
    explicit MPMVtkOutput(
        ModelPart& rModelPart,
        Parameters ThisParameters = Parameters(R"({})" )
        );

    ///@}
    ///@name Operations
    ///@{

    ///@}

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    static Parameters GetDefaultParameters();

    /// Turn back information as a string.
    std::string Info() const override
    {
        return " MPMVtkOutput object ";
    }

    /**
     * @brief Prints information about the class
     * @param rOStream ostream object where output is printed
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << " MPMVtkOutput object " << std::endl;
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

protected:

    /**
     * @brief Write the nodes in the rModelPart.
     * @param rModelPart modelpart which is beging output
     * @param rFileStream the file stream to which data is to be written.
     */
    void WriteNodesToFile(
        const ModelPart& rModelPart,
        std::ofstream& rFileStream
        ) const override;

    /**
     * @brief Write the element/condition WriteConnectivity provided the container they are in
     * @tparam TEntity Element/Condition
     * @param rContainer The container containing elements/conditions
     * @param rFileStream the file stream to which data is to be written.
     */
    template <typename TContainerType>
    void WriteConnectivity(
        const TContainerType& rContainer,
        std::ofstream& rFileStream
        ) const;

    /**
     * @brief Write the element/condition cell types provided the container they are in
     * @tparam TEntity Element/Condition
     * @param rContainer The container containing elements/conditions
     * @param rFileStream the file stream to which data is to be written.
     */
    template <typename TContainerType>
    void WriteCellType(
        const TContainerType& rContainer,
        std::ofstream& rFileStream
        ) const;

    /**
     * @brief Write the elements and conditions in rModelPart.
     *        IMPORTANT : Need to write them together because of the CELLS block in VTK format
     * @param rModelPart modelpart which is beging output
     * @param rFileStream the file stream to which data is to be written.
     */
    void WriteConditionsAndElementsToFile(
        const ModelPart& rModelPart,
        std::ofstream& rFileStream
        ) const override;

    /**
     * @brief Write the results on the nodes.
     * @param rModelPart modelpart which is beging output
     * @param rFileStream the file stream to which data is to be written.
     */
    void WriteNodalResultsToFile(
        const ModelPart& rModelPart,
        std::ofstream& rFileStream
        ) override;

};

} // namespace Kratos
