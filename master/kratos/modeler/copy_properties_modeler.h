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
//

#ifndef KRATOS_COPY_PROPERTIES_MODELER_H_INCLUDED
#define KRATOS_COPY_PROPERTIES_MODELER_H_INCLUDED


// System includes


// External includes


// Project includes
#include "modeler/modeler.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @author Miguel Maso Sotomayor
 * @class CopyPropertiesModeler
 * @brief Copy the properties from one model part to another.
 * @details The properties of the elements and conditions of the destination model part are replaced by the copies.
 */
class KRATOS_API(KRATOS_CORE) CopyPropertiesModeler : public Modeler
{
public:
    ///@name Pointer Definition
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(CopyPropertiesModeler);

    ///@}
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    ///@brief Default constructor.
    CopyPropertiesModeler() : Modeler() {}

    /**
     * @brief Constructor with Model and Parameters
     * @param rModel Reference of the Model
     * @param ModelerParameters Parameters of the Modeler
     */
    CopyPropertiesModeler(
        Model& rModel,
        Parameters ModelerParameters);

    /**
     * @brief Constructor with origin and destination model parts
     * @param rOriginModelPart Reference of the origin ModelPart
     * @param rDestinationModelPart Reference of the destination ModelPart
     */
    CopyPropertiesModeler(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart);

    ///@brief Copy constructor.
    CopyPropertiesModeler(CopyPropertiesModeler const& rOther) = delete;

    ///@brief Destructor.
    ~CopyPropertiesModeler() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@brief Assignment operator.
    CopyPropertiesModeler & operator=(CopyPropertiesModeler const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates the Modeler Pointer and returns a pointer to a new CopyPropertiesModeler
     * @param rModel Reference of the Model
     * @param ModelerParameters Parameters of the discretization
     * @return a Pointer to the new Modeler
     */
    Modeler::Pointer Create(
        Model& rModel,
        const Parameters ModelParameters) const override;

    /**
     * @brief Get the Default Parameters object
     * @return * This const 
     */
    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Stages
    ///@{

    /**
     * @brief Create a copy of the properties from the origin ModelPart to the destination ModelPart
     */
    void SetupModelPart() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "CopyPropertiesModeler";
    }

    ///@}

private:

    ///@name Member Variables
    ///@{

    Model* mpModel = nullptr;

    ///@}
    ///@name Operations
    ///@{

    void RecursivelyCopyProperties(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart);

    template<class TContainerType>
    void ReplaceProperties(
        TContainerType& rContainer,
        const ModelPart& rModelPart);

    ///@}
};

///@}
///@name Input and output
///@{

///@brief input stream function
inline std::istream& operator>>(std::istream& rIStream,
    CopyPropertiesModeler& rThis)
{
    return rIStream;
}

///@brief output stream function
inline std::ostream& operator<<(std::ostream& rOStream, const CopyPropertiesModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

} // namespace Kratos.

#endif //KRATOS_COPY_PROPERTIES_MODELER_H_INCLUDED defined
