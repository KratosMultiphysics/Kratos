#pragma once

// System includes

// External includes

// Project includes
#include "modeler/modeler.h"
#include "containers/model.h"
#include "includes/kratos_parameters.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @class CombineModelPartModeler
 * @ingroup KratosCore
 * @brief This modeler combines several model parts into a new one. All input model parts should share
 *   the same root model part.
 * @details ProcessInfo is shallow-copied from root origin model part.
 * @author Daniel Diez @ddiezrod
 */
class KRATOS_API(KRATOS_CORE) CombineModelPartModeler
    : public Modeler {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CombineModelPartModeler
    KRATOS_CLASS_POINTER_DEFINITION(CombineModelPartModeler);

    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor.
    CombineModelPartModeler() : Modeler()
    {
    }

    /**
     * @brief Constructor using a Model and Parameters
     * @param rModel Reference of the Model
     * @param ModelerParameters Parameters of the discretization
     */
    CombineModelPartModeler(
        Model& rModel,
        Parameters ModelerParameters
    );

    /// Destructor.
    virtual ~CombineModelPartModeler() = default;

    ///@}
    ///@name Operations
    ///@{

    const Parameters GetDefaultParameters() const override;

    /**
     * @brief Creates the Modeler Pointer and returns a pointer to a new
     * CombineModelPartModeler, created using the given input
     * @param rModel Reference of the Model
     * @param ModelerParameters Parameters of the discretization
     * @return a Pointer to the new Modeler
     */
    Modeler::Pointer Create(
        Model& rModel,
        const Parameters ModelParameters
        ) const override;

    /**
     * @brief Convert the geometry into an analysis suitable ModePart
     */
    void SetupModelPart() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "CombineModelPartModeler";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
protected:
    ///@name Protected Operations
    ///@{


    void DuplicateElements(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const Element& rReferenceElement) const;

    void DuplicateConditions(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const Condition& rReferenceBoundaryCondition) const;

    void DuplicateCommunicatorData(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart) const;

    void DuplicateSubModelParts(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart) const;

    void CopyCommonData(
        ModelPart& rCombinedModelPart) const;

    void DuplicateMesh() const;

    void CreateSubModelParts();

    void CreateCommunicators();

    void PopulateCommunicators();

    void PopulateLocalMesh(
        Communicator& rReferenceComm,
        Communicator& rDestinationComm,
        ModelPart& rDestinationModelPart
    ) const;

    void ResetModelPart(ModelPart& rCombinedModelPart) const;

    void CheckOriginModelPartsAndAssignRoot();

    ///@}
private:
    ///@name Member Variables
    ///@{

    Model* mpModel = nullptr; /// The model of the problem

    Parameters mParameters;   /// The configuration parameters

    ModelPart* mpOriginRootModelPart = nullptr; /// Root model part of all origin model parts

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
    }

    void load(Serializer& rSerializer)
    {
    }

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}

}; // Class CombineModelPartModeler

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator>>(std::istream& rIStream,
    CombineModelPartModeler& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream, const CombineModelPartModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

} // namespace Kratos.

