//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez
//

#pragma once

// Project includes
#include "includes/model_part.h"
#include "operations/operation.h"
#include "compressible_potential_flow_application_variables.h"

namespace Kratos
{
///@addtogroup CompressiblePotentialFlowApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class Define3DWakeOperation
 * @ingroup CompressiblePotentialFlowApplication
 * @brief This operation define the wake in 2d problems.
 * @details For a given model, this operation creates a line wake and
 *  uses it to select the wake and kutta elements. 
 * @authors Inigo Lopez and Marc Nunez 
 */
class KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) Define2DWakeOperation : public Operation
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Define2DWakeOperation
    KRATOS_CLASS_POINTER_DEFINITION(Define2DWakeOperation);

    typedef Node NodeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Define2DWakeOperation() : Operation() {};

    /// Constructor

    /**
     * @brief Constructor with Kratos parameters and Model container
     * @param rModel The Model container
     * @param rParameters Kratos parameters encapsulating the settings
     */
    Define2DWakeOperation(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    ~Define2DWakeOperation() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Define2DWakeOperation& operator=(Define2DWakeOperation const& rOther) = delete;

    /// Copy constructor.
    Define2DWakeOperation(Define2DWakeOperation const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    Operation::Pointer Create(
        Model& rModel,
        Parameters ThisParameters) const override
    {
        return Kratos::make_shared<Define2DWakeOperation>(rModel, ThisParameters);
    }

    /**
     * @brief Execute method is used to execute the Operation algorithms.
     */
    void Execute() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        return "Define2DWakeOperation";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Define2DWakeOperation";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
    }

    ///@}
private:
    ///@name Member Variables
    ///@{
        
    /// Registry current operation
    KRATOS_REGISTRY_ADD_PROTOTYPE("Operations.KratosMultiphysics.CompressiblePotentialFlowApplication", Operation, Define2DWakeOperation)
    KRATOS_REGISTRY_ADD_PROTOTYPE("Operations.All", Operation, Define2DWakeOperation)


    // Model parts
    ModelPart* mpBodyModelPart = nullptr;
    ModelPart* mpWakeModelPart = nullptr;
    ModelPart* mpRootModelPart = nullptr;

    NodeType* mpTrailingEdgeNode;

    // Tolerances 
    double mWakeDistanceTolerance;                  // To avoid nodes laying exactly on the wake 

    BoundedVector<double, 3> mWakeDirection;
    BoundedVector<double, 3> mWakeNormal;
    
    int mEchoLevel;

    // Vector to store the trailing edge elements ids
    std::vector<std::size_t> mTrailingEdgeElementsOrderedIds;

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method initializes nodal and elemental variables such as 
     * WAKE, WAKE_DISTANCE, saves the WAKE_NORMAL and computes the 
     * wake_direction and span direction.
     */
    void InitializeVariables();

    /**
     * @brief This method initializes the trailing edge submodelpart.
     */
    void InitializeTrailingEdgeSubModelpart() const;

    /**
     * @brief This method initializes the wake submodelpart.
     */
    void InitializeWakeSubModelpart() const;

    /**
     * @brief This method saves the trailing edge for further computations.
     */
    void SaveTrailingEdgeNode();

    /**
     * @brief This method checks which elements are cut and mark them as wake. 
     */
    void MarkWakeElements();

    /**
     * @brief This method checks if the element is touching the trailing edge.
     */
    void CheckIfTrailingEdgeElement(Element& rElement);

    /**
     * @brief This method selects the elements downstream the trailing edge as
     * potentially wake elements.
     */
    bool CheckIfPotentiallyWakeElement(const Element& rElement) const;

    /**
     * @brief This method computes the distance of the element nodes to the wake.
     */
    const BoundedVector<double, 3> ComputeNodalDistancesToWake(const Element& rElement) const;

    /**
     * @brief This method adds the trailing edge elements in the
     * trailing_edge_sub_model_part.
     */
    void AddTrailingEdgeAndWakeElements(std::vector<std::size_t>& rWakeElementsOrderedIds);

    /**
     * @brief This method marks the elements touching the trailing edge from below as kutta. 
     */
    void MarkKuttaElements() const;

    /**
     * @brief This methid marks the trailing edge element that is further downstream as wake.
     */
    void MarkWakeTrailingEdgeElement() const;

    /**
     * @brief This method checks if the element is cut by the wake.
     */
    bool CheckIfTrailingEdgeElementIsCutByWake(const Element& rElement) const;

    /**
     * @brief This method computes the distance from the trailing edge to an input point.
     */
    const BoundedVector<double, 3> ComputeDistanceFromTrailingEdgeToPoint(const Point& rInputPoint) const;
    ///@}

}; // Class Define2DWakeOperation

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  Define2DWakeOperation& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Define2DWakeOperation& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.