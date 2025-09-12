//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "geometries/geometry.h"

namespace Kratos
{
///@name Type Definitions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class TetrahedralMeshOrientationCheck
 * @ingroup KratosCore
 * @brief Check a triangular or tetrahedral mesh to ensure that local connectivities follow the expected convention.
 * @details This process checks all elements to verify that their Jacobian has positive determinant and face conditions to ensure that all face normals point outwards.
 * @note Note that, as a side result of the procedure used, nodal normals (not normalized) are computed and stored on olution step data NORMAL.
 * @author Pooyan Dadvand
 * @author Riccardo Rossi
 */
class KRATOS_API(KRATOS_CORE) TetrahedralMeshOrientationCheck
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    // DEFINITION OF FLAGS TO CONTROL THE BEHAVIOUR
    KRATOS_DEFINE_LOCAL_FLAG(ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS);
    KRATOS_DEFINE_LOCAL_FLAG(COMPUTE_NODAL_NORMALS);
    KRATOS_DEFINE_LOCAL_FLAG(COMPUTE_CONDITION_NORMALS);
    KRATOS_DEFINE_LOCAL_FLAG(MAKE_VOLUMES_POSITIVE);
    KRATOS_DEFINE_LOCAL_FLAG(ALLOW_REPEATED_CONDITIONS);

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(TetrahedralMeshOrientationCheck);

    /// The definition of the index type
    using IndexType = std::size_t;

    /// The definition of the size type
    using SizeType = std::size_t ;

    /// Definition of the geometry
    using GeometryType = Geometry<Node>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor for TetrahedralMeshOrientationCheck Process
     * @param rModelPart The model part to check.
     * @param ThrowErrors If true, an error will be thrown if the input model part contains malformed elements or conditions.
     * @param Options The flags to be set
     */
    TetrahedralMeshOrientationCheck(
        ModelPart& rModelPart,
        bool ThrowErrors,
        const Flags Options = COMPUTE_NODAL_NORMALS.AsFalse() | COMPUTE_CONDITION_NORMALS.AsFalse() | ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS.AsFalse() | ALLOW_REPEATED_CONDITIONS.AsFalse()
        ):  Process(),
            mrModelPart(rModelPart),
            mThrowErrors(ThrowErrors), //to be changed to a flag
            mrOptions(Options)

    {
    }

    /**
     * @brief Constructor for TetrahedralMeshOrientationCheck Process (simplified)
     * @param rModelPart The model part to check.
     * @param Options The flags to be set
     */
    TetrahedralMeshOrientationCheck(
        ModelPart& rModelPart,
        const Flags Options = COMPUTE_NODAL_NORMALS.AsFalse() | COMPUTE_CONDITION_NORMALS.AsFalse() | ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS.AsFalse() | ALLOW_REPEATED_CONDITIONS.AsFalse()
        ):  Process(),
            mrModelPart(rModelPart),
            mThrowErrors(false),
            mrOptions(Options)
    {
    }

    /// Destructor.
    ~TetrahedralMeshOrientationCheck() override {}

    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Checks elements for positive Jacobian and conditions for outward-pointing normals.
     * @details This method iterates through all tetrahedral elements in the `ModelPart` to verify that their Jacobian is positive,
     * which ensures they are not inverted. It also checks that the face normals of all conditions point outwards from the mesh.
     * If the `ThrowErrors` flag is set in the constructor, an exception will be thrown upon finding an invalid element or condition.
     */
    void Execute() override;

    /**
     * @brief Swaps the orientation of all tetrahedral elements in the mesh.
     * @details This method reverses the connectivity (node order) of every tetrahedral element and condition in the `ModelPart`.
     * This is useful for correcting the global orientation of a mesh.
     */
    void SwapAll();

    /**
     * @brief Swaps the orientation of elements with a negative Jacobian.
     * @details This method is specifically designed to correct inverted elements. It iterates through the elements, and if an element's
     * Jacobian is found to be negative, its connectivity is reversed to correct its orientation and ensure a positive Jacobian.
     */
    void SwapNegativeElements();

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "TetrahedralMeshOrientationCheck";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "TetrahedralMeshOrientationCheck";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        this->PrintInfo(rOStream);
    }

    ///@}
protected:
    ///@name Operations
    ///@{

    /**
     * @brief tells if the element is linear or higher order
     * @param rGeometry The element geometry
     * @return true if linear, false if not
     */
    bool LinearElement(const GeometryType& rGeometry);

    /**
     * @brief tells if the element is supported for the process
     * @param rGeometry The element geometry
     * @return true if supported, false if not
     */
    bool SupportedElement(const GeometryType& rGeometry);

    /**
     * @brief tells if the condition is supported for the process
     * @param rGeometry The condition geometry
     * @return true if supported, false if not
     */
    bool SupportedCondition(const GeometryType& rGeometry);

    /**
     * @brief Finds how many boundaries an element has
     * @param rGeometry The element geometry
     * @return the number of boundaries
     */
    SizeType BoundariesEntitiesNumber(const GeometryType& rGeometry);

    /**
     * @brief how many nodes each of the boundaries of the element contains
     * @param rGeometry The element geometry
     * @return number of nodes per boundary
     */
    SizeType NumberOfNodesInEachBoundary(const GeometryType& rGeometry);

    /**
     * @brief Finds the nodes of each of the boundaries of the element
     * @param rGeometry The element geometry
     * @param rNodesIds the retured matrix with the IDs
     */
    void NodesOfBoundaries(const GeometryType& rGeometry, DenseMatrix<int>& rNodesIds);

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    const bool mThrowErrors;
    Flags mrOptions;

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Checks the orientation of a given geometry and corrects it if necessary.
     * @param rGeom The geometry whose orientation is to be checked.
     * @return `true` if the geometry was originally correctly oriented or was successfully corrected; `false` otherwise.
     * @details This method verifies if the geometry has a positive volume (or "Jacobian" for elements). If the volume is found to be negative, the method
     * attempts to reverse the connectivity of the geometry (e.g., swapping node indices) to correct its orientation.
     */
    bool Orient(GeometryType& rGeom);

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    TetrahedralMeshOrientationCheck& operator=(TetrahedralMeshOrientationCheck const& rOther);

    /// Copy constructor.
    TetrahedralMeshOrientationCheck(TetrahedralMeshOrientationCheck const& rOther);

    ///@}

}; // Class Process

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  TetrahedralMeshOrientationCheck& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const TetrahedralMeshOrientationCheck& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos
