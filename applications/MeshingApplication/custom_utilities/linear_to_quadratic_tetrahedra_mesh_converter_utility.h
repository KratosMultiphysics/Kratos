// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ariadna Cortés
//

#if !defined(KRATOS_LINEAR_TO_QUADRATIC_TETRAHEDRA_MESH_CONVERTER_UTILITY)
#define  KRATOS_LINEAR_TO_QUADRATIC_TETRAHEDRA_MESH_CONVERTER_UTILITY

// System includes

/* Project includes */
#include "custom_utilities/local_refine_tetrahedra_mesh.hpp"
#include "geometries/tetrahedra_3d_10.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
///@name Kratos Globals
///@{
///@}
///@name Type Definitions
///@{
///@}
///@name  Enum's
///@{
///@}
///@name  Functions
///@{
///@}
///@name Kratos Classes
///@{
class KRATOS_API(MESHING_APPLICATION) LinearToQuadraticTetrahedraMeshConverter : public LocalRefineTetrahedraMesh
{
public:

    ///@name Type Definitions
    ///@{
    ///@}

    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef GeometryType::Pointer GeometryPtrType;

     /// Pointer definition of VoxelInsideVolume
    KRATOS_CLASS_POINTER_DEFINITION( LinearToQuadraticTetrahedraMeshConverter );

    ///@name Life Cycle
    ///@{

    /// Default constructors
    LinearToQuadraticTetrahedraMeshConverter(ModelPart& ModelPart) : LocalRefineTetrahedraMesh(ModelPart)
    {

    }

    /// Destructor
    ~LinearToQuadraticTetrahedraMeshConverter() 
    = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
        
    /**
    * Replaces the geometry of the mesh of Tetrahedra3D4 elements. Resulting is a mesh of Tetrahedra3D10 created
    * by adding intermediate nodes to each Tetrahedra3D4. Does the same for conditions, replacing the Triangle3D3 
    * by Triangle3D6
    * @param RefineOnReference: Boolean that defines if refine or not the mesh according to the reference
    * @param InterpolateInternalVariables: Boolean that defines if to interpolate or not the internal variables
    */
    void LocalConvertLinearToQuadraticTetrahedraMesh(
        bool RefineOnReference, 
        bool InterpolateInternalVariables);
  
private:
    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    
    /**
    * Creates a new tetrahedra3D10
    * @return The new tetrahedra
    */
    Tetrahedra3D10<Node> GenerateTetrahedra(
        ModelPart& rThisModelPart, 
        const std::vector<int>& rNodeIds);

    /**
    * Creates a new triangle3D6
    * @return The new triangle
    */
    Triangle3D6<Node> GenerateTriangle3D6(
        ModelPart& rThisModelPart, 
        const array_1d<int, 6>& rNodeIds);

        /**
    * It erases the old Tetrahedra3D4 elements and it creates the new Tetrahedra3D10 ones
    * @param Coord: The compressed matrix containing at (i,j) the id of the node created between nodes i,j
    * @param New_Elements: The new elements created
    * @param InterpolateInternalVariables: A boolean that defines if it is necessary to interpolate the internal variables
    * @return rThisModelPart: The model part of the model (it is the input too)
    */

    void EraseOldElementAndCreateNewElement(
            ModelPart& rThisModelPart,
            const compressed_matrix<int>& Coord,
            PointerVector< Element >& NewElements,
            bool InterpolateInternalVariables
    ) override;

    /**
    * It replaces the old Tetrahedra3D4 elements in its corresponding submodelpart by its new Tetrahedra3D10
    * @param rThisModelPart: The modelpart or submodelpart to replace the elements
    */
    void ReplaceElementsInSubModelPart(ModelPart& rThisModelPart);

    /**
    * Remove the old Trangle3D3 conditions and replace them by the new Trangle3D6 ones
    * @param Coord: The coordinates of the nodes of the geometry
    * @return rThisModelPart: The model part of the model (it is the input too)
    */

    void EraseOldConditionsAndCreateNew(
            ModelPart& rThisModelPart,
            const compressed_matrix<int>& Coord
            ) override;

    /**
    * It replaces the old Triangle3D3 elements in its corresponding submodelpart by its new Triangle3D6
    * @param rThisModelPart: The modelpart or submodelpart to replace the elements
    */
    void ReplaceConditionsInSubModelPart(ModelPart& rThisModelPart);

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Private LifeCycle
    ///@{
    ///@}

}; // Class LinearToQuadraticTetrahedraMeshConverter

} // namespace Kratos.

#endif // KRATOS_LINEAR_TO_QUADRATIC_TETRAHEDRA_MESH_CONVERTER_UTILITY  defined
