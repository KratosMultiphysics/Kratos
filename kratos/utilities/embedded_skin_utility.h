//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined(KRATOS_GENERATE_EMBEDDED_SKIN_UTILITY_H_INCLUDED )
#define  KRATOS_GENERATE_EMBEDDED_SKIN_UTILITY_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "geometries/geometry_data.h"
#include "modified_shape_functions/modified_shape_functions.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/divide_geometry.h"
#include "utilities/math_utils.h"


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

/// Utility to compute the skin representation from a distance function.
/** Provided either a continuous or discontinuous distance function, this
 *  utility reconstructs the skin representation coming from such distance
 *  function. This is done by computing the element intersections and saving
 *  them in an empty provided model part. Note that such skin representation
 *  is discontinuous even for a provided continuous distance field.
 */
template<std::size_t TDim>
class KRATOS_API(KRATOS_CORE) EmbeddedSkinUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EmbeddedSkinUtility
    KRATOS_CLASS_POINTER_DEFINITION(EmbeddedSkinUtility);

    typedef std::unordered_map<
        Node<3>::Pointer,
        std::tuple< const Element::Pointer, const unsigned int >,
        SharedPointerHasher<Node<3>::Pointer>,
        SharedPointerComparator<Node<3>::Pointer> > EdgeNodesMapType;

    ///@}
    ///@name  Enum's
    ///@{

    enum LevelSetTypeEnum
    {
        Continuous = 1,
        Discontinuous = 2
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    EmbeddedSkinUtility(
        ModelPart &rModelPart,
        ModelPart &rSkinModelPart,
        const std::string LevelSetType = "continuous",
        const std::vector<std::string> InterpolatedSkinVariables = {}) :
        mrModelPart(rModelPart),
        mrSkinModelPart(rSkinModelPart),
        mLevelSetType(LevelSetType == "continuous" ? Continuous : Discontinuous),
        mrConditionPrototype(KratosComponents<Condition>::Get(this->GetConditionType())),
        mInterpolatedSkinVariables(InterpolatedSkinVariables) {};

    /// Destructor.
    virtual ~EmbeddedSkinUtility() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Call to generate the embedded skin model part
     * This method collects all the operations required to generate
     * the embedded skin model part. The new geometries will be stored
     * inside the skin model part provided in the constructor.
     */
    void GenerateSkin();

    /**
     * @brief InterpolateMeshVariableToSkin double specialization
     * Double type specialization of the InterpolateMeshVariableToSkin method
     * @param rMeshVariable background mesh origin variable
     * @param rSkinVariable skin mesh destination variable
     */
    void InterpolateMeshVariableToSkin(
        const Variable<double> &rMeshVariable,
		const Variable<double> &rSkinVariable);

    /**
     * @brief InterpolateMeshVariableToSkin array specialization
     * Array type specialization of the InterpolateMeshVariableToSkin method
     * @param rMeshVariable background mesh origin variable
     * @param rSkinVariable skin mesh destination variable
     */
    void InterpolateMeshVariableToSkin(
        const Variable<array_1d<double,3>> &rMeshVariable,
		const Variable<array_1d<double,3>> &rSkinVariable);

    /**
     * @brief Discontinuous InterpolateMeshVariableToSkin double specialization
     * Double type specialization of the InterpolateMeshVariableToSkin method
     * for discontinuous level set type formulation
     * @param rMeshVariable background mesh origin variable
     * @param rSkinVariable skin mesh destination variable
     * @param rInterfaceSide interface side ("positive" or "negative")
     */
    void InterpolateDiscontinuousMeshVariableToSkin(
        const Variable<double> &rMeshVariable,
		const Variable<double> &rSkinVariable,
        const std::string &rInterfaceSide);

    /**
     * @brief Discontinuous InterpolateMeshVariableToSkin array specialization
     * Array type specialization of the InterpolateMeshVariableToSkin method
     * for discontinuous level set type formulation
     * @param rMeshVariable background mesh origin variable
     * @param rSkinVariable skin mesh destination variable
     * @param rInterfaceSide interface side ("positive" or "negative")
     */
    void InterpolateDiscontinuousMeshVariableToSkin(
        const Variable<array_1d<double,3>> &rMeshVariable,
		const Variable<array_1d<double,3>> &rSkinVariable,
        const std::string &rInterfaceSide);

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    ///@}
    ///@name Friends
    ///@{


    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart &mrModelPart;
    ModelPart &mrSkinModelPart;
    EdgeNodesMapType mEdgeNodesMap;
    const LevelSetTypeEnum mLevelSetType;
    const Condition &mrConditionPrototype;
    std::vector<std::string> mInterpolatedSkinVariables;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Geometry clear operation
     * This method is called before any construction of the skin.
     * It removes the existent nodes, elements and conditions in
     * the provided skin model part.
     */
    void Clear();

    /**
     * @brief Computes the skin entities of one element
     * For an intersected element, this method computes the skin
     * intersection geometries.
     * @param pElement element of interest
     * @param rNodalDistances vector containing the element node distances
     * @param rTempNodeId temporal id for the new nodes
     * @param rTempCondId temporal id for the new conditions
     * @param pCondProp pointer to the properties for the new skin conditions
     * @param rNewNodesVect vector to temporary store the new skin nodes
     * @param rNewCondsVect vector to temporary store the new skin conditions
     */
    void ComputeElementSkin(
        const Element::Pointer pElement,
        const Vector &rNodalDistances,
        unsigned int &rTempNodeId,
        unsigned int &rTempCondId,
        Properties::Pointer pCondProp,
        ModelPart::NodesContainerType &rNewNodesVect,
        ModelPart::ConditionsContainerType &rNewCondsVect);

    /**
     * @brief Checks if an element is split
     * This method checks if an element geometry is split
     * @param rGeometry geometry of the element of interest
     * @param rNodalDistances vector containing the element node distances
     * @return true if the element is split
     * @return false if the element is not split
     */
    bool inline ElementIsSplit(
        const Geometry<Node<3>> &rGeometry,
        const Vector &rNodalDistances);

    /**
     * @brief InterpolateMeshVariableToSkin specialization method
     * For a provided set of variables, this method interpolates the values
     * from the brackground fluid mesh to an embedded skin mesh. This can
     * be done for either the positive or negative sides of the interface.
     * @tparam TVarType variable type of the variable to be interpolated
     * @param rMeshVariable variable in the background mesh to interpolate from
     * @param rSkinVariable variable in the skin model part to interpolate to
     * @param rInterfaceSide interface side in where the shape functions
     * are to be computed. Must be either "positive" or "negative"
     */
    template<class TVarType>
    void InterpolateMeshVariableToSkinSpecialization(
        const Variable<TVarType> &rMeshVariable,
		const Variable<TVarType> &rSkinVariable,
        const std::string &rInterfaceSide = "positive")
    {
		// Check requested variables
		KRATOS_ERROR_IF((mrModelPart.NodesBegin())->SolutionStepsDataHas(rMeshVariable) == false)
			<< "Mesh model part solution step data missing variable: " << rMeshVariable << std::endl;
		KRATOS_ERROR_IF((mrSkinModelPart.NodesBegin())->SolutionStepsDataHas(rSkinVariable) == false)
			<< "Generated skin model part solution step data missing variable: " << rSkinVariable << std::endl;

        // Check that the mesh model part has elements
        KRATOS_ERROR_IF(mrModelPart.NumberOfElements() == 0) << "Mesh model part has no elements.";

        // Loop the edge intersection nodes to set their values
        unsigned int i_edge;
        Element::Pointer p_elem;
        #pragma omp parallel for private (i_edge, p_elem)
        for (int i_node = 0; i_node < static_cast<int>(mrSkinModelPart.NumberOfNodes()); ++i_node) {
            // Get the current node
            auto it_node = mrSkinModelPart.NodesBegin() + i_node;
            Node<3>::Pointer p_node = &(*it_node);

            // Search for the current node in the intersected edges map
            const auto i_node_info = mEdgeNodesMap.find(p_node);
            if (i_node_info != mEdgeNodesMap.end()){
                // Get the cut node info from the map tuple iterator
                std::tie(p_elem, i_edge) = std::get<1>(*i_node_info);

                // Set the modified shape functions for the parent element
                const auto p_elem_geom = p_elem->pGetGeometry();
                const auto elem_dist = this->SetDistancesVector(*p_elem);
                const auto p_mod_sh_func = pCreateModifiedShapeFunctions(p_elem_geom, elem_dist);

                // Get interface modified shape function values
                const auto edge_sh_func = this->GetModifiedShapeFunctionsValuesOnEdge(
                    p_mod_sh_func,
                    rInterfaceSide);

                // Compute the interpolation
                const auto edge_N = row(edge_sh_func, i_edge);
                const auto &r_elem_geom = p_elem->GetGeometry();
                auto &r_value = it_node->FastGetSolutionStepValue(rSkinVariable);
                r_value = rSkinVariable.Zero();
                for (unsigned int i_elem_node = 0; i_elem_node < r_elem_geom.PointsNumber(); ++i_elem_node) {
                    r_value += edge_N[i_elem_node] * r_elem_geom[i_elem_node].FastGetSolutionStepValue(rMeshVariable);
                }
            } else{
                KRATOS_ERROR << "Intersected edge node " << it_node->Id() << " not found in intersected edges nodes map" << std::endl;
            }
        }
    };

    /**
     * @brief Renumber and saves the new skin entities
     * This method renumbers the new skin geometrical entities (MPI compatible)
     * and add them to the provided skin model part.
     * @param rNewNodesVect vector that stores the new skin nodes
     * @param rNewCondsVect vector that stores the new skin conditions
     */
    void RenumberAndAddSkinEntities(
        const ModelPart::NodesContainerType &rNewNodesVect,
        const ModelPart::ConditionsContainerType &rNewCondsVect);

    /**
     * @brief Set the Distances Vector object
     * For a given element, this method sets the vector containing the
     * element node distances.
     * @param ItElem iterator to the element of interest
     * @return const Vector vector containing the element node distances
     */
    const Vector SetDistancesVector(const Element &rElement);

    /**
     * Sets the the divide geometry utility according to the geometry type.
     * @param pGeometry Pointer to the element geometry
     * @param rNodalDistances Vector containing the distance values
     * @return A pointer to the divide geometry utility
     */
    DivideGeometry::Pointer SetDivideGeometryUtility(
        const Geometry<Node<3>> &rGeometry,
        const Vector &rNodalDistances);

    /**
     * Creates the new interface condition geometry
     * @param rOriginGeometryType Interface subgeometry type
     * @param rNewNodesArray Nodes that conform the new interface geometry
     * @return A pointer to the new geometry
     */
    Geometry< Node<3> >::Pointer pCreateNewConditionGeometry(
        const GeometryData::KratosGeometryType &rOriginGeometryType,
        const Condition::NodesArrayType &rNewNodesArray);

    /**
     * @brief Creates a pointer to a new skin condition
     * From the split pattern of an intersected element, this method
     * creates and returns a pointer to a new skin condition.
     * @param rOriginGeometryType GeometryType of the condition to be created
     * @param rNewNodesArray array containing pointers to the nodes that will
     * conform the condition
     * @param rConditionId new condition identifier
     * @param pConditionProperties pointer to the new condition properties
     * @return Condition::Pointer pointer to a new skin condition
     */
    Condition::Pointer pCreateNewCondition(
        const GeometryData::KratosGeometryType &rOriginGeometryType,
        const Condition::NodesArrayType &rNewNodesArray,
        const unsigned int &rConditionId,
        const Properties::Pointer pConditionProperties);

    /**
     * @brief Get the condition type
     * Depending on the dimension template argument, this method returns
     * the condition type name (LineCondition2D2N or SurfaceCondition3D3N)
     * @return std::sting condition type name
     */
    static const std::string GetConditionType();

    /**
     * @brief Set the Skin Entities Properties
     * This method checks which is the last properties id. and
     * sets a new one accordingly to be used as skin conditions property
     * @return Properties::Pointer pointer to the new skin entities property
     */
    Properties::Pointer SetSkinEntitiesProperties();

    /**
     * @brief Creates a pointer to the Modified Shape Functions
     * For an intersected element, sets the modified shape functions utility
     * @param pGeometry Pointer to the intersected element geometry
     * @param rNodalDistances Vector containing the nodal distances
     * @return ModifiedShapeFunctions::UniquePointer Unique pointer
     * to the current element modified shape functions utility
     */
    ModifiedShapeFunctions::UniquePointer pCreateModifiedShapeFunctions(
        const Geometry<Node<3>>::Pointer pGeometry,
        const Vector& rNodalDistances);

    /**
     * @brief Get the Modified Shape Functions Values object
     * This method returns the shape function values in an element intersection
     * @param rpModifiedShapeFunctions pointer to the modified shape functions util
     * @param rInterfaceSide interface side in where the shape functions
     * are to be computed. Must be either "positive" or "negative"
     * @return Matrix matrix containing the split shape function values
     */
    Matrix GetModifiedShapeFunctionsValues(
        const ModifiedShapeFunctions::UniquePointer &rpModifiedShapeFunctions,
        const std::string &rInterfaceSide) const;

    /**
     * @brief Get the Modified Shape Functions Values On Edge object
     * This method returns the shape function values in the intersected edges
     * @param rpModifiedShapeFunctions pointer to the modified shape functions util
     * @param rInterfaceSide interface side in where the shape functions
     * are to be computed. Must be either "positive" or "negative"
     * @return Matrix matrix containing the split shape function values
     */
    Matrix GetModifiedShapeFunctionsValuesOnEdge(
        const ModifiedShapeFunctions::UniquePointer &rpModifiedShapeFunctions,
        const std::string &rInterfaceSide) const;

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    EmbeddedSkinUtility& operator=(EmbeddedSkinUtility const& rOther) = delete;

    /// Copy constructor.
    EmbeddedSkinUtility(EmbeddedSkinUtility const& rOther) = delete;

    ///@}
}; // Class EmbeddedSkinUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}
}  // namespace Kratos.

#endif // KRATOS_GENERATE_EMBEDDED_SKIN_UTILITY_H_INCLUDED  defined
