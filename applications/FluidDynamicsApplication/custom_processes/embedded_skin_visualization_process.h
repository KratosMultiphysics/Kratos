//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

#ifndef KRATOS_EMBEDDED_SKIN_VISUALIZATION_PROCESS_H
#define KRATOS_EMBEDDED_SKIN_VISUALIZATION_PROCESS_H

// System includes
#include <string>
#include <iostream>
#include <unordered_map>

// External includes
#include <boost/functional/hash.hpp> //TODO: remove this dependence when Kratos has en internal one
#include <boost/unordered_map.hpp>   //TODO: remove this dependence when Kratos has en internal one

// Project includes
#include "includes/define.h"
#include "includes/key_hash.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "utilities/divide_geometry.h"
#include "modified_shape_functions/modified_shape_functions.h"

// Application includes


namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

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

/**
 * @class EmbeddedSkinVisualizationProcess
 * @ingroup FluidDynamicsApplication
 * @brief This process saves the intersected elements in a different model part for its visualization.
 * @details For a given model part, this process checks if its elements are intersected. If they are,
 * calls the corresponding splitting utility to get the subgeometries that conform the splitting
 * pattern. Then, it saves that subgeometries in another model part for its visualization.
 * It has to be mentioned that all the origin model part nodes are kept. Then, the unique nodes
 * that are created are that ones in the intersection edge points.
 * Finally, the values in the visualization model part are computed using the corresponding
 * modify shape functions utility.
 * @author Ruben Zorrilla
 */
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) EmbeddedSkinVisualizationProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    struct Hash{
        std::size_t operator()(const std::pair<unsigned int,bool>& k) const{
            std::size_t h1 = std::hash<unsigned int>()(std::get<0>(k));
            std::size_t h2 = std::hash<bool>()(std::get<1>(k));
            return h1 ^ (h2 << 1);
        }
    };

    struct KeyEqual{
        bool operator()(const std::pair<unsigned int,bool>& lhs, const std::pair<unsigned int,bool>& rhs) const{
            return ((std::get<0>(lhs) == std::get<0>(rhs)) && (std::get<1>(lhs) == std::get<1>(rhs)));
        }
    };

    /// Pointer definition of EmbeddedSkinVisualizationProcess
    KRATOS_CLASS_POINTER_DEFINITION(EmbeddedSkinVisualizationProcess);

    typedef std::unordered_map<
        Node<3>::Pointer,
        std::tuple< const Node<3>::Pointer, const Node<3>::Pointer, const double, const double >,
        SharedPointerHasher<Node<3>::Pointer>,
        SharedPointerComparator<Node<3>::Pointer> > CutNodesMapType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.

    /**
     * @brief Default constructor
     * @param rModelPart The origin model part
     * @param rVisualizationModelPart The visualization model part to be filled
     * @param rVisualizationScalarVariables Scalar variables to be interpolated in the visualization model part
     * @param rVisualizationVectorVariables Vector variables to be interpolated in the visualization model part
     * @param rVisualizationComponentVariables Component variables to be interpolated in the visualization model part
     * @param rShapeFunctions Shape functions type. So far "standard" and "ausas" are implemented
     * @param ReformModelPartAtEachTimeStep Redo visualization model part at each time step flag
     */
    EmbeddedSkinVisualizationProcess(
        ModelPart& rModelPart,
        ModelPart& rVisualizationModelPart,
        const std::vector<Variable< double> >& rVisualizationScalarVariables,
        const std::vector<Variable< array_1d<double, 3> > >& rVisualizationVectorVariables,
        const std::vector<VariableComponent<VectorComponentAdaptor< array_1d< double, 3> > > >& rVisualizationComponentVariables,
        const std::string& rShapeFunctions = "standard",
        const bool ReformModelPartAtEachTimeStep = false);

    /**
     * @brief Constructor with Kratos parameters
     * @param rModelPart The origin model part
     * @param rVisualizationModelPart The visualization model part to be filled
     * @param rParameters Kratos parameters encapsulating the settings
     */
    EmbeddedSkinVisualizationProcess(
        ModelPart& rModelPart,
        ModelPart& rVisualizationModelPart,
        Parameters& rParameters);

    /// Destructor.
    ~EmbeddedSkinVisualizationProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    void ExecuteInitialize() override;

    void ExecuteInitializeSolutionStep() override;

    void ExecuteBeforeOutputStep() override;

    void ExecuteAfterOutputStep() override;

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "EmbeddedSkinVisualizationProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "EmbeddedSkinVisualizationProcess";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}

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

    ModelPart&                                                                          mrModelPart;
    ModelPart&                                                                          mrVisualizationModelPart;

    CutNodesMapType                                                                     mCutNodesMap;

    ModelPart::ElementsContainerType                                                    mNewElementsPointers;

    std::vector<Variable< double> >                                                     mVisualizationScalarVariables;
    std::vector<Variable< array_1d<double, 3> > >                                       mVisualizationVectorVariables;
    std::vector<VariableComponent<VectorComponentAdaptor< array_1d< double, 3> > > >    mVisualizationComponentVariables;

    std::string                                                                         mShapeFunctions;

    bool                                                                                mReformModelPartAtEachTimeStep = false;
    bool                                                                                mSetVisualizationMesh = true;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Create a Visualization Mesh object
     * Fills the visualization model part with the corresponding geometrical entities
     */
    void CreateVisualizationMesh();

    /**
     * Computes the interpolation in the new (interface) nodes
     */
    void ComputeNewNodesInterpolation();

    /**
     * Copies the non-interface nodes from the origin model part to the visualization one
     */
    void CopyOriginNodes();

    /**
     * Copies the non-interface nodes values from the origin model part to the visualization one
     */
    void CopyOriginNodalValues();

    /**
     * Creates the new geometrical entities (elements and conditions) in the visualization model part
     */
    void CreateVisualizationGeometries();

    /**
     * Checks wether the element is split or not
     * @param pGeometry Pointer to the element geometry
     * @param rNodalDistances Vector containing the distance values
     * @return True if it is split and false if not
     */
    bool ElementIsSplit(
        Geometry<Node<3>>::Pointer pGeometry,
        const Vector &rNodalDistances);

    /**
     * Checks wether the element is in the positive side or not
     * @param pGeometry Pointer to the element geometry
     * @param rNodalDistances Vector containing the distance values
     * @return True if it is split and false if not
     */
    bool ElementIsPositive(
        Geometry<Node<3>>::Pointer pGeometry,
        const Vector &rNodalDistances);

    /**
     * Sets the distance values. If Ausas shape functions are used,
     * it takes the ELEMENTAL_DISTANCES. Otherwise, the nodal ones
     * @param ItElem Element iterator
     * @return Vector containing the distance values
     */
    const Vector SetDistancesVector(ModelPart::ElementIterator ItElem);

    /**
     * Sets the the modified shape functions utility according to the
     * distance values.
     * @param pGeometry Pointer to the element geometry
     * @param rNodalDistances Vector containing the distance values
     * @return A pointer to the modified shape functions utility
     */
    ModifiedShapeFunctions::Pointer SetModifiedShapeFunctionsUtility(
        const Geometry<Node<3>>::Pointer pGeometry,
        const Vector &rNodalDistances);

    /**
     * Sets the new interface condition geometry
     * @param rOriginGeometryType Interface subgeometry type
     * @param rNewNodesArray Nodes that conform the new interface geometry
     * @return A pointer to the new geometry
     */
    Geometry< Node<3> >::Pointer SetNewConditionGeometry(
        const GeometryData::KratosGeometryType &rOriginGeometryType,
        const Condition::NodesArrayType &rNewNodesArray);

    /**
     * @brief Sets the visualization properties (one for the positive side and one for the negative)
     * Set the properties for the new elements depending if they are in the positive or negative side of the cut.
     * By doing this, two different layers will be created when printing this model part GiD.
     * @return std::tuple< Properties::Pointer , Properties::Pointer > Tuple containing two properties pointers
     */
    std::tuple< Properties::Pointer , Properties::Pointer > SetVisualizationProperties();

    /**
     * @brief Removes the visualization properties
     * When it is required, this function searchs for the visualization properties to remove them.
     */
    void RemoveVisualizationProperties();

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Default constructor.
    EmbeddedSkinVisualizationProcess() = delete;

    /// Assignment operator.
    EmbeddedSkinVisualizationProcess& operator=(EmbeddedSkinVisualizationProcess const& rOther) = delete;

    /// Copy constructor.
    EmbeddedSkinVisualizationProcess(EmbeddedSkinVisualizationProcess const& rOther) = delete;

    ///@}

}; // Class EmbeddedSkinVisualizationProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_EMBEDDED_SKIN_VISUALIZATION_PROCESS_H
