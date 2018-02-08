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

// struct NodeKeyComparor {
//     bool operator()(const Node<3>::Pointer& p_lhs, const Node<3>::Pointer& p_rhs) const {
//         if (p_lhs->Id() != p_rhs->Id()){
//             return false;
//         }

//         return true;
//     }
// };

// struct NodeKeyHasher {
//     std::size_t operator()(const Node<3>::Pointer& k) const {
//         return k->Id();
//     }
// };

///@}
///@name Kratos Classes
///@{

/// This process saves the intersected elements in a different model part for its visualization.
/** For a given model part, this process checks if its elements are intersected. If they are, 
 *  calls the corresponding splitting utility to get the subgeometries that conform the splitting
 *  pattern. Then, it saves that subgeometries in another model part for its visualization.
 * 
 *  It has to be mentioned that all the origin model part nodes are kept. Then, the unique nodes
 *  that are created are that ones in the intersection edge points.
 * 
 *  Finally, the values in the visualization model part are computed using the corresponding 
 *  modify shape functions utility. 
 */
class EmbeddedSkinVisualizationProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

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
    EmbeddedSkinVisualizationProcess(
        ModelPart& rModelPart,
        ModelPart& rVisualizationModelPart,
        const std::vector<Variable< double> >& rVisualizationScalarVariables,
        const std::vector<Variable< array_1d<double, 3> > >& rVisualizationVectorVariables,
        const std::vector<VariableComponent<VectorComponentAdaptor< array_1d< double, 3> > > >& rVisualizationComponentVariables,
        const std::string& rShapeFunctions = "standard",
        const bool ReformModelPartAtEachTimeStep = false);

    /// Constructor with Kratos parameters.
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

    void ExecuteBeforeSolutionLoop() override;

    void ExecuteInitializeSolutionStep() override;

    void ExecuteBeforeOutputStep() override;

    void ExecuteFinalizeSolutionStep() override;

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

    const bool ElementIsSplit(
        Geometry<Node<3>>::Pointer pGeometry,
        Vector& rNodalDistances);

    const bool ElementIsPositive(
        Geometry<Node<3>>::Pointer pGeometry);

    ModifiedShapeFunctions::Pointer SetModifiedShapeFunctionsUtility(
        const Geometry<Node<3>>::Pointer pGeometry,
        const Vector &rNodalDistances);

    Geometry< Node<3> >::Pointer SetNewConditionGeometry(
        const GeometryData::KratosGeometryType &rOriginGeometryType,
        const Condition::NodesArrayType &rNewNodesArray);

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
