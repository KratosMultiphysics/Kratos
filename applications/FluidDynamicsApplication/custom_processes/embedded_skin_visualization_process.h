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

// External includes
#include <boost/functional/hash.hpp> //TODO: remove this dependence when Kratos has en internal one
#include <boost/unordered_map.hpp>   //TODO: remove this dependence when Kratos has en internal one

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "utilities/divide_geometry.h"

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

struct KeyComparor {
    bool operator()(const vector<int> &lhs, const vector<int> &rhs) const {
        if (lhs.size() != rhs.size()){
            return false;
        }

        for (unsigned int i = 0; i < lhs.size(); i++){
            if (lhs[i] != rhs[i]){
                return false;
            }
        }

        return true;
    }
};

struct KeyHasher {
    std::size_t operator()(const vector<int> &k) const {
        return boost::hash_range(k.begin(), k.end());
    }
};

///@}
///@name Kratos Classes
///@{

/// This process saves the intersected elements in a different model part for its visualization.
/** For a given model part, this process checks if its elements are intersected. If they are, 
 *  calls the corresponding splitting utility to get the subgeometries that conform the splitting
 *  pattern. Then, it saves the such subgeometries in another model part for its visualization.
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

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    EmbeddedSkinVisualizationProcess(
        const ModelPart& rModelPart,
        ModelPart& rVisualizationModelPart);

    /// Constructor with Kratos parameters.
    EmbeddedSkinVisualizationProcess(
        const ModelPart& rModelPart,
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

    const ModelPart&        mrModelPart;
    ModelPart&              mrVisualizationModelPart;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    DivideGeometry::Pointer GetGeometrySplitUtility(
        const Geometry<Node<3>> &rGeometry,
        const Vector &rNodalDistances);

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
