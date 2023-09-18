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
//

#ifndef KRATOS_GAUSS_POINT_ERROR_UTILITY_H
#define KRATOS_GAUSS_POINT_ERROR_UTILITY_H

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
 * @class GaussPointErrorUtility
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
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) GaussPointErrorUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GaussPointErrorUtility
    KRATOS_CLASS_POINTER_DEFINITION(GaussPointErrorUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Get the Default Settings object
     * Static method to get the default settings inside the contructor
     * @return Parameters The parameters object containing the default settings
     */
    static Parameters GetDefaultSettings();

    /// Constructor.

    /**
     * @brief Constructor with Kratos parameters
     * @param rModelPart The origin model part
     * @param rVisualizationModelPart The visualization model part to be filled
     * @param rParameters Kratos parameters encapsulating the settings
     */
    GaussPointErrorUtility(
        ModelPart& rModelPart,
        Parameters rParameters);

    /**
     * @brief Constructor with Kratos parameters and Model container
     * @param rModel The Model container
     * @param rParameters Kratos parameters encapsulating the settings
     */
    GaussPointErrorUtility(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    ~GaussPointErrorUtility() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    double Execute();

    int Check();

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "GaussPointErrorUtility" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const {rOStream << "GaussPointErrorUtility";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const {}

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

    // Reference to the origin model part
    ModelPart& mrModelPart;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void CalculateExactSolution(
        const array_1d<double,3>& rCoords,
        array_1d<double,3>& rExactSolution);

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
    GaussPointErrorUtility() = delete;

    /// Assignment operator.
    GaussPointErrorUtility& operator=(GaussPointErrorUtility const& rOther) = delete;

    /// Copy constructor.
    GaussPointErrorUtility(GaussPointErrorUtility const& rOther) = delete;

    ///@}

}; // Class GaussPointErrorUtility

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_GAUSS_POINT_ERROR_UTILITY_H
