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
class KRATOS_API(CONVECTION_DIFFUSION_APPLICATION) GaussPointErrorUtility
{
public:
    ///@name Type Definitions
    ///@{

    /**
     * @brief Enum class with the available level set types
     * Auxiliary enum class to store the available level set types
     * Continuous: standard nodal based continuous level set function
     * Discontinuous: elementwise discontinuous level set function
     */
    enum class LevelSetType
    {
        Continuous,
        Discontinuous
    };

    /**
     * @brief Enum class with the available shape functions type
     * Auxiliary enum class to store the available shape functions types
     * Available shape functions are the standard linear FE shape functions (Standard)
     * and the Aausas FE space (see Ausas et. al. 2010 https://doi.org/10.1016/j.cma.2009.11.011)
     */
    enum class ShapeFunctionsType
    {
        Ausas,
        Standard
    };

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

    /**
     * @brief Check and return the level set type
     * This method checks and return the user provided level set type
     * @param rParameters Kratos parameters encapsulating the settings. These settings are assumed to be already validated.
     * @param rLevelSetType The validated level set type
     */
    static void CheckAndSetLevelSetType(
        const Parameters rParameters,
        LevelSetType& rLevelSetType);

    /**
     * @brief Check and return the shape functions
     * This method checks and return the user provided shape functions type
     * @param rParameters Kratos parameters encapsulating the settings. These settings are assumed to be already validated.
     * @param rShapeFunctionsType The validated shape functions type
     */
    static void CheckAndSetShapeFunctionsType(
        const Parameters rParameters,
        ShapeFunctionsType& rShapeFunctionsType);

    /**
     * @brief Checks and returns the distance variable name
     * This method checks the user provided distance variable name
     * If the variable is not prescribed, it returns a default one according to the selected level set type
     * @param rParameters Kratos parameters encapsulating the settings. These settings are assumed to be already validated.
     * @return std::string The distance variable name
     */
    static const std::string CheckAndReturnDistanceVariableName(
        const Parameters rParameters,
        const LevelSetType& rLevelSetType);

    /// Constructor.

    /**
     * @brief Default constructor
     * @param rModelPart The origin model part
     * @param rLevelSetType Level set type. So far "continuous" and "discontinuous" are implemented
     * @param rShapeFunctionsType Shape functions type. So far "standard" and "ausas" are implemented
     */
    GaussPointErrorUtility(
        ModelPart& rModelPart,
        const LevelSetType& rLevelSetType,
        const ShapeFunctionsType& rShapeFunctionsType);

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

    double ExecuteGradient();

    double ExecuteOnConditions(ModelPart& rModelPart);

    double ExecuteOnConditionsSolution(ModelPart& rModelPart);

    // std::array<double,2> Execute();

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

    // Level set type. Current available options continuous and discontinuous
    const LevelSetType mLevelSetType;

    // Shape functions type. Current available options ausas and standard
    const ShapeFunctionsType mShapeFunctionsType;

    // Pointer to the variable that stores the nodal level set function
    const Variable<double>* mpNodalDistanceVariable;

    // Pointer to the variable that stores the elemental level set function
    const Variable<Vector>* mpElementalDistanceVariable;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Auxiliary get double value method
     * This is an auxiliary method to get scalar values from the nodal databases
     * It is specialized for the historical and non-historical databases
     * @tparam IsHistorical Template argument to indicate the database. Historical (true) and non-historical (false)
     * @param rNode Node from which the values are retrieved
     * @param rVariable Scalar variable to be retrieved
     * @return double& Reference to the retrieved value
     */
    template<bool IsHistorical>
    double& AuxiliaryGetValue(
        Node<3>& rNode,
        const Variable<double>& rVariable);

    /**
     * @brief Auxiliary get vector value method
     * This is an auxiliary method to get vector values from the nodal databases
     * It is specialized for the historical and non-historical databases
     * @tparam IsHistorical Template argument to indicate the database. Historical (true) and non-historical (false)
     * @param rNode Node from which the values are retrieved
     * @param rVariable Vector variable to be retieved
     * @return array_1d<double,3>& Reference to the retrieved value
     */
    template<bool IsHistorical>
    array_1d<double,3>& AuxiliaryGetValue(
        Node<3>& rNode,
        const Variable<array_1d<double,3>>& rVariable);

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

    double CalculatePressureExactSolution(const array_1d<double,3>& rCoords);

    double CalculateTemperatureExactSolution(const array_1d<double,3>& rCoords);

    array_1d<double,3> CalculateTemperatureExactSolutionGradient(const array_1d<double,3>& rCoords);

    array_1d<double,3> CalculateVelocityExactSolution(const array_1d<double,3>& rCoords);

    double CalculateTemperatureFluxExactSolution(
        const array_1d<double,3>& rCoords,
        const array_1d<double,3>& rNormal);

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
