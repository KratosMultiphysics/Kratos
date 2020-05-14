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

    /**
     * @brief Get the Default Settings object
     * Static method to get the default settings inside the contructor
     * @return Parameters The parameters object containing the default settings
     */
    static Parameters GetDefaultSettings();

    /**
     * @brief Create a And Prepare Visualization Model Part object
     * This method creates the visualization model part and prepares it for the computation
     * @param rModel The model container
     * @param rParameters Kratos parameters encaptulating the settings. These settings are assumed to be already validated.
     * @return ModelPart Visualization model part
     */
    static ModelPart& CreateAndPrepareVisualizationModelPart(
        Model& rModel,
        const Parameters rParameters
    );

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

    /**
     * @brief Check and fill a visualization variable list
     * This method checks the user provided variable list and saves it in the corresponding vector list
     * @tparam TDataType The variable data type
     * @param rParameters Kratos parameters encapsulating the settings. These settings are assumed to be alreeeady validated.
     * @param rVariableList The filled variable data list
     */
    template<class TDataType>
    static void FillVariablesList(
        const Parameters rParameters,
        std::vector<const Variable<TDataType>*>& rVariablesList);

    /// Constructor.

    /**
     * @brief Default constructor
     * @param rModelPart The origin model part
     * @param rVisualizationModelPart The visualization model part to be filled
     * @param rVisualizationScalarVariables Scalar variables to be interpolated in the visualization model part
     * @param rVisualizationVectorVariables Vector variables to be interpolated in the visualization model part
     * @param rLevelSetType Level set type. So far "continuous" and "discontinuous" are implemented
     * @param rShapeFunctionsType Shape functions type. So far "standard" and "ausas" are implemented
     * @param ReformModelPartAtEachTimeStep Redo visualization model part at each time step flag
     */
    EmbeddedSkinVisualizationProcess(
        ModelPart& rModelPart,
        ModelPart& rVisualizationModelPart,
        const std::vector<const Variable< double>* >& rVisualizationScalarVariables,
        const std::vector<const Variable< array_1d<double, 3> >* >& rVisualizationVectorVariables,
        const std::vector<const Variable< double>* >& rVisualizationNonHistoricalScalarVariables,
        const std::vector<const Variable< array_1d<double, 3> >* >& rVisualizationNonHistoricalVectorVariables,
        const LevelSetType& rLevelSetType,
        const ShapeFunctionsType& rShapeFunctionsType,
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
        Parameters rParameters);

    /**
     * @brief Constructor with Kratos parameters and Model container
     * @param rModel The Model container
     * @param rParameters Kratos parameters encapsulating the settings
     */
    EmbeddedSkinVisualizationProcess(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    ~EmbeddedSkinVisualizationProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void ExecuteBeforeSolutionLoop() override;

    void ExecuteBeforeOutputStep() override;

    void ExecuteAfterOutputStep() override;

    int Check() override;

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

    // Unordered map to relate the newly created nodes and the origin mesh ones
    CutNodesMapType mCutNodesMap;

    // Container with the splitting newly created elements
    ModelPart::ElementsContainerType mNewElementsPointers;

    // Reference to the origin model part
    ModelPart& mrModelPart;

    // Reference to the visualization model part
    ModelPart& mrVisualizationModelPart;

    // Level set type. Current available options continuous and discontinuous
    const LevelSetType mLevelSetType;

    // Shape functions type. Current available options ausas and standard
    const ShapeFunctionsType mShapeFunctionsType;

    // If true, the visualization model part is created each time step (required in case the level set function is not constant)
    const bool mReformModelPartAtEachTimeStep;

    // Pointer to the variable that stores the nodal level set function
    const Variable<double>* mpNodalDistanceVariable;

    // Pointer to the variable that stores the elemental level set function
    const Variable<Vector>* mpElementalDistanceVariable;

    // Vector containing the scalar variables to be interpolated in the visualization mesh
    const std::vector<const Variable<double>*> mVisualizationScalarVariables;

    // Vector containing the vector variables to be interpolated in the visualization mesh
    const std::vector<const Variable<array_1d<double, 3>>*> mVisualizationVectorVariables;

    // Vector containing the non-historical scalar variables to be interpolated in the visualization mesh
    const std::vector<const Variable<double>*> mVisualizationNonHistoricalScalarVariables;

    // Vector containing the non-historical vector variables to be interpolated in the visualization mesh
    const std::vector<const Variable<array_1d<double, 3>>*> mVisualizationNonHistoricalVectorVariables;

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
     * @brief Create a Visualization Mesh object
     * Fills the visualization model part with the corresponding geometrical entities
     */
    void CreateVisualizationMesh();

    /**
     * Computes the interpolation in the new (interface) nodes
     */
    void ComputeNewNodesInterpolation();

    /**
     * @brief Auxiliary method to calculate the interpolation
     * For a given variables list, this method does the interpolation from to the intersected edges values
     * @tparam TDataType Variable list data type
     * @tparam IsHistorical Template argument to indicate the database. Historical (true) and non-historical (false)
     * @param rpNode Pointer to the visualization node in which the interpolation is calculated
     * @param rpNodeI Pointer to the I node of the intersected edge
     * @param rpNodeJ Pointer to the J node of the intersected edge
     * @param WeightI Weight of the I node value
     * @param WeightJ Weight of the J node value
     * @param rVariablesList List containing the variables to be interpolated
     */
    template<class TDataType, bool IsHistorical>
    void InterpolateVariablesListValues(
        const Node<3>::Pointer& rpNode,
        const Node<3>::Pointer& rpNodeI,
        const Node<3>::Pointer& rpNodeJ,
        const double WeightI,
        const double WeightJ,
        const std::vector<const Variable<TDataType>*>& rVariablesList);
    
    /**
     * Copies the non-interface nodes from the origin model part to the visualization one
     */
    void CopyOriginNodes();

    /**
     * Copies the non-interface nodes values from the origin model part to the visualization one
     */
    void CopyOriginNodalValues();

    /**
     * @brief Copy the values from the origin to the visualization mesh
     * For the nodes that are no created from an intersected edge (i.e. those already existent in the origin model part),
     * this method copies the values of the variables in the provided list from the origin to the visualization model part
     * @tparam TDataType Variable list data type
     * @tparam IsHistorical Template argument to indicate the database. Historical (true) and non-historical (false)
     * @param rItOriginNode Origin node in the origin model part
     * @param rItVisualizationNode Destination node in the visualization model part
     * @param rVariablesList List containing the variables whose values are to be copied
     */
    template<class TDataType, bool IsHistorical>
    void CopyVariablesListValues(
        const ModelPart::NodeIterator& rItOriginNode,
        ModelPart::NodeIterator& rItVisualizationNode,
        const std::vector<const Variable<TDataType>*>& rVariablesList);

    /**
     * Creates the new geometrical entities (elements and conditions) in the visualization model part
     */
    void CreateVisualizationGeometries();

    /**
     * @brief Initializes the non historical database in the visualization model part
     * This method initializes the non-historical variables in the provided variables list
     * @tparam TDataType The data type of the variables in the list
     * @param rNonHistoricalVariablesVector The vector containing the non-historical variables to be initialized
     */
    template<class TDataType>
    void InitializeNonHistoricalVariables(const std::vector<const Variable<TDataType>*>& rNonHistoricalVariablesVector);

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
