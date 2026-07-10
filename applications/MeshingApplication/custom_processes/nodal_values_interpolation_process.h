// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:         BSD License
//                   license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "meshing_application_variables.h"
#include "processes/process.h"

/* Several includes */
#include "includes/model_part.h"
#include "includes/key_hash.h"
#include "includes/kratos_parameters.h"

/* Utilities */
#include "utilities/binbased_fast_point_locator.h"

/* Tree structures */
// #include "spatial_containers/bounding_volume_tree.h" // k-DOP
#include "spatial_containers/spatial_containers.h" // kd-tree

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// Defining the integers
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

namespace NodalInterpolationFunctions
{
    /**
     * @brief This methoid creates the list of non-historical variables fro nodal interpolation
     * @param rModelPart The model part to extract the non-historical variables
     * @param rVariableList The list of non-historical variables
     */
    void KRATOS_API(MESHING_APPLICATION) GetListNonHistoricalVariables(
        const ModelPart& rModelPart,
        std::unordered_set<std::string>& rVariableList
        );
};

///@}
///@name Kratos Classes
///@{

/**
 * @ingroup MeshingApplication
 * @class PointBoundary
 * @brief Custom Point container to be used to look in the boundary skin
 * @details The main difference with this point and the base one is that it contains the pointer to condition where the center of the points belongs
 * @author Vicente Mataix Ferrandiz
 */
class PointBoundary
    : public Point
{
public:
    ///@name Type Definitions
    ///@{

    typedef Point BaseType;

    /// Counted pointer of PointBoundary
    KRATOS_CLASS_POINTER_DEFINITION( PointBoundary );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    PointBoundary():
        BaseType(),
        mpOriginCond(nullptr)
    {}

    PointBoundary(const array_1d<double, 3>& Coords)
        :BaseType(Coords),
         mpOriginCond(nullptr)
    {}

    PointBoundary(Condition::Pointer pCond):
        mpOriginCond(pCond)
    {
        UpdatePoint();
    }

    PointBoundary(
        const array_1d<double, 3>& Coords,
        Condition::Pointer pCond
    ):
        BaseType(Coords),
        mpOriginCond(pCond)
    {}

    ///Copy constructor  (not really required)
    PointBoundary(const PointBoundary& rhs):
        BaseType(rhs),
        mpOriginCond(rhs.mpOriginCond)
    {
    }

    /// Destructor.
    ~PointBoundary() override= default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Returns the point
     * @return The point
     */
    BaseType GetPoint()
    {
        BaseType Point(this->Coordinates());
        return Point;
    }

    /**
     * @brief Set the point
     * @param Point The point
     */
    void SetPoint(const BaseType Point)
    {
        this->Coordinates() = Point.Coordinates();
    }

    /**
     * @brief Sets the condition associated to the point
     * @param pCond The pointer to the condition
     */
    void SetCondition(Condition::Pointer pCond)
    {
        mpOriginCond = pCond;
    }

    /**
     * @brief Returns the condition associated to the point
     * @return mpOriginCond The pointer to the condition associated to the point
     */
    Condition::Pointer GetCondition()
    {
        KRATOS_DEBUG_ERROR_IF(mpOriginCond.get() == nullptr) << "Condition no initialized in the PointBoundary class" << std::endl;
        return mpOriginCond;
    }

    /**
     * @brief This method checks everything is right
     */
    void Check()
    {
        KRATOS_TRY;

        auto aux_coord = Kratos::make_shared<array_1d<double, 3>>(this->Coordinates());
        KRATOS_ERROR_IF(!aux_coord) << "Coordinates no initialized in the PointBoundary class" << std::endl;
        KRATOS_ERROR_IF(mpOriginCond.get() == nullptr) << "Condition no initialized in the PointBoundary class" << std::endl;

        KRATOS_CATCH("");
    }

    /**
     * @brief This function updates the database, using as base for the coordinates the condition center
     */
    void UpdatePoint()
    {
        noalias(this->Coordinates()) = mpOriginCond->GetGeometry().Center().Coordinates();
    }

private:
    ///@name Member Variables
    ///@{
    Condition::Pointer mpOriginCond; /// Condition pointer
    ///@}

}; // Class PointBoundary

/**
 * @class NodalValuesInterpolationProcess
 * @ingroup MeshingApplication
 * @brief This utility has as objective to interpolate the values inside elements (and conditions?) in a model part, using as input the original model part and the new one
 * @details The process employs the projection.h from MeshingApplication, which works internally using a kd-tree. Additionally if it can't found the node inside the reference mesh it will try to extrapolate from the skin (if the option is activated)
 * @author Vicente Mataix Ferrandiz
 */
template<SizeType TDim>
class KRATOS_API(MESHING_APPLICATION) NodalValuesInterpolationProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    // General type definitions
    typedef ModelPart::NodesContainerType                    NodesArrayType;
    typedef ModelPart::ElementsContainerType              ElementsArrayType;
    typedef ModelPart::ConditionsContainerType          ConditionsArrayType;
    typedef Node                                                NodeType;
    typedef Geometry<NodeType>                                 GeometryType;
    typedef Point                                                 PointType;
    typedef PointType::CoordinatesArrayType            CoordinatesArrayType;

    // Type definitions for the tree
    typedef PointBoundary                                 PointBoundaryType;
    typedef PointBoundaryType::Pointer                     PointTypePointer;
    typedef std::vector<PointTypePointer>                       PointVector;
    typedef PointVector::iterator                             PointIterator;
    typedef std::vector<double>                              DistanceVector;
    typedef DistanceVector::iterator                       DistanceIterator;

    // KDtree definitions
    typedef Bucket< 3ul, PointBoundaryType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTreeType;

    /// Pointer definition of NodalValuesInterpolationProcess
    KRATOS_CLASS_POINTER_DEFINITION( NodalValuesInterpolationProcess );

    ///@}
    ///@name  Enum's
    ///@{

    /**
     * @brief This enums allows to differentiate the working framework
     */
    enum class FrameworkEulerLagrange {EULERIAN = 0, LAGRANGIAN = 1, ALE = 2};

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief The constructor of the search utility uses the following inputs:
     * @param rOriginMainModelPart The model part from where interpolate values
     * @param rDestinationMainModelPart The model part where we want to interpolate the values
     * @param ThisParameters The parameters containing all the information needed
     */
    NodalValuesInterpolationProcess(
        ModelPart& rOriginMainModelPart,
        ModelPart& rDestinationMainModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Destructor
    ~NodalValuesInterpolationProcess() override= default;;

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief We execute the search relative to the old and new model part
     */
    void Execute() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /************************************ GET INFO *************************************/
    /***********************************************************************************/

    std::string Info() const override
    {
        return "NodalValuesInterpolationProcess";
    }

    /************************************ PRINT INFO ***********************************/
    /***********************************************************************************/

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

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

    ModelPart& mrOriginMainModelPart;               /// The origin model part
    ModelPart& mrDestinationMainModelPart;          /// The destination model part
    Parameters mThisParameters;                     /// Here the configuration parameters are stored
    std::unordered_set<std::string> mListVariables; /// List of non-historical variables

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This converts the framework string to an enum
     * @param Str The string
     * @return FrameworkEulerLagrange: The equivalent enum
     */
    static inline FrameworkEulerLagrange ConvertFramework(const std::string& Str)
    {
        if(Str == "Lagrangian" || Str == "LAGRANGIAN")
            return FrameworkEulerLagrange::LAGRANGIAN;
        else if(Str == "Eulerian" || Str == "EULERIAN")
            return FrameworkEulerLagrange::EULERIAN;
        else if(Str == "ALE")
            return FrameworkEulerLagrange::ALE;
        else
            return FrameworkEulerLagrange::EULERIAN;
    }

    /**
     * @brief It calculates the data (DataContainer) interpolated to the node
     * @param pNode The node pointer
     * @param pEntity The element pointer
     * @param rShapeFunctions The shape functions
     */
    template<class TEntity>
    void CalculateData(
        NodeType::Pointer pNode,
        const typename TEntity::Pointer& pEntity,
        const Vector& rShapeFunctions
        )
    {
        // The nodal data (non-historical) of each node of the original mesh
        GeometryType& r_geometry = pEntity->GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        // Now we interpolate the values of each node
        double aux_coeff = 0.0;
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            aux_coeff += rShapeFunctions[i];
        }
        for (auto& r_variable_name : mListVariables) {
            if (KratosComponents<Variable<double>>::Has(r_variable_name)) {
                const Variable<double>& r_variable = KratosComponents<Variable<double>>::Get(r_variable_name);
                if (std::abs(aux_coeff) > std::numeric_limits<double>::epsilon()) {
                    aux_coeff = 1.0/aux_coeff;
                    double aux_value = 0.0;
                    for (IndexType i = 0; i < number_of_nodes; ++i) {
                        if (r_geometry[i].Has(r_variable)) {
                            aux_value += rShapeFunctions[i] * r_geometry[i].GetValue(r_variable);
                        }
                    }
                    pNode->SetValue(r_variable, aux_coeff * aux_value);
                }
            } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(r_variable_name)) {
                const Variable<array_1d<double, 3>>& r_variable = KratosComponents<Variable<array_1d<double, 3>>>::Get(r_variable_name);
                if (std::abs(aux_coeff) > std::numeric_limits<double>::epsilon()) {
                    aux_coeff = 1.0/aux_coeff;
                    array_1d<double, 3> aux_value = ZeroVector(3);
                    for (IndexType i = 0; i < number_of_nodes; ++i) {
                        if (r_geometry[i].Has(r_variable)) {
                            aux_value += rShapeFunctions[i] * r_geometry[i].GetValue(r_variable);
                        }
                    }
                    pNode->SetValue(r_variable, aux_coeff * aux_value);
                }
            } else if (KratosComponents<Variable<Vector>>::Has(r_variable_name)) {
                const Variable<Vector>& r_variable = KratosComponents<Variable<Vector>>::Get(r_variable_name);
                if (std::abs(aux_coeff) > std::numeric_limits<double>::epsilon()) {
                    aux_coeff = 1.0/aux_coeff;
                    Vector aux_value = ZeroVector(r_geometry[0].GetValue(r_variable).size());
                    for (IndexType i = 0; i < number_of_nodes; ++i) {
                        if (r_geometry[i].Has(r_variable)) {
                            aux_value += rShapeFunctions[i] * r_geometry[i].GetValue(r_variable);
                        }
                    }
                    pNode->SetValue(r_variable, aux_coeff * aux_value);
                }
            } else if (KratosComponents<Variable<Matrix>>::Has(r_variable_name)) {
                const Variable<Matrix>& r_variable = KratosComponents<Variable<Matrix>>::Get(r_variable_name);
                if (std::abs(aux_coeff) > std::numeric_limits<double>::epsilon()) {
                    aux_coeff = 1.0/aux_coeff;
                    Matrix aux_value = ZeroMatrix(r_geometry[0].GetValue(r_variable).size1(), r_geometry[0].GetValue(r_variable).size2());
                    for (IndexType i = 0; i < number_of_nodes; ++i) {
                        if (r_geometry[i].Has(r_variable)) {
                            aux_value += rShapeFunctions[i] * r_geometry[i].GetValue(r_variable);
                        }
                    }
                    pNode->SetValue(r_variable, aux_coeff * aux_value);
                }
            }
        }
    }

    /**
     * @brief It calculates the Step data interpolated to the node
     * @param pNode The node pointer
     * @param pEntity The element pointer
     * @param rShapeFunctions The shape functions
     * @param Step The current time step
     */
    template<class TEntity>
    void CalculateStepData(
        NodeType::Pointer pNode,
        const typename TEntity::Pointer& pEntity,
        const Vector& rShapeFunctions,
        const IndexType Step
        )
    {
        // The nodal data (historical)
        double* step_data = pNode->SolutionStepData().Data(Step);
        for (int j = 0; j < mThisParameters["step_data_size"].GetInt(); ++j)
            step_data[j] = 0;

        // The nodal data (historical) of each node of the original mesh
        GeometryType& geom = pEntity->GetGeometry();
        const SizeType number_of_nodes = geom.size();
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const double* nodal_data = geom[i].SolutionStepData().Data(Step);
            // Now we interpolate the values of each node
            for (int j = 0; j < mThisParameters["step_data_size"].GetInt(); ++j) {
                step_data[j] += rShapeFunctions[i] * nodal_data[j];
            }
        }
    }

    /**
     * @brief This methoid a boundary model part in the reference and target model part
     * @param rAuxiliarNameModelPart The name of the model part to be created
     */
    void GenerateBoundary(const std::string& rAuxiliarNameModelPart);

    /**
     * @brief This methoid a boundary model part in the reference and target model part
     * @param rModelPart The model part to compute
     * @param rAuxiliarNameModelPart The name of the model part to be created
     */
    void GenerateBoundaryFromElements(
        ModelPart& rModelPart,
        const std::string& rAuxiliarNameModelPart
        );

    /**
     * @brief This methoid a boundary model part in the reference and target model part
     * @param rAuxiliarNameModelPart The name of the model part to be created
     * @param rToExtrapolateNodes The list of nodes to extrapolate
     */
    void ExtrapolateValues(
        const std::string& rAuxiliarNameModelPart,
        std::vector<NodeType::Pointer>& rToExtrapolateNodes
        );

    /**
     * @brief This method computes the normal in a skin model part (saved in non historical variables)
     * @param rModelPart The input model part with the skin
     */
    void ComputeNormalSkin(ModelPart& rModelPart);

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class NodalValuesInterpolationProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/****************************** INPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

// template<class TPointType, class TPointerType>
// inline std::istream& operator >> (std::istream& rIStream,
//                                   NodalValuesInterpolationProcess& rThis);

/***************************** OUTPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

// template<class TPointType, class TPointerType>
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const NodalValuesInterpolationProcess& rThis)
// {
//     return rOStream;
// }

///@}

}  // namespace Kratos.
