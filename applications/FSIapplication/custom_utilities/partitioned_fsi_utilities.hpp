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

#if !defined( KRATOS_PARTITIONED_FSI_UTILITIES )
#define  KRATOS_PARTITIONED_FSI_UTILITIES


/* System includes */
#include <set>
#include <typeinfo>

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/mesh_moving_variables.h"
#include "includes/fsi_variables.h"
#include "containers/array_1d.h"
#include "includes/model_part.h"
#include "includes/communicator.h"
#include "includes/ublas_interface.h"
#include "utilities/math_utils.h"
#include "utilities/openmp_utils.h"
#include "utilities/variable_utils.h"


namespace Kratos
{

/**@name Kratos Globals */
/*@{ */

/*@} */
/**@name Type Definitions */
/*@{ */

/*@} */
/**@name  Enum's */
/*@{ */

/*@} */
/**@name  Functions */
/*@{ */

/*@} */
/**@name Kratos Classes */
/*@{ */

/** Short class definition.
Detail class definition.
*/

template<class TSpace, class TValueType, unsigned int TDim>
class PartitionedFSIUtilities
{

public:

    /** Type Definitions
    */

    /*@{ */
    typedef typename TSpace::VectorType                     VectorType;
    typedef typename TSpace::MatrixType                     MatrixType;

    typedef typename TSpace::VectorPointerType              VectorPointerType;
    typedef typename TSpace::MatrixPointerType              MatrixPointerType;

    //~ /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( PartitionedFSIUtilities );
    /*@} */

    /** Constructor.
     */

    /**
    * Empty constructor
    */
    PartitionedFSIUtilities(){}

    /** Copy constructor.
    */

    /*@{ */
    PartitionedFSIUtilities(const PartitionedFSIUtilities& Other);
    /*@{ */

    /** Destructor.
     */

    virtual ~PartitionedFSIUtilities(){}

    /*@} */
    /**@name Public Operators*/
    /*@{ */

    /**
     * @brief Creates an element based skin
     * For a modelpart defining the skin using conditions, this method
     * copies such skin to the elements of an auxiliar model part. Note
     * that the same geometry is used so the nodes of the auxiliar geometry
     * are actually the ones in the origin modelpart.
     * @param rOriginInterfaceModelPart
     * @param rDestinationInterfaceModelPart
     */
    void CopySkinToElements(
        const ModelPart& rOriginInterfaceModelPart,
        ModelPart& rDestinationInterfaceModelPart)
    {
        // Add the origin interface nodes to the destination interface model part
        rDestinationInterfaceModelPart.AddNodes(
            rOriginInterfaceModelPart.NodesBegin(),
            rOriginInterfaceModelPart.NodesEnd());

        // Create new elements emulating the condition based interface
        ModelPart::ElementsContainerType new_elems_vect;
        for (int i_cond = 0; i_cond < rOriginInterfaceModelPart.NumberOfConditions(); ++i_cond) {
            const auto &it_cond = rOriginInterfaceModelPart.ConditionsBegin() + i_cond;
            auto p_elem = Kratos::make_shared<Element>(it_cond->Id(), it_cond->pGetGeometry());
            rDestinationInterfaceModelPart.AddElement(p_elem);
        }
    }

    /**
     * This function computes the interface residual size as the
     * number of interface nodes times the problem domain size.
     * @return the model part residual size
     */
    int GetInterfaceResidualSize(ModelPart& rInterfaceModelPart)
    {
        // Check the block size. 1 for double coupling variables an TDim for array type coupling variables
        double A; // Fake double to check the type from
        unsigned int block_size = typeid(TValueType).hash_code() == typeid(A).hash_code() ? 1 : TDim;
        int local_number_of_nodes = (rInterfaceModelPart.GetCommunicator().LocalMesh().NumberOfNodes()) * block_size;
        rInterfaceModelPart.GetCommunicator().SumAll(local_number_of_nodes);
        return local_number_of_nodes;
    }

    /**
     * This function returns the interface length in 2D or the interface area in 3D.
     * @param rInterfaceModelPart: interface modelpart in where the are is computed
     * @return the given modelpart interface length
     */
    double GetInterfaceArea(ModelPart& rInterfaceModelPart)
    {
        double interface_area = 0.0;

        auto& rLocalMesh = rInterfaceModelPart.GetCommunicator().LocalMesh();
        ModelPart::ConditionIterator local_mesh_conditions_begin = rLocalMesh.ConditionsBegin();
        #pragma omp parallel for firstprivate(local_mesh_conditions_begin) reduction(+:interface_area)
        for(int k=0; k < static_cast<int>(rLocalMesh.NumberOfConditions()); ++k) {
            const ModelPart::ConditionIterator it_cond = local_mesh_conditions_begin+k;
            const Condition::GeometryType& rGeom = it_cond->GetGeometry();
            interface_area += rGeom.Length();
        }

        rInterfaceModelPart.GetCommunicator().SumAll(interface_area);

        return interface_area;
    }

    /**
     * This function resizes and sets to zero an interface vector (length equal to the
     * residual size).
     * @param rInterfaceModelPart: interface modelpart in where the residual is computed
     * @return pointer to the vector that has been created
     */
    virtual VectorPointerType SetUpInterfaceVector(ModelPart& rInterfaceModelPart)
    {
        VectorPointerType p_int_vector = TSpace::CreateEmptyVectorPointer();
        const unsigned int residual_size = this->GetInterfaceResidualSize(rInterfaceModelPart);
        if (TSpace::Size(*p_int_vector) != residual_size){
            TSpace::Resize(p_int_vector, residual_size);
        }
        TSpace::SetToZero(*p_int_vector);
        return p_int_vector;
    }

    /**
     * @brief Compute array variable interface residual vector
     * This function computes (and stores in a vector) the residual of a vector variable over the fluid interface.
     * The residual is defined as the OriginalVariable value minus the ModifiedVariable value.
     * The nodal values of the residual are stored in the FSI_INTERFACE_RESIDUAL variable.
     * Besides, the norm of the residual vector is stored in the ProcessInfo.
     * @param rInterfaceModelPart interface modelpart in where the residual is computed
     * @param rOriginalVariable origin variable to compute the residual
     * @param rModifiedVariable end variable to compute the residual
     * @param rInterfaceResidual reference to the residual vector
     * @param ResidualType residual computation type (nodal or consistent)
     * @param rArrayResidualVariable variable to save the residual values
     * @param rResidualNormVariable variable to save the residual norm
     */
    virtual void ComputeInterfaceResidualVector(
        ModelPart &rInterfaceModelPart,
        const Variable<TValueType> &rOriginalVariable,
        const Variable<TValueType> &rModifiedVariable,
        const Variable<TValueType> &rResidualVariable,
        VectorType &rInterfaceResidual,
        const std::string ResidualType = "nodal",
        const Variable<double> &rResidualNormVariable = FSI_INTERFACE_RESIDUAL_NORM)
    {
        TSpace::SetToZero(rInterfaceResidual);

        // Call compute the provided variables residual
        if (ResidualType == "nodal") {
            // Compute node-by-node residual
            this->ComputeNodeByNodeResidual(rInterfaceModelPart, rOriginalVariable, rModifiedVariable, rResidualVariable);
        } else if (ResidualType == "consistent") {
            // Compute consitent residual
            this->ComputeConsistentResidual(rInterfaceModelPart, rOriginalVariable, rModifiedVariable, rResidualVariable);
        } else {
            KRATOS_ERROR << "Provided interface residual type " << ResidualType << " is not available. Available options are \"nodal\" and \"consistent\"" << std::endl;
        }

        // Assemble the final consistent residual values
        auto& rLocalMesh = rInterfaceModelPart.GetCommunicator().LocalMesh();
        ModelPart::NodeIterator local_mesh_nodes_begin = rLocalMesh.NodesBegin();
        #pragma omp parallel for firstprivate(local_mesh_nodes_begin)
        for(int k = 0; k < static_cast<int>(rLocalMesh.NumberOfNodes()); ++k) {
            const ModelPart::NodeIterator it_node = local_mesh_nodes_begin+k;
            const auto &r_res_value = it_node->FastGetSolutionStepValue(rResidualVariable);
            this->AuxSetLocalValue(rInterfaceResidual, r_res_value, k);
        }

        // Store the L2 norm of the error in the fluid process info
        rInterfaceModelPart.GetProcessInfo().GetValue(rResidualNormVariable) = TSpace::TwoNorm(rInterfaceResidual);
    }

    /**
     * @brief Computes the interface residual norm
     * This method computes the interface residual norm. To do that it firstly
     * computes the residual vector as is done in ComputeInterfaceResidualVector.
     * @param rInterfaceModelPart interface modelpart in where the residual is computed
     * @param rOriginalVariable origin variable to compute the residual
     * @param rModifiedVariable end variable to compute the residual
     * @param rInterfaceResidual reference to the residual vector
     * @param ResidualType residual computation type (nodal or consistent)
     * @return double interface residual norm
     */
    double ComputeInterfaceResidualNorm(
        ModelPart &rInterfaceModelPart,
        const Variable<TValueType> &rOriginalVariable,
        const Variable<TValueType> &rModifiedVariable,
        const Variable<TValueType> &rResidualVariable,
        const std::string ResidualType = "nodal")
    {
        // Set an auxiliar interface vector
        VectorPointerType p_interface_residual = this->SetUpInterfaceVector(rInterfaceModelPart);

        // Call compute the provided variables residual
        if (ResidualType == "nodal") {
            // Compute node-by-node residual
            this->ComputeNodeByNodeResidual(rInterfaceModelPart, rOriginalVariable, rModifiedVariable, rResidualVariable);
        } else if (ResidualType == "consistent") {
            // Compute consitent residual
            this->ComputeConsistentResidual(rInterfaceModelPart, rOriginalVariable, rModifiedVariable, rResidualVariable);
        } else {
            KRATOS_ERROR << "Provided interface residual type " << ResidualType << " is not available. Available options are \"nodal\" and \"consistent\"" << std::endl;
        }

        // Assemble the final consistent residual values
        auto &rLocalMesh = rInterfaceModelPart.GetCommunicator().LocalMesh();
        ModelPart::NodeIterator local_mesh_nodes_begin = rLocalMesh.NodesBegin();
        #pragma omp parallel for firstprivate(local_mesh_nodes_begin)
        for(int k=0; k<static_cast<int>(rLocalMesh.NumberOfNodes()); ++k) {
            const ModelPart::NodeIterator it_node = local_mesh_nodes_begin + k;
            const auto &r_res_value = it_node->FastGetSolutionStepValue(rResidualVariable);
            this->AuxSetLocalValue(*p_interface_residual, r_res_value, k);
        }

        // Return the L2 norm of the error in the fluid process info
        return TSpace::TwoNorm(*p_interface_residual);
    }

    /**
     * @brief Auxiliar call to SetLocalValue for double type
     * This method serves as an auxiliar call to the SetLocalValue function for double variables
     * @param rInterfaceResidual Interface residual vector in where the value is set
     * @param rResidualValue Double residual value
     * @param AuxPosition Position in rInterfaceResidual where the value is to be set
     */
    virtual void AuxSetLocalValue(
        VectorType &rInterfaceResidual,
        const double &rResidualValue,
        const int AuxPosition)
    {
        this->SetLocalValue(rInterfaceResidual, AuxPosition, rResidualValue);
    }

    /**
     * @brief Auxiliar call to SetLocalValue for array type
     * This method serves as an auxiliar call to the SetLocalValue function for array variables
     * @param rInterfaceResidual Interface residual vector in where the value is set
     * @param rResidualValue Double residual value
     * @param AuxPosition Position in rInterfaceResidual where the value is to be set. Note that
     * the final position is computed using the TDim template value, since each entry in the
     * residual vector contains a component of the residual variable value.
     */
    virtual void AuxSetLocalValue(
        VectorType &rInterfaceResidual,
        const array_1d<double,3> &rResidualValue,
        const int AuxPosition)
    {
        const unsigned int base_i = AuxPosition * TDim;
        for (unsigned int jj = 0; jj < TDim; ++jj) {
            this->SetLocalValue(rInterfaceResidual, base_i + jj, rResidualValue[jj]);
        }
    }

    /**
     * Sets the values in the corrected guess vector inside the selected nodal array variable.
     * @param rInterfaceModelPart: interface modelpart in where the residual is computed
     * @param rSolutionVariable: variable in where the corrected solution is to be stored
     * @param rCorrectedGuess: vector containing the interface corrected values
     */
    virtual void UpdateInterfaceValues(
        ModelPart &rInterfaceModelPart,
        const Variable<TValueType> &rSolutionVariable,
        const VectorType &rCorrectedGuess)
    {
        auto& rLocalMesh = rInterfaceModelPart.GetCommunicator().LocalMesh();
        ModelPart::NodeIterator local_mesh_nodes_begin = rLocalMesh.NodesBegin();
        #pragma omp parallel for firstprivate(local_mesh_nodes_begin)
        for(int k = 0; k < static_cast<int>(rLocalMesh.NumberOfNodes()); ++k){
            const ModelPart::NodeIterator it_node = local_mesh_nodes_begin + k;
            TValueType &r_updated_value = it_node->FastGetSolutionStepValue(rSolutionVariable);
            this->UpdateInterfaceLocalValue(rCorrectedGuess, r_updated_value, k);
        }

        rInterfaceModelPart.GetCommunicator().SynchronizeVariable(rSolutionVariable);
    }

    /**
     * @brief Corrected guess local double value update
     * This auxiliar method helps to update the nodal database using the values from
     * a corrected guess vector values.
     * @param rCorrectedGuess Corrected coupling guess vector
     * @param rValueToUpdate Reference in where the updated value is to be stored
     * @param AuxPosition Position in rCorrectedGuess from where the value is retrieved
     */
    void UpdateInterfaceLocalValue(
        const VectorType &rCorrectedGuess,
        double &rValueToUpdate,
        const int AuxPosition)
    {
        rValueToUpdate = this->GetLocalValue(rCorrectedGuess, AuxPosition);
    }

    /**
     * @brief Corrected guess local array value update
     * This auxiliar method helps to update the nodal database using the values from
     * a corrected guess vector values. Note that it performs the array components loop.
     * @param rCorrectedGuess Corrected coupling guess vector
     * @param rValueToUpdate Reference in where the updated value is to be stored
     * @param AuxPosition Position in rCorrectedGuess from where the value is retrieved
     */
    void UpdateInterfaceLocalValue(
        const VectorType &rCorrectedGuess,
        array_1d<double, 3> &rValueToUpdate,
        const int AuxPosition)
    {
        const int base_i = AuxPosition * TDim;
        for (unsigned int jj = 0; jj < TDim; ++jj){
            rValueToUpdate[jj] = this->GetLocalValue(rCorrectedGuess, base_i + jj);
        }
    }

    /**
     * Computes and prints the fluid interface residual norms for debugging purposes
     * @param rInterfaceModelPart: interface modelpart in where the residual is computed
     */
    virtual void ComputeAndPrintFluidInterfaceNorms(ModelPart& rInterfaceModelPart)
    {
        double p_norm = 0.0;
        double vx_norm = 0.0;
        double vy_norm = 0.0;
        double vz_norm = 0.0;
        double rx_norm = 0.0;
        double ry_norm = 0.0;
        double rz_norm = 0.0;
        double ux_mesh_norm = 0.0;
        double uy_mesh_norm = 0.0;
        double uz_mesh_norm = 0.0;

        auto& rLocalMesh = rInterfaceModelPart.GetCommunicator().LocalMesh();
        ModelPart::NodeIterator local_mesh_nodes_begin = rLocalMesh.NodesBegin();
        #pragma omp parallel for firstprivate(local_mesh_nodes_begin) reduction(+ : p_norm, vx_norm, vy_norm, vz_norm, rx_norm, ry_norm, rz_norm, ux_mesh_norm, uy_mesh_norm, uz_mesh_norm)
        for(int k=0; k<static_cast<int>(rLocalMesh.NumberOfNodes()); ++k)
        {
            const ModelPart::NodeIterator it_node = local_mesh_nodes_begin+k;

            p_norm += std::pow(it_node->FastGetSolutionStepValue(PRESSURE), 2);
            vx_norm += std::pow(it_node->FastGetSolutionStepValue(VELOCITY_X), 2);
            vy_norm += std::pow(it_node->FastGetSolutionStepValue(VELOCITY_Y), 2);
            vz_norm += std::pow(it_node->FastGetSolutionStepValue(VELOCITY_Z), 2);
            rx_norm += std::pow(it_node->FastGetSolutionStepValue(REACTION_X), 2);
            ry_norm += std::pow(it_node->FastGetSolutionStepValue(REACTION_Y), 2);
            rz_norm += std::pow(it_node->FastGetSolutionStepValue(REACTION_Z), 2);
            ux_mesh_norm += std::pow(it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT_X), 2);
            uy_mesh_norm += std::pow(it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT_Y), 2);
            uz_mesh_norm += std::pow(it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT_Z), 2);
        }

        rInterfaceModelPart.GetCommunicator().SumAll(p_norm);
        rInterfaceModelPart.GetCommunicator().SumAll(vx_norm);
        rInterfaceModelPart.GetCommunicator().SumAll(vy_norm);
        rInterfaceModelPart.GetCommunicator().SumAll(vz_norm);
        rInterfaceModelPart.GetCommunicator().SumAll(rx_norm);
        rInterfaceModelPart.GetCommunicator().SumAll(ry_norm);
        rInterfaceModelPart.GetCommunicator().SumAll(rz_norm);
        rInterfaceModelPart.GetCommunicator().SumAll(ux_mesh_norm);
        rInterfaceModelPart.GetCommunicator().SumAll(uy_mesh_norm);
        rInterfaceModelPart.GetCommunicator().SumAll(uz_mesh_norm);

        if (rInterfaceModelPart.GetCommunicator().MyPID() == 0)
        {
            std::cout << " " << std::endl;
            std::cout << "|p_norm| = " << std::sqrt(p_norm) << std::endl;
            std::cout << "|vx_norm| = " << std::sqrt(vx_norm) << std::endl;
            std::cout << "|vy_norm| = " << std::sqrt(vy_norm) << std::endl;
            std::cout << "|vz_norm| = " << std::sqrt(vz_norm) << std::endl;
            std::cout << "|rx_norm| = " << std::sqrt(rx_norm) << std::endl;
            std::cout << "|ry_norm| = " << std::sqrt(ry_norm) << std::endl;
            std::cout << "|rz_norm| = " << std::sqrt(rz_norm) << std::endl;
            std::cout << "|ux_mesh_norm| = " << std::sqrt(ux_mesh_norm) << std::endl;
            std::cout << "|uy_mesh_norm| = " << std::sqrt(uy_mesh_norm) << std::endl;
            std::cout << "|uz_mesh_norm| = " << std::sqrt(uz_mesh_norm) << std::endl;
            std::cout << " " << std::endl;
        }
    }

    /**
     * Computes and prints the structure interface residual norms for debugging purposes
     * @param rInterfaceModelPart: interface modelpart in where the residual is computed
     */
    virtual void ComputeAndPrintStructureInterfaceNorms(ModelPart& rInterfaceModelPart)
    {
        double ux_norm = 0.0;
        double uy_norm = 0.0;
        double uz_norm = 0.0;

        auto& rLocalMesh = rInterfaceModelPart.GetCommunicator().LocalMesh();
        ModelPart::NodeIterator local_mesh_nodes_begin = rLocalMesh.NodesBegin();
        #pragma omp parallel for firstprivate(local_mesh_nodes_begin) reduction(+ : ux_norm, uy_norm, uz_norm)
        for(int k=0; k<static_cast<int>(rLocalMesh.NumberOfNodes()); ++k)
        {
            const ModelPart::NodeIterator it_node = local_mesh_nodes_begin+k;
            const array_1d<double, 3>& disp = it_node->FastGetSolutionStepValue(DISPLACEMENT);

            ux_norm += std::pow(disp[0], 2);
            uy_norm += std::pow(disp[1], 2);
            uz_norm += std::pow(disp[2], 2);
        }

        rInterfaceModelPart.GetCommunicator().SumAll(ux_norm);
        rInterfaceModelPart.GetCommunicator().SumAll(uy_norm);
        rInterfaceModelPart.GetCommunicator().SumAll(uz_norm);

        if (rInterfaceModelPart.GetCommunicator().MyPID() == 0)
        {
            std::cout << " " << std::endl;
            std::cout << "|ux_norm| = " << std::sqrt(ux_norm) << std::endl;
            std::cout << "|uy_norm| = " << std::sqrt(uy_norm) << std::endl;
            std::cout << "|uz_norm| = " << std::sqrt(uz_norm) << std::endl;
            std::cout << " " << std::endl;
        }
    }

    /**
     * Checks if X = X0 + deltaX
     * @param rModelPart: reference to the model part in where the check has to be performed
     */
    virtual void CheckCurrentCoordinatesFluid(ModelPart& rModelPart, const double tolerance)
    {
        auto& rLocalMesh = rModelPart.GetCommunicator().LocalMesh();
        ModelPart::NodeIterator local_mesh_nodes_begin = rLocalMesh.NodesBegin();
        #pragma omp parallel for firstprivate(local_mesh_nodes_begin)
        for(int k=0; k<static_cast<int>(rLocalMesh.NumberOfNodes()); ++k)
        {
            const ModelPart::NodeIterator it_node = local_mesh_nodes_begin+k;
            const array_1d<double, 3>& disp = it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT);

            if (std::fabs(it_node->X() - (it_node->X0() + disp[0])) > tolerance)
            {
                KRATOS_ERROR << "Node " << it_node->Id() << " X != X0 + deltaX";
            }
            if (std::fabs(it_node->Y() - (it_node->Y0() + disp[1])) > tolerance)
            {
                KRATOS_ERROR << "Node " << it_node->Id() << " Y != Y0 + deltaY";
            }
            if (std::fabs(it_node->Z() - (it_node->Z0() + disp[2])) > tolerance)
            {
                KRATOS_ERROR << "Node " << it_node->Id() << " Z != Z0 + deltaZ";
            }
        }
    }

    /**
     * Checks if X = X0 + deltaX
     * @param rModelPart: reference to the model part in where the check has to be performed
     */
    virtual void CheckCurrentCoordinatesStructure(ModelPart& rModelPart, const double tolerance)
    {
        auto& rLocalMesh = rModelPart.GetCommunicator().LocalMesh();
        ModelPart::NodeIterator local_mesh_nodes_begin = rLocalMesh.NodesBegin();
        #pragma omp parallel for firstprivate(local_mesh_nodes_begin)
        for(int k=0; k<static_cast<int>(rLocalMesh.NumberOfNodes()); ++k)
        {
            const ModelPart::NodeIterator it_node = local_mesh_nodes_begin+k;
            const array_1d<double, 3>& disp = it_node->FastGetSolutionStepValue(DISPLACEMENT);

            if (std::fabs(it_node->X() - (it_node->X0() + disp[0])) > tolerance)
            {
                KRATOS_ERROR << "Node " << it_node->Id() << " X != X0 + deltaX";
            }
            if (std::fabs(it_node->Y() - (it_node->Y0() + disp[1])) > tolerance)
            {
                KRATOS_ERROR << "Node " << it_node->Id() << " Y != Y0 + deltaY";
            }
            if (std::fabs(it_node->Z() - (it_node->Z0() + disp[2])) > tolerance)
            {
                KRATOS_ERROR << "Node " << it_node->Id() << " Z != Z0 + deltaZ";
            }
        }
    }

    /**
     * This function computes the nodal error of a vector magnitude in a consistent manner.
     * The error is defined as the integral over the interface of a tests function times
     * the difference between rOriginalVariable and rModifiedVariable.
     * @param rInterfaceModelPart: interface modelpart in where the residual is computed
     * @param rOriginalVariable: variable with the reference value
     * @param rModifiedVariable: variable with the computed vvalue
     * @param rErrorStorageVariable: variable to store the error nodal value
     */
    void ComputeConsistentResidual(
        ModelPart& rInterfaceModelPart,
        const Variable<TValueType>& rOriginalVariable,
        const Variable<TValueType>& rModifiedVariable,
        const Variable<TValueType>& rErrorStorageVariable)
    {
        // Initialize the interface residual variable
        VariableUtils().SetVariable<TValueType>(rErrorStorageVariable, rErrorStorageVariable.Zero(), rInterfaceModelPart.Nodes());

        #pragma omp parallel for
        for(int i_cond = 0; i_cond < static_cast<int>(rInterfaceModelPart.NumberOfConditions()); ++i_cond) {

            auto it_cond = rInterfaceModelPart.ConditionsBegin() + i_cond;

            auto& rGeom = it_cond->GetGeometry();
            // const unsigned int BlockSize = typeid(TValueType).name() == "double" ? 1 : TDim;
            const unsigned int n_nodes = rGeom.PointsNumber();

            // Auxiliar array to save the computed consisted residual
            // Vector cons_res_vect = ZeroVector(BlockSize*n_nodes);

            // Initialize auxiliar array to save the condition nodes consistent residual
            std::vector<TValueType> cons_res_vect(n_nodes);
            for (unsigned int i = 0; i < n_nodes; ++i) {
                cons_res_vect[i] = rErrorStorageVariable.Zero();
            }

            // Compute the consistent residual
            const auto &r_int_pts = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
            const auto N_container = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
            Vector jac_gauss;
            rGeom.DeterminantOfJacobian(jac_gauss, GeometryData::GI_GAUSS_2);

            for (unsigned int i_gauss = 0; i_gauss < r_int_pts.size(); ++i_gauss) {
                // Compute condition Gauss pt. data
                const Vector N_gauss = row(N_container, i_gauss);
                const double w_gauss = jac_gauss[i_gauss] * r_int_pts[i_gauss].Weight();

                // Add the current Gauss pt. residual contribution
                for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
                    for (unsigned int j_node = 0; j_node < n_nodes; ++j_node) {
                        const double aux_val = w_gauss * N_gauss[i_node] * N_gauss[j_node];
                        const TValueType value = rGeom[j_node].FastGetSolutionStepValue(rOriginalVariable);
                        const TValueType value_projected = rGeom[j_node].FastGetSolutionStepValue(rModifiedVariable);
                        cons_res_vect[i_node] += aux_val * (value - value_projected);
                    }
                }
            }

            // Save the computed consistent residual values
            for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
                rGeom[i_node].SetLock(); // So it is safe to write in the condition node in OpenMP
                rGeom[i_node].FastGetSolutionStepValue(rErrorStorageVariable) += cons_res_vect[i_node];
                rGeom[i_node].UnSetLock(); // Free the condition node for other threads
            }

            // Synchronize the computed error values
            rInterfaceModelPart.GetCommunicator().AssembleCurrentData(rErrorStorageVariable);
        }
    }

    /*@} */
protected:
    /**@name Protected static Member Variables */
    /*@{ */

    /*@} */
    /**@name Protected member Variables */
    /*@{ */

    /*@} */
    /**@name Protected Operators*/
    /*@{ */


    /*@} */
    /**@name Protected Operations*/
    /*@{ */

    /**
     * This function computes the nodal error of a vector magnitude. The error is defined
     * as OriginalVariable minus ModifiedVariable.
     * @param rInterfaceModelPart: interface modelpart in where the residual is computed
     * @param rOriginalVariable: variable with the reference value
     * @param rModifiedVariable: variable with the computed vvalue
     * @param rErrorStorageVariable: variable to store the error nodal value
     */
    void ComputeNodeByNodeResidual(
        ModelPart& rInterfaceModelPart,
        const Variable<TValueType>& rOriginalVariable,
        const Variable<TValueType>& rModifiedVariable,
        const Variable<TValueType>& rErrorStorageVariable)
    {
        auto& rLocalMesh = rInterfaceModelPart.GetCommunicator().LocalMesh();
        ModelPart::NodeIterator local_mesh_nodes_begin = rLocalMesh.NodesBegin();
        #pragma omp parallel for firstprivate(local_mesh_nodes_begin)
        for(int k = 0; k < static_cast<int>(rLocalMesh.NumberOfNodes()); ++k) {
            ModelPart::NodeIterator it_node = local_mesh_nodes_begin+k;
            auto &r_error_storage = it_node->FastGetSolutionStepValue(rErrorStorageVariable);
            const auto &value_origin = it_node->FastGetSolutionStepValue(rOriginalVariable);
            const auto &value_modified = it_node->FastGetSolutionStepValue(rModifiedVariable);
            r_error_storage = value_origin - value_modified;
        }
    }

    virtual void SetLocalValue(VectorType& rVector, int LocalRow, double Value) const
    {
        TSpace::SetValue(rVector,LocalRow,Value);
    }

    virtual double GetLocalValue(const VectorType& rVector, int LocalRow) const
    {
        return TSpace::GetValue(rVector,LocalRow);
    }

    /**@} */
    /**@name Protected  Access */
    /*@{ */


    /*@} */
    /**@name Protected Inquiry */
    /*@{ */


    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */

private:

    /*@} */
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */


    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */


    /*@} */
    /**@name Private  Access */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */

    /*@} */

}; /* Class PartitionedFSIUtilities */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_PARTITIONED_FSI_UTILITIES  defined */
