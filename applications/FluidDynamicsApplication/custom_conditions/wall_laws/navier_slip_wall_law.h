//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes


// External includes


// Project includes
#include "includes/cfd_variables.h"
#include "includes/condition.h"
#include "includes/define.h"
#include "includes/process_info.h"

// Application includes
#include "fluid_dynamics_application_variables.h"


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
 * @brief Navier-slip law LHS and RHS contribution implementation
 * This class implements the LHS and RHS Gauss point contributions of the Navier-slip wall model
 * This class should be used in combination with an incompressible Navier-Stokes (or Stokes) condition
 * implementing the remaining terms (see @NavierStokesWallCondition).
 * More information about can be found in https://onlinelibrary.wiley.com/doi/abs/10.1002/fld.663
 * @tparam TDim Number of dimensions
 * @tparam TNumNodes Number of nodes
 */
template<std::size_t TDim, std::size_t TNumNodes>
class NavierSlipWallLaw
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NavierSlipWallLaw
    KRATOS_CLASS_POINTER_DEFINITION(NavierSlipWallLaw);

    static constexpr std::size_t BlockSize = TDim+1;

    static constexpr std::size_t LocalSize = TNumNodes*BlockSize;

    using SizeType = Condition::SizeType;

    using IndexType = Condition::IndexType;

    using VectorType = array_1d<double,LocalSize>;

    using MatrixType = BoundedMatrix<double,LocalSize,LocalSize>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NavierSlipWallLaw() = delete;

    /// Copy constructor.
    NavierSlipWallLaw(NavierSlipWallLaw const& rOther) = delete;

    /// Destructor.
    ~NavierSlipWallLaw() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Add the LHS and RHS Navier-slip contributions
     * This method adds the Navier-slip LHS and RHS Gauss point contributions to the provided matrices
     * @tparam TConditionDataContainer Condition data container
     * @param rLeftHandSideMatrix Reference to the output LHS matrix
     * @param rRightHandSideVector Reference to the output RHS matrix
     * @param pCondition Pointer to the current condition
     * @param rCurrentProcessInfo Reference to current ProcessInfo container
     * @param rConditionData Reference to the condition data container
     */
    template<class TConditionDataContainer>
    static void AddLocalSystemGaussPointContribution(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const Condition* pCondition,
        const ProcessInfo& rCurrentProcessInfo,
        const TConditionDataContainer& rConditionData)
    {
        // Create and fill the wall law auxiliar data container
        WallLawDataContainer wall_law_data;
        wall_law_data.Initialize(*pCondition, rConditionData);

        // Set the tangential projection matrix
        BoundedMatrix<double,TDim,TDim> tang_proj_mat;
        SetTangentialProjectionMatrix(rConditionData.Normal, tang_proj_mat);

        // Calculate the Navier-slip required magnitudes
        const double beta = wall_law_data.DynamicViscosity / wall_law_data.SlipLength;
        const double aux_val = beta * rConditionData.wGauss;
        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            for (IndexType j_node = 0; j_node < TNumNodes; ++j_node) {
                const auto& r_v_j = wall_law_data.NodalWallVelocities[j_node];
                for (IndexType d1 = 0; d1 < TDim; ++d1) {
                    for (IndexType d2 = 0; d2 < TDim; ++d2) {
                        rRightHandSideVector[i_node*BlockSize + d2] += aux_val * rConditionData.N[i_node] * rConditionData.N[j_node] * tang_proj_mat(d1,d2) * r_v_j[d1];
                        rLeftHandSideMatrix(i_node*BlockSize + d1, j_node*BlockSize + d2) -= aux_val * rConditionData.N[i_node] * rConditionData.N[j_node] * tang_proj_mat(d1,d2);
                    }
                }
            }
        }
    }

    /**
     * @brief Add the LHS Navier-slip contribution
     * This method adds the Navier-slip LHS Gauss point contribution to the provided matrix
     * @tparam TConditionDataContainer Condition data container
     * @param rLeftHandSideMatrix Reference to the output LHS matrix
     * @param pCondition Pointer to the current condition
     * @param rCurrentProcessInfo Reference to current ProcessInfo container
     * @param rConditionData Reference to the condition data container
     */
    template<class TConditionDataContainer>
    static void AddLeftHandSideGaussPointContribution(
        MatrixType& rLeftHandSideMatrix,
        const Condition* pCondition,
        const ProcessInfo& rCurrentProcessInfo,
        const TConditionDataContainer& rConditionData)
    {
        // Create and fill the wall law auxiliar data container
        WallLawDataContainer wall_law_data;
        wall_law_data.Initialize(*pCondition, rConditionData);

        // Set the tangential projection matrix
        BoundedMatrix<double,TDim,TDim> tang_proj_mat;
        SetTangentialProjectionMatrix(rConditionData.Normal, tang_proj_mat);

        // Calculate the Navier-slip required magnitudes
        const double beta = wall_law_data.DynamicViscosity / wall_law_data.SlipLength;
        const double aux_val = beta * rConditionData.wGauss;
        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            for (IndexType j_node = 0; j_node < TNumNodes; ++j_node) {
                for (IndexType d1 = 0; d1 < TDim; ++d1) {
                    for (IndexType d2 = 0; d2 < TDim; ++d2) {
                        rLeftHandSideMatrix(i_node*BlockSize + d1, j_node*BlockSize + d2) -= aux_val * rConditionData.N[i_node] * rConditionData.N[j_node] * tang_proj_mat(d1,d2);
                    }
                }
            }
        }
    }

    /**
     * @brief Add the RHS Navier-slip contribution
     * This method adds the Navier-slip RHS Gauss point contribution to the provided vector
     * @tparam TConditionDataContainer Condition data container
     * @param rRightHandSideVector Reference to the output RHS vector
     * @param pCondition Pointer to the current condition
     * @param rCurrentProcessInfo Reference to current ProcessInfo container
     * @param rConditionData Reference to the condition data container
     */
    template<class TConditionDataContainer>
    static void AddRightHandSideGaussPointContribution(
        VectorType& rRightHandSideVector,
        const Condition* pCondition,
        const ProcessInfo& rCurrentProcessInfo,
        const TConditionDataContainer& rConditionData)
    {
        // Create and fill the wall law auxiliar data container
        WallLawDataContainer wall_law_data;
        wall_law_data.Initialize(*pCondition, rConditionData);

        // Set the tangential projection matrix
        BoundedMatrix<double,TDim,TDim> tang_proj_mat;
        SetTangentialProjectionMatrix(rConditionData.Normal, tang_proj_mat);

        // Calculate the Navier-slip required magnitudes
        const double beta = wall_law_data.DynamicViscosity / wall_law_data.SlipLength;
        const double aux_val = beta * rConditionData.wGauss;
        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            for (IndexType j_node = 0; j_node < TNumNodes; ++j_node) {
                const auto& r_v_j = wall_law_data.NodalWallVelocities[j_node];
                for (IndexType d1 = 0; d1 < TDim; ++d1) {
                    for (IndexType d2 = 0; d2 < TDim; ++d2) {
                        rRightHandSideVector[i_node*BlockSize + d2] += aux_val * rConditionData.N[i_node] * rConditionData.N[j_node] * tang_proj_mat(d1,d2) * r_v_j[d1];
                    }
                }
            }
        }
    }

    /**
     * @brief Check function
     * This function checks the current wall law input parameters
     * @param pCondition Pointer to the current condition
     * @param rCurrentProcessInfo Reference to the ProcessInfo data container
     * @return int Check status
     */
    static int Check(
        const Condition* pCondition,
        const ProcessInfo& rCurrentProcessInfo)
    {
        return 0;
    }

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
        buffer << "NavierSlipWallLaw";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "NavierSlipWallLaw";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const {}

    ///@}
    ///@name Friends
    ///@{


    ///@}
private:
    ///@name Private operations
    ///@{

    /**
     * @brief Auxiliary data container
     * This class serves as a temporary data container for the wall-related data
     */
    struct WallLawDataContainer
    {
        double SlipLength;
        double DynamicViscosity;
        array_1d<array_1d<double,3>,TNumNodes> NodalWallVelocities;

        template<class TConditionDataContainer>
        void Initialize(
            const Condition& rCondition,
            const TConditionDataContainer& rConditionData)
        {
            // Get dynamic viscosity from parent element properties
            // Note that this assumes constant viscosity within the element
            const auto& r_parent_element = rCondition.GetValue(NEIGHBOUR_ELEMENTS)[0];
            DynamicViscosity = r_parent_element.GetProperties().GetValue(DYNAMIC_VISCOSITY);

            // Save the nodal velocities and interpolate the slip length at current integration point
            SlipLength = 0.0;
            const auto& r_geom = rCondition.GetGeometry();
            for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node) {
                const auto& r_node = r_geom[i_node];
                SlipLength += rConditionData.N[i_node] * r_node.GetValue(SLIP_LENGTH);
                noalias(NodalWallVelocities[i_node]) = r_node.FastGetSolutionStepValue(VELOCITY) - r_node.FastGetSolutionStepValue(MESH_VELOCITY);
            }
            KRATOS_ERROR_IF(SlipLength < 1.0e-12) << "Negative or zero 'SLIP_LENGTH' in condition " << rCondition.Id() << "." << std::endl;
        }
    };

    //TODO: Remove this by a common implementation (FluidElementUtilities cannot be used yet)
    static void SetTangentialProjectionMatrix(
        const array_1d<double,3>& rUnitNormal,
        BoundedMatrix<double,TDim,TDim>& rTangProjMat)
    {
        noalias(rTangProjMat) = IdentityMatrix(TDim,TDim);
        for (std::size_t d1 = 0; d1 < TDim; ++d1) {
            for (std::size_t d2 = 0; d2 < TDim; ++d2) {
                rTangProjMat(d1,d2) -= rUnitNormal[d1]*rUnitNormal[d2];
            }
        }
    }

    ///@}
}; // Class NavierSlipWallLaw

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}
///@} addtogroup block

}  // namespace Kratos.
