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
 * This class implements the LHS and RHS contributions of the Navier-slip wall model
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

    using VectorType = Condition::VectorType;

    using MatrixType = Condition::MatrixType;

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
     * This method adds the Navier-slip LHS and RHS contributions to the provided matrices
     * @param rLeftHandSideMatrix Reference to the output LHS matrix
     * @param rRightHandSideVector Reference to the output RHS matrix
     * @param pCondition Pointer to the current condition
     * @param rCurrentProcessInfo Reference to current ProcessInfo con
     */
    static void AddWallModelLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const Condition* pCondition,
        const ProcessInfo& rCurrentProcessInfo)
    {
        // Create and fill the wall law auxiliar data container
        WallLawDataContainer wall_law_data;
        wall_law_data.Initialize(*pCondition);

        // Set the tangential projection matrix
        BoundedMatrix<double,TDim,TDim> tang_proj_mat;
        SetTangentialProjectionMatrix(wall_law_data.Normal, tang_proj_mat);

        // Calculate the Navier-slip contribution
        const SizeType n_gauss = wall_law_data.GaussPtsWeights.size();
        for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
            // Get current Gauss point data
            const double w_gauss = wall_law_data.GaussPtsWeights[i_gauss];
            const auto& N_gauss = row(wall_law_data.ShapeFunctionsContainer, i_gauss);
            
            // Interpolate slip length at current Gauss point
            double gp_slip_length = 0.0;
            for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
                gp_slip_length += N_gauss[i_node] * wall_law_data.NodalSlipLength[i_node];
            }

            // Assemble RHS and LHS contributions
            const double aux_val = w_gauss * wall_law_data.DynamicViscosity / gp_slip_length;
            for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
                for (IndexType j_node = 0; j_node < TNumNodes; ++j_node) {
                    const auto& r_v_j = wall_law_data.NodalWallVelocities[j_node];
                    for (IndexType d1 = 0; d1 < TDim; ++d1) {
                        for (IndexType d2 = 0; d2 < TDim; ++d2) {
                            rRightHandSideVector[i_node*BlockSize + d2] += aux_val * N_gauss[i_node] * N_gauss[j_node] * tang_proj_mat(d1,d2) * r_v_j[d1];
                            rLeftHandSideMatrix(i_node*BlockSize + d1, j_node*BlockSize + d2) -= aux_val * N_gauss[i_node] * N_gauss[j_node] * tang_proj_mat(d1,d2);
                        }
                    }
                }
            }
        }
    }

    /**
     * @brief Add the LHS Navier-slip contribution
     * This method adds the Navier-slip LHS contribution to the provided matrix
     * @param rLeftHandSideMatrix Reference to the output LHS matrix
     * @param pCondition Pointer to the current condition
     * @param rCurrentProcessInfo Reference to current ProcessInfo container
     */
    static void AddWallModelLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const Condition* pCondition,
        const ProcessInfo& rCurrentProcessInfo)
    {
        // Create and fill the wall law auxiliar data container
        WallLawDataContainer wall_law_data;
        wall_law_data.Initialize(*pCondition);

        // Set the tangential projection matrix
        BoundedMatrix<double,TDim,TDim> tang_proj_mat;
        SetTangentialProjectionMatrix(wall_law_data.Normal, tang_proj_mat);

        // Calculate the Navier-slip contribution
        const SizeType n_gauss = wall_law_data.GaussPtsWeights.size();
        for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
            // Get current Gauss point data
            const double w_gauss = wall_law_data.GaussPtsWeights[i_gauss];
            const auto& N_gauss = row(wall_law_data.ShapeFunctionsContainer, i_gauss);
            
            // Interpolate slip length at current Gauss point
            double gp_slip_length = 0.0;
            for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
                gp_slip_length += N_gauss[i_node] * wall_law_data.NodalSlipLength[i_node];
            }

            // Assemble RHS and LHS contributions
            const double aux_val = w_gauss * wall_law_data.DynamicViscosity / gp_slip_length;
            for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
                for (IndexType j_node = 0; j_node < TNumNodes; ++j_node) {
                    for (IndexType d1 = 0; d1 < TDim; ++d1) {
                        for (IndexType d2 = 0; d2 < TDim; ++d2) {
                            rLeftHandSideMatrix(i_node*BlockSize + d1, j_node*BlockSize + d2) -= aux_val * N_gauss[i_node] * N_gauss[j_node] * tang_proj_mat(d1,d2);
                        }
                    }
                }
            }
        }
    }

    /**
     * @brief Add the RHS Navier-slip contribution
     * This method adds the Navier-slip RHS contribution to the provided vector
     * @param rRightHandSideVector Reference to the output RHS vector
     * @param pCondition Pointer to the current condition
     * @param rCurrentProcessInfo Reference to current ProcessInfo container
     */
    static void AddWallModelRightHandSide(
        VectorType& rRightHandSideVector,
        const Condition* pCondition,
        const ProcessInfo& rCurrentProcessInfo)
    {
        // Create and fill the wall law auxiliar data container
        WallLawDataContainer wall_law_data;
        wall_law_data.Initialize(*pCondition);

        // Set the tangential projection matrix
        BoundedMatrix<double,TDim,TDim> tang_proj_mat;
        SetTangentialProjectionMatrix(wall_law_data.Normal, tang_proj_mat);

        // Calculate the Navier-slip contribution
        const SizeType n_gauss = wall_law_data.GaussPtsWeights.size();
        for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
            // Get current Gauss point data
            const double w_gauss = wall_law_data.GaussPtsWeights[i_gauss];
            const auto& N_gauss = row(wall_law_data.ShapeFunctionsContainer, i_gauss);
            
            // Interpolate slip length at current Gauss point
            double gp_slip_length = 0.0;
            for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
                gp_slip_length += N_gauss[i_node] * wall_law_data.NodalSlipLength[i_node];
            }

            // Assemble RHS and LHS contributions
            const double aux_val = w_gauss * wall_law_data.DynamicViscosity / gp_slip_length;
            for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
                for (IndexType j_node = 0; j_node < TNumNodes; ++j_node) {
                    const auto& r_v_j = wall_law_data.NodalWallVelocities[j_node];
                    for (IndexType d1 = 0; d1 < TDim; ++d1) {
                        for (IndexType d2 = 0; d2 < TDim; ++d2) {
                            rRightHandSideVector[i_node*BlockSize + d2] += aux_val * N_gauss[i_node] * N_gauss[j_node] * tang_proj_mat(d1,d2) * r_v_j[d1];
                        }
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
        double DynamicViscosity;
        array_1d<double,3> Normal;
        Vector GaussPtsWeights;
        Matrix ShapeFunctionsContainer;
        array_1d<double, TNumNodes> NodalSlipLength;
        array_1d<array_1d<double,3>,TNumNodes> NodalWallVelocities;

        void Initialize(const Condition& rCondition)
        {
            // Get dynamic viscosity from parent element properties
            // Note that this assumes constant viscosity within the element
            const auto& r_parent_element = rCondition.GetValue(NEIGHBOUR_ELEMENTS)[0];
            DynamicViscosity = r_parent_element.GetProperties().GetValue(DYNAMIC_VISCOSITY);

            // Compute condition unit normal vector
            // Note that in here we are assuming a constant normal over the entire condition
            const auto& r_geom = rCondition.GetGeometry();
            Normal = r_geom.UnitNormal(0, GeometryData::IntegrationMethod::GI_GAUSS_1);

            // Calculate the required Gauss points integration data
            const auto& r_integration_pts = r_geom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_2);
            const SizeType n_gauss = r_integration_pts.size();
            r_geom.DeterminantOfJacobian(GaussPtsWeights, GeometryData::IntegrationMethod::GI_GAUSS_2);
            for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
                GaussPtsWeights[i_gauss] = GaussPtsWeights[i_gauss] * r_integration_pts[i_gauss].Weight();
            }
            ShapeFunctionsContainer = r_geom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_2);

            // Save the nodal velocities and interpolate the slip length at current integration point
            for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node) {
                const auto& r_node = r_geom[i_node];
                const double aux_slip_length = r_node.GetValue(SLIP_LENGTH);
                KRATOS_ERROR_IF(aux_slip_length < 1.0e-12) << "Negative or zero 'SLIP_LENGTH' at node " << r_node.Id() << "." << std::endl;
                NodalSlipLength[i_node] = aux_slip_length;
                noalias(NodalWallVelocities[i_node]) = r_node.FastGetSolutionStepValue(VELOCITY) - r_node.FastGetSolutionStepValue(MESH_VELOCITY);
            }
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
