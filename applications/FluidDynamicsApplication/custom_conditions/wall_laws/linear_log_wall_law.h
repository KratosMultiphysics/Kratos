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
//                   Jordi Cotela
//                   Riccardo Rossi
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

///@name Kratos Classes
///@{

/**
 * @brief Linear-log law LHS and RHS contribution implementation
 * This class implements the LHS and RHS Gauss point contributions of the linear-log wall model
 * This class should be used in combination with an incompressible Navier-Stokes (or Stokes) condition
 * implementing the remaining terms (see @NavierStokesWallCondition).
 * @tparam TDim Number of dimensions
 * @tparam TNumNodes Number of nodes
 */
template<std::size_t TDim, std::size_t TNumNodes>
class LinearLogWallLaw
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LinearLogWallLaw
    KRATOS_CLASS_POINTER_DEFINITION(LinearLogWallLaw);

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
    LinearLogWallLaw() = delete;

    /// Copy constructor.
    LinearLogWallLaw(LinearLogWallLaw const& rOther) = delete;

    /// Destructor.
    ~LinearLogWallLaw() = default;

    ///@}
    ///@name Operators
    ///@name Operations
    ///@{

    /**
     * @brief Add the LHS and RHS Navier-slip contributions
     * This method adds the Navier-slip LHS and RHS Gauss point contributions to the provided matrices
     * @param rLeftHandSideMatrix Reference to the output LHS matrix
     * @param rRightHandSideVector Reference to the output RHS matrix
     * @param pCondition Pointer to the current condition
     * @param rCurrentProcessInfo Reference to current ProcessInfo container
     * @param rConditionData Reference to the condition data container
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

        // Calculate the linear-log law contribution
        // Note that we use a second order Gauss-Lobatto quadrature in this case in order to obtain a diagonal matrix contribution
        // This is of special interest to ease convergence by avoiding the off-diagonal terms in presence of large friction velocity
        const auto& r_geom = pCondition->GetGeometry();
        const double w_gauss_lobatto = r_geom.DomainSize() / TNumNodes;
        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            // Calculate current Gauss-Lobatto point (node) wall velocity norm
            const auto& r_aux_v = wall_law_data.NodalWallVelocities[i_node];
            const double wall_vel_norm = norm_2(r_aux_v);

            // Do nothing in case of zero wall velocity
            if (wall_vel_norm > ZeroTol) {
                // Calculate the friction (shear) velocity
                const double u_tau = wall_law_data.CalculateFrictionVelocity(wall_vel_norm, wall_law_data.NodalYWallValues[i_node]);

                // Assembly the current contribution
                const double tmp = w_gauss_lobatto * std::pow(u_tau,2) * wall_law_data.Density / wall_vel_norm;
                for (IndexType d = 0; d < TDim; ++d) {
                    rRightHandSideVector(i_node*BlockSize + d) -= tmp * r_aux_v[d];
                    rLeftHandSideMatrix(i_node*BlockSize + d, i_node*BlockSize + d) += tmp;
                }
            }
        }
    }

    /**
     * @brief Add the LHS Navier-slip contribution
     * This method adds the Navier-slip LHS Gauss point contribution to the provided matrix
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

        // Calculate the linear-log law contribution
        // Note that we use a second order Gauss-Lobatto quadrature in this case in order to obtain a diagonal matrix contribution
        // This is of special interest to ease convergence by avoiding the off-diagonal terms in presence of large friction velocity
        const auto& r_geom = pCondition->GetGeometry();
        const double w_gauss_lobatto = r_geom.DomainSize() / TNumNodes;
        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            // Calculate current Gauss-Lobatto point (node) wall velocity norm
            const auto& r_aux_v = wall_law_data.NodalWallVelocities[i_node];
            const double wall_vel_norm = norm_2(r_aux_v);

            // Do nothing in case of zero wall velocity
            if (wall_vel_norm > ZeroTol) {
                // Calculate the friction (shear) velocity
                const double u_tau = wall_law_data.CalculateFrictionVelocity(wall_vel_norm, wall_law_data.NodalYWallValues[i_node]);

                // Assembly the current contribution
                const double tmp = w_gauss_lobatto * std::pow(u_tau,2) * wall_law_data.Density / wall_vel_norm;
                for (IndexType d = 0; d < TDim; ++d) {
                    rLeftHandSideMatrix(i_node*BlockSize + d, i_node*BlockSize + d) += tmp;
                }
            }
        }
    }

    /**
     * @brief Add the RHS Navier-slip contribution
     * This method adds the Navier-slip RHS Gauss point contribution to the provided vector
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

        // Calculate the linear-log law contribution
        // Note that we use a second order Gauss-Lobatto quadrature in this case in order to obtain a diagonal matrix contribution
        // This is of special interest to ease convergence by avoiding the off-diagonal terms in presence of large friction velocity
        const auto& r_geom = pCondition->GetGeometry();
        const double w_gauss_lobatto = r_geom.DomainSize() / TNumNodes;
        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            // Calculate current Gauss-Lobatto point (node) wall velocity norm
            const auto& r_aux_v = wall_law_data.NodalWallVelocities[i_node];
            const double wall_vel_norm = norm_2(r_aux_v);

            // Do nothing in case of zero wall velocity
            if (wall_vel_norm > ZeroTol) {
                // Calculate the friction (shear) velocity
                const double u_tau = wall_law_data.CalculateFrictionVelocity(wall_vel_norm, wall_law_data.NodalYWallValues[i_node]);

                // Assembly the current contribution
                const double tmp = w_gauss_lobatto * std::pow(u_tau,2) * wall_law_data.Density / wall_vel_norm;
                for (IndexType d = 0; d < TDim; ++d) {
                    rRightHandSideVector(i_node*BlockSize + d) -= tmp * r_aux_v[d];
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
        buffer << "LinearLogWallLaw";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "LinearLogWallLaw";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const {}

    ///@}
    ///@name Friends
    ///@{


    ///@}
private:
    ///@name Private variables
    ///@{
    
    static constexpr double BPlus = 5.2; // Dimensionless velocity u+ constant
    static constexpr double InvKappa = 1.0/0.41; // Inverse of Von Karman's kappa
    static constexpr double YPlusLimit = 10.9931899; // Limit between linear and log regions
    static constexpr double ZeroTol = 1.0e-12; // Auxiliary tolerance to check zero values

    ///@name Private operations
    ///@{

    /**
     * @brief Auxiliary data container
     * This class serves as a temporary data container for the wall-related data
     */
    struct WallLawDataContainer
    {
        double Density;
        double KinematicViscosity;
        array_1d<double,TNumNodes> NodalYWallValues;
        array_1d<array_1d<double,3>,TNumNodes> NodalWallVelocities;

        void Initialize(const Condition& rCondition)
        {
            // Get fluid parameters from parent element properties
            // Note that this assumes constant viscosity within the element
            const auto& r_parent_element = rCondition.GetValue(NEIGHBOUR_ELEMENTS)[0];
            Density = r_parent_element.GetProperties().GetValue(DENSITY);
            KinematicViscosity = r_parent_element.GetProperties().GetValue(DYNAMIC_VISCOSITY) / Density;

            // Save the nodal velocities and interpolate the slip length at current integration point
            const auto& r_geom = rCondition.GetGeometry();
            for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node) {
                const auto& r_node = r_geom[i_node];
                const double aux_y_wall = r_node.GetValue(Y_WALL);
                KRATOS_ERROR_IF(aux_y_wall < ZeroTol) << "Negative or zero 'Y_WALL' in condition " << rCondition.Id() << "." << std::endl;
                NodalYWallValues[i_node] = aux_y_wall;
                noalias(NodalWallVelocities[i_node]) = r_node.FastGetSolutionStepValue(VELOCITY) - r_node.FastGetSolutionStepValue(MESH_VELOCITY);
            }
        }

        double CalculateFrictionVelocity(
            const double WallVelocityNorm,
            const double YWall)
        {
            // If wall velocity is non-zero add the wall contribution
            double u_tau = 0.0;
            if (WallVelocityNorm > ZeroTol) {
                // Linear region
                u_tau = sqrt(WallVelocityNorm * KinematicViscosity / YWall); // Friction velocity (shear velocity)
                double y_plus = YWall * u_tau / KinematicViscosity; // Wall coordinate

                // Log region (local Newton-Raphson iterative problem)
                if (y_plus > YPlusLimit) {
                    // Initialize the iterative procedure
                    double dx = 1e10;
                    const double rel_tol = 1e-6;
                    const std::size_t max_it = 100;
                    double u_plus = InvKappa * log(y_plus) + BPlus; // Dimensionless velocity

                    std::size_t it = 0;
                    while (it < 100 && std::abs(dx) > rel_tol * u_tau) {
                        // Newton-Raphson iteration
                        double f = u_tau * u_plus - WallVelocityNorm;
                        double df = u_plus + InvKappa;
                        dx = f/df;

                        // Update variables
                        u_tau -= dx;
                        y_plus = YWall * u_tau / KinematicViscosity;
                        u_plus = InvKappa * log(y_plus) + BPlus;
                        ++it;
                    }
                    KRATOS_WARNING_IF("LinearLogWallLaw", it == max_it) << "Wall condition Newton-Raphson did not converge. Residual is " << dx << "." << std::endl;
                }
            }

            return u_tau;
        }
    };

    ///@}
}; // Class LinearLogWallLaw

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}
///@} addtogroup block

}  // namespace Kratos.
