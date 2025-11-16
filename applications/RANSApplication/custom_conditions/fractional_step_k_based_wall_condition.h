//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#ifndef KRATOS_FS_HIGH_RE_K_WALL_CONDITION_H
#define KRATOS_FS_HIGH_RE_K_WALL_CONDITION_H

// System includes
#include <iostream>
#include <string>

#include "includes/deprecated_variables.h"
#include "includes/kratos_flags.h"

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/condition.h"
#include "includes/define.h"
#include "includes/process_info.h"
#include "includes/serializer.h"
#include "includes/variables.h"

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

namespace Kratos
{
///@addtogroup RANSApplication
///@name Kratos Classes
///@{

/// Implements a wall condition for the monolithic formulation.
/**
  It is intended to be used in combination with ASGS and VMS elements or their derived classes
  and the ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent time scheme, which supports
  slip conditions.
  This condition will add a wall stress term to all nodes identified with SLIP==true (in the
  nodal flags). This stress term is determined according to the wall distance provided as Y_WALL.
  @see ASGS2D,ASGS3D,VMS,ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent
 */
template <unsigned int TDim, unsigned int TNumNodes = TDim>
class FractionalStepKBasedWallCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FractionalStepKBasedWallCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(FractionalStepKBasedWallCondition);

    using NodeType = Node;

    using PropertiesType = Properties;

    using GeometryType = Geometry<NodeType>;

    using NodesArrayType = Geometry<NodeType>::PointsArrayType;

    using VectorType = Vector;

    using MatrixType = Matrix;

    using IndexType = std::size_t;

    using SizeType = std::size_t;

    using EquationIdVectorType = std::vector<std::size_t>;

    using DofsVectorType = std::vector<Dof<double>::Pointer>;

    using DofsArrayType = PointerVectorSet<Dof<double>, IndexedObject>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
      @param NewId Index for the new condition
      */
    FractionalStepKBasedWallCondition(
        IndexType NewId = 0)
    : Condition(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    FractionalStepKBasedWallCondition(
        IndexType NewId,
        const NodesArrayType& ThisNodes)
    : Condition(NewId, ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    FractionalStepKBasedWallCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    FractionalStepKBasedWallCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
    {
    }

    /// Copy constructor.
    FractionalStepKBasedWallCondition(
        FractionalStepKBasedWallCondition const& rOther)
    : Condition(rOther)
    {
    }

    /// Destructor.
    ~FractionalStepKBasedWallCondition() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Copy constructor
    FractionalStepKBasedWallCondition& operator=(
        FractionalStepKBasedWallCondition const& rOther)
    {
        Condition::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create a new FractionalStepKBasedWallCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<FractionalStepKBasedWallCondition>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    Condition::Pointer Create(
        IndexType NewId,
        Condition::GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<FractionalStepKBasedWallCondition>(
            NewId, pGeom, pProperties);
    }

    /**
     * Clones the selected element variables, creating a new one
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */

    Condition::Pointer Clone(
        IndexType NewId,
        NodesArrayType const& rThisNodes) const override
    {
        Condition::Pointer p_new_condition =
            Create(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

        p_new_condition->SetData(this->GetData());
        p_new_condition->SetFlags(this->GetFlags());

        return p_new_condition;
    }

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        VectorType RHS;
        this->CalculateLocalSystem(rLeftHandSideMatrix, RHS, rCurrentProcessInfo);
    }

    /// Calculate wall stress term for all nodes with SLIP set.
    /**
      @param rDampingMatrix Left-hand side matrix
      @param rRightHandSideVector Right-hand side vector
      @param rCurrentProcessInfo ProcessInfo instance (unused)
      */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const unsigned int step = rCurrentProcessInfo[FRACTIONAL_STEP];
        if (step == 1) {
            // Initialize local contributions
            const SizeType local_size = TDim * TNumNodes;

            if (rLeftHandSideMatrix.size1() != local_size) {
                rLeftHandSideMatrix.resize(local_size, local_size, false);
            }

            if (rRightHandSideVector.size() != local_size) {
                rRightHandSideVector.resize(local_size, false);
            }

            noalias(rLeftHandSideMatrix) = ZeroMatrix(local_size, local_size);
            noalias(rRightHandSideVector) = ZeroVector(local_size);

            this->ApplyNeumannCondition(rLeftHandSideMatrix, rRightHandSideVector);

            this->ApplyWallLaw(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
        } else {
            if (this->Is(INTERFACE) && step == 5) {
                // add here a mass matrix in the form Dt/rho_equivalent_structure to the lhs alone
                const double N = 1.0 / static_cast<double>(TNumNodes);
                array_1d<double, 3> r_normal;
                this->CalculateNormal(r_normal); // this already contains the area
                const double Area = norm_2(r_normal);

                if (rLeftHandSideMatrix.size1() != TNumNodes) {
                    rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
                }

                if (rRightHandSideVector.size() != TNumNodes) {
                    rRightHandSideVector.resize(TNumNodes, false);
                }

                noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);
                noalias(rRightHandSideVector) = ZeroVector(TNumNodes);

                const double dt = rCurrentProcessInfo[DELTA_TIME];
                const double equivalent_structural_density = rCurrentProcessInfo[DENSITY];
                const double diag_term = dt * Area * N / (equivalent_structural_density);

                for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode) {
                    rLeftHandSideMatrix(iNode, iNode) = diag_term;
                }
            } else {
                if (rLeftHandSideMatrix.size1() != 0) {
                    rLeftHandSideMatrix.resize(0, 0, false);
                }

                if (rRightHandSideVector.size() != 0) {
                    rRightHandSideVector.resize(0, false);
                }
            }
        }
    }

    /// Check that all data required by this condition is available and reasonable
    int Check(const ProcessInfo& rCurrentProcessInfo) const override
    {
        KRATOS_TRY;

        int Check = Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area > 0

        if (Check != 0) {
            return Check;
        } else {
            // Checks on nodes

            // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
            for (unsigned int i = 0; i < this->GetGeometry().size(); ++i) {
                if (this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
                    KRATOS_THROW_ERROR(std::invalid_argument,
                                       "missing VELOCITY variable on solution "
                                       "step data for node ",
                                       this->GetGeometry()[i].Id());
                if (this->GetGeometry()[i].SolutionStepsDataHas(MESH_VELOCITY) == false)
                    KRATOS_THROW_ERROR(std::invalid_argument,
                                       "missing MESH_VELOCITY variable on "
                                       "solution step data for node ",
                                       this->GetGeometry()[i].Id());
                if (this->GetGeometry()[i].SolutionStepsDataHas(NORMAL) == false)
                    KRATOS_THROW_ERROR(std::invalid_argument,
                                       "missing NORMAL variable on solution "
                                       "step data for node ",
                                       this->GetGeometry()[i].Id());
                if (this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
                    this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
                    this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
                    KRATOS_THROW_ERROR(
                        std::invalid_argument,
                        "missing VELOCITY component degree of freedom on node ",
                        this->GetGeometry()[i].Id());
            }

            return Check;
        }

        KRATOS_CATCH("");
    }

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        if (RansCalculationUtilities::IsWallFunctionActive(this->GetGeometry())) {
            this->SetValue(GAUSS_RANS_Y_PLUS, Vector(this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod())));
        }

        KRATOS_CATCH("");
    }

    /// Provides the global indices for each one of this element's local rows.
    /** This determines the elemental equation ID vector for all elemental DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo the current process info object (unused)
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(
        DofsVectorType& ConditionDofList,
        const ProcessInfo& CurrentProcessInfo) const override;

    /// Returns VELOCITY_X, VELOCITY_Y, (VELOCITY_Z,) for each node
    /**
     * @param Values Vector of nodal unknowns
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetValuesVector(
        Vector& Values,
        int Step = 0) const override
    {
        const SizeType local_size = TDim * TNumNodes;
        unsigned int LocalIndex = 0;

        if (Values.size() != local_size) {
            Values.resize(local_size, false);
        }

        for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode) {
            const array_1d<double, 3>& rVelocity =
                this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, Step);
            for (unsigned int d = 0; d < TDim; ++d) {
                Values[LocalIndex++] = rVelocity[d];
            }
        }
    }

    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3>>& rVariable,
        std::vector<array_1d<double, 3>>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        rValues.resize(1);
        if (rVariable == NORMAL) {
            this->CalculateNormal(rValues[0]);
        } else {
            /*
             The cast is done to avoid modification of the element's data. Data
             modification would happen if rVariable is not stored now (would
             initialize a pointer to &rVariable with associated value of 0.0).
             This is catastrophic if the variable referenced goes out of scope.
             */
            const FractionalStepKBasedWallCondition* const_this =
                static_cast<const FractionalStepKBasedWallCondition*>(this);
            rValues[0] = const_this->GetValue(rVariable);
        }
    }

    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        rValues.resize(1);
        /*
         The cast is done to avoid modification of the element's data. Data
         modification would happen if rVariable is not stored now (would
         initialize a pointer to &rVariable with associated value of 0.0). This
         is catastrophic if the variable referenced goes out of scope.
         */
        const FractionalStepKBasedWallCondition* const_this =
            static_cast<const FractionalStepKBasedWallCondition*>(this);
        rValues[0] = const_this->GetValue(rVariable);
    }

    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 6>>& rVariable,
        std::vector<array_1d<double, 6>>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        rValues.resize(1);
        const FractionalStepKBasedWallCondition* const_this =
            static_cast<const FractionalStepKBasedWallCondition*>(this);
        rValues[0] = const_this->GetValue(rVariable);
    }

    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        rValues.resize(1);
        const FractionalStepKBasedWallCondition* const_this =
            static_cast<const FractionalStepKBasedWallCondition*>(this);
        rValues[0] = const_this->GetValue(rVariable);
    }

    void CalculateOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        rValues.resize(1);
        const FractionalStepKBasedWallCondition* const_this =
            static_cast<const FractionalStepKBasedWallCondition*>(this);
        rValues[0] = const_this->GetValue(rVariable);
    }

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        this->PrintInfo(buffer);
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FractionalStepKBasedWallCondition" << TDim << "D #" << this->Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        this->pGetGeometry()->PrintData(rOStream);
    }

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    void CalculateNormal(array_1d<double, 3>& An) const;

    /// Commpute the wall stress and add corresponding terms to the system contributions.
    /**
      @param rLocalMatrix Local system matrix
      @param rLocalVector Local right hand side
      */
    void ApplyWallLaw(
        MatrixType& rLocalMatrix,
        VectorType& rLocalVector,
        const ProcessInfo& rCurrentProcessInfo)
    {
        using namespace RansCalculationUtilities;

        if (IsWallFunctionActive(this->GetGeometry())) {
            const auto& r_geometry = this->GetGeometry();
            // Get Shape function data
            Vector gauss_weights;
            Matrix shape_functions;
            CalculateConditionGeometryData(r_geometry, this->GetIntegrationMethod(),
                                           gauss_weights, shape_functions);
            const IndexType num_gauss_points = gauss_weights.size();
            const double eps = std::numeric_limits<double>::epsilon();
            const double c_mu_25 = std::pow(rCurrentProcessInfo[TURBULENCE_RANS_C_MU], 0.25);
            const double kappa = rCurrentProcessInfo[VON_KARMAN];

            // get parent element
            auto& r_parent_element = this->GetValue(NEIGHBOUR_ELEMENTS)[0];

            // get fluid properties from parent element
            const auto& r_elem_properties = r_parent_element.GetProperties();
            const double rho = r_elem_properties[DENSITY];
            const double nu = r_elem_properties[DYNAMIC_VISCOSITY] / rho;

            // get surface properties from condition
            const PropertiesType& r_cond_properties = this->GetProperties();
            const double beta = r_cond_properties.GetValue(WALL_SMOOTHNESS_BETA);
            const double y_plus_limit = r_cond_properties.GetValue(RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT);
            const double inv_kappa = 1.0 / kappa;

            double tke;
            array_1d<double, 3> wall_velocity, fluid_velocity, mesh_velocity;
            const double wall_height = this->GetValue(DISTANCE);

            for (size_t g = 0; g < num_gauss_points; ++g)
            {
                const Vector& gauss_shape_functions = row(shape_functions, g);

                FluidCalculationUtilities::EvaluateInPoint(
                    r_geometry, gauss_shape_functions,
                    std::tie(tke, TURBULENT_KINETIC_ENERGY),
                    std::tie(fluid_velocity, VELOCITY),
                    std::tie(mesh_velocity, MESH_VELOCITY));

                noalias(wall_velocity) = fluid_velocity - mesh_velocity;

                const double wall_velocity_magnitude = norm_2(wall_velocity);

                double y_plus{0.0}, u_tau{0.0};
                CalculateYPlusAndUtau(y_plus, u_tau, wall_velocity_magnitude,
                                      wall_height, nu, kappa, beta);
                y_plus = std::max(y_plus, y_plus_limit);

                u_tau = std::max(c_mu_25 * std::sqrt(std::max(tke, 0.0)),
                                 wall_velocity_magnitude /
                                     (inv_kappa * std::log(y_plus) + beta));

                if (wall_velocity_magnitude > eps) {
                    const double value = rho * std::pow(u_tau, 2) *
                                         gauss_weights[g] / wall_velocity_magnitude;
                    for (size_t a = 0; a < r_geometry.PointsNumber(); ++a) {
                        for (size_t dim = 0; dim < TDim; ++dim) {
                            for (size_t b = 0; b < r_geometry.PointsNumber(); ++b) {
                                rLocalMatrix(a * TDim + dim, b * TDim + dim) +=
                                    gauss_shape_functions[a] *
                                    gauss_shape_functions[b] * value;
                            }
                            rLocalVector[a * TDim + dim] -=
                                gauss_shape_functions[a] * value * wall_velocity[dim];
                        }
                    }
                }
            }
        }
    }

    /// Apply boundary terms to allow imposing a pressure (normal stress), with a correction to prevent inflow.
    /** This correction should prevent numerical problems arising from inflow in
     * outflow areas, typically due to vortices. exiting the domain.
     * @param rLocalMatrix Local LHS matrix
     * @param rLocalVector Local RHS vector
     */
    void ApplyNeumannCondition(
        MatrixType& rLocalMatrix,
        VectorType& rLocalVector);

    ///@}

private:
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_TRY

        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);

        KRATOS_CATCH("");
    }
    void load(Serializer& rSerializer) override
    {
        KRATOS_TRY

        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);

        KRATOS_CATCH("");
    }

    ///@}

}; // Class FractionalStepKBasedWallCondition

///@}
///@name Input and output
///@{

/// input stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::istream& operator>>(
    std::istream& rIStream,
    FractionalStepKBasedWallCondition<TDim, TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const FractionalStepKBasedWallCondition<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_FS_HIGH_RE_K_WALL_CONDITION_H
