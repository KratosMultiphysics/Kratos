/*
==============================================================================
Kratos Fluid Dynamics Application
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

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

    typedef Node<3> NodeType;

    typedef Properties PropertiesType;

    typedef Geometry<NodeType> GeometryType;

    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector<Dof<double>::Pointer> DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
      @param NewId Index for the new condition
      */
    FractionalStepKBasedWallCondition(IndexType NewId = 0) : Condition(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    FractionalStepKBasedWallCondition(IndexType NewId, const NodesArrayType& ThisNodes)
        : Condition(NewId, ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    FractionalStepKBasedWallCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    FractionalStepKBasedWallCondition(IndexType NewId,
                                      GeometryType::Pointer pGeometry,
                                      PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {
    }

    /// Copy constructor.
    FractionalStepKBasedWallCondition(FractionalStepKBasedWallCondition const& rOther)
        : Condition(rOther)
    {
    }

    /// Destructor.
    ~FractionalStepKBasedWallCondition() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /// Copy constructor
    FractionalStepKBasedWallCondition& operator=(FractionalStepKBasedWallCondition const& rOther)
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
    Condition::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<FractionalStepKBasedWallCondition>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    Condition::Pointer Create(IndexType NewId,
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

    Condition::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const override
    {
        Condition::Pointer pNewCondition =
            Create(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

        pNewCondition->SetData(this->GetData());
        pNewCondition->SetFlags(this->GetFlags());

        return pNewCondition;
    }

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
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
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override
    {
        const ProcessInfo& r_process_info = rCurrentProcessInfo;
        unsigned int step = r_process_info[FRACTIONAL_STEP];
        if (step == 1)
        {
            // Initialize local contributions
            const SizeType LocalSize = TDim * TNumNodes;

            if (rLeftHandSideMatrix.size1() != LocalSize)
                rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
            if (rRightHandSideVector.size() != LocalSize)
                rRightHandSideVector.resize(LocalSize, false);

            noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
            noalias(rRightHandSideVector) = ZeroVector(LocalSize);

            this->ApplyNeumannCondition(rLeftHandSideMatrix, rRightHandSideVector);

            this->ApplyWallLaw(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
        }
        else
        {
            if (this->Is(INTERFACE) && step == 5)
            {
                // add here a mass matrix in the form Dt/rho_equivalent_structure to the lhs alone
                const double N = 1.0 / static_cast<double>(TNumNodes);
                array_1d<double, 3> rNormal;
                this->CalculateNormal(rNormal); // this already contains the area
                const double Area = norm_2(rNormal);

                if (rLeftHandSideMatrix.size1() != TNumNodes)
                    rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
                if (rRightHandSideVector.size() != TNumNodes)
                    rRightHandSideVector.resize(TNumNodes, false);

                noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);
                noalias(rRightHandSideVector) = ZeroVector(TNumNodes);

                const double dt = r_process_info[DELTA_TIME];
                const double equivalent_structural_density = r_process_info[DENSITY];
                const double diag_term = dt * Area * N / (equivalent_structural_density);

                for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
                {
                    rLeftHandSideMatrix(iNode, iNode) = diag_term;
                }
            }
            else
            {
                if (rLeftHandSideMatrix.size1() != 0)
                    rLeftHandSideMatrix.resize(0, 0, false);

                if (rRightHandSideVector.size() != 0)
                    rRightHandSideVector.resize(0, false);
            }
        }
    }

    /// Check that all data required by this condition is available and reasonable
    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        int Check = Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area > 0

        if (Check != 0)
        {
            return Check;
        }
        else
        {
            // Check that all required variables have been registered
            if (VELOCITY.Key() == 0)
                KRATOS_THROW_ERROR(std::invalid_argument,
                                   "VELOCITY Key is 0. Check if the "
                                   "application was correctly registered.",
                                   "");
            if (MESH_VELOCITY.Key() == 0)
                KRATOS_THROW_ERROR(std::invalid_argument,
                                   "MESH_VELOCITY Key is 0. Check if the "
                                   "application was correctly registered.",
                                   "");
            if (NORMAL.Key() == 0)
                KRATOS_THROW_ERROR(std::invalid_argument,
                                   "NORMAL Key is 0. Check if the application "
                                   "was correctly registered.",
                                   "")
            if (Y_WALL.Key() == 0)
                KRATOS_THROW_ERROR(std::invalid_argument,
                                   "Y_WALL Key is 0. Check if the application "
                                   "was correctly registered.",
                                   "")

            // Checks on nodes

            // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
            for (unsigned int i = 0; i < this->GetGeometry().size(); ++i)
            {
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

    void Initialize() override
    {
        KRATOS_TRY;

        if (RansCalculationUtilities::IsWallFunctionActive(*this))
        {
            const array_1d<double, 3>& rNormal = this->GetValue(NORMAL);
            KRATOS_ERROR_IF(norm_2(rNormal) == 0.0)
                << "NORMAL must be calculated before using this "
                << this->Info() << "\n";

            KRATOS_ERROR_IF(this->GetValue(NEIGHBOUR_ELEMENTS).size() == 0)
                << this->Info() << " cannot find parent element\n";

            const double nu = RansCalculationUtilities::EvaluateInParentCenter(
                KINEMATIC_VISCOSITY, *this);
            KRATOS_ERROR_IF(nu == 0.0) << "KINEMATIC_VISCOSITY is not defined "
                                          "in the parent element of "
                                       << this->Info() << "\n.";

            mWallHeight = RansCalculationUtilities::CalculateWallHeight(*this, rNormal);
        }

        KRATOS_CATCH("");
    }

    /// Provides the global indices for each one of this element's local rows.
    /** This determines the elemental equation ID vector for all elemental DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo the current process info object (unused)
     */
    void EquationIdVector(EquationIdVectorType& rResult,
                          const ProcessInfo& rCurrentProcessInfo) const override;

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(DofsVectorType& ConditionDofList,
                    const ProcessInfo& CurrentProcessInfo) const override;

    /// Returns VELOCITY_X, VELOCITY_Y, (VELOCITY_Z,) for each node
    /**
     * @param Values Vector of nodal unknowns
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetValuesVector(Vector& Values, int Step = 0) const override
    {
        const SizeType LocalSize = TDim * TNumNodes;
        unsigned int LocalIndex = 0;

        if (Values.size() != LocalSize)
            Values.resize(LocalSize, false);

        for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
        {
            const array_1d<double, 3>& rVelocity =
                this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, Step);
            for (unsigned int d = 0; d < TDim; ++d)
                Values[LocalIndex++] = rVelocity[d];
        }
    }

    void GetValueOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                     std::vector<array_1d<double, 3>>& rValues,
                                     const ProcessInfo& rCurrentProcessInfo) override
    {
        rValues.resize(1);
        if (rVariable == NORMAL)
        {
            this->CalculateNormal(rValues[0]);
        }
        else
        {
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

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
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

    void GetValueOnIntegrationPoints(const Variable<array_1d<double, 6>>& rVariable,
                                     std::vector<array_1d<double, 6>>& rValues,
                                     const ProcessInfo& rCurrentProcessInfo) override
    {
        rValues.resize(1);
        const FractionalStepKBasedWallCondition* const_this =
            static_cast<const FractionalStepKBasedWallCondition*>(this);
        rValues[0] = const_this->GetValue(rVariable);
    }

    void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                                     std::vector<Vector>& rValues,
                                     const ProcessInfo& rCurrentProcessInfo) override
    {
        rValues.resize(1);
        const FractionalStepKBasedWallCondition* const_this =
            static_cast<const FractionalStepKBasedWallCondition*>(this);
        rValues[0] = const_this->GetValue(rVariable);
    }

    void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                     std::vector<Matrix>& rValues,
                                     const ProcessInfo& rCurrentProcessInfo) override
    {
        rValues.resize(1);
        const FractionalStepKBasedWallCondition* const_this =
            static_cast<const FractionalStepKBasedWallCondition*>(this);
        rValues[0] = const_this->GetValue(rVariable);
    }

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

    void CalculateNormal(array_1d<double, 3>& An);

    /// Commpute the wall stress and add corresponding terms to the system contributions.
    /**
      @param rLocalMatrix Local system matrix
      @param rLocalVector Local right hand side
      */
    void ApplyWallLaw(MatrixType& rLocalMatrix, VectorType& rLocalVector, const ProcessInfo& rCurrentProcessInfo)
    {
        if (RansCalculationUtilities::IsWallFunctionActive(*this))
        {
            const GeometryType& r_geometry = this->GetGeometry();
            // Get Shape function data
            Vector gauss_weights;
            Matrix shape_functions;
            RansCalculationUtilities::CalculateConditionGeometryData(
                r_geometry, this->GetIntegrationMethod(), gauss_weights, shape_functions);
            const IndexType num_gauss_points = gauss_weights.size();

            const double c_mu_25 =
                std::pow(rCurrentProcessInfo[TURBULENCE_RANS_C_MU], 0.25);
            const double kappa = rCurrentProcessInfo[WALL_VON_KARMAN];
            const double inv_kappa = 1.0 / kappa;
            const double beta = rCurrentProcessInfo[WALL_SMOOTHNESS_BETA];
            const double y_plus_limit =
                rCurrentProcessInfo[RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT];

            const double eps = std::numeric_limits<double>::epsilon();

            for (size_t g = 0; g < num_gauss_points; ++g)
            {
                const Vector& gauss_shape_functions = row(shape_functions, g);

                const array_1d<double, 3>& r_wall_velocity =
                    RansCalculationUtilities::EvaluateInPoint(
                        r_geometry, VELOCITY, gauss_shape_functions);
                const double wall_velocity_magnitude = norm_2(r_wall_velocity);

                const double tke = RansCalculationUtilities::EvaluateInPoint(
                    r_geometry, TURBULENT_KINETIC_ENERGY, gauss_shape_functions);
                const double rho = RansCalculationUtilities::EvaluateInPoint(
                    r_geometry, DENSITY, gauss_shape_functions);
                const double nu = RansCalculationUtilities::EvaluateInPoint(
                    r_geometry, KINEMATIC_VISCOSITY, gauss_shape_functions);

                double y_plus{0.0}, u_tau{0.0};
                RansCalculationUtilities::CalculateYPlusAndUtau(
                    y_plus, u_tau, wall_velocity_magnitude, mWallHeight, nu, kappa, beta);
                y_plus = std::max(y_plus, y_plus_limit);

                u_tau = std::max(c_mu_25 * std::sqrt(std::max(tke, 0.0)),
                                 wall_velocity_magnitude /
                                     (inv_kappa * std::log(y_plus) + beta));

                if (wall_velocity_magnitude > eps)
                {
                    const double value = rho * std::pow(u_tau, 2) *
                                         gauss_weights[g] / wall_velocity_magnitude;
                    for (size_t a = 0; a < r_geometry.PointsNumber(); ++a)
                    {
                        for (size_t dim = 0; dim < TDim; ++dim)
                        {
                            for (size_t b = 0; b < r_geometry.PointsNumber(); ++b)
                            {
                                rLocalMatrix(a * TDim + dim, b * TDim + dim) +=
                                    gauss_shape_functions[a] *
                                    gauss_shape_functions[b] * value;
                            }
                            rLocalVector[a * TDim + dim] -=
                                gauss_shape_functions[a] * value * r_wall_velocity[dim];
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
    void ApplyNeumannCondition(MatrixType& rLocalMatrix, VectorType& rLocalVector);

    ///@}

private:
    ///@name Member Variables
    ///@{

    double mWallHeight;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
    }

    ///@}

}; // Class FractionalStepKBasedWallCondition

///@}
///@name Input and output
///@{

/// input stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::istream& operator>>(std::istream& rIStream,
                                FractionalStepKBasedWallCondition<TDim, TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::ostream& operator<<(std::ostream& rOStream,
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
