//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#ifndef KRATOS_RANS_EVM_EPSILON_WALL_CONDITION_H
#define KRATOS_RANS_EVM_EPSILON_WALL_CONDITION_H

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/condition.h"
#include "includes/define.h"
#include "includes/process_info.h"
#include "includes/serializer.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "includes/cfd_variables.h"
#include "rans_modelling_application_variables.h"

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
 * @brief Epsilon Neumann wall condition
 *
 * This is a Neumann wall condition for $\epsilon$ equation in $k-\epsilon$ formulation of RANS
 * based on eddy viscosity model formulation.
 *
 * @tparam TNumNodes Number of nodes in the wall condition
 */

template <unsigned int TNumNodes>
class RansEvmEpsilonWallCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansEvmEpsilonWallCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(RansEvmEpsilonWallCondition);

    typedef Node<3> NodeType;

    typedef Properties PropertiesType;

    typedef Geometry<NodeType> GeometryType;

    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector<Dof<double>::Pointer> DofsVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
      @param NewId Index for the new condition
      */
    explicit RansEvmEpsilonWallCondition(IndexType NewId = 0) : Condition(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    RansEvmEpsilonWallCondition(IndexType NewId, const NodesArrayType& ThisNodes)
        : Condition(NewId, ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    RansEvmEpsilonWallCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    RansEvmEpsilonWallCondition(IndexType NewId,
                            GeometryType::Pointer pGeometry,
                            PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {
    }

    /// Copy constructor.
    RansEvmEpsilonWallCondition(RansEvmEpsilonWallCondition const& rOther)
        : Condition(rOther)
    {
    }

    /// Destructor.
    ~RansEvmEpsilonWallCondition() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    RansEvmEpsilonWallCondition& operator=(RansEvmEpsilonWallCondition const& rOther)
    {
        Condition::operator=(rOther);

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create a new RansEvmEpsilonWallCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<RansEvmEpsilonWallCondition>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    Condition::Pointer Create(IndexType NewId,
                              GeometryType::Pointer pGeom,
                              PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<RansEvmEpsilonWallCondition>(NewId, pGeom, pProperties);
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

    /// Return local contributions of the correct size, filled with zeros (for compatibility with time schemes).
    /** The actual local contributions are computed in the Damping functions
      @see CalculateLocalVelocityContribution
      */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override
    {
        if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes)
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes);

        if (rRightHandSideVector.size() != TNumNodes)
            rRightHandSideVector.resize(TNumNodes);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);
        noalias(rRightHandSideVector) = ZeroVector(TNumNodes);
    }

    /// Return a matrix of the correct size, filled with zeros (for compatibility with time schemes).
    /** The actual local contributions are computed in the Damping functions
      @see DampingMatrix
      */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               ProcessInfo& rCurrentProcessInfo) override
    {
        if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes)
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);
    }

    /// Return local right hand side of the correct size, filled with zeros (for compatibility with time schemes).
    /** The actual local contributions are computed in the Damping functions
      @see CalculateLocalVelocityContribution
      */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo) override
    {
        if (rRightHandSideVector.size() != TNumNodes)
            rRightHandSideVector.resize(TNumNodes);

        noalias(rRightHandSideVector) = ZeroVector(TNumNodes);
    }

    void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) override
    {
        VectorType RHS;
        this->CalculateLocalVelocityContribution(rDampingMatrix, RHS, rCurrentProcessInfo);
    }

    /// Calculate wall stress term for all nodes with STRUCTURE == true
    /**
      @param rDampingMatrix Left-hand side matrix
      @param rRightHandSideVector Right-hand side vector
      @param rCurrentProcessInfo ProcessInfo instance (unused)
      */
    void CalculateLocalVelocityContribution(MatrixType& rDampingMatrix,
                                            VectorType& rRightHandSideVector,
                                            ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rRightHandSideVector.size() != TNumNodes)
            rRightHandSideVector.resize(TNumNodes);

        if (rDampingMatrix.size1() != TNumNodes || rDampingMatrix.size2() != TNumNodes)
            rDampingMatrix.resize(TNumNodes, TNumNodes);

        rRightHandSideVector.clear();
        rDampingMatrix.clear();

        const GeometryType& r_geometry = this->GetGeometry();

        if (!this->Is(STRUCTURE))
            return;

        for (size_t i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node)
            if (!r_geometry[i_node].GetDof(TURBULENT_ENERGY_DISSIPATION_RATE).IsFree())
                return;

        RansCalculationUtilities rans_calculation_utilities;

        // Get Shape function data
        const GeometryType::IntegrationPointsArrayType& integration_points =
            r_geometry.IntegrationPoints(GeometryData::GI_GAUSS_2);
        const std::size_t num_gauss_points = integration_points.size();
        MatrixType shape_functions =
            r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

        const double area = r_geometry.DomainSize();

        // CAUTION: "Jacobian" is 2.0*A for triangles but 0.5*A for lines
        double J = (TNumNodes == 2) ? 0.5 * area : 2.0 * area;

        const double epsilon_sigma =
            rCurrentProcessInfo[TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA];
        const double c_mu_25 = std::pow(rCurrentProcessInfo[TURBULENCE_RANS_C_MU], 0.25);
        const double eps = std::numeric_limits<double>::epsilon();
        for (unsigned int g = 0; g < num_gauss_points; g++)
        {
            const Vector& gauss_shape_functions = row(shape_functions, g);
            const double weight = J * integration_points[g].Weight();

            const double nu = rans_calculation_utilities.EvaluateInPoint(
                r_geometry, KINEMATIC_VISCOSITY, gauss_shape_functions);
            const double nu_t = rans_calculation_utilities.EvaluateInPoint(
                r_geometry, TURBULENT_VISCOSITY, gauss_shape_functions);
            const double tke = rans_calculation_utilities.EvaluateInPoint(
                r_geometry, TURBULENT_KINETIC_ENERGY, gauss_shape_functions);
            const double epsilon = rans_calculation_utilities.EvaluateInPoint(
                r_geometry, TURBULENT_ENERGY_DISSIPATION_RATE, gauss_shape_functions);
            const double y_plus = rans_calculation_utilities.EvaluateInPoint(
                r_geometry, RANS_Y_PLUS, gauss_shape_functions);

            const double u_tau = c_mu_25 * std::sqrt(std::max(tke, 0.0));

            if (y_plus > eps)
            {
                const double value =
                    weight * (nu + nu_t / epsilon_sigma) * u_tau / (y_plus * nu);

                for (unsigned int a = 0; a < TNumNodes; ++a)
                {
                    for (unsigned int b = 0; b < TNumNodes; ++b)
                    {
                        rDampingMatrix(a, b) -= gauss_shape_functions[a] *
                                                gauss_shape_functions[b] * value;
                    }
                    rRightHandSideVector[a] += value * gauss_shape_functions[a] * epsilon;
                }
            }
        }

        KRATOS_CATCH("");
    }

    /// Check that all data required by this condition is available and reasonable
    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        int Check = Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area > 0

        KRATOS_CHECK_VARIABLE_KEY(TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA);
        KRATOS_CHECK_VARIABLE_KEY(TURBULENCE_RANS_C_MU);
        KRATOS_CHECK_VARIABLE_KEY(KINEMATIC_VISCOSITY);
        KRATOS_CHECK_VARIABLE_KEY(TURBULENT_VISCOSITY);
        KRATOS_CHECK_VARIABLE_KEY(TURBULENT_KINETIC_ENERGY);
        KRATOS_CHECK_VARIABLE_KEY(TURBULENT_ENERGY_DISSIPATION_RATE);

        const GeometryType& r_geometry = this->GetGeometry();

        for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
        {
            const NodeType& r_node = r_geometry[i_node];

            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(KINEMATIC_VISCOSITY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_VISCOSITY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_Y_PLUS, r_node);

            KRATOS_CHECK_DOF_IN_NODE(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
        }

        return Check;

        KRATOS_CATCH("");
    }

    /// Provides the global indices for each one of this element's local rows.
    /** This determines the elemental equation ID vector for all elemental DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo the current process info object (unused)
     */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override
    {
        if (rResult.size() != TNumNodes)
            rResult.resize(TNumNodes, false);

        for (unsigned int i = 0; i < TNumNodes; i++)
            rResult[i] = Condition::GetGeometry()[i]
                             .GetDof(TURBULENT_ENERGY_DISSIPATION_RATE)
                             .EquationId();
    }

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(DofsVectorType& ConditionDofList, ProcessInfo& CurrentProcessInfo) override
    {
        if (ConditionDofList.size() != TNumNodes)
            ConditionDofList.resize(TNumNodes);

        for (unsigned int i = 0; i < TNumNodes; i++)
            ConditionDofList[i] =
                Condition::GetGeometry()[i].pGetDof(TURBULENT_ENERGY_DISSIPATION_RATE);
    }

    void GetValuesVector(VectorType& rValues, int Step = 0) override
    {
        this->GetFirstDerivativesVector(rValues, Step);
    }

    /// Returns TURBULENT_ENERGY_DISSIPATION_RATE for each node
    /**
     * @param Values Vector of nodal unknowns
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) override
    {
        if (rValues.size() != TNumNodes)
            rValues.resize(TNumNodes, false);

        GeometryType& rGeom = this->GetGeometry();
        IndexType LocalIndex = 0;
        for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
        {
            rValues[LocalIndex++] = rGeom[iNode].FastGetSolutionStepValue(
                TURBULENT_ENERGY_DISSIPATION_RATE, Step);
        }
    }

    /// Returns TURBULENT_ENERGY_DISSIPATION_RATE_2 0 for each node
    /**
     * @param Values Vector of nodal second derivatives
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) override
    {
        if (rValues.size() != TNumNodes)
            rValues.resize(TNumNodes, false);

        GeometryType& rGeom = this->GetGeometry();
        IndexType LocalIndex = 0;
        for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
        {
            rValues[LocalIndex++] = rGeom[iNode].FastGetSolutionStepValue(
                TURBULENT_ENERGY_DISSIPATION_RATE_2, Step);
        }
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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "RansEvmEpsilonWallCondition" << TNumNodes << "N";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "RansEvmEpsilonWallCondition";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

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
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

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

}; // Class RansEvmEpsilonWallCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <unsigned int TNumNodes>
inline std::istream& operator>>(std::istream& rIStream,
                                RansEvmEpsilonWallCondition<TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template <unsigned int TNumNodes>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansEvmEpsilonWallCondition<TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_EVM_EPSILON_WALL_CONDITION_H
