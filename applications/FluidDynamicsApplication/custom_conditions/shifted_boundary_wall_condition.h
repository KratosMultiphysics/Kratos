//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Franziska Wahl
//

#pragma once

// System includes
#include <cstddef>
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/properties.h"
#include "includes/condition.h"
#include "geometries/geometry.h"
#include "includes/variables.h"

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

template<std::size_t TDim>
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) ShiftedBoundaryWallCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ShiftedBoundaryWallCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ShiftedBoundaryWallCondition);

    static constexpr std::size_t VoigtSize = 3 * (TDim-1);  // 3 in 2D, 6 in 3D - StrainSize of FluidElementData
    static constexpr std::size_t BlockSize = TDim + 1;

    typedef Node NodeType;

    typedef Properties PropertiesType;

    typedef Geometry<NodeType> GeometryType;

    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    ShiftedBoundaryWallCondition(
        IndexType NewId,
        const NodesArrayType& ThisNodes)
        : Condition(NewId, ThisNodes) {}

    ShiftedBoundaryWallCondition(
        IndexType NewId,
        Geometry<Node>::Pointer pGeometry)
        : Condition(NewId, pGeometry) {}

    ShiftedBoundaryWallCondition(
        IndexType NewId,
        Geometry<Node>::Pointer pGeometry,
        Properties::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties) {}

    /// Copy constructor.
    ShiftedBoundaryWallCondition(ShiftedBoundaryWallCondition const& rOther) = delete;

    /// Destructor.
    ~ShiftedBoundaryWallCondition() override = default;


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    ShiftedBoundaryWallCondition& operator=(ShiftedBoundaryWallCondition const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<ShiftedBoundaryWallCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<ShiftedBoundaryWallCondition>(NewId, pGeom, pProperties);
    }

    /**
     * @brief This function calculates the LHS matrix and RHS vector of the local system.
     * All terms should be added which are necessary for the imposition of a shifted-boundary wall condition
     * for the given geometry of an integration point to the system.
     * The geometry should include all nodes that contribute to the calculation of the value at the integration point.
     * AddNitscheImposition is called in order to use Nitsche imposition of a Navier-slip boundary condition at the integration point.
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rRightHandSideVector reference to the RHS vector
     * @param rCurrentProcessInfo reference to the ProcessInfo
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This function calculates the LHS matrix of the local system.
     * Inside the function CalculateLocalSystem is called.
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rRightHandSideVector reference to the RHS vector
     * @param rCurrentProcessInfo reference to the ProcessInfo
     */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This function calculates the RHS vector of the local system.
     * Inside the function CalculateLocalSystem is called.
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rRightHandSideVector reference to the RHS vector
     * @param rCurrentProcessInfo reference to the ProcessInfo
     */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This function collects the Equation IDs of all DOFs of all nodes of the condition's geometry.
     * The geometry should include all nodes that contribute to the calculation of the value at the integration point.
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rRightHandSideVector reference to the RHS vector
     * @param rCurrentProcessInfo reference to the ProcessInfo (unused)
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * @brief This function collects the DOFs of all nodes of the condition's geometry.
     * The geometry should include all nodes that contribute to the calculation of the value at the integration point.
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rRightHandSideVector reference to the RHS vector
     * @param rCurrentProcessInfo reference to the ProcessInfo (unused)
     */
    void GetDofList(
        DofsVectorType& ConditionalDofList,
        const ProcessInfo& CurrentProcessInfo) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

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
        buffer << "ShiftedBoundaryWallCondition #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ShiftedBoundaryWallCondition #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "ShiftedBoundaryWallCondition #" << Id() << std::endl;
        this->GetGeometry().PrintData(rOStream);
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

    /**
     * @brief This function adds the terms for imposing a Navier-slip boundary condition at an integration point to the system.
     * NOTE that the condition does NOT account for mesh motion so far.
     * The stabilized Nitsche imposition of the Navier-slip boundary condition (Robin-type BC) consists of a
     * no penetration constraint in wall normal direction and a shear force imposition in tangential direction.
     * It behaves as a linear wall-law using the slip length parameter epsilon (no-slip for epsilon towards zero,
     * slip for epsilon towards infinity).
     * Reference: M. Winter, B. Schott, A. Massing, W. Wall, A nitsche cut finite element method for the oseen problem with general
     * navier boundary conditions, Comput. Methods Appl. Mech. Engrg. 330 (2018) 220–252, http://dx.doi.org/10.1016/j.cma.2017.10.023.
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rRightHandSideVector reference to the RHS vector
     * @param rCurrentProcessInfo reference to the ProcessInfo
     */
    void AddNitscheImposition(
        MatrixType& rLHS,
        VectorType& rRHS,
        const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief This function penalizes a violation of a Dirichlet boundary condition at an integration point using a penalty constant.
     * PENALTY_COEFFICIENT is taken as penalty constant from rCurrentProcessInfo.
     * Penalization is added to LHS and RHS as violation of a zero velocity constraint.
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rRightHandSideVector reference to the RHS vector
     * @param rCurrentProcessInfo reference to the ProcessInfo
     */
    void AddDirichletPenalization(
        MatrixType& rLHS,
        VectorType& rRHS,
        const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief This function builds the strain matrix from the shape function derivatives utilizing Voigt notation.
     * @param rDN_DX matrix of shape function derivatives at all cloud points
     * @param NumNodes number of nodes of the geometry (cloud points)
     * @param rB B_matrix/ strain matrix
     */
    void CalculateStrainMatrix(
        const Matrix& rDN_DX,
        const std::size_t NumNodes,
        Matrix& rB);

    /**
     * This function computes the penalty coefficient for the Nitsche normal imposition (penalization and stabilization)
     * @param rN the current Gauss pt. shape functions vector
     * @param DeltaTime time step
     * @param Gamma Nitsche penalty coefficient (gamma)
     * @param ParentSize size/ volume of the parent element
     * @param EffectiveViscosity effective viscosity
     * @return double The normal penalty coefficient value
     */
    double ComputeSlipNormalPenaltyCoefficient(
        const Vector& rN,
        const double DeltaTime,
        const double Gamma,
        const double ParentSize,
        const double EffectiveViscosity) const;

    /**
     * This function computes the penalty coefficients for the Nitsche tangential imposition
     * @param SlipLength slip length for Navier-slip (zero for no-slip)
     * @param Gamma Nitsche penalty coefficient (gamma)
     * @param GammaShear Nitsche penalty coefficient for slip (gamma)
     * @param ParentSize size/ volume of the parent element
     * @param EffectiveViscosity effective viscosity
     * @return a pair of double containing the two coefficients
     */
    std::pair<const double, const double> ComputeSlipTangentialPenaltyCoefficients(
        const double SlipLength,
        const double Gamma,
        const double GammaShear,
        const double ParentSize,
        const double EffectiveViscosity) const;

    /**
     * This function computes the Nitsche coefficients for the Nitsche tangential imposition
     * @param SlipLength slip length for Navier-slip (zero for no-slip)
     * @param GammaShear Nitsche penalty coefficient (gamma)
     * @param CharactLength Characteristic length of the problem
     * @param EffectiveViscosity effective viscosity
     * @return a pair of double containing the two coefficients
     */
    std::pair<const double, const double> ComputeSlipTangentialNitscheCoefficients(
        const double SlipLength,
        const double GammaShear,
        const double CharactLength,
        const double EffectiveViscosity) const;

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{

    // Internal default constructor for serialization
    ShiftedBoundaryWallCondition() = default;

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
}; // Class ShiftedBoundaryWallCondition

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template< std::size_t TDim>
inline std::istream& operator >> (std::istream& rIStream, ShiftedBoundaryWallCondition<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template< std::size_t TDim>
inline std::ostream& operator << (std::ostream& rOStream, const ShiftedBoundaryWallCondition<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
///@} addtogroup block
}  // namespace Kratos.
