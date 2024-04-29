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
#include <string>
#include <iostream>

// External includes

// Project includes
#include "geometries/geometry_data.h"
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/model_part.h"
#include "includes/serializer.h"
#include "includes/process_info.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "navier_stokes_wall_condition.h"

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
 * @brief Implements a wall condition for the Navier-Stokes (and Stokes) monolithic formulations
 * This condition is intended to be used in combination with Navier-Stokes (or Stokes) P2-P1 (continuous
 * pressure monolithic formulations. It supports the Neumann BC contribution as well as the addition of
 * a wall law contribution through the TWallModel template argument. Such TWallModel must be a class
 * implementing the wall model RHS and LHS Gauss point contributions (as example see @NavierSlipWallLaw).
 * Current condition also has optional features that help numerical stability such as the outlet
 * inflow energy correction or the spurious tangential velocity correction for pure slip boundaries.
 * TODO: Implement the drag calculation (this requires implementing the viscous stress from the parent)
 * TODO: Implement the slip boundaries spurious velocity correction (this requires implementing the viscous stress from the parent)
 * TODO: Implement the wall contributions (this requires Lobatto quadratures and also I/O considerations if current ApplyWallLawProcess is used)
 * @tparam TDim Number of dimensions
 * @tparam TWallModel Optional class implementing a LHS and RHS wall contribution
 */
template<unsigned int TDim, class... TWallModel>
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) NavierStokesP2P1ContinuousWallCondition : public Condition
{
public:

    static_assert(sizeof...(TWallModel) == 0, "Wall models are not supported in 'NavierStokesP2P1ContinuousWallCondition' yet.");

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(NavierStokesP2P1ContinuousWallCondition);

    static constexpr std::size_t VoigtSize = 3*(TDim-1);

    static constexpr std::size_t VelocityNumNodes = TDim == 2 ? 3 : 6;

    static constexpr std::size_t PressureNumNodes = TDim == 2 ? 2 : 3;

    static constexpr std::size_t LocalSize = VelocityNumNodes*TDim + PressureNumNodes;

    static constexpr GeometryData::IntegrationMethod IntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_3;

    using BaseType = Condition;

    using SizeType = typename BaseType::SizeType;

    using IndexType = typename BaseType::IndexType;

    using GeometryType = typename BaseType::GeometryType;

    using NodesArrayType = typename BaseType::NodesArrayType;

    using VectorType = typename BaseType::VectorType;

    using MatrixType = typename BaseType::MatrixType;

    using EquationIdVectorType = typename BaseType::EquationIdVectorType;

    using DofsVectorType = typename BaseType::DofsVectorType;

    struct ConditionDataStruct
    {
        double Weight;                                  // Gauss point weight
        array_1d<double, TDim> UnitNormal;              // Condition unit normal
        array_1d<double, VelocityNumNodes> N_v;         // Gauss point velocity shape functions values
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NavierStokesP2P1ContinuousWallCondition(IndexType NewId = 0)
        : BaseType(NewId)
    {
    }

    /// Constructor using an array of nodes
    NavierStokesP2P1ContinuousWallCondition(
        IndexType NewId,
        const NodesArrayType& ThisNodes)
        : BaseType(NewId, ThisNodes)
    {
    }

    /// Constructor using Geometry
    NavierStokesP2P1ContinuousWallCondition(
        IndexType NewId,
        typename GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {
    }

    /// Constructor using Properties
    NavierStokesP2P1ContinuousWallCondition(
        IndexType NewId,
        typename GeometryType::Pointer pGeometry,
        Properties::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /// Copy constructor.
    NavierStokesP2P1ContinuousWallCondition(NavierStokesP2P1ContinuousWallCondition const& rOther)
        : BaseType(rOther)
    {
    }

    /// Destructor.
    ~NavierStokesP2P1ContinuousWallCondition() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    NavierStokesP2P1ContinuousWallCondition& operator=(NavierStokesP2P1ContinuousWallCondition const& rOther)
    {
        Condition::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        typename Properties::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<NavierStokesP2P1ContinuousWallCondition>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    Condition::Pointer Create(
        IndexType NewId,
        typename GeometryType::Pointer pGeom,
        typename Properties::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<NavierStokesP2P1ContinuousWallCondition>(NewId, pGeom, pProperties);
    }

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(
        DofsVectorType& rConditionDofList,
        const ProcessInfo& rCurrentProcessInfo) const override;

    void Calculate(
        const Variable< array_1d<double,3>>& rVariable,
        array_1d<double,3>& Output,
        const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "NavierStokesP2P1ContinuousWallCondition" << TDim << "D";
        return buffer.str();
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NavierStokesP2P1ContinuousWallCondition";
    }

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

    /**
     * @brief Calculates the Gauss point RHS contribution
     * This method calculates the current Gauss point RHS contribution and saves it
     * in the provided array. Note that the input data container is expected to
     * already contain the data at the Gauss point of interest.
     * @param rLHS Reference to the RHS output vector
     * @param rData Condition data container
     * @param rProcessInfo Reference to the ProcessInfo container
     */
    void AddGaussPointRHSContribution(
        VectorType& rRHS,
        const ConditionDataStruct& rData,
        const ProcessInfo& rProcessInfo);

    /**
     * @brief Calculates the RHS Neumann BC contribution
     * This method calculates the Neumann BC pressure flux contribution
     * Note that the Neumann BC value is expected to be stored in the historical
     * database within the EXTERNAL_PRESSURE variable.
     * @param rRHS Reference to the RHS output vector
     * @param data Condition data container
     */
    void ComputeRHSNeumannContribution(
        VectorType& rRHS,
        const ConditionDataStruct& data);

    /**
     * @brief Calculates and adds the RHS outlet inflow prevention contribution
     * This method calculates and adds an extra numerical contribution to the RHS in order
     * to prevent uncontrolled system energy growth coming from inflow in free-boundaries.
     * More information can be found in Dong et al. 2014 (https://doi.org/10.1016/j.jcp.2013.12.042).
     * @param rRHS Reference to RHS vector
     * @param rData Condition data container
     * @param rProcessInfo Reference to the ProcessInfo container
     */
    void ComputeRHSOutletInflowContribution(
        VectorType& rRHS,
        const ConditionDataStruct& rData,
        const ProcessInfo& rProcessInfo);

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

    /**
     * @brief Calculate the condition unit normal
     * This method calculates the current condition unit normal
     * @param rUnitNormal Reference to the current condition unit normal
     */
    void CalculateUnitNormal(array_1d<double, TDim>& rUnitNormal);

    template<typename TWallModelType>
    int WallModelCheckCall(const ProcessInfo& rProcessInfo) const
    {
        return TWallModelType::Check(this, rProcessInfo);
    }

    template<typename TWallModelType>
    void AddWallModelRightHandSideCall(
        VectorType& rRHS,
        const ProcessInfo& rProcessInfo)
    {
        TWallModelType::AddWallModelRightHandSide(rRHS, this, rProcessInfo);
    }

    template<typename TWallModelType>
    void AddWallModelLeftHandSideCall(
        MatrixType& rLHS,
        const ProcessInfo& rProcessInfo)
    {
        TWallModelType::AddWallModelLeftHandSide(rLHS, this, rProcessInfo);
    }

    template<typename TWallModelType>
    void AddWallModelLocalSystemCall(
        MatrixType& rLHS,
        VectorType& rRHS,
        const ProcessInfo& rProcessInfo)
    {
        TWallModelType::AddWallModelLocalSystem(rLHS, rRHS, this, rProcessInfo);
    }

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
}; // Class NavierStokesP2P1ContinuousWallCondition

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim, unsigned int TNumNodes, class TWallModel >
inline std::istream& operator >> (std::istream& rIStream, NavierStokesP2P1ContinuousWallCondition<TDim,TWallModel>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim, unsigned int TNumNodes, class TWallModel >
inline std::ostream& operator << (std::ostream& rOStream, const NavierStokesP2P1ContinuousWallCondition<TDim,TWallModel>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block


}  // namespace Kratos.
