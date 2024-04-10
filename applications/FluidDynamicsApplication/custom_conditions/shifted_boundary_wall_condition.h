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

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) ShiftedBoundaryWallCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ShiftedBoundaryWallCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ShiftedBoundaryWallCondition);

    /*typedef NavierStokesWallCondition<TDim, TNumNodes, TWallModel...> BaseType;

    using BaseType::LocalSize;

    typedef typename BaseType::ConditionDataStruct ConditionDataStruct;

    typedef Node NodeType;

    typedef Properties PropertiesType;

    typedef Geometry<NodeType> GeometryType;

    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;*/

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

    //void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

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

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

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
     * This function computes the penalty coefficient for the Nitsche normal imposition (penalization and stabilization)
     * @param rN the current Gauss pt. shape functions vector
     * @param TODO
     * @return double The normal penalty coefficient value
     */
    double ComputeSlipNormalPenaltyCoefficient(
        const Vector& rN,
        const double DeltaTime,
        const double Penalty,
        const double ParentSize) const;

    /**
     * This function computes the penalty coefficients for the Nitsche tangential imposition
     * @param TODO
     * @return a pair of double containing the two coefficients
     */
    std::pair<const double, const double> ComputeSlipTangentialPenaltyCoefficients(
        const double SlipLength,
        const double Penalty,
        const double ParentSize) const;

    /**
     * This function computes the Nitsche coefficients for the Nitsche tangential imposition
     * @param TODO
     * @return a pair of double containing the two coefficients
     */
    std::pair<const double, const double> ComputeSlipTangentialNitscheCoefficients(
        const double SlipLength,
        const double Penalty,
        const double ParentSize) const;

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
inline std::istream& operator >> (std::istream& rIStream, ShiftedBoundaryWallCondition& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const ShiftedBoundaryWallCondition& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
///@} addtogroup block
}  // namespace Kratos.
