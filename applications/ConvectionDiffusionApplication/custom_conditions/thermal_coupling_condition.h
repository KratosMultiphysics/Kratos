// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined(KRATOS_THERMAL_COUPLING_CONDITION_H_INCLUDED)
#define  KRATOS_THERMAL_COUPLING_CONDITION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/condition.h"
#include "includes/define.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/properties.h"
#include "includes/variables.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @brief Thermal coupling condition
 * This class implements a thermal coupling condition between two bodies
 */
template<std::size_t TDim, std::size_t TNumNodes>
class KRATOS_API(CONVECTION_DIFFUSION_APPLICATION) ThermalCouplingCondition: public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ThermalCouplingCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ThermalCouplingCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    ThermalCouplingCondition(
        IndexType NewId,
        Geometry< Node<3> >::Pointer pGeometry);

    ThermalCouplingCondition(
        IndexType NewId,
        Geometry< Node<3> >::Pointer pGeometry,
        Properties::Pointer pProperties);

    /// Destructor.
    ~ThermalCouplingCondition() override;

    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        Properties::Pointer pProperties) const override;

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        Properties::Pointer pProperties) const override;

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
        DofsVectorType& rConditionDofList,
        const ProcessInfo& CurrentProcessInfo) const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}

protected:

    ///@name Protected Life Cycle
    ///@{

    // Internal default constructor for serialization
    ThermalCouplingCondition();

    ///@}
    ///@name Protected Operations
    ///@{

    double CalculateNodalIntegrationWeight() const;

    void FillLeftHandSideMatrix(
        const double LeftHandSideCoefficient,
        MatrixType& rLeftHandSideMatrix) const;

    void FillRightHandSideVector(
        const double RightHandSideCoefficient,
        const Variable<double>& rUnknownVariable,
        VectorType& rRightHandSideVector) const;

    ///@}
private:
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ThermalCouplingCondition& operator=(ThermalCouplingCondition const& rOther);

    /// Copy constructor.
    ThermalCouplingCondition(ThermalCouplingCondition const& rOther);

    ///@}
}; // Class ThermalCouplingCondition

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_THERMAL_COUPLING_CONDITION_H_INCLUDED  defined
