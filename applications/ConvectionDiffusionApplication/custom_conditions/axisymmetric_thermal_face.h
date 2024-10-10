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
#include "thermal_face.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @brief Axisymmetric thermal face condition
 * Extension of the base ThermalFace class to be used in axisymmetric problems
 */
class KRATOS_API(CONVECTION_DIFFUSION_APPLICATION) AxisymmetricThermalFace : public ThermalFace
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AxisymmetricThermalFace
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(AxisymmetricThermalFace);

    /// Base thermal face condition type
    using BaseType = ThermalFace;

    /// Geometry type
    using GeometryType = typename BaseType::GeometryType;

    /// Properties type
    using PropertiesType = typename BaseType::PropertiesType;

    /// Index type
    using IndexType = typename BaseType::IndexType;

    /// Size type
    using SizeType = typename BaseType::SizeType;

    /// Vector type
    using VectorType = typename BaseType::VectorType;

    /// Matrix type
    using MatrixType = typename BaseType::MatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    AxisymmetricThermalFace(
        IndexType NewId,
        typename GeometryType::Pointer pGeometry);

    AxisymmetricThermalFace(
        IndexType NewId,
        typename GeometryType::Pointer pGeometry,
        typename PropertiesType::Pointer pProperties);

    /// Destructor.
    ~AxisymmetricThermalFace() override = default;

    /// Copy constructor.
    AxisymmetricThermalFace(AxisymmetricThermalFace const &rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    AxisymmetricThermalFace &operator=(AxisymmetricThermalFace const &rOther) = delete;

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
    AxisymmetricThermalFace();

    ///@}
    ///@name Protected Operations
    ///@{

    void SetIntegrationWeight(
        const IndexType IntegrationPointIndex,
        const typename GeometryType::IntegrationPointsArrayType &rIntegrationPoints,
        const Vector &rJacobianDeterminantsVector,
        ConditionDataStruct &rData) override;

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


    ///@}
}; // Class AxisymmetricThermalFace

///@}

///@} addtogroup block

}  // namespace Kratos.
