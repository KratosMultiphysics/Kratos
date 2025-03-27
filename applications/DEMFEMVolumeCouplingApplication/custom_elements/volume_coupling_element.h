// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    





#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/solid_elements/small_displacement.h"


namespace Kratos
{

class KRATOS_API(DEMFEM_VOLUME_COUPLING_APPLICATION) VolumeCouplingElement
    : public SmallDisplacement
{
public:
    // Type Definitions
    typedef Element::Pointer Pointer;

    // Constructors
    VolumeCouplingElement(IndexType NewId, GeometryType::Pointer pGeometry);
    VolumeCouplingElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);
    VolumeCouplingElement(); // Default constructor needed for serialization

    // Create and Clone methods
    Element::Pointer Create(IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties) const override;
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const override;

    // Destructor
    virtual ~VolumeCouplingElement();

    // Other Public Methods
    virtual double GetIntegrationWeight(
        const GeometryType::IntegrationPointsArrayType& rThisIntegrationPoints,
        const IndexType PointNumber,
        const double detJ
        ) const override;

    void CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    ) override;
    
    // Other declarations can be added here as needed

};

} // namespace Kratos.


// #pragma once

// // System includes

// // External includes

// // Project includes
// #include "includes/define.h"
// // #include "includes/element.h"
// // #include "utilities/integration_utilities.h"
// // #include "structural_mechanics_application_variables.h"
// // #include "utilities/geometrical_sensitivity_utility.h"
// // #include "custom_utilities/structural_mechanics_element_utilities.h"
// #include "custom_elements/small_displacement.h"

// namespace Kratos
// {

// class KRATOS_API(DEMFEM_VOLUME_COUPLING_APPLICATION) VolumeCouplingElement
//     : public SmallDisplacement
// {


// public:
// VolumeCouplingElement(IndexType NewId, GeometryType::Pointer pGeometry);
// VolumeCouplingElement();// Default constructor needed for serialization
   
//     virtual double GetIntegrationWeight(
//         const GeometryType::IntegrationPointsArrayType& rThisIntegrationPoints,
//         const IndexType PointNumber,
//         const double detJ
//         ) const override;


// };
// } // namespace Kratos.
