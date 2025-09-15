//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "elements/mesh_element.h"

namespace Kratos
{
MeshElement::MeshElement(IndexType NewId)
    : BaseType(NewId)
{
}

/***********************************************************************************/
/***********************************************************************************/

MeshElement::MeshElement(
    IndexType NewId,
    const NodesArrayType& rThisNodes
    ) : BaseType(NewId, rThisNodes)
{
}

/***********************************************************************************/
/***********************************************************************************/

MeshElement::MeshElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry
    ) : BaseType(NewId, pGeometry)
{
}

/***********************************************************************************/
/***********************************************************************************/

MeshElement::MeshElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties
    ) : BaseType(NewId,pGeometry, pProperties)
{
}

/***********************************************************************************/
/***********************************************************************************/

MeshElement::MeshElement(MeshElement const& rOther)
    : BaseType(rOther)
{
}

/***********************************************************************************/
/***********************************************************************************/

MeshElement::~MeshElement()
{
}

/***********************************************************************************/
/***********************************************************************************/

MeshElement& MeshElement::operator=(MeshElement const& rOther)
{
    //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

    Element::operator=(rOther);

    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer MeshElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<MeshElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer MeshElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<MeshElement>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer MeshElement::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Element::Pointer p_new_elem = Kratos::make_intrusive<MeshElement>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));
    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void MeshElement::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    const Variable<double >& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

void MeshElement::AddExplicitContribution(
    const VectorType& rRHS,
    const Variable<VectorType>& rRHSVariable,
    const Variable<array_1d<double,3> >& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

void MeshElement::AddExplicitContribution(
    const MatrixType& rLHSMatrix,
    const Variable<MatrixType>& rLHSVariable,
    const Variable<Matrix>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

void MeshElement::CalculateOnIntegrationPoints(
    const Variable<bool>& rVariable,
    std::vector<bool>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_geometry = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(r_geometry.GetDefaultIntegrationMethod());
    const SizeType integration_points_number = r_integration_points.size();
    if (rOutput.size() != integration_points_number) {
        rOutput.resize(integration_points_number, rVariable.Zero());
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MeshElement::CalculateOnIntegrationPoints(
    const Variable<int>& rVariable,
    std::vector<int>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_geometry = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(r_geometry.GetDefaultIntegrationMethod());
    const SizeType integration_points_number = r_integration_points.size();
    if (rOutput.size() != integration_points_number) {
        rOutput.resize(integration_points_number, rVariable.Zero());
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MeshElement::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_geometry = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(r_geometry.GetDefaultIntegrationMethod());
    const SizeType integration_points_number = r_integration_points.size();
    if (rOutput.size() != integration_points_number) {
        rOutput.resize(integration_points_number, rVariable.Zero());
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MeshElement::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector< array_1d<double, 3 > >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_geometry = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(r_geometry.GetDefaultIntegrationMethod());
    const SizeType integration_points_number = r_integration_points.size();
    if (rOutput.size() != integration_points_number) {
        rOutput.resize(integration_points_number, rVariable.Zero());
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MeshElement::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 4 > >& rVariable,
    std::vector< array_1d<double, 4 > >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_geometry = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(r_geometry.GetDefaultIntegrationMethod());
    const SizeType integration_points_number = r_integration_points.size();
    if (rOutput.size() != integration_points_number) {
        rOutput.resize(integration_points_number, rVariable.Zero());
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MeshElement::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 6 > >& rVariable,
    std::vector< array_1d<double, 6 > >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_geometry = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(r_geometry.GetDefaultIntegrationMethod());
    const SizeType integration_points_number = r_integration_points.size();
    if (rOutput.size() != integration_points_number) {
        rOutput.resize(integration_points_number, rVariable.Zero());
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MeshElement::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 9 > >& rVariable,
    std::vector< array_1d<double, 9 > >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_geometry = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(r_geometry.GetDefaultIntegrationMethod());
    const SizeType integration_points_number = r_integration_points.size();
    if (rOutput.size() != integration_points_number) {
        rOutput.resize(integration_points_number, rVariable.Zero());
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MeshElement::CalculateOnIntegrationPoints(
    const Variable<Vector >& rVariable,
    std::vector< Vector >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_geometry = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(r_geometry.GetDefaultIntegrationMethod());
    const SizeType integration_points_number = r_integration_points.size();
    if (rOutput.size() != integration_points_number) {
        rOutput.resize(integration_points_number, rVariable.Zero());
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MeshElement::CalculateOnIntegrationPoints(
    const Variable<Matrix >& rVariable,
    std::vector< Matrix >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_geometry = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(r_geometry.GetDefaultIntegrationMethod());
    const SizeType integration_points_number = r_integration_points.size();
    if (rOutput.size() != integration_points_number) {
        rOutput.resize(integration_points_number, rVariable.Zero());
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MeshElement::CalculateOnIntegrationPoints(
    const Variable<ConstitutiveLaw::Pointer>& rVariable,
    std::vector<ConstitutiveLaw::Pointer>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_geometry = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(r_geometry.GetDefaultIntegrationMethod());
    const SizeType integration_points_number = r_integration_points.size();
    if (rOutput.size() != integration_points_number) {
        rOutput.resize(integration_points_number, nullptr);
    }
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters MeshElement::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : [],
        "framework"                  : "lagrangian",
        "symmetric_lhs"              : true,
        "positive_definite_lhs"      : true,
        "output"                     : {
            "gauss_point"            : [],
            "nodal_historical"       : [],
            "nodal_non_historical"   : [],
            "entity"                 : []
        },
        "required_variables"         : [],
        "required_dofs"              : [],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Line2D2", "Triangle2D3", "Triangle2D6", "Triangle2D10", "Triangle2D15", "Quadrilateral2D4", "Quadrilateral2D8", "Quadrilateral2D9", "Line3D2", "Triangle3D3", "Tetrahedra3D4", "Prism3D6", "Prism3D15", "Hexahedra3D8", "Hexahedra3D20", "Hexahedra3D27", "Tetrahedra3D10"],
        "element_integrates_in_time" : true,
        "compatible_constitutive_laws": {
            "type"        : [],
            "dimension"   : [],
            "strain_size" : []
        },
        "required_polynomial_degree_of_geometry" : -1,
        "documentation"   : "This is a pure geometric element, no computation"
    })");
    return specifications;
}

/***********************************************************************************/
/***********************************************************************************/

void MeshElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
}

void MeshElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
}

} // Namespace Kratos


