//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Riccardo Tosi
//

#include "symbolic_eulerian_convection_diffusion_explicit.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

SymbolicEulerianConvectionDiffusionExplicit::SymbolicEulerianConvectionDiffusionExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

SymbolicEulerianConvectionDiffusionExplicit::SymbolicEulerianConvectionDiffusionExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

SymbolicEulerianConvectionDiffusionExplicit::~SymbolicEulerianConvectionDiffusionExplicit() {}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

Element::Pointer SymbolicEulerianConvectionDiffusionExplicit::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicEulerianConvectionDiffusionExplicit>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}


Element::Pointer SymbolicEulerianConvectionDiffusionExplicit::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicEulerianConvectionDiffusionExplicit>(NewId, pGeom, pProperties);
}

void SymbolicEulerianConvectionDiffusionExplicit::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const unsigned int LocalSize = r_geometry.size();
    const unsigned int Dimension = r_geometry.WorkingSpaceDimension();

    // Resize and intialize output
    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize, false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    // Element variables
    ElementVariables rVariables;
    this->InitializeEulerianElement(rVariables,rCurrentProcessInfo);

    // Reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    // Define local variables
    Element::GeometryType::JacobiansType J0;
    Matrix DN(LocalSize,Dimension);
    Matrix InvJ0(Dimension,Dimension);
    Vector temp(LocalSize);
    double DetJ0;
    MatrixType lhs;
    VectorType rhs;

    // Compute Jacobian
    r_geometry.Jacobian(J0,this->GetIntegrationMethod());

    // Iterate over integration points to evaluate local contribution
    for (unsigned int g = 0; g < integration_points.size(); g++) {

        const double k = rVariables.diffusivity;
        const Vector f = rVariables.forcing;
        const Vector phi = rVariables.unknown;

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(J0[g],InvJ0,DetJ0);
        // Calculate the cartesian derivatives on integration point "g"
        noalias(DN) = prod(DN_De[g],InvJ0);
        // Caluclate N on the gauss point "g"
        auto N = row(N_gausspoint,g);
        // Compute weight
        const double IntToReferenceWeight = integration_points[g].Weight() * DetJ0;

        if (Dimension == 2){

            //substitute_lhs_2D

            //substitute_rhs_2D
        }
        else if (Dimension == 3){

            //substitute_lhs_3D

            //substitute_rhs_3D
        }

        noalias(rLeftHandSideMatrix) += lhs * IntToReferenceWeight;
        noalias(rRightHandSideVector) += rhs * IntToReferenceWeight;

    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

int SymbolicEulerianConvectionDiffusionExplicit::Check(const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;
    int out = Element::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;

    return 0;

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public I/O


///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected operations

void SymbolicEulerianConvectionDiffusionExplicit::InitializeEulerianElement(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const ProcessInfo& r_process_info = rCurrentProcessInfo;
        ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
        auto& r_settings = *p_settings;

        const auto& r_geometry = GetGeometry();
        const unsigned int LocalSize = r_geometry.size();

        rVariables.diffusivity = r_settings.GetDiffusionVariable();

        for(unsigned int node_element = 0; node_element<LocalSize; node_element++)
    {
        rVariables.unknown[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetUnknownVariable());
        rVariables.forcing[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVolumeSourceVariable());
    }

        KRATOS_CATCH( "" )
    }

Element::IntegrationMethod SymbolicEulerianConvectionDiffusionExplicit::GetIntegrationMethod() const
{
    return GeometryData::GI_GAUSS_1;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private serialization

void SymbolicEulerianConvectionDiffusionExplicit::save(Serializer& rSerializer) const
{
    using BaseType = Element;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType );
}


void SymbolicEulerianConvectionDiffusionExplicit::load(Serializer& rSerializer)
{
    using BaseType = Element;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

}