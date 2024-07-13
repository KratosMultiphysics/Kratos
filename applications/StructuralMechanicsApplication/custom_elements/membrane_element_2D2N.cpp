// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "utilities/atomic_utilities.h"

// Application includes
#include "custom_elements/membrane_element_2D2N.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{

void MembraneElement2D2N::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        if(this->UseGeometryIntegrationMethod()) {
            if( GetProperties().Has(INTEGRATION_ORDER) ) {
                const SizeType integration_order = GetProperties()[INTEGRATION_ORDER];
                switch (integration_order)
                {
                case 1:
                    mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
                    break;
                case 2:
                    mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
                    break;
                case 3:
                    mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_3;
                    break;
                case 4:
                    mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_4;
                    break;
                case 5:
                    mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
                    break;
                default:
                    KRATOS_WARNING("MembraneElement2D2N") << "Integration order " << integration_order << " is not available, using default integration order for the geometry" << std::endl;
                    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
                }
            } else {
                mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
            }
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
}

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::InitializeNonLinearIteration( const ProcessInfo& rCurrentProcessInfo )
{
}

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::FinalizeNonLinearIteration( const ProcessInfo& rCurrentProcessInfo )
{
}

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::FinalizeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
{
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer MembraneElement2D2N::Create(
    IndexType NewId,
    NodesArrayType const &ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<MembraneElement2D2N>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer MembraneElement2D2N::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<MembraneElement2D2N>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer MembraneElement2D2N::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes) const
{
    KRATOS_TRY

    MembraneElement2D2N::Pointer p_new_elem = Kratos::make_intrusive<MembraneElement2D2N>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(mThisIntegrationMethod);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const SizeType n_nodes = r_geom.PointsNumber();
    if (rResult.size() != n_nodes) {
        rResult.resize(n_nodes, false);
    }

    rResult[0] = r_geom[0].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[1] = r_geom[1].GetDof(DISPLACEMENT_Y).EquationId();

    KRATOS_CATCH("")
};

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const SizeType n_nodes = r_geom.PointsNumber();
    if (rElementalDofList.size() != n_nodes) {
        rElementalDofList.resize(n_nodes);
    }

    rElementalDofList[0] = r_geom[0].pGetDof(DISPLACEMENT_Y);
    rElementalDofList[1] = r_geom[1].pGetDof(DISPLACEMENT_Y);

    KRATOS_CATCH("")
};

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::GetFirstDerivativesVector(
    Vector& rValues,
    int Step) const
{
    const auto& r_geom = GetGeometry();
    const SizeType n_nodes = r_geom.PointsNumber();
    if (rValues.size() != n_nodes) {
        rValues.resize(n_nodes, false);
    }

    rValues[0] = r_geom[0].FastGetSolutionStepValue(VELOCITY_Y, Step);
    rValues[1] = r_geom[1].FastGetSolutionStepValue(VELOCITY_Y, Step);
}

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::GetSecondDerivativesVector(
    Vector& rValues,
    int Step) const
{
    const auto& r_geom = GetGeometry();
    const SizeType n_nodes = r_geom.PointsNumber();
    if (rValues.size() != n_nodes) {
        rValues.resize(n_nodes, false);
    }

    rValues[0] = r_geom[0].FastGetSolutionStepValue(ACCELERATION_Y, Step);
    rValues[1] = r_geom[1].FastGetSolutionStepValue(ACCELERATION_Y, Step);
}

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Get element data
    const auto& r_geom = GetGeometry();
    const auto& r_prop = GetProperties();
    const SizeType n_nodes = r_geom.PointsNumber();

    // Get membrane data
    const double h = r_prop.GetValue(THICKNESS);
    const double young = r_prop.GetValue(YOUNG_MODULUS);
    const double stress_0 = 1.0e3; //TODO: Properly get the prestress

    // Calculate linearised tensile stress (note that we are assuming it constant in the element)
    const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength2D2N(*this);
    const double l_0 = StructuralMechanicsElementUtilities::CalculateReferenceLength2D2N(*this);
    const double strain = (l - l_0) / l_0;
    const double stress = (young * strain + stress_0) * h;

    // Check RHS size
    if (rRightHandSideVector.size() != n_nodes) {
        rRightHandSideVector.resize(n_nodes, false);
    }

    // Check LHS size
    if (rLeftHandSideMatrix.size1() != n_nodes || rLeftHandSideMatrix.size2() != n_nodes) {
        rLeftHandSideMatrix.resize(n_nodes, n_nodes, false);
    }

    // Calculate the RHS and LHS contributions
    rLeftHandSideMatrix.clear();
    rRightHandSideVector.clear();

    Vector N(n_nodes);
    Matrix DN_De(n_nodes, 1);
    const SizeType n_gauss = r_geom.IntegrationPointsNumber(GetIntegrationMethod());
    const auto& r_integration_points = r_geom.IntegrationPoints(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {

        // Calculate kinematics
        r_geom.ShapeFunctionsValues(N, r_integration_points[i_gauss].Coordinates());
        r_geom.ShapeFunctionsLocalGradients(DN_De, r_integration_points[i_gauss].Coordinates());
        const double detJ = r_geom.DeterminantOfJacobian(r_integration_points[i_gauss].Coordinates());
        const double w_gauss = detJ * r_integration_points[i_gauss].Weight();

        // Calculate the linearised membrane stress
        double lin_rot = 0.0;
        for (IndexType i = 0; i < n_nodes; ++i) {
            lin_rot += (1.0/detJ) * DN_De(i,0) * r_geom[i].FastGetSolutionStepValue(DISPLACEMENT_Y);
        }
        const double lin_stress_term = stress /  std::sqrt(std::pow(1.0 + std::pow(lin_rot, 2), 3));

        // Assemble nodal contributions
        for (IndexType i = 0; i < n_nodes; ++i) {
            for (IndexType j = 0; j < n_nodes; ++j) {
                // Pressure contribution
                const double p = r_geom[j].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE) - r_geom[j].FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE);
                rRightHandSideVector[i] += w_gauss * N[i] * N[j] * p; //TODO: Extend to also consider PRESSURE from the element database.

                // Linearised stress contribution
                const double aux = w_gauss * (DN_De(i,0) / detJ) * lin_stress_term * (DN_De(j,0) / detJ);
                rLeftHandSideMatrix(i,j) += aux;
                rRightHandSideVector[i] -= aux * r_geom[j].FastGetSolutionStepValue(DISPLACEMENT_Y);
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Get element data
    const auto& r_geom = GetGeometry();
    const auto& r_prop = GetProperties();
    const SizeType n_nodes = r_geom.PointsNumber();

    // Get membrane data
    const double h = r_prop.GetValue(THICKNESS);
    const double young = r_prop.GetValue(YOUNG_MODULUS);
    const double stress_0 = 0.0; //TODO: Properly get the prestress

    // Calculate linearised tensile stress (note that we are assuming it constant in the element)
    const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength2D2N(*this);
    const double l_0 = StructuralMechanicsElementUtilities::CalculateReferenceLength2D2N(*this);
    const double strain = (l - l_0) / l_0;
    const double stress = (young * strain + stress_0) * h;;

    // Check LHS size
    if (rLeftHandSideMatrix.size1() != n_nodes || rLeftHandSideMatrix.size2() != n_nodes) {
        rLeftHandSideMatrix.resize(n_nodes, n_nodes, false);
    }

    // Calculate the LHS contribution
    rLeftHandSideMatrix.clear();

    Vector N(n_nodes);
    Matrix DN_De(n_nodes, 1);
    const SizeType n_gauss = r_geom.IntegrationPointsNumber(GetIntegrationMethod());
    const auto& r_integration_points = r_geom.IntegrationPoints(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {

        // Calculate kinematics
        r_geom.ShapeFunctionsValues(N, r_integration_points[i_gauss].Coordinates());
        r_geom.ShapeFunctionsLocalGradients(DN_De, r_integration_points[i_gauss].Coordinates());
        const double detJ = r_geom.DeterminantOfJacobian(r_integration_points[i_gauss].Coordinates());
        const double w_gauss = detJ * r_integration_points[i_gauss].Weight();

        // Calculate the linearised membrane stress
        double lin_rot = 0.0;
        for (IndexType i = 0; i < n_nodes; ++i) {
            lin_rot += (1.0/detJ) * DN_De(i,0) * r_geom[i].FastGetSolutionStepValue(DISPLACEMENT_Y);
        }
        const double lin_stress_term = stress /  std::sqrt(std::pow(1.0 + std::pow(lin_rot, 2), 3));

        // Assemble nodal contributions
        for (IndexType i = 0; i < n_nodes; ++i) {
            for (IndexType j = 0; j < n_nodes; ++j) {
                // Linearised stress contribution
                rLeftHandSideMatrix(i,j) += w_gauss * (DN_De(i,0) / detJ) * lin_stress_term * (DN_De(j,0) / detJ);
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Get element data
    const auto& r_geom = GetGeometry();
    const auto& r_prop = GetProperties();
    const SizeType n_nodes = r_geom.PointsNumber();

    // Get membrane data
    const double h = r_prop.GetValue(THICKNESS);
    const double young = r_prop.GetValue(YOUNG_MODULUS);
    const double stress_0 = 0.0; //TODO: Properly get the prestress

    // Calculate linearised tensile stress (note that we are assuming it constant in the element)
    const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength2D2N(*this);
    const double l_0 = StructuralMechanicsElementUtilities::CalculateReferenceLength2D2N(*this);
    const double strain = (l - l_0) / l_0;
    const double stress = (young * strain + stress_0) * h;;

    // Check RHS size
    if (rRightHandSideVector.size() != n_nodes) {
        rRightHandSideVector.resize(n_nodes, false);
    }

    // Calculate the RHS contribution
    rRightHandSideVector.clear();

    Vector N(n_nodes);
    Matrix DN_De(n_nodes, 1);
    const SizeType n_gauss = r_geom.IntegrationPointsNumber(GetIntegrationMethod());
    const auto& r_integration_points = r_geom.IntegrationPoints(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {

        // Calculate kinematics
        r_geom.ShapeFunctionsValues(N, r_integration_points[i_gauss].Coordinates());
        r_geom.ShapeFunctionsLocalGradients(DN_De, r_integration_points[i_gauss].Coordinates());
        const double detJ = r_geom.DeterminantOfJacobian(r_integration_points[i_gauss].Coordinates());
        const double w_gauss = detJ * r_integration_points[i_gauss].Weight();

        // Calculate the linearised membrane stress
        double lin_rot = 0.0;
        for (IndexType i = 0; i < n_nodes; ++i) {
            lin_rot += (1.0/detJ) * DN_De(i,0) * r_geom[i].FastGetSolutionStepValue(DISPLACEMENT_Y);
        }
        const double lin_stress_term = stress /  std::sqrt(std::pow(1.0 + std::pow(lin_rot, 2), 3));

        // Assemble nodal contributions
        for (IndexType i = 0; i < n_nodes; ++i) {
            for (IndexType j = 0; j < n_nodes; ++j) {
                // Pressure contribution
                const double p = r_geom[j].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE) - r_geom[j].FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE);
                rRightHandSideVector[i] += w_gauss * N[i] * N[j] * p; //TODO: Extend to also consider PRESSURE from the element database.

                // Linearised stress contribution
                rRightHandSideVector[i] -=  w_gauss * (DN_De(i,0) / detJ) * lin_stress_term * (DN_De(j,0) / detJ) * r_geom[j].FastGetSolutionStepValue(DISPLACEMENT_Y);
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Get element data
    const auto &r_geom = GetGeometry();
    const auto &r_prop = GetProperties();
    const SizeType n_nodes = r_geom.PointsNumber();

    // Get membrane data
    const double h = r_prop.GetValue(THICKNESS);
    const double rho = r_prop.GetValue(DENSITY);

    // Check the mass matrix size
    if (rMassMatrix.size1() != n_nodes || rMassMatrix.size2() != n_nodes) {
        rMassMatrix.resize(n_nodes, n_nodes, false);
    }

    // Calculate the mass matrixcontributions
    rMassMatrix.clear();

    Vector N(n_nodes);
    const SizeType n_gauss = r_geom.IntegrationPointsNumber(GetIntegrationMethod());
    const auto& r_integration_points = r_geom.IntegrationPoints(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {

        // Calculate kinematics
        r_geom.ShapeFunctionsValues(N, r_integration_points[i_gauss].Coordinates());
        const double detJ = r_geom.DeterminantOfJacobian(r_integration_points[i_gauss].Coordinates());
        const double w_gauss = detJ * r_integration_points[i_gauss].Weight();

        // Assemble nodal contributions
        const double aux = w_gauss * h * rho;
        for (IndexType i = 0; i < n_nodes; ++i) {
            for (IndexType j = 0; j < n_nodes; ++j) {
                rMassMatrix(i,j) += aux * N[i] * N[j];
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    StructuralMechanicsElementUtilities::CalculateRayleighDampingMatrix(
        *this,
        rDampingMatrix,
        rCurrentProcessInfo,
        GetGeometry().PointsNumber());

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

int  MembraneElement2D2N::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    int check = Element::Check(rCurrentProcessInfo);

    return check;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters MembraneElement2D2N::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["implicit"],
        "framework"                  : "lagrangian",
        "symmetric_lhs"              : true,
        "positive_definite_lhs"      : true,
        "output"                     : {
            "gauss_point"            : [],
            "nodal_historical"       : ["DISPLACEMENT","VELOCITY","ACCELERATION"],
            "nodal_non_historical"   : [],
            "entity"                 : []
        },
        "required_variables"         : ["DISPLACEMENT","VELOCITY","ACCELERATION","POSITIVE_FACE_PRESSURE","NEGATIVE_FACE_PRESSURE"],
        "required_dofs"              : ["DISPLACEMENT_Y"],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Line2D2N"],
        "element_integrates_in_time" : false,
        "compatible_constitutive_laws": {
            "type"        : [],
            "dimension"   : [],
            "strain_size" : []
        },
        "required_polynomial_degree_of_geometry" : -1,
        "documentation"   : "This element implements a simplified linear elastic two-dimensional membrane model."
    })");

    return specifications;
}

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
    int IntMethod = int(this->GetIntegrationMethod());
    rSerializer.save("IntegrationMethod",IntMethod);
}

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
}

} // Namespace Kratos
