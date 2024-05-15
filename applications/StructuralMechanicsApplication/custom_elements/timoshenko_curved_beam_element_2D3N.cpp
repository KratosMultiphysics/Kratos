// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Alejandro Cornejo
//
//

// System includes

// External includes

// Project includes

// Application includes
#include "custom_elements/timoshenko_curved_beam_element_2D3N.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

void LinearTimoshenkoCurvedBeamElement2D3N::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        if (this->UseGeometryIntegrationMethod()) {
            if (GetProperties().Has(INTEGRATION_ORDER) ) {
                mThisIntegrationMethod = static_cast<GeometryData::IntegrationMethod>(GetProperties()[INTEGRATION_ORDER] - 1);
            } else {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
            }
        }

        const auto& r_integration_points = this->IntegrationPoints(mThisIntegrationMethod);

        // Constitutive Law initialisation
        if (mConstitutiveLawVector.size() != r_integration_points.size())
            mConstitutiveLawVector.resize(r_integration_points.size());
        InitializeMaterial();
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::InitializeMaterial()
{
    KRATOS_TRY

    if (GetProperties()[CONSTITUTIVE_LAW] != nullptr) {
        const auto& r_geometry   = GetGeometry();
        const auto& r_properties = GetProperties();
        auto N_values            = Vector();
        for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
            mConstitutiveLawVector[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[point_number]->InitializeMaterial(r_properties, r_geometry, N_values);
        }
    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer LinearTimoshenkoCurvedBeamElement2D3N::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    LinearTimoshenkoCurvedBeamElement2D3N::Pointer p_new_elem = Kratos::make_intrusive<LinearTimoshenkoCurvedBeamElement2D3N>
        (NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    const auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dofs_per_node = DoFperNode; // u, v, theta

    IndexType local_index = 0;

    if (rResult.size() != dofs_per_node * number_of_nodes)
        rResult.resize(dofs_per_node * number_of_nodes, false);

    const IndexType xpos    = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);
    const IndexType rot_pos = this->GetGeometry()[0].GetDofPosition(ROTATION_Z);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_X, xpos    ).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_Y, xpos + 1).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(ROTATION_Z    , rot_pos ).EquationId();
    }
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType dofs_per_node = DoFperNode; // u, v, theta
    rElementalDofList.resize(dofs_per_node * number_of_nodes);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const SizeType index = i * dofs_per_node;
        rElementalDofList[index]     = r_geom[i].pGetDof(DISPLACEMENT_X);
        rElementalDofList[index + 1] = r_geom[i].pGetDof(DISPLACEMENT_Y);
        rElementalDofList[index + 2] = r_geom[i].pGetDof(ROTATION_Z    );
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

const double LinearTimoshenkoCurvedBeamElement2D3N::GetJacobian(const double xi)
{
    GlobalSizeVector r_dN_dxi;
    GetLocalFirstDerivativesNu0ShapeFunctionsValues(r_dN_dxi, xi);
    const auto& r_geom = GetGeometry();

    double dx_dxi = 0.0;
    double dy_dxi = 0.0;

    for (IndexType i = 0; i < NumberOfNodes; ++i) {
        const IndexType u_coord = DoFperNode * i;
        const auto &r_coords_node = r_geom[i].GetInitialPosition();
        dx_dxi += r_coords_node[0] * r_dN_dxi[u_coord];
        dy_dxi += r_coords_node[1] * r_dN_dxi[u_coord];
    }
    return std::sqrt(std::pow(dx_dxi, 2) + std::pow(dy_dxi, 2));
}

/***********************************************************************************/
/***********************************************************************************/

const double LinearTimoshenkoCurvedBeamElement2D3N::GetGeometryCurvature(
    const double J,
    const double xi
    )
{
    GlobalSizeVector dN_dxi, d2N_dxi2;
    GetLocalFirstDerivativesNu0ShapeFunctionsValues (dN_dxi,   xi);
    GetLocalSecondDerivativesNu0ShapeFunctionsValues(d2N_dxi2, xi);
    const auto& r_geom = GetGeometry();

    double dx_dxi = 0.0;
    double dy_dxi = 0.0;

    double d2x_dxi2 = 0.0;
    double d2y_dxi2 = 0.0;

    for (IndexType i = 0; i < NumberOfNodes; ++i) {
        const IndexType u_coord = DoFperNode * i;
        const auto &r_coords_node = r_geom[i].GetInitialPosition();
        dx_dxi += r_coords_node[0] * dN_dxi[u_coord];
        dy_dxi += r_coords_node[1] * dN_dxi[u_coord];

        d2x_dxi2 += r_coords_node[0] * d2N_dxi2[u_coord];
        d2y_dxi2 += r_coords_node[1] * d2N_dxi2[u_coord];
    }
    return (dx_dxi * d2y_dxi2 - dy_dxi * d2x_dxi2) / std::pow(J, 3);
}

/***********************************************************************************/
/***********************************************************************************/

const double LinearTimoshenkoCurvedBeamElement2D3N::GetBendingShearStiffnessRatio()
{
    const auto &r_props = GetProperties();
    const double E  = r_props[YOUNG_MODULUS];
    const double I  = r_props[I33];
    const double As = r_props[AREA_EFFECTIVE_Y];
    const double G  = ConstitutiveLawUtilities<3>::CalculateShearModulus(r_props);

    if (As == 0.0) // If effective area is null -> Euler Bernouilli case
        return 0.0;
    else
        return E * I / (G * As);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetShapeFunctionsValues(
    GlobalSizeVector& rN,
    const double J,
    const double xi
    )
{
    // Nodal values of the Jacobian
    const double J1 = GetJacobian(-1.0);
    const double J2 = GetJacobian( 1.0);
    const double J3 = GetJacobian( 0.0);

    // Nodal values of the geometric curvatures
    const double k01 = GetGeometryCurvature(J1, -1.0);
    const double k02 = GetGeometryCurvature(J2,  1.0);
    const double k03 = GetGeometryCurvature(J3,  0.0);

    const double k_s = GetBendingShearStiffnessRatio();

    const double cr_N0 = std::pow(J2, 2);
    const double cr_N1 = 6.0*k_s;
    const double cr_N2 = cr_N0*cr_N1;
    const double cr_N3 = std::pow(xi, 3);
    const double cr_N4 = std::pow(k_s, 2);
    const double cr_N5 = 180.0*cr_N4;
    const double cr_N6 = cr_N3*cr_N5;
    const double cr_N7 = cr_N5*xi;
    const double cr_N8 = std::pow(J3, 2);
    const double cr_N9 = cr_N0*cr_N8;
    const double cr_N10 = cr_N9*xi;
    const double cr_N11 = std::pow(xi, 4);
    const double cr_N12 = cr_N11*cr_N9;
    const double cr_N13 = std::pow(xi, 2);
    const double cr_N14 = cr_N13*cr_N9;
    const double cr_N15 = cr_N3*cr_N9;
    const double cr_N16 = 12.0*k_s;
    const double cr_N17 = cr_N0*cr_N3;
    const double cr_N18 = k_s*xi;
    const double cr_N19 = 27.0*cr_N8;
    const double cr_N20 = cr_N18*cr_N19;
    const double cr_N21 = cr_N0*xi;
    const double cr_N22 = cr_N3*k_s;
    const double cr_N23 = cr_N19*cr_N22;
    const double cr_N24 = 72.0*cr_N4;
    const double cr_N25 = cr_N11*cr_N24;
    const double cr_N26 = cr_N16*cr_N8;
    const double cr_N27 = cr_N11*cr_N26 - cr_N13*cr_N26 + cr_N24 - cr_N25;
    const double cr_N28 = std::pow(k_s, 3);
    const double cr_N29 = 8640.0*cr_N28;
    const double cr_N30 = std::pow(J1, 2);
    const double cr_N31 = 648.0*cr_N4;
    const double cr_N32 = cr_N4*cr_N8;
    const double cr_N33 = 1296.0*cr_N32;
    const double cr_N34 = cr_N0*cr_N30;
    const double cr_N35 = cr_N34*k_s;
    const double cr_N36 = 48.0*cr_N35;
    const double cr_N37 = cr_N30*cr_N8;
    const double cr_N38 = 78.0*k_s;
    const double cr_N39 = cr_N30*cr_N9;
    const double cr_N40 = xi/(-cr_N0*cr_N31 - cr_N29 - cr_N30*cr_N31 + cr_N33 - cr_N36 + cr_N37*cr_N38 + cr_N38*cr_N9 + 4.0*cr_N39);
    const double cr_N41 = std::pow(J1, 3)*cr_N40*(cr_N10 - cr_N11*cr_N2 + cr_N12 - cr_N14 - cr_N15 + cr_N16*cr_N17 - cr_N16*cr_N21 + cr_N2 + cr_N20 - cr_N23 + cr_N27 + cr_N6 - cr_N7);
    const double cr_N42 = 576.0*cr_N4;
    const double cr_N43 = 864.0*cr_N4;
    const double cr_N44 = cr_N29*xi;
    const double cr_N45 = cr_N33*xi;
    const double cr_N46 = cr_N30*cr_N4;
    const double cr_N47 = 216.0*cr_N11;
    const double cr_N48 = cr_N21*cr_N4;
    const double cr_N49 = cr_N46*xi;
    const double cr_N50 = 8.0*cr_N10*cr_N30;
    const double cr_N51 = cr_N3*cr_N36;
    const double cr_N52 = 18.0*k_s;
    const double cr_N53 = cr_N15*k_s;
    const double cr_N54 = cr_N11*cr_N37;
    const double cr_N55 = 66.0*k_s;
    const double cr_N56 = 57.0*k_s;
    const double cr_N57 = cr_N37*xi;
    const double cr_N58 = 207.0*k_s;
    const double cr_N59 = 4.0*cr_N30;
    const double cr_N60 = cr_N15*cr_N59;
    const double cr_N61 = 96.0*cr_N35;
    const double cr_N62 = cr_N61*xi;
    const double cr_N63 = 96.0*k_s;
    const double cr_N64 = cr_N13*cr_N37;
    const double cr_N65 = 144.0*k_s;
    const double cr_N66 = cr_N3*cr_N37;
    const double cr_N67 = cr_N66*k_s;
    const double cr_N68 = cr_N12*cr_N30;
    const double cr_N69 = 144.0*cr_N11*cr_N32 - 12.0*cr_N11*cr_N35 - 1440.0*cr_N13*cr_N32 - 10.0*cr_N14*cr_N30 + cr_N29 + 60.0*cr_N35 + 6.0*cr_N68;
    const double cr_N70 = 1296.0*cr_N4;
    const double cr_N71 = 156.0*k_s;
    const double cr_N72 = xi/(-cr_N0*cr_N70 - 17280.0*cr_N28 - cr_N30*cr_N70 + 2592.0*cr_N32 + cr_N37*cr_N71 + 8.0*cr_N39 - cr_N61 + cr_N71*cr_N9);
    const double cr_N73 = cr_N1*cr_N30;
    const double cr_N74 = cr_N16*cr_N30;
    const double cr_N75 = std::pow(J2, 3)*cr_N40*(-cr_N11*cr_N73 - cr_N20 + cr_N23 + cr_N27 - cr_N3*cr_N74 + cr_N54 - cr_N57 - cr_N6 - cr_N64 + cr_N66 + cr_N7 + cr_N73 + cr_N74*xi);
    const double cr_N76 = 4320.0*cr_N28;
    const double cr_N77 = 324.0*cr_N4;
    const double cr_N78 = cr_N31*cr_N8;
    const double cr_N79 = 24.0*cr_N35;
    const double cr_N80 = 39.0*k_s;
    const double cr_N81 = -cr_N0*cr_N77 - cr_N30*cr_N77 + cr_N37*cr_N80 + 2.0*cr_N39 - cr_N76 + cr_N78 - cr_N79 + cr_N80*cr_N9;
    const double cr_N82 = 1.0/cr_N81;
    const double cr_N83 = 24.0*cr_N18;
    const double cr_N84 = 2.0*cr_N34;
    const double cr_N85 = cr_N13*cr_N4;
    const double cr_N86 = cr_N30*k_s;
    const double cr_N87 = 54.0*cr_N13;
    const double cr_N88 = 15.0*cr_N11;
    const double cr_N89 = cr_N0*k_s;
    const double cr_N90 = std::pow(J3, 3)*cr_N82*xi*(cr_N0*cr_N80 - cr_N0*cr_N83 + cr_N11*cr_N84 - 4.0*cr_N13*cr_N34 + 24.0*cr_N17*k_s - 24.0*cr_N22*cr_N30 + cr_N25 + cr_N30*cr_N80 + cr_N30*cr_N83 + cr_N31 + cr_N84 - 720.0*cr_N85 - cr_N86*cr_N87 + cr_N86*cr_N88 - cr_N87*cr_N89 + cr_N88*cr_N89);
    const double cr_N91 = std::pow(xi, 5);
    const double cr_N92 = cr_N24*cr_N30;
    const double cr_N93 = 504.0*cr_N85;
    const double cr_N94 = cr_N11*cr_N5;
    const double cr_N95 = cr_N16*cr_N91;
    const double cr_N96 = 27.0*k_s;
    rN[0]=-cr_N41*k01;
    rN[1]=cr_N72*(cr_N0*cr_N25 + cr_N0*cr_N42 + cr_N0*cr_N6 + cr_N10*cr_N56 + cr_N12*cr_N52 - cr_N14*cr_N63 + 540.0*cr_N3*cr_N46 + cr_N30*cr_N43 - cr_N44 + cr_N45 - cr_N46*cr_N47 - 828.0*cr_N48 - 1188.0*cr_N49 + cr_N50 + cr_N51 + 21.0*cr_N53 + cr_N54*cr_N55 + cr_N57*cr_N58 - cr_N60 - cr_N62 - cr_N64*cr_N65 - 129.0*cr_N67 + cr_N69);
    rN[2]=cr_N41;
    rN[3]=-cr_N75*k02;
    rN[4]=-cr_N72*(-cr_N0*cr_N4*cr_N47 + cr_N0*cr_N43 - cr_N10*cr_N58 + cr_N12*cr_N55 - cr_N14*cr_N65 - 540.0*cr_N17*cr_N4 + cr_N25*cr_N30 + cr_N30*cr_N42 - cr_N30*cr_N6 + cr_N44 - cr_N45 + 1188.0*cr_N48 + 828.0*cr_N49 - cr_N50 - cr_N51 + cr_N52*cr_N54 + 129.0*cr_N53 - cr_N56*cr_N57 + cr_N60 + cr_N62 - cr_N63*cr_N64 - 21.0*cr_N67 + cr_N69);
    rN[5]=cr_N75;
    rN[6]=-cr_N90*k03;
    rN[7]=cr_N82*(-cr_N0*cr_N24*cr_N91 + cr_N0*cr_N93 - cr_N0*cr_N94 - cr_N11*cr_N79 + cr_N12*cr_N96 + cr_N13*cr_N36 + cr_N13*cr_N76 - cr_N13*cr_N78 - cr_N14*cr_N55 - cr_N14*cr_N59 - cr_N15*cr_N16 + cr_N16*cr_N66 + cr_N21*cr_N24 + cr_N30*cr_N93 - cr_N30*cr_N94 - cr_N37*cr_N95 + cr_N54*cr_N96 - cr_N55*cr_N64 + 2.0*cr_N68 + cr_N81 + cr_N9*cr_N95 + cr_N91*cr_N92 - cr_N92*xi);
    rN[8]=cr_N90;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetFirstDerivativesShapeFunctionsValues(
    GlobalSizeVector& rdN,
    const double J,
    const double xi
    )
{
    // Nodal values of the Jacobian
    const double J1 = GetJacobian(-1.0);
    const double J2 = GetJacobian( 1.0);
    const double J3 = GetJacobian( 0.0);

    // Nodal values of the geometric curvatures
    const double k01 = GetGeometryCurvature(J1, -1.0);
    const double k02 = GetGeometryCurvature(J2,  1.0);
    const double k03 = GetGeometryCurvature(J3,  0.0);

    const double k_s = GetBendingShearStiffnessRatio();

    const double crdN0 = std::pow(k_s, 3);
    const double crdN1 = 8640.0*crdN0;
    const double crdN2 = std::pow(J1, 2);
    const double crdN3 = std::pow(k_s, 2);
    const double crdN4 = 648.0*crdN3;
    const double crdN5 = std::pow(J2, 2);
    const double crdN6 = std::pow(J3, 2);
    const double crdN7 = crdN3*crdN6;
    const double crdN8 = 1296.0*crdN7;
    const double crdN9 = crdN2*crdN5;
    const double crdN10 = 48.0*k_s;
    const double crdN11 = crdN6*k_s;
    const double crdN12 = 78.0*crdN11;
    const double crdN13 = crdN6*crdN9;
    const double crdN14 = 1.0/(-crdN1 - crdN10*crdN9 + crdN12*crdN2 + crdN12*crdN5 + 4.0*crdN13 - crdN2*crdN4 - crdN4*crdN5 + crdN8);
    const double crdN15 = 6.0*k_s;
    const double crdN16 = std::pow(xi, 3);
    const double crdN17 = crdN16*crdN3;
    const double crdN18 = 720.0*crdN17;
    const double crdN19 = 360.0*crdN3;
    const double crdN20 = crdN19*xi;
    const double crdN21 = crdN5*crdN6;
    const double crdN22 = 4*crdN16;
    const double crdN23 = std::pow(xi, 2);
    const double crdN24 = 3*crdN23;
    const double crdN25 = 2*xi;
    const double crdN26 = std::pow(xi, 4);
    const double crdN27 = 5*crdN26;
    const double crdN28 = crdN10*crdN16;
    const double crdN29 = crdN11*xi;
    const double crdN30 = 54.0*crdN29;
    const double crdN31 = crdN5*k_s;
    const double crdN32 = 24.0*xi;
    const double crdN33 = 30.0*crdN26;
    const double crdN34 = crdN11*crdN16;
    const double crdN35 = 108.0*crdN34;
    const double crdN36 = 72.0*crdN3;
    const double crdN37 = crdN19*crdN26;
    const double crdN38 = crdN11*crdN26;
    const double crdN39 = 60.0*crdN38;
    const double crdN40 = crdN11*crdN23;
    const double crdN41 = 36.0*crdN40;
    const double crdN42 = crdN36 - crdN37 + crdN39 - crdN41;
    const double crdN43 = std::pow(J1, 3)*crdN14*(crdN15*crdN5 + crdN18 - crdN20 - crdN21*crdN22 - crdN21*crdN24 + crdN21*crdN25 + crdN21*crdN27 + crdN28*crdN5 + crdN30 - crdN31*crdN32 - crdN31*crdN33 - crdN35 + crdN42);
    const double crdN44 = 17280.0*crdN0;
    const double crdN45 = 1296.0*crdN3;
    const double crdN46 = 2592.0*crdN7;
    const double crdN47 = crdN9*k_s;
    const double crdN48 = 96.0*crdN47;
    const double crdN49 = 156.0*crdN11;
    const double crdN50 = 8.0*crdN13;
    const double crdN51 = 1.0/(-crdN2*crdN45 + crdN2*crdN49 - crdN44 - crdN45*crdN5 + crdN46 - crdN48 + crdN49*crdN5 + crdN50);
    const double crdN52 = 576.0*crdN3;
    const double crdN53 = 864.0*crdN3;
    const double crdN54 = crdN44*xi;
    const double crdN55 = crdN37*crdN5;
    const double crdN56 = crdN18*crdN5;
    const double crdN57 = crdN46*xi;
    const double crdN58 = 2160.0*crdN17;
    const double crdN59 = crdN2*crdN3;
    const double crdN60 = 1080.0*crdN26;
    const double crdN61 = 1656.0*xi;
    const double crdN62 = crdN3*crdN5;
    const double crdN63 = 2376.0*xi;
    const double crdN64 = 16.0*crdN13;
    const double crdN65 = crdN64*xi;
    const double crdN66 = 192.0*crdN47;
    const double crdN67 = crdN16*crdN66;
    const double crdN68 = crdN34*crdN5;
    const double crdN69 = 90.0*crdN38;
    const double crdN70 = 114.0*crdN29;
    const double crdN71 = 330.0*crdN38;
    const double crdN72 = 414.0*crdN29;
    const double crdN73 = crdN16*crdN64;
    const double crdN74 = crdN66*xi;
    const double crdN75 = 288.0*crdN40;
    const double crdN76 = 432.0*crdN40;
    const double crdN77 = crdN2*crdN34;
    const double crdN78 = 60.0*crdN47;
    const double crdN79 = crdN1 - 30.0*crdN13*crdN23 + crdN13*crdN33 - 4320.0*crdN23*crdN7 + 720.0*crdN26*crdN7 - crdN26*crdN78 + crdN78;
    const double crdN80 = crdN2*crdN6;
    const double crdN81 = crdN2*k_s;
    const double crdN82 = std::pow(J2, 3)*crdN14*(crdN15*crdN2 - crdN18 - crdN2*crdN28 + crdN20 + crdN22*crdN80 - crdN24*crdN80 - crdN25*crdN80 + crdN27*crdN80 - crdN30 + crdN32*crdN81 - crdN33*crdN81 + crdN35 + crdN42);
    const double crdN83 = -crdN18*crdN2 + crdN2*crdN37;
    const double crdN84 = 324.0*crdN3;
    const double crdN85 = 39.0*k_s;
    const double crdN86 = crdN2*crdN85;
    const double crdN87 = crdN5*crdN85;
    const double crdN88 = 2.0*crdN9;
    const double crdN89 = 1.0/(-4320.0*crdN0 - crdN2*crdN84 + crdN4*crdN6 - 24.0*crdN47 - crdN5*crdN84 + crdN6*crdN86 + crdN6*crdN87 + crdN6*crdN88);
    const double crdN90 = crdN10*xi;
    const double crdN91 = 162.0*crdN23;
    const double crdN92 = 96.0*crdN16;
    const double crdN93 = 75.0*crdN26;
    const double crdN94 = std::pow(J3, 3)*crdN89*(crdN2*crdN90 - 2160.0*crdN23*crdN3 - 12.0*crdN23*crdN9 + 10.0*crdN26*crdN9 - crdN31*crdN91 + crdN31*crdN92 + crdN31*crdN93 + crdN37 + crdN4 - crdN5*crdN90 - crdN81*crdN91 - crdN81*crdN92 + crdN81*crdN93 + crdN86 + crdN87 + crdN88);
    const double crdN95 = 1008.0*xi;
    const double crdN96 = 132.0*crdN29;

    rdN[0]=-crdN43*k01;
    rdN[1]=crdN51*(crdN2*crdN53 + crdN2*crdN58 + crdN2*crdN71 + crdN2*crdN72 - crdN2*crdN76 + crdN5*crdN52 + crdN5*crdN69 + crdN5*crdN70 - crdN5*crdN75 - crdN54 + crdN55 + crdN56 + crdN57 - crdN59*crdN60 - crdN59*crdN63 - crdN61*crdN62 + crdN65 + crdN67 + 84.0*crdN68 - crdN73 - crdN74 - 516.0*crdN77 + crdN79);
    rdN[2]=crdN43;
    rdN[3]=-crdN82*k02;
    rdN[4]=-crdN51*(crdN2*crdN52 + crdN2*crdN69 - crdN2*crdN70 - crdN2*crdN75 + crdN5*crdN53 - crdN5*crdN58 + crdN5*crdN71 - crdN5*crdN72 - crdN5*crdN76 + crdN54 - crdN57 + crdN59*crdN61 - crdN60*crdN62 + crdN62*crdN63 - crdN65 - crdN67 + 516.0*crdN68 + crdN73 + crdN74 - 84.0*crdN77 + crdN79 + crdN83);
    rdN[5]=crdN82;
    rdN[6]=-crdN94*k03;
    rdN[7]=crdN89*(crdN1*xi - crdN16*crdN48 + crdN16*crdN50 + crdN2*crdN35 - crdN2*crdN36 - crdN2*crdN39 + crdN2*crdN41 - crdN2*crdN96 + crdN35*crdN5 + crdN36*crdN5 + crdN39*crdN5 - crdN41*crdN5 + crdN48*xi - crdN5*crdN96 - crdN50*xi - crdN55 - crdN56 + crdN59*crdN95 + crdN62*crdN95 - crdN8*xi + crdN83);
    rdN[8]=crdN94;

    rdN /= J;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetSecondDerivativesShapeFunctionsValues(
    GlobalSizeVector& rd2N,
    const double J,
    const double xi
    )
{
    // Nodal values of the Jacobian
    const double J1 = GetJacobian(-1.0);
    const double J2 = GetJacobian( 1.0);
    const double J3 = GetJacobian( 0.0);

    // Nodal values of the geometric curvatures
    const double k01 = GetGeometryCurvature(J1, -1.0);
    const double k02 = GetGeometryCurvature(J2,  1.0);
    const double k03 = GetGeometryCurvature(J3,  0.0);

    const double k_s = GetBendingShearStiffnessRatio();

    const double crd2N0 = std::pow(k_s, 3);
    const double crd2N1 = 8640.0*crd2N0;
    const double crd2N2 = std::pow(J1, 2);
    const double crd2N3 = std::pow(k_s, 2);
    const double crd2N4 = 648.0*crd2N3;
    const double crd2N5 = std::pow(J2, 2);
    const double crd2N6 = std::pow(J3, 2);
    const double crd2N7 = crd2N3*crd2N6;
    const double crd2N8 = 1296.0*crd2N7;
    const double crd2N9 = 48.0*k_s;
    const double crd2N10 = crd2N2*crd2N9;
    const double crd2N11 = crd2N6*k_s;
    const double crd2N12 = 78.0*crd2N11;
    const double crd2N13 = crd2N5*crd2N6;
    const double crd2N14 = crd2N13*crd2N2;
    const double crd2N15 = 1.0/(-crd2N1 - crd2N10*crd2N5 + crd2N12*crd2N2 + crd2N12*crd2N5 + 4.0*crd2N14 - crd2N2*crd2N4 - crd2N4*crd2N5 + crd2N8);
    const double crd2N16 = 360.0*crd2N3;
    const double crd2N17 = 2*crd2N6;
    const double crd2N18 = 24.0*k_s;
    const double crd2N19 = crd2N18*crd2N5;
    const double crd2N20 = 54.0*crd2N11;
    const double crd2N21 = std::pow(xi, 2);
    const double crd2N22 = crd2N21*crd2N3;
    const double crd2N23 = 2160.0*crd2N22;
    const double crd2N24 = std::pow(xi, 3);
    const double crd2N25 = 20*crd2N24;
    const double crd2N26 = 6*xi;
    const double crd2N27 = 12*crd2N21;
    const double crd2N28 = crd2N5*k_s;
    const double crd2N29 = 120.0*crd2N24;
    const double crd2N30 = crd2N11*crd2N21;
    const double crd2N31 = 324.0*crd2N30;
    const double crd2N32 = 144.0*crd2N21;
    const double crd2N33 = crd2N24*crd2N3;
    const double crd2N34 = 1440.0*crd2N33;
    const double crd2N35 = crd2N11*xi;
    const double crd2N36 = 72.0*crd2N35;
    const double crd2N37 = crd2N11*crd2N24;
    const double crd2N38 = 240.0*crd2N37;
    const double crd2N39 = crd2N34 + crd2N36 - crd2N38;
    const double crd2N40 = std::pow(J1, 3)*crd2N15*(-crd2N13*crd2N25 + crd2N13*crd2N26 + crd2N13*crd2N27 + crd2N16 - crd2N17*crd2N5 + crd2N19 - crd2N20 - crd2N23 + crd2N28*crd2N29 - crd2N28*crd2N32 + crd2N31 + crd2N39);
    const double crd2N41 = 1296.0*crd2N3;
    const double crd2N42 = 8.0*crd2N14;
    const double crd2N43 = 156.0*crd2N11;
    const double crd2N44 = crd2N2*crd2N28;
    const double crd2N45 = 96.0*crd2N44;
    const double crd2N46 = -17280.0*crd2N0 + 2592.0*crd2N7;
    const double crd2N47 = 1.0/(-crd2N2*crd2N41 + crd2N2*crd2N43 - crd2N41*crd2N5 + crd2N42 + crd2N43*crd2N5 - crd2N45 + crd2N46);
    const double crd2N48 = 1656.0*crd2N3;
    const double crd2N49 = 2376.0*crd2N3;
    const double crd2N50 = crd2N34*crd2N5;
    const double crd2N51 = 2880.0*crd2N24*crd2N7;
    const double crd2N52 = 114.0*crd2N11;
    const double crd2N53 = crd2N23*crd2N5;
    const double crd2N54 = 414.0*crd2N11;
    const double crd2N55 = 6480.0*crd2N22;
    const double crd2N56 = 4320.0*crd2N33;
    const double crd2N57 = 8640.0*crd2N7*xi;
    const double crd2N58 = crd2N14*crd2N29;
    const double crd2N59 = 360.0*crd2N37;
    const double crd2N60 = 252.0*crd2N30;
    const double crd2N61 = 1320.0*crd2N37;
    const double crd2N62 = 576.0*crd2N35;
    const double crd2N63 = 60.0*crd2N14*xi;
    const double crd2N64 = 240.0*crd2N24*crd2N44;
    const double crd2N65 = 864.0*crd2N35;
    const double crd2N66 = 1548.0*crd2N30;
    const double crd2N67 = crd2N14*crd2N21;
    const double crd2N68 = 16.0*crd2N14 + 576.0*crd2N21*crd2N44 - 192.0*crd2N44 + crd2N46 - 48.0*crd2N67;
    const double crd2N69 = crd2N2*crd2N6;
    const double crd2N70 = crd2N2*k_s;
    const double crd2N71 = std::pow(J2, 3)*crd2N15*(-crd2N16 + crd2N17*crd2N2 - crd2N18*crd2N2 + crd2N20 + crd2N23 - crd2N25*crd2N69 + crd2N26*crd2N69 - crd2N27*crd2N69 + crd2N29*crd2N70 - crd2N31 + crd2N32*crd2N70 + crd2N39);
    const double crd2N72 = crd2N2*crd2N23;
    const double crd2N73 = crd2N2*crd2N34;
    const double crd2N74 = 324.0*crd2N3;
    const double crd2N75 = 39.0*crd2N11;
    const double crd2N76 = 1.0/(-4320.0*crd2N0 + 2.0*crd2N14 - crd2N19*crd2N2 - crd2N2*crd2N74 + crd2N2*crd2N75 - crd2N5*crd2N74 + crd2N5*crd2N75 + 648.0*crd2N7);
    const double crd2N77 = 324.0*xi;
    const double crd2N78 = crd2N2*crd2N5;
    const double crd2N79 = 288.0*crd2N21;
    const double crd2N80 = 300.0*crd2N24;
    const double crd2N81 = std::pow(J3, 3)*crd2N76*(crd2N10 + 40.0*crd2N24*crd2N78 - crd2N28*crd2N77 + crd2N28*crd2N79 + crd2N28*crd2N80 - 4320.0*crd2N3*xi + crd2N34 - crd2N5*crd2N9 - crd2N70*crd2N77 - crd2N70*crd2N79 + crd2N70*crd2N80 - 24.0*crd2N78*xi);
    const double crd2N82 = 1008.0*crd2N3;
    const double crd2N83 = 132.0*crd2N11;
    rd2N[0]=crd2N40*k01;
    rd2N[1]=crd2N47*(-crd2N2*crd2N49 + crd2N2*crd2N54 + crd2N2*crd2N55 - crd2N2*crd2N56 + crd2N2*crd2N61 - crd2N2*crd2N65 - crd2N2*crd2N66 - crd2N48*crd2N5 + crd2N5*crd2N52 + crd2N5*crd2N59 + crd2N5*crd2N60 - crd2N5*crd2N62 + crd2N50 + crd2N51 + crd2N53 - crd2N57 + crd2N58 - crd2N63 - crd2N64 + crd2N68);
    rd2N[2]=-crd2N40;
    rd2N[3]=crd2N71*k02;
    rd2N[4]=crd2N47*(-crd2N2*crd2N48 + crd2N2*crd2N52 - crd2N2*crd2N59 + crd2N2*crd2N60 + crd2N2*crd2N62 - crd2N49*crd2N5 + crd2N5*crd2N54 + crd2N5*crd2N55 + crd2N5*crd2N56 - crd2N5*crd2N61 + crd2N5*crd2N65 - crd2N5*crd2N66 - crd2N51 + crd2N57 - crd2N58 + crd2N63 + crd2N64 + crd2N68 + crd2N72 - crd2N73);
    rd2N[5]=-crd2N71;
    rd2N[6]=-crd2N81*k03;
    rd2N[7]=crd2N76*(crd2N1 + crd2N2*crd2N31 + crd2N2*crd2N36 - crd2N2*crd2N38 + crd2N2*crd2N82 - crd2N2*crd2N83 + crd2N31*crd2N5 - crd2N36*crd2N5 + crd2N38*crd2N5 - crd2N42 - crd2N44*crd2N79 + crd2N45 + crd2N5*crd2N82 - crd2N5*crd2N83 - crd2N50 - crd2N53 + 24.0*crd2N67 - crd2N72 + crd2N73 - crd2N8);
    rd2N[8]=crd2N81;

    rd2N /= std::pow(J, 2);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetThirdDerivativesShapeFunctionsValues(
    GlobalSizeVector& rd3N,
    const double J,
    const double xi
    )
{
    // Nodal values of the Jacobian
    const double J1 = GetJacobian(-1.0);
    const double J2 = GetJacobian( 1.0);
    const double J3 = GetJacobian( 0.0);

    // Nodal values of the geometric curvatures
    const double k01 = GetGeometryCurvature(J1, -1.0);
    const double k02 = GetGeometryCurvature(J2,  1.0);
    const double k03 = GetGeometryCurvature(J3,  0.0);

    const double k_s = GetBendingShearStiffnessRatio();

    const double crd3N0 = std::pow(k_s, 3);
    const double crd3N1 = std::pow(J1, 2);
    const double crd3N2 = std::pow(k_s, 2);
    const double crd3N3 = 648.0*crd3N2;
    const double crd3N4 = std::pow(J2, 2);
    const double crd3N5 = std::pow(J3, 2);
    const double crd3N6 = crd3N2*crd3N5;
    const double crd3N7 = crd3N1*crd3N4;
    const double crd3N8 = crd3N7*k_s;
    const double crd3N9 = crd3N5*k_s;
    const double crd3N10 = 78.0*crd3N9;
    const double crd3N11 = crd3N5*crd3N7;
    const double crd3N12 = 1.0/(-8640.0*crd3N0 + crd3N1*crd3N10 - crd3N1*crd3N3 + crd3N10*crd3N4 + 4.0*crd3N11 - crd3N3*crd3N4 + 1296.0*crd3N6 - 48.0*crd3N8);
    const double crd3N13 = 6*crd3N5;
    const double crd3N14 = 4320.0*crd3N2;
    const double crd3N15 = crd3N14*xi;
    const double crd3N16 = crd3N4*crd3N5;
    const double crd3N17 = std::pow(xi, 2);
    const double crd3N18 = 60*crd3N17;
    const double crd3N19 = 24*xi;
    const double crd3N20 = crd3N4*k_s;
    const double crd3N21 = 360.0*crd3N17;
    const double crd3N22 = crd3N9*xi;
    const double crd3N23 = 648.0*crd3N22;
    const double crd3N24 = 288.0*xi;
    const double crd3N25 = 72.0*crd3N9;
    const double crd3N26 = crd3N14*crd3N17;
    const double crd3N27 = crd3N17*crd3N9;
    const double crd3N28 = 720.0*crd3N27;
    const double crd3N29 = crd3N25 + crd3N26 - crd3N28;
    const double crd3N30 = std::pow(J1, 3)*crd3N12*(crd3N13*crd3N4 - crd3N15 - crd3N16*crd3N18 + crd3N16*crd3N19 + crd3N20*crd3N21 - crd3N20*crd3N24 + crd3N23 + crd3N29);
    const double crd3N31 = 1296.0*crd3N2;
    const double crd3N32 = 156.0*crd3N9;
    const double crd3N33 = 1.0/(-17280.0*crd3N0 - crd3N1*crd3N31 + crd3N1*crd3N32 + 8.0*crd3N11 - crd3N31*crd3N4 + crd3N32*crd3N4 + 2592.0*crd3N6 - 96.0*crd3N8);
    const double crd3N34 = 8640.0*crd3N6;
    const double crd3N35 = crd3N15*crd3N4;
    const double crd3N36 = crd3N26*crd3N4;
    const double crd3N37 = crd3N17*crd3N34;
    const double crd3N38 = 12960.0*crd3N2;
    const double crd3N39 = crd3N1*xi;
    const double crd3N40 = 576.0*crd3N9;
    const double crd3N41 = 60.0*crd3N11;
    const double crd3N42 = 864.0*crd3N9;
    const double crd3N43 = crd3N17*crd3N38;
    const double crd3N44 = crd3N11*crd3N21;
    const double crd3N45 = 504.0*crd3N22;
    const double crd3N46 = 1080.0*crd3N27;
    const double crd3N47 = 3960.0*crd3N27;
    const double crd3N48 = 720.0*crd3N17*crd3N8;
    const double crd3N49 = 3096.0*crd3N22;
    const double crd3N50 = crd3N11*xi;
    const double crd3N51 = -96.0*crd3N50 + 1152.0*crd3N8*xi;
    const double crd3N52 = crd3N1*crd3N5;
    const double crd3N53 = crd3N1*k_s;
    const double crd3N54 = std::pow(J2, 3)*crd3N12*(crd3N1*crd3N13 + crd3N15 - crd3N18*crd3N52 - crd3N19*crd3N52 + crd3N21*crd3N53 - crd3N23 + crd3N24*crd3N53 + crd3N29);
    const double crd3N55 = crd3N1*crd3N15;
    const double crd3N56 = crd3N1*crd3N26;
    const double crd3N57 = 324.0*crd3N2;
    const double crd3N58 = 24.0*crd3N7;
    const double crd3N59 = 39.0*crd3N9;
    const double crd3N60 = 1.0/(-4320.0*crd3N0 - crd3N1*crd3N57 + crd3N1*crd3N59 + 2.0*crd3N11 - crd3N4*crd3N57 + crd3N4*crd3N59 - crd3N58*k_s + 648.0*crd3N6);
    const double crd3N61 = 324.0*k_s;
    const double crd3N62 = 576.0*xi;
    const double crd3N63 = 900.0*crd3N17;
    const double crd3N64 = std::pow(J3, 3)*crd3N60*(-crd3N1*crd3N61 - crd3N14 + 120.0*crd3N17*crd3N7 + crd3N20*crd3N62 + crd3N20*crd3N63 + crd3N26 - 576.0*crd3N39*k_s - crd3N4*crd3N61 + crd3N53*crd3N63 - crd3N58);
    rd3N[0]=crd3N30*k01;
    rd3N[1]=crd3N33*(-crd3N1*crd3N42 - crd3N1*crd3N43 + crd3N1*crd3N47 - crd3N1*crd3N49 - crd3N34 + crd3N35 + crd3N36 + crd3N37 + crd3N38*crd3N39 - crd3N4*crd3N40 + crd3N4*crd3N45 + crd3N4*crd3N46 - crd3N41 + crd3N44 - crd3N48 + crd3N51);
    rd3N[2]=-crd3N30;
    rd3N[3]=crd3N54*k02;
    rd3N[4]=crd3N33*(crd3N1*crd3N40 + crd3N1*crd3N45 - crd3N1*crd3N46 + crd3N34 - crd3N37 + crd3N38*crd3N4*xi + crd3N4*crd3N42 + crd3N4*crd3N43 - crd3N4*crd3N47 - crd3N4*crd3N49 + crd3N41 - crd3N44 + crd3N48 + crd3N51 + crd3N55 - crd3N56);
    rd3N[5]=-crd3N54;
    rd3N[6]=-crd3N64*k03;
    rd3N[7]=crd3N60*(crd3N1*crd3N23 + crd3N1*crd3N25 - crd3N1*crd3N28 + crd3N23*crd3N4 - crd3N25*crd3N4 + crd3N28*crd3N4 - crd3N35 - crd3N36 + 48.0*crd3N50 - crd3N55 + crd3N56 - crd3N62*crd3N8);
    rd3N[8]=crd3N64;

    rd3N /= std::pow(J, 3);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetFourthDerivativesShapeFunctionsValues(
    GlobalSizeVector& rd4N,
    const double J,
    const double xi
    )
{
    // Nodal values of the Jacobian
    const double J1 = GetJacobian(-1.0);
    const double J2 = GetJacobian( 1.0);
    const double J3 = GetJacobian( 0.0);

    // Nodal values of the geometric curvatures
    const double k01 = GetGeometryCurvature(J1, -1.0);
    const double k02 = GetGeometryCurvature(J2,  1.0);
    const double k03 = GetGeometryCurvature(J3,  0.0);

    const double k_s = GetBendingShearStiffnessRatio();

    const double crd4N0 = std::pow(k_s, 2);
    const double crd4N1 = 4320.0*crd4N0;
    const double crd4N2 = std::pow(J2, 2);
    const double crd4N3 = std::pow(J3, 2);
    const double crd4N4 = 24*crd4N3;
    const double crd4N5 = 288.0*k_s;
    const double crd4N6 = crd4N3*k_s;
    const double crd4N7 = 648.0*crd4N6;
    const double crd4N8 = crd4N2*xi;
    const double crd4N9 = 120*crd4N3;
    const double crd4N10 = 720.0*k_s;
    const double crd4N11 = crd4N0*xi;
    const double crd4N12 = 8640.0*crd4N11;
    const double crd4N13 = crd4N6*xi;
    const double crd4N14 = 1440.0*crd4N13;
    const double crd4N15 = -crd4N12 + crd4N14;
    const double crd4N16 = std::pow(k_s, 3);
    const double crd4N17 = std::pow(J1, 2);
    const double crd4N18 = 648.0*crd4N0;
    const double crd4N19 = crd4N0*crd4N3;
    const double crd4N20 = crd4N17*crd4N2;
    const double crd4N21 = crd4N20*k_s;
    const double crd4N22 = 78.0*crd4N6;
    const double crd4N23 = crd4N20*crd4N3;
    const double crd4N24 = 1.0/(-8640.0*crd4N16 - crd4N17*crd4N18 + crd4N17*crd4N22 - crd4N18*crd4N2 + 1296.0*crd4N19 + crd4N2*crd4N22 - 48.0*crd4N21 + 4.0*crd4N23);
    const double crd4N25 = std::pow(J1, 3)*crd4N24*(crd4N1 - crd4N10*crd4N8 + crd4N15 - crd4N2*crd4N4 + crd4N2*crd4N5 - crd4N7 + crd4N8*crd4N9);
    const double crd4N26 = 1296.0*crd4N0;
    const double crd4N27 = 156.0*crd4N6;
    const double crd4N28 = 1.0/(-17280.0*crd4N16 - crd4N17*crd4N26 + crd4N17*crd4N27 + 2592.0*crd4N19 - crd4N2*crd4N26 + crd4N2*crd4N27 - 96.0*crd4N21 + 8.0*crd4N23);
    const double crd4N29 = crd4N1*crd4N2;
    const double crd4N30 = 12960.0*crd4N0;
    const double crd4N31 = 1152.0*crd4N21;
    const double crd4N32 = crd4N2*crd4N6;
    const double crd4N33 = crd4N12*crd4N2;
    const double crd4N34 = 96.0*crd4N23;
    const double crd4N35 = crd4N17*crd4N6;
    const double crd4N36 = 25920.0*crd4N11;
    const double crd4N37 = crd4N6*crd4N8;
    const double crd4N38 = crd4N13*crd4N17;
    const double crd4N39 = 17280.0*crd4N19*xi - 1440.0*crd4N21*xi + 720.0*crd4N23*xi;
    const double crd4N40 = crd4N17*xi;
    const double crd4N41 = std::pow(J2, 3)*crd4N24*(-crd4N1 - crd4N10*crd4N40 + crd4N15 + crd4N17*crd4N4 - crd4N17*crd4N5 + crd4N40*crd4N9 + crd4N7);
    const double crd4N42 = -crd4N1*crd4N17 + crd4N12*crd4N17;
    const double crd4N43 = 576.0*k_s;
    const double crd4N44 = crd4N2*crd4N43;
    const double crd4N45 = 1800.0*k_s;
    const double crd4N46 = 324.0*crd4N0;
    const double crd4N47 = 1.0/(-4320.0*crd4N16 - crd4N17*crd4N46 + 648.0*crd4N19 - crd4N2*crd4N46 - 24.0*crd4N21 + 2.0*crd4N23 + 39.0*crd4N32 + 39.0*crd4N35);
    const double crd4N48 = std::pow(J3, 3)*crd4N47*(crd4N12 - crd4N17*crd4N43 + 240.0*crd4N20*xi + crd4N40*crd4N45 + crd4N44 + crd4N45*crd4N8);
    rd4N[0]=-crd4N25*k01;
    rd4N[1]=crd4N28*(crd4N17*crd4N30 - crd4N17*crd4N36 + crd4N29 + crd4N31 + 504.0*crd4N32 + crd4N33 - crd4N34 - 3096.0*crd4N35 + 2160.0*crd4N37 + 7920.0*crd4N38 + crd4N39);
    rd4N[2]=crd4N25;
    rd4N[3]=-crd4N41*k02;
    rd4N[4]=-crd4N28*(-crd4N2*crd4N30 - crd4N2*crd4N36 - crd4N31 + 3096.0*crd4N32 + crd4N34 - 504.0*crd4N35 + 7920.0*crd4N37 + 2160.0*crd4N38 + crd4N39 + crd4N42);
    rd4N[5]=crd4N41;
    rd4N[6]=-crd4N48*k03;
    rd4N[7]=crd4N47*(-crd4N14*crd4N17 - crd4N17*crd4N44 + crd4N17*crd4N7 + crd4N2*crd4N7 + 48.0*crd4N23 - crd4N29 - crd4N33 + 1440.0*crd4N37 + crd4N42);
    rd4N[8]=crd4N48;

    rd4N /= std::pow(J, 4);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetNThetaShapeFunctionsValues(
    GlobalSizeVector& rN,
    const double J,
    const double xi
    )
{
    GlobalSizeVector dN, d3N, Nu;
    GetFirstDerivativesShapeFunctionsValues(dN,  J, xi);
    GetThirdDerivativesShapeFunctionsValues(d3N, J, xi);
    GetNu0ShapeFunctionsValues(Nu, xi);

    const double k0 = GetGeometryCurvature(J, xi);
    const double k_s = GetBendingShearStiffnessRatio();
    // v' + ks * v''' + k0 * u
    noalias(rN) = dN + k_s * d3N + k0 * Nu;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetFirstDerivativesNThetaShapeFunctionsValues(
    GlobalSizeVector& rN,
    const double J,
    const double xi
    )
{
    GlobalSizeVector dN_dxi, d2N_dxi2;
    GetLocalFirstDerivativesNu0ShapeFunctionsValues (dN_dxi,   xi);
    GetLocalSecondDerivativesNu0ShapeFunctionsValues(d2N_dxi2, xi);
    const auto& r_geom = GetGeometry();

    double dx_dxi = 0.0;
    double dy_dxi = 0.0;

    double d2x_dxi2 = 0.0;
    double d2y_dxi2 = 0.0;

    for (IndexType i = 0; i < NumberOfNodes; ++i) {
        const IndexType u_coord = DoFperNode * i;
        const auto &r_coords_node = r_geom[i].GetInitialPosition();
        dx_dxi += r_coords_node[0] * dN_dxi[u_coord];
        dy_dxi += r_coords_node[1] * dN_dxi[u_coord];

        d2x_dxi2 += r_coords_node[0] * d2N_dxi2[u_coord];
        d2y_dxi2 += r_coords_node[1] * d2N_dxi2[u_coord];
    }

    const double dk0_ds = -3.0 / std::pow(J, 6) * (dx_dxi * d2y_dxi2 - dy_dxi * d2x_dxi2) * (dx_dxi * d2x_dxi2 + dy_dxi * d2y_dxi2);

    GlobalSizeVector d2N, d4N, dNu, Nu;
    GetSecondDerivativesShapeFunctionsValues (d2N,  J, xi);
    GetFourthDerivativesShapeFunctionsValues (d4N, J, xi);
    GetNu0ShapeFunctionsValues(Nu, xi);
    GetFirstDerivativesNu0ShapeFunctionsValues(dNu, J, xi);

    const double k0 = GetGeometryCurvature(J, xi);
    const double k_s = GetBendingShearStiffnessRatio();

    // v'' + ks * v'''' + k0 * u' + dk0*u
    noalias(rN) = d2N + k_s * d4N + k0 * dNu + dk0_ds * Nu;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetNu0ShapeFunctionsValues(
    GlobalSizeVector& rNu,
    const double xi
    )
{
    rNu.clear();
    rNu[0] = 0.5 * xi * (xi - 1.0);
    rNu[3] = 0.5 * xi * (xi + 1.0);
    rNu[6] = 1.0 - std::pow(xi, 2);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetFirstDerivativesNu0ShapeFunctionsValues(
    GlobalSizeVector& rNu,
    const double J,
    const double xi
    )
{
    GetLocalFirstDerivativesNu0ShapeFunctionsValues(rNu, xi);
    rNu /= J;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetLocalFirstDerivativesNu0ShapeFunctionsValues(
    GlobalSizeVector& rNu,
    const double xi
    )
{
    rNu.clear();
    rNu[0] = xi - 0.5;
    rNu[3] = xi + 0.5;
    rNu[6] = -2.0 * xi;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetLocalSecondDerivativesNu0ShapeFunctionsValues(
    GlobalSizeVector& rNu,
    const double xi
    )
{
    rNu.clear();
    rNu[0] = 1.0;
    rNu[3] = 1.0;
    rNu[6] = -2.0;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetSecondDerivativesNu0ShapeFunctionsValues(
    GlobalSizeVector& rNu,
    const double J,
    const double xi
    )
{
    GetLocalSecondDerivativesNu0ShapeFunctionsValues(rNu, xi);
    rNu /= std::pow(J, 2);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetNodalValuesVector(
    GlobalSizeVector& rNodalValues,
    const double angle1,
    const double angle2,
    const double angle3
    )
{
    const auto &r_geom = GetGeometry();
    BoundedMatrix<double, 3, 3> T1, T2, T3;
    GlobalSizeVector global_values;
    BoundedMatrix<double, 9, 9> global_size_T;
    StructuralMechanicsElementUtilities::BuildRotationMatrixForBeam(T1, angle1);
    StructuralMechanicsElementUtilities::BuildRotationMatrixForBeam(T2, angle2);
    StructuralMechanicsElementUtilities::BuildRotationMatrixForBeam(T3, angle3);

    global_size_T.clear();

    global_size_T(0, 0) = T1(0, 0);
    global_size_T(0, 1) = T1(0, 1);
    global_size_T(0, 2) = T1(0, 2);

    global_size_T(1, 0) = T1(1, 0);
    global_size_T(1, 1) = T1(1, 1);
    global_size_T(1, 2) = T1(1, 2);

    global_size_T(2, 0) = T1(2, 0);
    global_size_T(2, 1) = T1(2, 1);
    global_size_T(2, 2) = T1(2, 2);

    global_size_T(3, 3) = T2(0, 0);
    global_size_T(3, 4) = T2(0, 1);
    global_size_T(3, 5) = T2(0, 2);

    global_size_T(4, 3) = T2(1, 0);
    global_size_T(4, 4) = T2(1, 1);
    global_size_T(4, 5) = T2(1, 2);

    global_size_T(5, 3) = T2(2, 0);
    global_size_T(5, 4) = T2(2, 1);
    global_size_T(5, 5) = T2(2, 2);

    global_size_T(6, 6) = T3(0, 0);
    global_size_T(6, 7) = T3(0, 1);
    global_size_T(6, 8) = T3(0, 2);

    global_size_T(7, 6) = T3(1, 0);
    global_size_T(7, 7) = T3(1, 1);
    global_size_T(7, 8) = T3(1, 2);

    global_size_T(8, 6) = T3(2, 0);
    global_size_T(8, 7) = T3(2, 1);
    global_size_T(8, 8) = T3(2, 2);

    const auto& r_displ_0 = r_geom[0].FastGetSolutionStepValue(DISPLACEMENT);
    const auto& r_displ_1 = r_geom[1].FastGetSolutionStepValue(DISPLACEMENT);
    const auto& r_displ_2 = r_geom[2].FastGetSolutionStepValue(DISPLACEMENT);

    global_values[0] = r_displ_0[0];
    global_values[1] = r_displ_0[1];
    global_values[2] = r_geom[0].FastGetSolutionStepValue(ROTATION_Z);

    global_values[3] = r_displ_1[0];
    global_values[4] = r_displ_1[1];
    global_values[5] = r_geom[1].FastGetSolutionStepValue(ROTATION_Z);

    global_values[6] = r_displ_2[0];
    global_values[7] = r_displ_2[1];
    global_values[8] = r_geom[2].FastGetSolutionStepValue(ROTATION_Z);

    // We rotate to local tangent s,y axes
    noalias(rNodalValues) = prod(trans(global_size_T), global_values);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::CalculateGeneralizedStrainsVector(
    VectorType& rStrain,
    const double J,
    const double xi,
    const GlobalSizeVector &rNodalValues
    )
{
    rStrain[0] = CalculateAxialStrain     (J, xi, rNodalValues); // El
    rStrain[1] = CalculateBendingCurvature(J, xi, rNodalValues); // Kappa
    rStrain[2] = CalculateShearStrain     (J, xi, rNodalValues); // Gamma_xy
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoCurvedBeamElement2D3N::CalculateAxialStrain(
    const double J,
    const double xi,
    const GlobalSizeVector& rNodalValues
    )
{
    GlobalSizeVector dNu, N, dN;
    GetFirstDerivativesNu0ShapeFunctionsValues(dNu, J, xi);
    GetShapeFunctionsValues(N, J, xi);
    GetFirstDerivativesShapeFunctionsValues(dN, J, xi);
    const double k0 = GetGeometryCurvature(J, xi);
    return inner_prod(dNu - k0 * N, rNodalValues) + 0.5 * std::pow(inner_prod(dNu, rNodalValues), 2) + 0.5 * std::pow(inner_prod(dN, rNodalValues), 2);
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoCurvedBeamElement2D3N::CalculateShearStrain(
    const double J,
    const double xi,
    const GlobalSizeVector& rNodalValues
    )
{
    GlobalSizeVector dN, Nu, Ntheta;
    GetFirstDerivativesShapeFunctionsValues(dN, J, xi);
    GetNu0ShapeFunctionsValues(Nu, xi);
    GetNThetaShapeFunctionsValues(Ntheta, J, xi);
    const double k0 = GetGeometryCurvature(J, xi);

    return inner_prod(dN + k0 * Nu - Ntheta, rNodalValues);
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoCurvedBeamElement2D3N::CalculateBendingCurvature(
    const double J,
    const double xi,
    const GlobalSizeVector& rNodalValues
    )
{
    GlobalSizeVector dNtheta;
    GetFirstDerivativesNThetaShapeFunctionsValues(dNtheta, J, xi);
    return inner_prod(dNtheta, rNodalValues);
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3> LinearTimoshenkoCurvedBeamElement2D3N::GetLocalAxesBodyForce(
    const Element &rElement,
    const GeometryType::IntegrationPointsArrayType &rIntegrationPoints,
    const IndexType PointNumber,
    const double angle
    )
{
    const auto body_force = StructuralMechanicsElementUtilities::GetBodyForce(*this, rIntegrationPoints, PointNumber);

    const double c = std::cos(angle);
    const double s = std::sin(angle);
    array_1d<double, 3> local_body_force = ZeroVector(3);
    local_body_force[0] = c * body_force[0] + s * body_force[1];
    local_body_force[1] = -s * body_force[0] + c * body_force[1];
    return local_body_force;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::CalculateLocalSystem(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY;
    const auto &r_props = GetProperties();
    const auto &r_geometry = GetGeometry();

    if (rLHS.size1() != SystemSize || rLHS.size2() != SystemSize) {
        rLHS.resize(SystemSize, SystemSize, false);
    }
    noalias(rLHS) = ZeroMatrix(SystemSize, SystemSize);

    if (rRHS.size() != SystemSize) {
        rRHS.resize(SystemSize, false);
    }
    noalias(rRHS) = ZeroVector(SystemSize);

    const auto& integration_points = IntegrationPoints(GetIntegrationMethod());

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    const double area = r_props[CROSS_AREA];

    // Let's initialize the cl values
    VectorType strain_vector(StrainSize), stress_vector(StrainSize);
    MatrixType constitutive_matrix(StrainSize, StrainSize);
    strain_vector.clear();
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);
    GlobalSizeVector nodal_values, aux_array, d_el_du;

    const double angle1 = GetAngle(-1.0);
    const double angle2 = GetAngle(1.0);
    const double angle3 = GetAngle(0.0);

    GlobalSizeVector dNu, dN_theta, N_shape, Nu, N_s, N_theta, dN_shape;

    // Loop over the integration points
    for (SizeType IP = 0; IP < integration_points.size(); ++IP) {

        const double xi     = integration_points[IP].X();
        const double weight = integration_points[IP].Weight();
        const double J  = GetJacobian(xi);
        const double k0 = GetGeometryCurvature(J, xi);
        const double jacobian_weight = weight * J;

        GetNodalValuesVector(nodal_values, angle1, angle2, angle3);

        // We fill the strain vector with El, kappa and Gamma_sy
        CalculateGeneralizedStrainsVector(strain_vector, J, xi, nodal_values);

        // We fill the resulting stresses
        mConstitutiveLawVector[IP]->CalculateMaterialResponseCauchy(cl_values);
        const Vector &r_generalized_stresses = cl_values.GetStressVector();
        const double N = r_generalized_stresses[0];
        const double M = r_generalized_stresses[1];
        const double V = r_generalized_stresses[2];

        // And its derivatives with respect to the gen. strains
        const MatrixType& r_constitutive_matrix = cl_values.GetConstitutiveMatrix();
        const double dN_dEl    = r_constitutive_matrix(0, 0);
        const double dM_dkappa = r_constitutive_matrix(1, 1);
        const double dV_dgamma = r_constitutive_matrix(2, 2);

        // Let's compute the required arrays
        GetFirstDerivativesNu0ShapeFunctionsValues(dNu, J, xi);
        GetFirstDerivativesNThetaShapeFunctionsValues(dN_theta, J, xi);
        GetNThetaShapeFunctionsValues(N_theta, J, xi);
        GetFirstDerivativesShapeFunctionsValues(dN_shape, J, xi);
        GetShapeFunctionsValues(N_shape, J, xi);
        GetNu0ShapeFunctionsValues(Nu, xi);
        noalias(N_s) = dN_shape - N_theta + k0 * Nu;
        noalias(aux_array) = dNu - k0 * N_shape;

        const double du = inner_prod(dNu, nodal_values);
        const double dv = inner_prod(dN_shape, nodal_values);
        noalias(d_el_du) = aux_array + du * dNu + dv * dN_shape;

        // Axial contributions
        noalias(rRHS) -= d_el_du * N * jacobian_weight;
        noalias(rLHS) += outer_prod(d_el_du, d_el_du) * dN_dEl * jacobian_weight;
        noalias(rLHS) += outer_prod(dNu, dNu) * N * jacobian_weight;
        noalias(rLHS) += outer_prod(dN_shape, dN_shape) * N * jacobian_weight;

        // Bending contributions
        noalias(rLHS) += outer_prod(dN_theta, dN_theta) * dM_dkappa * jacobian_weight;
        noalias(rRHS) -= dN_theta * M * jacobian_weight;

        // Shear contributions
        noalias(rLHS) += outer_prod(N_s, N_s) * dV_dgamma * jacobian_weight;
        noalias(rRHS) -= N_s * V * jacobian_weight;

        // Now we add the body forces contributions
        auto local_body_forces = GetLocalAxesBodyForce(*this, integration_points, IP, GetAngle(xi));
        noalias(rRHS) += Nu      * local_body_forces[0] * jacobian_weight * area;
        noalias(rRHS) += N_shape * local_body_forces[1] * jacobian_weight * area;

    } // IP loop
    RotateAll(rLHS, rRHS, angle1, angle2, angle3);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::CalculateLeftHandSide(
    MatrixType& rLHS,
    const ProcessInfo& rProcessInfo
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::CalculateRightHandSide(
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY;
    const auto &r_props = GetProperties();
    const auto &r_geometry = GetGeometry();

    if (rRHS.size() != SystemSize) {
        rRHS.resize(SystemSize, false);
    }
    noalias(rRHS) = ZeroVector(SystemSize);

    const auto& integration_points = IntegrationPoints(GetIntegrationMethod());

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    const double area = r_props[CROSS_AREA];

    // Let's initialize the cl values
    VectorType strain_vector(StrainSize), stress_vector(StrainSize);
    MatrixType constitutive_matrix(StrainSize, StrainSize);
    strain_vector.clear();
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);
    GlobalSizeVector nodal_values, aux_array, d_el_du;

    GlobalSizeVector dNu, dN_theta, N_shape, Nu, N_s, N_theta, dN_shape;

    const double angle1 = GetAngle(-1.0);
    const double angle2 = GetAngle(1.0);
    const double angle3 = GetAngle(0.0);

    // Loop over the integration points
    for (SizeType IP = 0; IP < integration_points.size(); ++IP) {

        const double xi     = integration_points[IP].X();
        const double weight = integration_points[IP].Weight();
        const double J  = GetJacobian(xi);
        const double k0 = GetGeometryCurvature(J, xi);
        const double jacobian_weight = weight * J;
        // const double angle = GetAngle(xi);

        GetNodalValuesVector(nodal_values, angle1, angle2, angle3);

        // We fill the strain vector with El, kappa and Gamma_sy
        CalculateGeneralizedStrainsVector(strain_vector, J, xi, nodal_values);

        // We fill the resulting stresses
        mConstitutiveLawVector[IP]->CalculateMaterialResponseCauchy(cl_values);
        const Vector &r_generalized_stresses = cl_values.GetStressVector();
        const double N = r_generalized_stresses[0];
        const double M = r_generalized_stresses[1];
        const double V = r_generalized_stresses[2];

        // And its derivatives with respect to the gen. strains
        const MatrixType& r_constitutive_matrix = cl_values.GetConstitutiveMatrix();
        const double dN_dEl    = r_constitutive_matrix(0, 0);
        const double dM_dkappa = r_constitutive_matrix(1, 1);
        const double dV_dgamma = r_constitutive_matrix(2, 2);

        // Let's compute the required arrays
        GetFirstDerivativesNu0ShapeFunctionsValues(dNu, J, xi);
        GetFirstDerivativesNThetaShapeFunctionsValues(dN_theta, J, xi);
        GetNThetaShapeFunctionsValues(N_theta, J, xi);
        GetFirstDerivativesShapeFunctionsValues(dN_shape, J, xi);
        GetShapeFunctionsValues(N_shape, J, xi);
        GetNu0ShapeFunctionsValues(Nu, xi);
        noalias(N_s) = dN_shape - N_theta + k0 * Nu;
        noalias(aux_array) = dNu - k0 * N_shape;

        const double du = inner_prod(dNu, nodal_values);
        const double dv = inner_prod(dN_shape, nodal_values);
        noalias(d_el_du) = aux_array + du * dNu + dv * dN_shape;

        // Axial contributions
        noalias(rRHS) -= d_el_du * N * jacobian_weight;

        // Bending contributions
        noalias(rRHS) -= dN_theta * M * jacobian_weight;

        // Shear contributions
        noalias(rRHS) -= N_s * V * jacobian_weight;

        // Now we add the body forces contributions
        auto local_body_forces = GetLocalAxesBodyForce(*this, integration_points, IP, GetAngle(xi));
        noalias(rRHS) += Nu      * local_body_forces[0] * jacobian_weight * area;
        noalias(rRHS) += N_shape * local_body_forces[1] * jacobian_weight * area;

    } // IP loop
    RotateRHS(rRHS, angle1, angle2, angle3);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoCurvedBeamElement2D3N::GetAngle(
    const double xi
    )
{
    GlobalSizeVector dN_dxi;
    GetLocalFirstDerivativesNu0ShapeFunctionsValues(dN_dxi, xi);
    const auto r_geom = GetGeometry();

    double dx_dxi = 0.0;
    double dy_dxi = 0.0;

    for (IndexType i = 0; i < NumberOfNodes; ++i) {
        const IndexType u_coord = DoFperNode * i;
        const auto &r_coords_node = r_geom[i].GetInitialPosition();
        dx_dxi += r_coords_node[0] * dN_dxi[u_coord];
        dy_dxi += r_coords_node[1] * dN_dxi[u_coord];
    }
    return std::atan2(dy_dxi, dx_dxi);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::RotateLHS(
    MatrixType& rLHS,
    const double angle1
)
{
    // BoundedMatrix<double, 3, 3> T, Tt;
    // BoundedMatrix<double, 9, 9> global_size_T, aux_product;
    // StructuralMechanicsElementUtilities::BuildRotationMatrixForBeam(T, angle);
    // StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D3NBeam(T, global_size_T);
    // noalias(aux_product) = prod(rLHS, trans(global_size_T));
    // noalias(rLHS) = prod(global_size_T, aux_product);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::RotateAll(
    MatrixType& rLHS,
    VectorType& rRHS,
    const double angle1,
    const double angle2,
    const double angle3
)
{
    BoundedMatrix<double, 3, 3> T1, T2, T3;
    BoundedMatrix<double, 9, 9> global_size_T, aux_product;
    BoundedVector<double, 9> local_rhs;
    StructuralMechanicsElementUtilities::BuildRotationMatrixForBeam(T1, angle1);
    StructuralMechanicsElementUtilities::BuildRotationMatrixForBeam(T2, angle2);
    StructuralMechanicsElementUtilities::BuildRotationMatrixForBeam(T3, angle3);
    // StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D3NBeam(T, global_size_T);

    global_size_T.clear();

    global_size_T(0, 0) = T1(0, 0);
    global_size_T(0, 1) = T1(0, 1);
    global_size_T(0, 2) = T1(0, 2);

    global_size_T(1, 0) = T1(1, 0);
    global_size_T(1, 1) = T1(1, 1);
    global_size_T(1, 2) = T1(1, 2);

    global_size_T(2, 0) = T1(2, 0);
    global_size_T(2, 1) = T1(2, 1);
    global_size_T(2, 2) = T1(2, 2);

    global_size_T(3, 3) = T2(0, 0);
    global_size_T(3, 4) = T2(0, 1);
    global_size_T(3, 5) = T2(0, 2);

    global_size_T(4, 3) = T2(1, 0);
    global_size_T(4, 4) = T2(1, 1);
    global_size_T(4, 5) = T2(1, 2);

    global_size_T(5, 3) = T2(2, 0);
    global_size_T(5, 4) = T2(2, 1);
    global_size_T(5, 5) = T2(2, 2);

    global_size_T(6, 6) = T3(0, 0);
    global_size_T(6, 7) = T3(0, 1);
    global_size_T(6, 8) = T3(0, 2);

    global_size_T(7, 6) = T3(1, 0);
    global_size_T(7, 7) = T3(1, 1);
    global_size_T(7, 8) = T3(1, 2);

    global_size_T(8, 6) = T3(2, 0);
    global_size_T(8, 7) = T3(2, 1);
    global_size_T(8, 8) = T3(2, 2);

    noalias(local_rhs) = rRHS;
    noalias(rRHS) = prod(global_size_T, local_rhs);

    noalias(aux_product) = prod(rLHS, trans(global_size_T));
    noalias(rLHS) = prod(global_size_T, aux_product);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::RotateRHS(
    VectorType& rRHS,
    const double angle1,
    const double angle2,
    const double angle3
)
{
    BoundedMatrix<double, 3, 3> T1, T2, T3;
    BoundedMatrix<double, 9, 9> global_size_T, aux_product;
    BoundedVector<double, 9> local_rhs;
    StructuralMechanicsElementUtilities::BuildRotationMatrixForBeam(T1, angle1);
    StructuralMechanicsElementUtilities::BuildRotationMatrixForBeam(T2, angle2);
    StructuralMechanicsElementUtilities::BuildRotationMatrixForBeam(T3, angle3);
    // StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D3NBeam(T, global_size_T);

    global_size_T.clear();

    global_size_T(0, 0) = T1(0, 0);
    global_size_T(0, 1) = T1(0, 1);
    global_size_T(0, 2) = T1(0, 2);

    global_size_T(1, 0) = T1(1, 0);
    global_size_T(1, 1) = T1(1, 1);
    global_size_T(1, 2) = T1(1, 2);

    global_size_T(2, 0) = T1(2, 0);
    global_size_T(2, 1) = T1(2, 1);
    global_size_T(2, 2) = T1(2, 2);

    global_size_T(3, 3) = T2(0, 0);
    global_size_T(3, 4) = T2(0, 1);
    global_size_T(3, 5) = T2(0, 2);

    global_size_T(4, 3) = T2(1, 0);
    global_size_T(4, 4) = T2(1, 1);
    global_size_T(4, 5) = T2(1, 2);

    global_size_T(5, 3) = T2(2, 0);
    global_size_T(5, 4) = T2(2, 1);
    global_size_T(5, 5) = T2(2, 2);

    global_size_T(6, 6) = T3(0, 0);
    global_size_T(6, 7) = T3(0, 1);
    global_size_T(6, 8) = T3(0, 2);

    global_size_T(7, 6) = T3(1, 0);
    global_size_T(7, 7) = T3(1, 1);
    global_size_T(7, 8) = T3(1, 2);

    global_size_T(8, 6) = T3(2, 0);
    global_size_T(8, 7) = T3(2, 1);
    global_size_T(8, 8) = T3(2, 2);


    noalias(local_rhs) = rRHS;
    noalias(rRHS) = prod(global_size_T, local_rhs);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rProcessInfo
    )
{
    const auto& integration_points = IntegrationPoints(GetIntegrationMethod());
    rOutput.resize(integration_points.size());
    const auto &r_props = GetProperties();

    if (rVariable == AXIAL_FORCE ||
        rVariable == BENDING_MOMENT ||
        rVariable == SHEAR_FORCE ||
        rVariable == INITIAL_GEOMETRIC_CURVATURE)
    {
        const auto &r_geometry = GetGeometry();

        ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
        auto &r_cl_options = cl_values.GetOptions();
        r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
        r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        // Let's initialize the cl values
        VectorType strain_vector(StrainSize), stress_vector(StrainSize);
        strain_vector.clear();
        cl_values.SetStrainVector(strain_vector);
        cl_values.SetStressVector(stress_vector);
        GlobalSizeVector nodal_values(SystemSize);

        const double angle1 = GetAngle(-1.0);
        const double angle2 = GetAngle(1.0);
        const double angle3 = GetAngle(0.0);

        // Loop over the integration points
        for (SizeType IP = 0; IP < integration_points.size(); ++IP) {
            const double xi = integration_points[IP].X();
            // const double angle = GetAngle(xi);
            GetNodalValuesVector(nodal_values, angle1, angle2, angle3);
            const double J = GetJacobian(xi);

            CalculateGeneralizedStrainsVector(strain_vector, J, xi, nodal_values);

            mConstitutiveLawVector[IP]->CalculateMaterialResponseCauchy(cl_values);
            const Vector &r_generalized_stresses = cl_values.GetStressVector();
            if (rVariable == AXIAL_FORCE)
                rOutput[IP] = r_generalized_stresses[0];
            else if (rVariable == BENDING_MOMENT)
                rOutput[IP] = r_generalized_stresses[1];
            else if (rVariable == SHEAR_FORCE)
                rOutput[IP] = r_generalized_stresses[2];
            else
                rOutput[IP] = GetGeometryCurvature(J, xi);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::CalculateOnIntegrationPoints(
    const Variable<ConstitutiveLaw::Pointer>& rVariable,
    std::vector<ConstitutiveLaw::Pointer>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rVariable == CONSTITUTIVE_LAW) {
        const SizeType integration_points_number = mConstitutiveLawVector.size();
        if (rValues.size() != integration_points_number) {
            rValues.resize(integration_points_number);
        }
        for (IndexType point_number = 0; point_number < integration_points_number; ++point_number) {
            rValues[point_number] = mConstitutiveLawVector[point_number];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

int LinearTimoshenkoCurvedBeamElement2D3N::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    return mConstitutiveLawVector[0]->Check(GetProperties(), GetGeometry(), rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    int IntMethod = int(this->GetIntegrationMethod());
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos
