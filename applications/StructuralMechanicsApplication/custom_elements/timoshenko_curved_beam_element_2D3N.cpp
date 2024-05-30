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
#include "timoshenko_curved_beam_element_2D3N.h"
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
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
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

    rN.clear();

    // Nodal values of the Jacobian
    const double J1 = GetJacobian(-1.0);
    const double J2 = GetJacobian( 1.0);
    const double J3 = GetJacobian( 0.0);
    const double k_s = GetBendingShearStiffnessRatio();

    const double cr_N0 = std::pow(J2, 2);
    const double cr_N1 = std::pow(k_s, 2);
    const double cr_N2 = 576.0*cr_N1;
    const double cr_N3 = std::pow(J1, 2);
    const double cr_N4 = 864.0*cr_N1;
    const double cr_N5 = std::pow(k_s, 3);
    const double cr_N6 = 8640.0*cr_N5;
    const double cr_N7 = cr_N6*xi;
    const double cr_N8 = std::pow(xi, 4);
    const double cr_N9 = 72.0*cr_N1;
    const double cr_N10 = cr_N8*cr_N9;
    const double cr_N11 = std::pow(xi, 3);
    const double cr_N12 = 180.0*cr_N1;
    const double cr_N13 = cr_N11*cr_N12;
    const double cr_N14 = std::pow(J3, 2);
    const double cr_N15 = cr_N1*cr_N14;
    const double cr_N16 = 1296.0*cr_N15;
    const double cr_N17 = cr_N16*xi;
    const double cr_N18 = cr_N1*cr_N3;
    const double cr_N19 = 540.0*cr_N11;
    const double cr_N20 = 216.0*cr_N8;
    const double cr_N21 = cr_N0*xi;
    const double cr_N22 = cr_N1*cr_N21;
    const double cr_N23 = cr_N18*xi;
    const double cr_N24 = cr_N0*cr_N14;
    const double cr_N25 = cr_N24*xi;
    const double cr_N26 = 8.0*cr_N25*cr_N3;
    const double cr_N27 = cr_N0*cr_N3;
    const double cr_N28 = cr_N27*k_s;
    const double cr_N29 = 48.0*cr_N28;
    const double cr_N30 = cr_N11*cr_N29;
    const double cr_N31 = cr_N24*cr_N8;
    const double cr_N32 = 18.0*k_s;
    const double cr_N33 = cr_N11*cr_N24;
    const double cr_N34 = cr_N33*k_s;
    const double cr_N35 = cr_N14*cr_N3;
    const double cr_N36 = cr_N35*cr_N8;
    const double cr_N37 = 66.0*k_s;
    const double cr_N38 = 57.0*k_s;
    const double cr_N39 = cr_N35*xi;
    const double cr_N40 = 207.0*k_s;
    const double cr_N41 = 4.0*cr_N3;
    const double cr_N42 = cr_N33*cr_N41;
    const double cr_N43 = 96.0*cr_N28;
    const double cr_N44 = cr_N43*xi;
    const double cr_N45 = std::pow(xi, 2);
    const double cr_N46 = cr_N24*cr_N45;
    const double cr_N47 = 96.0*k_s;
    const double cr_N48 = cr_N35*cr_N45;
    const double cr_N49 = 144.0*k_s;
    const double cr_N50 = cr_N11*cr_N35;
    const double cr_N51 = cr_N50*k_s;
    const double cr_N52 = cr_N3*cr_N31;
    const double cr_N53 = -1440.0*cr_N15*cr_N45 + 144.0*cr_N15*cr_N8 - 12.0*cr_N28*cr_N8 + 60.0*cr_N28 - 10.0*cr_N3*cr_N46 + 6.0*cr_N52 + cr_N6;
    const double cr_N54 = 1296.0*cr_N1;
    const double cr_N55 = 156.0*k_s;
    const double cr_N56 = cr_N24*cr_N3;
    const double cr_N57 = xi/(-cr_N0*cr_N54 + 2592.0*cr_N15 + cr_N24*cr_N55 - cr_N3*cr_N54 + cr_N35*cr_N55 - cr_N43 - 17280.0*cr_N5 + 8.0*cr_N56);
    const double cr_N58 = 6.0*k_s;
    const double cr_N59 = cr_N0*cr_N58;
    const double cr_N60 = cr_N12*xi;
    const double cr_N61 = 12.0*k_s;
    const double cr_N62 = cr_N11*cr_N61;
    const double cr_N63 = k_s*xi;
    const double cr_N64 = 27.0*cr_N14;
    const double cr_N65 = cr_N63*cr_N64;
    const double cr_N66 = cr_N11*k_s;
    const double cr_N67 = cr_N64*cr_N66;
    const double cr_N68 = cr_N14*cr_N61;
    const double cr_N69 = -cr_N10 - cr_N45*cr_N68 + cr_N68*cr_N8 + cr_N9;
    const double cr_N70 = 648.0*cr_N1;
    const double cr_N71 = 78.0*k_s;
    const double cr_N72 = xi/(-cr_N0*cr_N70 + cr_N16 + cr_N24*cr_N71 - cr_N29 - cr_N3*cr_N70 + cr_N35*cr_N71 + 4.0*cr_N56 - cr_N6);
    const double cr_N73 = cr_N0*cr_N1;
    const double cr_N74 = cr_N3*cr_N58;
    const double cr_N75 = cr_N3*xi;
    const double cr_N76 = 4320.0*cr_N5;
    const double cr_N77 = 324.0*cr_N1;
    const double cr_N78 = cr_N14*cr_N70;
    const double cr_N79 = 24.0*cr_N28;
    const double cr_N80 = 39.0*k_s;
    const double cr_N81 = -cr_N0*cr_N77 + cr_N24*cr_N80 - cr_N3*cr_N77 + cr_N35*cr_N80 + 2.0*cr_N56 - cr_N76 + cr_N78 - cr_N79;
    const double cr_N82 = 1.0/cr_N81;
    const double cr_N83 = std::pow(xi, 5);
    const double cr_N84 = cr_N83*cr_N9;
    const double cr_N85 = cr_N1*cr_N45;
    const double cr_N86 = 504.0*cr_N85;
    const double cr_N87 = cr_N12*cr_N8;
    const double cr_N88 = cr_N61*cr_N83;
    const double cr_N89 = 27.0*k_s;
    const double cr_N90 = 24.0*cr_N3;
    const double cr_N91 = 2.0*cr_N27;
    const double cr_N92 = cr_N3*k_s;
    const double cr_N93 = 54.0*cr_N45;
    const double cr_N94 = 15.0*cr_N8;
    const double cr_N95 = cr_N0*k_s;

    rN[1]=cr_N57*(cr_N0*cr_N10 + cr_N0*cr_N13 + cr_N0*cr_N2 + cr_N17 + cr_N18*cr_N19 - cr_N18*cr_N20 - 828.0*cr_N22 - 1188.0*cr_N23 + cr_N25*cr_N38 + cr_N26 + cr_N3*cr_N4 + cr_N30 + cr_N31*cr_N32 + 21.0*cr_N34 + cr_N36*cr_N37 + cr_N39*cr_N40 - cr_N42 - cr_N44 - cr_N46*cr_N47 - cr_N48*cr_N49 - 129.0*cr_N51 + cr_N53 - cr_N7);
    rN[2]=std::pow(J1, 3)*cr_N72*(cr_N0*cr_N62 + cr_N13 - cr_N21*cr_N61 + cr_N25 + cr_N31 - cr_N33 - cr_N46 - cr_N59*cr_N8 + cr_N59 - cr_N60 + cr_N65 - cr_N67 + cr_N69);

    rN[4]=-cr_N57*(cr_N0*cr_N4 + cr_N10*cr_N3 - cr_N13*cr_N3 - cr_N17 - cr_N19*cr_N73 + cr_N2*cr_N3 - cr_N20*cr_N73 + 1188.0*cr_N22 + 828.0*cr_N23 - cr_N25*cr_N40 - cr_N26 - cr_N30 + cr_N31*cr_N37 + cr_N32*cr_N36 + 129.0*cr_N34 - cr_N38*cr_N39 + cr_N42 + cr_N44 - cr_N46*cr_N49 - cr_N47*cr_N48 - 21.0*cr_N51 + cr_N53 + cr_N7);
    rN[5]=std::pow(J2, 3)*cr_N72*(-cr_N13 - cr_N3*cr_N62 + cr_N36 - cr_N39 - cr_N48 + cr_N50 + cr_N60 + cr_N61*cr_N75 - cr_N65 + cr_N67 + cr_N69 - cr_N74*cr_N8 + cr_N74);

    rN[7]=cr_N82*(-cr_N0*cr_N84 + cr_N0*cr_N86 - cr_N0*cr_N87 + cr_N21*cr_N9 + cr_N24*cr_N88 + cr_N29*cr_N45 + cr_N3*cr_N84 + cr_N3*cr_N86 - cr_N3*cr_N87 + cr_N31*cr_N89 - cr_N33*cr_N61 - cr_N35*cr_N88 + cr_N36*cr_N89 - cr_N37*cr_N46 - cr_N37*cr_N48 - cr_N41*cr_N46 + cr_N45*cr_N76 - cr_N45*cr_N78 + cr_N50*cr_N61 + 2.0*cr_N52 - cr_N75*cr_N9 - cr_N79*cr_N8 + cr_N81);
    rN[8]=std::pow(J3, 3)*cr_N82*xi*(24.0*cr_N0*cr_N66 + cr_N0*cr_N80 + cr_N10 - 24.0*cr_N21*k_s - 4.0*cr_N27*cr_N45 + cr_N3*cr_N80 + cr_N63*cr_N90 - cr_N66*cr_N90 + cr_N70 + cr_N8*cr_N91 - 720.0*cr_N85 + cr_N91 - cr_N92*cr_N93 + cr_N92*cr_N94 - cr_N93*cr_N95 + cr_N94*cr_N95);
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
    const double k_s = GetBendingShearStiffnessRatio();

    rdN.clear();

    const double crdN0 = std::pow(k_s, 3);
    const double crdN1 = 17280.0*crdN0;
    const double crdN2 = std::pow(J1, 2);
    const double crdN3 = std::pow(k_s, 2);
    const double crdN4 = 1296.0*crdN3;
    const double crdN5 = std::pow(J2, 2);
    const double crdN6 = std::pow(J3, 2);
    const double crdN7 = crdN3*crdN6;
    const double crdN8 = 2592.0*crdN7;
    const double crdN9 = crdN2*crdN5;
    const double crdN10 = crdN9*k_s;
    const double crdN11 = 96.0*crdN10;
    const double crdN12 = crdN6*k_s;
    const double crdN13 = 156.0*crdN12;
    const double crdN14 = crdN6*crdN9;
    const double crdN15 = 8.0*crdN14;
    const double crdN16 = 1.0/(-crdN1 - crdN11 + crdN13*crdN2 + crdN13*crdN5 + crdN15 - crdN2*crdN4 - crdN4*crdN5 + crdN8);
    const double crdN17 = 576.0*crdN3;
    const double crdN18 = 864.0*crdN3;
    const double crdN19 = crdN1*xi;
    const double crdN20 = std::pow(xi, 4);
    const double crdN21 = 360.0*crdN3;
    const double crdN22 = crdN20*crdN21;
    const double crdN23 = crdN22*crdN5;
    const double crdN24 = std::pow(xi, 3);
    const double crdN25 = crdN24*crdN3;
    const double crdN26 = 720.0*crdN25;
    const double crdN27 = crdN26*crdN5;
    const double crdN28 = crdN8*xi;
    const double crdN29 = 2160.0*crdN25;
    const double crdN30 = crdN2*crdN3;
    const double crdN31 = 1080.0*crdN20;
    const double crdN32 = crdN5*xi;
    const double crdN33 = crdN3*crdN32;
    const double crdN34 = crdN30*xi;
    const double crdN35 = 16.0*crdN14;
    const double crdN36 = crdN35*xi;
    const double crdN37 = 192.0*crdN10;
    const double crdN38 = crdN24*crdN37;
    const double crdN39 = crdN12*crdN24;
    const double crdN40 = crdN39*crdN5;
    const double crdN41 = crdN12*crdN20;
    const double crdN42 = 90.0*crdN41;
    const double crdN43 = crdN12*crdN32;
    const double crdN44 = 330.0*crdN41;
    const double crdN45 = crdN12*xi;
    const double crdN46 = crdN2*crdN45;
    const double crdN47 = crdN24*crdN35;
    const double crdN48 = crdN37*xi;
    const double crdN49 = std::pow(xi, 2);
    const double crdN50 = crdN12*crdN49;
    const double crdN51 = 288.0*crdN50;
    const double crdN52 = 432.0*crdN50;
    const double crdN53 = crdN2*crdN39;
    const double crdN54 = 8640.0*crdN0;
    const double crdN55 = 60.0*crdN10;
    const double crdN56 = 30.0*crdN20;
    const double crdN57 = -30.0*crdN14*crdN49 + crdN14*crdN56 - crdN20*crdN55 + 720.0*crdN20*crdN7 - 4320.0*crdN49*crdN7 + crdN54 + crdN55;
    const double crdN58 = 648.0*crdN3;
    const double crdN59 = 1296.0*crdN7;
    const double crdN60 = 78.0*crdN12;
    const double crdN61 = 1.0/(-48.0*crdN10 + 4.0*crdN14 - crdN2*crdN58 + crdN2*crdN60 - crdN5*crdN58 + crdN5*crdN60 - crdN54 + crdN59);
    const double crdN62 = 6.0*k_s;
    const double crdN63 = crdN21*xi;
    const double crdN64 = crdN5*crdN6;
    const double crdN65 = 4*crdN24;
    const double crdN66 = 3*crdN49;
    const double crdN67 = 2*crdN6;
    const double crdN68 = 5*crdN20;
    const double crdN69 = crdN5*k_s;
    const double crdN70 = 48.0*crdN24;
    const double crdN71 = 54.0*crdN45;
    const double crdN72 = 24.0*k_s;
    const double crdN73 = 108.0*crdN39;
    const double crdN74 = 72.0*crdN3;
    const double crdN75 = 60.0*crdN41;
    const double crdN76 = 36.0*crdN50;
    const double crdN77 = -crdN22 + crdN74 + crdN75 - crdN76;
    const double crdN78 = crdN2*crdN22 - crdN2*crdN26;
    const double crdN79 = crdN2*crdN6;
    const double crdN80 = crdN2*xi;
    const double crdN81 = crdN2*k_s;
    const double crdN82 = 324.0*crdN3;
    const double crdN83 = 39.0*k_s;
    const double crdN84 = crdN2*crdN83;
    const double crdN85 = crdN5*crdN83;
    const double crdN86 = 2.0*crdN9;
    const double crdN87 = 1.0/(-4320.0*crdN0 - 24.0*crdN10 - crdN2*crdN82 - crdN5*crdN82 + crdN58*crdN6 + crdN6*crdN84 + crdN6*crdN85 + crdN6*crdN86);
    const double crdN88 = 48.0*k_s;
    const double crdN89 = 162.0*crdN49;
    const double crdN90 = 96.0*crdN24;
    const double crdN91 = 75.0*crdN20;

    rdN[1]=crdN16*(crdN17*crdN5 + crdN18*crdN2 - crdN19 + crdN2*crdN29 + crdN2*crdN44 - crdN2*crdN52 + crdN23 + crdN27 + crdN28 - crdN30*crdN31 - 1656.0*crdN33 - 2376.0*crdN34 + crdN36 + crdN38 + 84.0*crdN40 + crdN42*crdN5 + 114.0*crdN43 + 414.0*crdN46 - crdN47 - crdN48 - crdN5*crdN51 - 516.0*crdN53 + crdN57);
    rdN[2]=std::pow(J1, 3)*crdN61*(crdN26 + crdN32*crdN67 - crdN32*crdN72 + crdN5*crdN62 - crdN56*crdN69 - crdN63 - crdN64*crdN65 - crdN64*crdN66 + crdN64*crdN68 + crdN69*crdN70 + crdN71 - crdN73 + crdN77);

    rdN[4]=-crdN16*(crdN17*crdN2 + crdN18*crdN5 + crdN19 + crdN2*crdN42 - crdN2*crdN51 - crdN28 - crdN29*crdN5 - crdN3*crdN31*crdN5 + 2376.0*crdN33 + 1656.0*crdN34 - crdN36 - crdN38 + 516.0*crdN40 - 414.0*crdN43 + crdN44*crdN5 - 114.0*crdN46 + crdN47 + crdN48 - crdN5*crdN52 - 84.0*crdN53 + crdN57 + crdN78);
    rdN[5]=std::pow(J2, 3)*crdN61*(crdN2*crdN62 - crdN26 - crdN56*crdN81 + crdN63 + crdN65*crdN79 - crdN66*crdN79 - crdN67*crdN80 + crdN68*crdN79 - crdN70*crdN81 - crdN71 + crdN72*crdN80 + crdN73 + crdN77);

    rdN[7]=crdN87*(-crdN11*crdN24 + crdN11*xi + crdN15*crdN24 - crdN15*xi + crdN2*crdN73 - crdN2*crdN74 - crdN2*crdN75 + crdN2*crdN76 - crdN23 - crdN27 + 1008.0*crdN33 + 1008.0*crdN34 - 132.0*crdN43 - 132.0*crdN46 + crdN5*crdN73 + crdN5*crdN74 + crdN5*crdN75 - crdN5*crdN76 + crdN54*xi - crdN59*xi + crdN78);
    rdN[8]=std::pow(J3, 3)*crdN87*(10.0*crdN20*crdN9 + crdN22 - 2160.0*crdN3*crdN49 - crdN32*crdN88 - 12.0*crdN49*crdN9 + crdN58 - crdN69*crdN89 + crdN69*crdN90 + crdN69*crdN91 + crdN80*crdN88 - crdN81*crdN89 - crdN81*crdN90 + crdN81*crdN91 + crdN84 + crdN85 + crdN86);

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

    const double k_s = GetBendingShearStiffnessRatio();
    rd2N.clear();

	const double crd2N0 = std::pow(J1, 2);
	const double crd2N1 = std::pow(k_s, 2);
	const double crd2N2 = 1296.0*crd2N1;
	const double crd2N3 = std::pow(J2, 2);
	const double crd2N4 = std::pow(J3, 2);
	const double crd2N5 = crd2N0*crd2N3;
	const double crd2N6 = crd2N4*crd2N5;
	const double crd2N7 = 8.0*crd2N6;
	const double crd2N8 = crd2N4*k_s;
	const double crd2N9 = 156.0*crd2N8;
	const double crd2N10 = crd2N5*k_s;
	const double crd2N11 = 96.0*crd2N10;
	const double crd2N12 = std::pow(k_s, 3);
	const double crd2N13 = crd2N1*crd2N4;
	const double crd2N14 = -17280.0*crd2N12 + 2592.0*crd2N13;
	const double crd2N15 = 1.0/(-crd2N0*crd2N2 + crd2N0*crd2N9 - crd2N11 + crd2N14 - crd2N2*crd2N3 + crd2N3*crd2N9 + crd2N7);
	const double crd2N16 = crd2N1*crd2N3;
	const double crd2N17 = crd2N0*crd2N1;
	const double crd2N18 = std::pow(xi, 3);
	const double crd2N19 = 1440.0*crd2N18;
	const double crd2N20 = crd2N16*crd2N19;
	const double crd2N21 = 2880.0*crd2N13*crd2N18;
	const double crd2N22 = crd2N3*crd2N8;
	const double crd2N23 = std::pow(xi, 2);
	const double crd2N24 = 2160.0*crd2N23;
	const double crd2N25 = crd2N16*crd2N24;
	const double crd2N26 = crd2N0*crd2N8;
	const double crd2N27 = 6480.0*crd2N23;
	const double crd2N28 = 4320.0*crd2N18;
	const double crd2N29 = 8640.0*crd2N13*xi;
	const double crd2N30 = 120.0*crd2N18;
	const double crd2N31 = crd2N30*crd2N6;
	const double crd2N32 = 360.0*crd2N18;
	const double crd2N33 = 252.0*crd2N23;
	const double crd2N34 = 1320.0*crd2N18;
	const double crd2N35 = 576.0*xi;
	const double crd2N36 = 60.0*crd2N6*xi;
	const double crd2N37 = 240.0*crd2N18;
	const double crd2N38 = crd2N10*crd2N37;
	const double crd2N39 = 864.0*xi;
	const double crd2N40 = 1548.0*crd2N23;
	const double crd2N41 = crd2N23*crd2N6;
	const double crd2N42 = 576.0*crd2N10*crd2N23 - 192.0*crd2N10 + crd2N14 - 48.0*crd2N41 + 16.0*crd2N6;
	const double crd2N43 = 8640.0*crd2N12;
	const double crd2N44 = 1296.0*crd2N13;
	const double crd2N45 = 48.0*k_s;
	const double crd2N46 = crd2N0*crd2N45;
	const double crd2N47 = 1.0/(-648.0*crd2N16 - 648.0*crd2N17 + 78.0*crd2N22 + 78.0*crd2N26 - crd2N3*crd2N46 - crd2N43 + crd2N44 + 4.0*crd2N6);
	const double crd2N48 = 360.0*crd2N1;
	const double crd2N49 = 2*crd2N4;
	const double crd2N50 = 24.0*k_s;
	const double crd2N51 = crd2N3*crd2N50;
	const double crd2N52 = 54.0*crd2N8;
	const double crd2N53 = crd2N1*crd2N24;
	const double crd2N54 = crd2N3*crd2N4;
	const double crd2N55 = 20*crd2N18;
	const double crd2N56 = 6*xi;
	const double crd2N57 = 12*crd2N23;
	const double crd2N58 = crd2N3*k_s;
	const double crd2N59 = 324.0*crd2N23;
	const double crd2N60 = crd2N59*crd2N8;
	const double crd2N61 = 144.0*crd2N23;
	const double crd2N62 = crd2N1*crd2N19;
	const double crd2N63 = 72.0*xi;
	const double crd2N64 = -crd2N37*crd2N8 + crd2N62 + crd2N63*crd2N8;
	const double crd2N65 = crd2N17*crd2N24;
	const double crd2N66 = crd2N17*crd2N19;
	const double crd2N67 = crd2N0*crd2N4;
	const double crd2N68 = crd2N0*k_s;
	const double crd2N69 = 1.0/(-crd2N0*crd2N51 - 4320.0*crd2N12 + 648.0*crd2N13 - 324.0*crd2N16 - 324.0*crd2N17 + 39.0*crd2N22 + 39.0*crd2N26 + 2.0*crd2N6);
	const double crd2N70 = 288.0*crd2N23;
	const double crd2N71 = 324.0*xi;
	const double crd2N72 = 300.0*crd2N18;

	rd2N[1]=crd2N15*(-1656.0*crd2N16 + crd2N17*crd2N27 - crd2N17*crd2N28 - 2376.0*crd2N17 + crd2N20 + crd2N21 + crd2N22*crd2N32 + crd2N22*crd2N33 - crd2N22*crd2N35 + 114.0*crd2N22 + crd2N25 + crd2N26*crd2N34 - crd2N26*crd2N39 - crd2N26*crd2N40 + 414.0*crd2N26 - crd2N29 + crd2N31 - crd2N36 - crd2N38 + crd2N42);
	rd2N[2]=-std::pow(J1, 3)*crd2N47*(-crd2N3*crd2N49 + crd2N30*crd2N58 + crd2N48 + crd2N51 - crd2N52 - crd2N53 - crd2N54*crd2N55 + crd2N54*crd2N56 + crd2N54*crd2N57 - crd2N58*crd2N61 + crd2N60 + crd2N64);

	rd2N[4]=crd2N15*(crd2N16*crd2N27 + crd2N16*crd2N28 - 2376.0*crd2N16 - 1656.0*crd2N17 - crd2N21 - crd2N22*crd2N34 + crd2N22*crd2N39 - crd2N22*crd2N40 + 414.0*crd2N22 - crd2N26*crd2N32 + crd2N26*crd2N33 + crd2N26*crd2N35 + 114.0*crd2N26 + crd2N29 - crd2N31 + crd2N36 + crd2N38 + crd2N42 + crd2N65 - crd2N66);
	rd2N[5]=-std::pow(J2, 3)*crd2N47*(crd2N0*crd2N49 - crd2N0*crd2N50 + crd2N30*crd2N68 - crd2N48 + crd2N52 + crd2N53 - crd2N55*crd2N67 + crd2N56*crd2N67 - crd2N57*crd2N67 - crd2N60 + crd2N61*crd2N68 + crd2N64);

	rd2N[7]=crd2N69*(-crd2N10*crd2N70 + crd2N11 + 1008.0*crd2N16 + 1008.0*crd2N17 - crd2N20 + crd2N22*crd2N37 + crd2N22*crd2N59 - crd2N22*crd2N63 - 132.0*crd2N22 - crd2N25 - crd2N26*crd2N37 + crd2N26*crd2N59 + crd2N26*crd2N63 - 132.0*crd2N26 + 24.0*crd2N41 + crd2N43 - crd2N44 - crd2N65 + crd2N66 - crd2N7);
	rd2N[8]=std::pow(J3, 3)*crd2N69*(-4320.0*crd2N1*xi + 40.0*crd2N18*crd2N5 - crd2N3*crd2N45 + crd2N46 - 24.0*crd2N5*xi + crd2N58*crd2N70 - crd2N58*crd2N71 + crd2N58*crd2N72 + crd2N62 - crd2N68*crd2N70 - crd2N68*crd2N71 + crd2N68*crd2N72);

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

    const double k_s = GetBendingShearStiffnessRatio();

    rd3N.clear();

	const double crd3N0 = std::pow(k_s, 3);
	const double crd3N1 = std::pow(J1, 2);
	const double crd3N2 = std::pow(k_s, 2);
	const double crd3N3 = 1296.0*crd3N2;
	const double crd3N4 = std::pow(J2, 2);
	const double crd3N5 = std::pow(J3, 2);
	const double crd3N6 = crd3N2*crd3N5;
	const double crd3N7 = crd3N1*crd3N4;
	const double crd3N8 = crd3N7*k_s;
	const double crd3N9 = crd3N5*k_s;
	const double crd3N10 = 156.0*crd3N9;
	const double crd3N11 = crd3N5*crd3N7;
	const double crd3N12 = 1.0/(-17280.0*crd3N0 + crd3N1*crd3N10 - crd3N1*crd3N3 + crd3N10*crd3N4 + 8.0*crd3N11 - crd3N3*crd3N4 + 2592.0*crd3N6 - 96.0*crd3N8);
	const double crd3N13 = 8640.0*crd3N6;
	const double crd3N14 = 4320.0*crd3N2;
	const double crd3N15 = crd3N14*xi;
	const double crd3N16 = crd3N15*crd3N4;
	const double crd3N17 = std::pow(xi, 2);
	const double crd3N18 = crd3N14*crd3N17;
	const double crd3N19 = crd3N18*crd3N4;
	const double crd3N20 = crd3N13*crd3N17;
	const double crd3N21 = 12960.0*crd3N2;
	const double crd3N22 = crd3N1*xi;
	const double crd3N23 = crd3N4*crd3N9;
	const double crd3N24 = 60.0*crd3N11;
	const double crd3N25 = crd3N1*crd3N9;
	const double crd3N26 = crd3N17*crd3N21;
	const double crd3N27 = 360.0*crd3N17;
	const double crd3N28 = crd3N11*crd3N27;
	const double crd3N29 = crd3N23*xi;
	const double crd3N30 = 1080.0*crd3N17;
	const double crd3N31 = 3960.0*crd3N17;
	const double crd3N32 = 720.0*crd3N17;
	const double crd3N33 = crd3N32*crd3N8;
	const double crd3N34 = crd3N22*crd3N9;
	const double crd3N35 = crd3N8*xi;
	const double crd3N36 = crd3N11*xi;
	const double crd3N37 = 1152.0*crd3N35 - 96.0*crd3N36;
	const double crd3N38 = 648.0*crd3N2;
	const double crd3N39 = 1.0/(-8640.0*crd3N0 - crd3N1*crd3N38 + 4.0*crd3N11 + 78.0*crd3N23 + 78.0*crd3N25 - crd3N38*crd3N4 + 1296.0*crd3N6 - 48.0*crd3N8);
	const double crd3N40 = 6*crd3N5;
	const double crd3N41 = crd3N4*crd3N5;
	const double crd3N42 = 60*crd3N17;
	const double crd3N43 = crd3N4*k_s;
	const double crd3N44 = 648.0*crd3N9;
	const double crd3N45 = crd3N44*xi;
	const double crd3N46 = crd3N43*xi;
	const double crd3N47 = 72.0*crd3N9;
	const double crd3N48 = crd3N18 - crd3N32*crd3N9 + crd3N47;
	const double crd3N49 = crd3N1*crd3N15;
	const double crd3N50 = crd3N1*crd3N18;
	const double crd3N51 = crd3N22*k_s;
	const double crd3N52 = crd3N1*k_s;
	const double crd3N53 = 324.0*crd3N2;
	const double crd3N54 = 24.0*crd3N7;
	const double crd3N55 = 1.0/(-4320.0*crd3N0 - crd3N1*crd3N53 + 2.0*crd3N11 + 39.0*crd3N23 + 39.0*crd3N25 - crd3N4*crd3N53 - crd3N54*k_s + 648.0*crd3N6);
	const double crd3N56 = 324.0*k_s;
	const double crd3N57 = 900.0*crd3N17;

	rd3N[1]=crd3N12*(-crd3N1*crd3N26 - crd3N13 + crd3N16 + crd3N19 + crd3N20 + crd3N21*crd3N22 + crd3N23*crd3N30 - 576.0*crd3N23 - crd3N24 + crd3N25*crd3N31 - 864.0*crd3N25 + crd3N28 + 504.0*crd3N29 - crd3N33 - 3096.0*crd3N34 + crd3N37);
	rd3N[2]=-std::pow(J1, 3)*crd3N39*(-crd3N15 + crd3N27*crd3N43 + crd3N4*crd3N40 - crd3N41*crd3N42 + 24*crd3N41*xi + crd3N45 - 288.0*crd3N46 + crd3N48);

	rd3N[4]=crd3N12*(crd3N13 - crd3N20 + crd3N21*crd3N4*xi - crd3N23*crd3N31 + 864.0*crd3N23 + crd3N24 - crd3N25*crd3N30 + 576.0*crd3N25 + crd3N26*crd3N4 - crd3N28 - 3096.0*crd3N29 + crd3N33 + 504.0*crd3N34 + crd3N37 + crd3N49 - crd3N50);
	rd3N[5]=-std::pow(J2, 3)*crd3N39*(crd3N1*crd3N40 - crd3N1*crd3N42*crd3N5 + crd3N15 - 24*crd3N22*crd3N5 + crd3N27*crd3N52 - crd3N45 + crd3N48 + 288.0*crd3N51);

	rd3N[7]=crd3N55*(crd3N1*crd3N47 - crd3N16 - crd3N19 + crd3N22*crd3N44 + crd3N23*crd3N32 - crd3N25*crd3N32 + 648.0*crd3N29 - 576.0*crd3N35 + 48.0*crd3N36 - crd3N4*crd3N47 - crd3N49 + crd3N50);
	rd3N[8]=std::pow(J3, 3)*crd3N55*(-crd3N1*crd3N56 - crd3N14 + 120.0*crd3N17*crd3N7 + crd3N18 - crd3N4*crd3N56 + crd3N43*crd3N57 + 576.0*crd3N46 - 576.0*crd3N51 + crd3N52*crd3N57 - crd3N54);

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

    const double k_s = GetBendingShearStiffnessRatio();

    rd4N.clear();

	const double crd4N0 = std::pow(k_s, 3);
	const double crd4N1 = std::pow(J1, 2);
	const double crd4N2 = std::pow(k_s, 2);
	const double crd4N3 = 1296.0*crd4N2;
	const double crd4N4 = std::pow(J2, 2);
	const double crd4N5 = std::pow(J3, 2);
	const double crd4N6 = crd4N2*crd4N5;
	const double crd4N7 = crd4N1*crd4N4;
	const double crd4N8 = crd4N7*k_s;
	const double crd4N9 = crd4N5*k_s;
	const double crd4N10 = 156.0*crd4N9;
	const double crd4N11 = crd4N5*crd4N7;
	const double crd4N12 = 1.0/(-17280.0*crd4N0 + crd4N1*crd4N10 - crd4N1*crd4N3 + crd4N10*crd4N4 + 8.0*crd4N11 - crd4N3*crd4N4 + 2592.0*crd4N6 - 96.0*crd4N8);
	const double crd4N13 = 4320.0*crd4N2;
	const double crd4N14 = crd4N13*crd4N4;
	const double crd4N15 = crd4N1*crd4N2;
	const double crd4N16 = 1152.0*crd4N8;
	const double crd4N17 = crd4N4*crd4N9;
	const double crd4N18 = crd4N2*crd4N4;
	const double crd4N19 = 8640.0*xi;
	const double crd4N20 = crd4N18*crd4N19;
	const double crd4N21 = 96.0*crd4N11;
	const double crd4N22 = crd4N1*crd4N9;
	const double crd4N23 = 25920.0*xi;
	const double crd4N24 = 2160.0*xi;
	const double crd4N25 = 7920.0*xi;
	const double crd4N26 = 1440.0*xi;
	const double crd4N27 = 720.0*crd4N11*xi - crd4N26*crd4N8 + 17280.0*crd4N6*xi;
	const double crd4N28 = 24*crd4N5;
	const double crd4N29 = 288.0*k_s;
	const double crd4N30 = 648.0*crd4N9;
	const double crd4N31 = crd4N4*xi;
	const double crd4N32 = 120*crd4N5;
	const double crd4N33 = 720.0*k_s;
	const double crd4N34 = crd4N19*crd4N2;
	const double crd4N35 = crd4N26*crd4N9 - crd4N34;
	const double crd4N36 = 1.0/(-8640.0*crd4N0 + 4.0*crd4N11 - 648.0*crd4N15 + 78.0*crd4N17 - 648.0*crd4N18 + 78.0*crd4N22 + 1296.0*crd4N6 - 48.0*crd4N8);
	const double crd4N37 = -crd4N1*crd4N13 + crd4N15*crd4N19;
	const double crd4N38 = crd4N1*xi;
	const double crd4N39 = 1.0/(-4320.0*crd4N0 + 2.0*crd4N11 - 324.0*crd4N15 + 39.0*crd4N17 - 324.0*crd4N18 + 39.0*crd4N22 + 648.0*crd4N6 - 24.0*crd4N8);
	const double crd4N40 = 576.0*k_s;
	const double crd4N41 = crd4N4*crd4N40;
	const double crd4N42 = 1800.0*k_s;

	rd4N[1]=crd4N12*(crd4N14 - crd4N15*crd4N23 + 12960.0*crd4N15 + crd4N16 + crd4N17*crd4N24 + 504.0*crd4N17 + crd4N20 - crd4N21 + crd4N22*crd4N25 - 3096.0*crd4N22 + crd4N27);
	rd4N[2]=std::pow(J1, 3)*crd4N36*(crd4N13 - crd4N28*crd4N4 + crd4N29*crd4N4 - crd4N30 + crd4N31*crd4N32 - crd4N31*crd4N33 + crd4N35);

	rd4N[4]=-crd4N12*(-crd4N16 + crd4N17*crd4N25 + 3096.0*crd4N17 - crd4N18*crd4N23 - 12960.0*crd4N18 + crd4N21 + crd4N22*crd4N24 - 504.0*crd4N22 + crd4N27 + crd4N37);
	rd4N[5]=std::pow(J2, 3)*crd4N36*(crd4N1*crd4N28 - crd4N1*crd4N29 - crd4N13 + crd4N30 + crd4N32*crd4N38 - crd4N33*crd4N38 + crd4N35);

	rd4N[7]=crd4N39*(crd4N1*crd4N30 - crd4N1*crd4N41 + 48.0*crd4N11 - crd4N14 + crd4N17*crd4N26 - crd4N20 - crd4N22*crd4N26 + crd4N30*crd4N4 + crd4N37);
	rd4N[8]=std::pow(J3, 3)*crd4N39*(-crd4N1*crd4N40 + crd4N31*crd4N42 + crd4N34 + crd4N38*crd4N42 + crd4N41 + 240.0*crd4N7*xi);

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
    GlobalSizeVector dN, d3N;
    GetFirstDerivativesShapeFunctionsValues(dN,  J, xi);
    GetThirdDerivativesShapeFunctionsValues(d3N, J, xi);

    const double k_s = GetBendingShearStiffnessRatio();
    // v' + ks * v'''
    noalias(rN) = dN + k_s * d3N;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetFirstDerivativesNThetaShapeFunctionsValues(
    GlobalSizeVector& rN,
    const double J,
    const double xi
    )
{

    GlobalSizeVector d2N, d4N;
    GetSecondDerivativesShapeFunctionsValues (d2N,  J, xi);
    GetFourthDerivativesShapeFunctionsValues (d4N, J, xi);

    const double k_s = GetBendingShearStiffnessRatio();

    // v'' + ks * v''''
    noalias(rN) = d2N + k_s * d4N;
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
    GlobalSizeVector& rNodalValues
    )
{
    const auto &r_geom = GetGeometry();

    const auto& r_displ_0 = r_geom[0].FastGetSolutionStepValue(DISPLACEMENT);
    const auto& r_displ_1 = r_geom[1].FastGetSolutionStepValue(DISPLACEMENT);
    const auto& r_displ_2 = r_geom[2].FastGetSolutionStepValue(DISPLACEMENT);

    rNodalValues[0] = r_displ_0[0];
    rNodalValues[1] = r_displ_0[1];
    rNodalValues[2] = r_geom[0].FastGetSolutionStepValue(ROTATION_Z);

    rNodalValues[3] = r_displ_1[0];
    rNodalValues[4] = r_displ_1[1];
    rNodalValues[5] = r_geom[1].FastGetSolutionStepValue(ROTATION_Z);

    rNodalValues[6] = r_displ_2[0];
    rNodalValues[7] = r_displ_2[1];
    rNodalValues[8] = r_geom[2].FastGetSolutionStepValue(ROTATION_Z);
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
    GlobalSizeVector dNu, dN;
    GetFirstDerivativesNu0ShapeFunctionsValues(dNu, J, xi);
    GetFirstDerivativesShapeFunctionsValues(dN, J, xi);
    return inner_prod(dNu, rNodalValues) + 0.5 * std::pow(inner_prod(dNu, rNodalValues), 2) +
        0.5 * std::pow(inner_prod(dN, rNodalValues), 2);
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoCurvedBeamElement2D3N::CalculateShearStrain(
    const double J,
    const double xi,
    const GlobalSizeVector& rNodalValues
    )
{
    GlobalSizeVector dN, Ntheta;
    GetFirstDerivativesShapeFunctionsValues(dN, J, xi);
    GetNThetaShapeFunctionsValues(Ntheta, J, xi);

    return inner_prod(dN - Ntheta, rNodalValues);
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

    GlobalSizeVector dNu, d2Nu, dN_theta, N_shape, Nu, N_s, N_theta, dN_shape;

    // Loop over the integration points
    for (SizeType IP = 0; IP < integration_points.size(); ++IP) {

        const double xi     = integration_points[IP].X();
        const double weight = integration_points[IP].Weight();
        const double J  = GetJacobian(xi);
        const double jacobian_weight = weight * J;

        GetNodalValuesVector(nodal_values);

        const double E    = r_props[YOUNG_MODULUS];
        const double A    = r_props[CROSS_AREA];
        const double I    = r_props[I33];
        const double G    = ConstitutiveLawUtilities<3>::CalculateShearModulus(r_props);
        const double A_s  = r_props[AREA_EFFECTIVE_Y];

        const double N1 = 0.5 * xi * (xi - 1.0);
        const double N2 = 0.5 * xi * (xi + 1.0);
        const double N3 = 1.0 - std::pow(xi, 2);

        const double dN1 = (xi - 0.5) / J;
        const double dN2 = (xi + 0.5) / J;
        const double dN3 = (-2.0 * xi) / J;

        // deflection v
        N_shape.clear();
        dN_shape.clear();
        N_shape[1] = N1;
        N_shape[4] = N2;
        N_shape[7] = N3;
        dN_shape[1] = dN1;
        dN_shape[4] = dN2;
        dN_shape[7] = dN3;

        // axial u
        dNu.clear();
        Nu.clear();
        Nu[0] = N1;
        Nu[3] = N2;
        Nu[6] = N3;
        dNu[0] = dN1;
        dNu[3] = dN2;
        dNu[6] = dN3;

        // rotation
        dN_theta.clear();
        N_theta.clear();
        N_theta[2] = N1;
        N_theta[5] = N2;
        N_theta[8] = N3;
        dN_theta[2] = dN1;
        dN_theta[5] = dN2;
        dN_theta[8] = dN3;

        // Initialize matrices and vectors...
        BoundedMatrix<double, 2, 2> C_gamma, frenet_serret;
        C_gamma.clear();
        C_gamma(0, 0) = E * A;
        C_gamma(1, 1) = G * A_s;
        frenet_serret.clear();
        BoundedVector<double, 2> N, Gamma;
        N.clear();
        Gamma.clear();
        BoundedMatrix<double, 2, 9> B_s, aux_B_s;
        B_s.clear();
        aux_B_s.clear();
        GlobalSizeVector B_b;

        VectorType t, n;
        GetTangentandTransverseUnitVectors(xi, t, n);
        noalias(frenet_serret) = GetFrenetSerretMatrix(xi, t, n);
        noalias(B_b) =  dN_theta;

        // we fill aux_B_s
        for (IndexType i = 0; i < SystemSize; ++i) {
            aux_B_s(0, i) = dNu[i] + t[1] * N_theta[i];
            aux_B_s(1, i) = dN_shape[i] - t[0] * N_theta[i];
        }
        noalias(B_s) = prod(frenet_serret, aux_B_s);

        noalias(Gamma) = prod(B_s, nodal_values);
        noalias(N) = prod(C_gamma, Gamma);

        const double curvature = inner_prod(B_b, nodal_values);
        const double bending_moment = E * I * curvature;

        noalias(rRHS) -= jacobian_weight * (prod(trans(B_s), N) + bending_moment * B_b);
        noalias(rLHS) += jacobian_weight * (prod(trans(B_s), Matrix(prod(C_gamma, B_s))) + E * I * outer_prod(B_b, B_b));

    } // IP loop
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::CalculateLeftHandSide(
    MatrixType& rLHS,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_ERROR << "do not enter here" << std::endl;
    // KRATOS_TRY;
    // KRATOS_CATCH("");
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

    VectorType local_rhs = rRHS;

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

    GlobalSizeVector dNu, d2Nu, dN_theta, N_shape, Nu, N_s, N_theta, dN_shape;

    // Loop over the integration points
    for (SizeType IP = 0; IP < integration_points.size(); ++IP) {
        const double xi     = integration_points[IP].X();
        const double weight = integration_points[IP].Weight();
        const double J  = GetJacobian(xi);
        const double jacobian_weight = weight * J;

        GetNodalValuesVector(nodal_values);

        const double E    = r_props[YOUNG_MODULUS];
        const double A    = r_props[CROSS_AREA];
        const double I    = r_props[I33];
        const double G    = ConstitutiveLawUtilities<3>::CalculateShearModulus(r_props);
        const double A_s  = r_props[AREA_EFFECTIVE_Y];

        const double N1 = 0.5 * xi * (xi - 1.0);
        const double N2 = 0.5 * xi * (xi + 1.0);
        const double N3 = 1.0 - std::pow(xi, 2);

        const double dN1 = (xi - 0.5) / J;
        const double dN2 = (xi + 0.5) / J;
        const double dN3 = (-2.0 * xi) / J;

        // deflection v
        N_shape.clear();
        dN_shape.clear();
        N_shape[1] = N1;
        N_shape[4] = N2;
        N_shape[7] = N3;
        dN_shape[1] = dN1;
        dN_shape[4] = dN2;
        dN_shape[7] = dN3;

        // axial u
        dNu.clear();
        Nu.clear();
        Nu[0] = N1;
        Nu[3] = N2;
        Nu[6] = N3;
        dNu[0] = dN1;
        dNu[3] = dN2;
        dNu[6] = dN3;

        // rotation
        dN_theta.clear();
        N_theta.clear();
        N_theta[2] = N1;
        N_theta[5] = N2;
        N_theta[8] = N3;
        dN_theta[2] = dN1;
        dN_theta[5] = dN2;
        dN_theta[8] = dN3;

        // Initialize matrices and vectors...
        BoundedMatrix<double, 2, 2> C_gamma, frenet_serret;
        C_gamma.clear();
        C_gamma(0, 0) = E * A;
        C_gamma(1, 1) = G * A_s;
        frenet_serret.clear();
        BoundedVector<double, 2> N, Gamma;
        N.clear();
        Gamma.clear();
        BoundedMatrix<double, 2, 9> B_s, aux_B_s;
        B_s.clear();
        aux_B_s.clear();
        GlobalSizeVector B_b;

        VectorType t, n;
        GetTangentandTransverseUnitVectors(xi, t, n);
        noalias(frenet_serret) = GetFrenetSerretMatrix(xi, t, n);

        noalias(B_b) = dN_theta;

        // we fill aux_B_s
        for (IndexType i = 0; i < SystemSize; ++i) {
            aux_B_s(0, i) = dNu[i] + t[1] * N_theta[i];
            aux_B_s(1, i) = dN_shape[i] - t[0] * N_theta[i];
        }
        noalias(B_s) = prod(frenet_serret, aux_B_s);

        noalias(Gamma) = prod(B_s, nodal_values);
        noalias(N) = prod(C_gamma, Gamma);

        const double curvature = inner_prod(B_b, nodal_values);
        const double bending_moment = E * I * curvature;

        noalias(rRHS) -= jacobian_weight * (prod(trans(B_s), N) + bending_moment * B_b);

    } // IP loop
    KRATOS_CATCH("");
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
        rVariable == INITIAL_GEOMETRIC_CURVATURE ||
        rVariable == SHEAR_ANGLE)
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

        // Loop over the integration points
        for (SizeType IP = 0; IP < integration_points.size(); ++IP) {
            const double xi = integration_points[IP].X();
            GetNodalValuesVector(nodal_values);
            const double J = GetJacobian(xi);

            GlobalSizeVector dNu, dN_theta, N_shape, Nu, N_s, N_theta, dN_shape;

            // const double k0 = GetGeometryCurvature(J, xi);

            // BoundedMatrix<double, 9, 9> global_T;
            // noalias(global_T) = GetGlobalSizeRotationMatrixGlobalToLocalAxes(xi);

            // we rotate to local axes
            // nodal_values = prod(global_T, nodal_values);

            const double E    = r_props[YOUNG_MODULUS];
            const double A    = r_props[CROSS_AREA];
            const double I    = r_props[I33];
            const double G    = ConstitutiveLawUtilities<3>::CalculateShearModulus(r_props);
            const double A_s  = r_props[AREA_EFFECTIVE_Y];

            const double N1 = 0.5 * xi * (xi - 1.0);
            const double N2 = 0.5 * xi * (xi + 1.0);
            const double N3 = 1.0 - std::pow(xi, 2);

            const double dN1 = (xi - 0.5) / J;
            const double dN2 = (xi + 0.5) / J;
            const double dN3 = (-2.0 * xi) / J;

            // deflection v
            N_shape.clear();
            dN_shape.clear();
            N_shape[1] = N1;
            N_shape[4] = N2;
            N_shape[7] = N3;
            dN_shape[1] = dN1;
            dN_shape[4] = dN2;
            dN_shape[7] = dN3;

            // axial u
            dNu.clear();
            Nu.clear();
            Nu[0] = N1;
            Nu[3] = N2;
            Nu[6] = N3;
            dNu[0] = dN1;
            dNu[3] = dN2;
            dNu[6] = dN3;

            // rotation
            dN_theta.clear();
            N_theta.clear();
            N_theta[2] = N1;
            N_theta[5] = N2;
            N_theta[8] = N3;
            dN_theta[2] = dN1;
            dN_theta[5] = dN2;
            dN_theta[8] = dN3;


            // strain_vector[0] =  (inner_prod(dNu, nodal_values) - k0 * inner_prod(N_shape, nodal_values) + 0.5 * (std::pow(inner_prod(dNu, nodal_values), 2)) + 0.5 * (std::pow(inner_prod(dN_shape, nodal_values), 2)));
            strain_vector[0] =  (inner_prod(dNu, nodal_values) + 0.5 * (std::pow(inner_prod(dNu, nodal_values), 2)) + 0.5 * (std::pow(inner_prod(dN_shape, nodal_values), 2)));
            strain_vector[1] =  inner_prod(dN_theta, nodal_values);
            // strain_vector[2] =  (inner_prod(dN_shape, nodal_values) + k0 * inner_prod(Nu, nodal_values) - inner_prod(N_theta, nodal_values));
            strain_vector[2] =  (inner_prod(dN_shape, nodal_values)  - inner_prod(N_theta, nodal_values));

            Vector &r_generalized_stresses = cl_values.GetStressVector();
            r_generalized_stresses[0] = E * A * strain_vector[0];
            r_generalized_stresses[1] = E * I * strain_vector[1];
            r_generalized_stresses[2] = G * A_s * strain_vector[2];

            if (rVariable == AXIAL_FORCE)
                rOutput[IP] = r_generalized_stresses[0];
            else if (rVariable == BENDING_MOMENT)
                rOutput[IP] = r_generalized_stresses[1];
            else if (rVariable == SHEAR_FORCE)
                rOutput[IP] = r_generalized_stresses[2];
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

void LinearTimoshenkoCurvedBeamElement2D3N::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto &integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    const SizeType number_of_integration_points = integration_points.size();
    if (rOutput.size() != number_of_integration_points)
        rOutput.resize( number_of_integration_points );

    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

    // if (rVariable == LOCAL_AXIS_1 || rVariable == LOCAL_AXIS_2) {
    //     for (IndexType IP = 0; IP < number_of_integration_points; ++IP) {
    //         rOutput[IP].clear();
    //         const double xi = integration_points[IP].X();
    //         const auto T = GetFrenetSerretMatrix(xi);
    //         if (rVariable == LOCAL_AXIS_1) {
    //             rOutput[IP][0] = T(0, 0);
    //             rOutput[IP][1] = T(0, 1);
    //         } else {
    //             rOutput[IP][0] = T(1, 0);
    //             rOutput[IP][1] = T(1, 1);
    //         }
    //     }
    // }

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
