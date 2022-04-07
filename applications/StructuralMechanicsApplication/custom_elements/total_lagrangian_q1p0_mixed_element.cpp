// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

// System includes

// External includes


// Project includes
#include "custom_elements/total_lagrangian_q1p0_mixed_element.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

TotalLagrangianQ1P0MixedElement::TotalLagrangianQ1P0MixedElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : TotalLagrangian(NewId, pGeometry)
{
    // DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

TotalLagrangianQ1P0MixedElement::TotalLagrangianQ1P0MixedElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : TotalLagrangian( NewId, pGeometry, pProperties )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer TotalLagrangianQ1P0MixedElement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<TotalLagrangianQ1P0MixedElement>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

//************************************************************************************
//************************************************************************************

Element::Pointer TotalLagrangianQ1P0MixedElement::Create( IndexType NewId,  GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<TotalLagrangianQ1P0MixedElement>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

TotalLagrangianQ1P0MixedElement::~TotalLagrangianQ1P0MixedElement()
{
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer TotalLagrangianQ1P0MixedElement::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    TotalLagrangianQ1P0MixedElement::Pointer p_new_elem = Kratos::make_intrusive<TotalLagrangianQ1P0MixedElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(BaseType::mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(BaseType::mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangianQ1P0MixedElement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY;

    const auto &r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const auto strain_size = GetStrainSize();
    const auto &r_props = GetProperties();

    KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
    ConstitutiveVariables this_constitutive_variables(strain_size);

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * dimension;

    if ( CalculateStiffnessMatrixFlag == true ) { // Calculation of the matrix is required
        if ( rLeftHandSideMatrix.size1() != mat_size )
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
    }

    // Resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) { // Calculation of the matrix is required
        if ( rRightHandSideVector.size() != mat_size )
            rRightHandSideVector.resize( mat_size, false );

        noalias(rRightHandSideVector) = ZeroVector( mat_size ); //resetting RHS
    }

    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    ConstitutiveLaw::Parameters Values(r_geometry,r_props,rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    if ( CalculateStiffnessMatrixFlag ) {
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    } else {
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    }

    // If strain has to be computed inside of the constitutive law with PK2
    Values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter

    // Some declarations
    array_1d<double, 3> body_force;
    double int_to_reference_weight;

    const double E = r_props[YOUNG_MODULUS];
    const double nu = r_props[POISSON_RATIO];
    const double mu = E / (2.0 * (1.0 + nu));
    const double C1 = 0.5 * mu;
    const double bulk_modulus = E / (3.0 * (1.0 - 2.0 * nu));

    double Kpp = 0.0;
    double Fp = 0.0;
    Vector Kup(mat_size);
    noalias(Kup) = ZeroVector(mat_size);

    // Computing in all integrations points
    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number) {

        // Contribution to external forces
        noalias(body_force) = this->GetBodyForce(integration_points, point_number);

        // Compute element kinematics B, F, DN_DX ...
        this->CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

        // Compute material reponse
        // this->CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, this->GetStressMeasure());

        // let's try to calculate split NeoHookean
        const Matrix C = prod(trans(this_kinematic_variables.F), this_kinematic_variables.F);
        CalculateNeoHookeanStressAndTangent(C, mPressure, C1, this_constitutive_variables.StressVector, this_constitutive_variables.D);

        Matrix inv_C;
        double det;
        MathUtils<double>::InvertMatrix3(C, inv_C, det);
        Vector inv_c_voigt = MathUtils<double>::StrainTensorToVector(inv_C, GetStrainSize());

        // Calculating weights for integration on the reference configuration
        int_to_reference_weight = GetIntegrationWeight(integration_points, point_number, this_kinematic_variables.detJ0);

        if (dimension == 2 && r_props.Has(THICKNESS))
            int_to_reference_weight *= r_props[THICKNESS];

        // we compute u-p entities
        inv_c_voigt[3] /= 2.0;
        inv_c_voigt[4] /= 2.0;
        inv_c_voigt[5] /= 2.0;
        noalias(Kup) -= int_to_reference_weight * this_kinematic_variables.detF * prod(trans(this_kinematic_variables.B), inv_c_voigt);
        Kpp          -= int_to_reference_weight / bulk_modulus;
        Fp           -= int_to_reference_weight * ((this_kinematic_variables.detF - 1.0) + (mPressure / bulk_modulus));

        if (CalculateStiffnessMatrixFlag) { // Calculation of the matrix is required
            // Contributions to stiffness matrix calculated on the reference config
            /* Material stiffness matrix */
            this->CalculateAndAddKm(rLeftHandSideMatrix, this_kinematic_variables.B, this_constitutive_variables.D, int_to_reference_weight);

            /* Geometric stiffness matrix */
            this->CalculateAndAddKg(rLeftHandSideMatrix, this_kinematic_variables.DN_DX, this_constitutive_variables.StressVector, int_to_reference_weight);
        }

        if ( CalculateResidualVectorFlag ) { // Calculation of the matrix is required
            this->CalculateAndAddResidualVector(rRightHandSideVector, this_kinematic_variables, rCurrentProcessInfo, body_force, this_constitutive_variables.StressVector, int_to_reference_weight);
        }
    } // IP loop
    if (CalculateStiffnessMatrixFlag)
        noalias(rLeftHandSideMatrix) -= outer_prod(Kup, Kup) / Kpp;
    if (CalculateResidualVectorFlag)
        noalias(rRightHandSideVector) += Kup * Fp / Kpp;

    Vector displ, displ_old;
    GetValuesVector(displ, 0);
    GetValuesVector(displ_old, 1);
    mPressure -= (Fp + inner_prod(Kup, displ - displ_old)) / Kpp;

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangianQ1P0MixedElement::InitializeMaterial()
{
    KRATOS_TRY

    BaseType::InitializeMaterial();

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangianQ1P0MixedElement::CalculateNeoHookeanStressAndTangent(
    const Matrix &rC,
    const double Pressure,
    const double C1,
    Vector &rStress,
    Matrix &rTangentTensor
    )
{
    double v[308];
    v[1]=rC(0,0);
    v[2]=rC(0,1);
    v[3]=rC(0,2);
    v[4]=rC(1,0);
    v[5]=rC(1,1);
    v[58]=-(v[2]*v[4])+v[1]*v[5];
    v[6]=rC(1,2);
    v[167]=-(v[3]*v[5])+v[2]*v[6];
    v[154]=v[3]*v[4]-v[1]*v[6];
    v[7]=rC(2,0);
    v[8]=rC(2,1);
    v[82]=-(v[5]*v[7])+v[4]*v[8];
    v[74]=v[2]*v[7]-v[1]*v[8];
    v[9]=rC(2,2);
    v[135]=v[3]*v[8]-v[2]*v[9];
    v[66]=v[6]*v[7]-v[4]*v[9];
    v[46]=-(v[3]*v[7])+v[1]*v[9];
    v[38]=-(v[6]*v[8])+v[5]*v[9];
    v[14]=v[1]*v[38]+v[2]*v[66]+v[3]*v[82];
    v[200]=0.5e0/std::pow(v[14],0.5e0);
    v[201]=1e0*v[200];
    v[84]=v[200]*v[82];
    v[76]=v[200]*v[74];
    v[68]=v[201]*v[66];
    v[60]=v[201]*v[58];
    v[50]=1e0*v[201]*v[46];
    v[47]=1e0/(v[14]*v[14]);
    v[83]=-(v[47]*v[82]);
    v[75]=-(v[47]*v[74]);
    v[152]=v[38]*v[75];
    v[67]=-(v[47]*v[66]);
    v[165]=v[58]*v[67];
    v[59]=-(v[47]*v[58]);
    v[178]=v[38]*v[59];
    v[176]=v[135]*v[59];
    v[163]=v[46]*v[59];
    v[48]=-(v[46]*v[47]);
    v[146]=v[38]*v[48];
    v[40]=1e0*v[201]*v[38];
    v[39]=-(v[38]*v[47]);
    v[174]=v[154]*v[39];
    v[99]=v[1]/v[14]+v[163];
    v[63]=v[178]+v[5]/v[14];
    v[54]=v[146]+v[9]/v[14];
    v[10]=C1;
    v[11]=Pressure;
    v[226]=v[11]*v[40];
    v[218]=v[11]*v[60];
    v[217]=-(v[11]*v[84]);
    v[216]=v[11]*v[76];
    v[214]=v[11]*v[68];
    v[212]=v[11]*v[50];
    v[12]=std::pow(v[14],0.5e0);
    v[202]=-(v[11]*v[12]);
    v[203]=((-4e0/3e0)*v[10])/std::pow(v[12],0.16666666666666669e1);
    v[86]=v[203]*v[84];
    v[78]=v[203]*v[76];
    v[70]=v[203]*v[68];
    v[62]=v[203]*v[60];
    v[53]=v[203]*v[50];
    v[42]=v[203]*v[40];
    v[31]=(2e0*v[10])/std::pow(v[12],0.6666666666666666e0);
    v[138]=-0.3333333333333333e0*v[31];
    v[13]=v[1]+v[5]+v[9];
    v[207]=1e0*v[13];
    v[204]=-0.3333333333333333e0*v[13];
    v[205]=1e0*v[204];
    v[225]=v[217]+v[204]*v[86];
    v[148]=v[204]*v[78];
    v[229]=v[148]-v[216];
    v[145]=v[205]*v[70];
    v[224]=v[145]-v[214];
    v[142]=v[138]+v[205]*v[62];
    v[227]=v[142]+v[218];
    v[139]=v[138]+1e0*v[205]*v[53];
    v[223]=v[139]+v[212];
    v[134]=v[138]+1e0*v[205]*v[42];
    v[220]=v[134]+v[226];
    v[117]=1e0*v[205]*v[99];
    v[114]=1e0*v[205]*v[63];
    v[93]=1e0*v[205]*v[54];
    v[34]=1e0*v[205]*v[31];
    v[222]=v[34]*v[48];
    v[219]=v[34]*v[58];
    v[215]=v[34]*v[46];
    v[211]=v[34]*v[38];
    v[228]=v[34]*v[59];
    v[221]=v[34]*v[39];
    v[208]=v[202]+v[34];
    v[15]=v[38]/v[14];
    v[206]=-(v[11]*v[15]);
    v[210]=v[206]*v[60]+v[202]*v[63];
    v[209]=v[206]*v[50]+v[202]*v[54];
    v[55]=-0.3333333333333333e0*v[15];
    v[44]=1e0+v[207]*v[55];
    v[17]=v[135]/v[14];
    v[18]=v[167]/v[14];
    v[20]=v[46]/v[14];
    v[213]=-(v[20]*v[218])+v[202]*v[99];
    v[96]=-0.3333333333333333e0*v[20];
    v[91]=1e0+v[207]*v[96];
    v[21]=v[154]/v[14];
    v[24]=v[58]/v[14];
    v[116]=-0.3333333333333333e0*v[24];
    v[112]=1e0+1e0*v[116]*v[207];
    rStress[0]=v[12]*v[206]+v[31]*v[44];
    rStress[1]=v[20]*v[202]+v[31]*v[91];
    rStress[2]=v[202]*v[24]+v[112]*v[31];
    rStress[3]=v[17]*v[208];
    rStress[4]=v[208]*v[21];
    rStress[5]=v[18]*v[208];
    rTangentTensor(0,0)=2e0*(-(v[206]*v[40])+v[42]*v[44]+v[31]*(v[205]*v[38]*v[39]+v[55]));
    rTangentTensor(0,1)=2e0*(v[209]+v[44]*v[53]+v[31]*(v[55]+v[93]));
    rTangentTensor(0,2)=2e0*(v[210]+v[31]*(v[114]+v[55])+v[44]*v[62]);
    rTangentTensor(0,3)=v[211]*v[67]-v[206]*v[68]+v[44]*v[70];
    rTangentTensor(0,4)=v[206]*v[76]+v[44]*v[78]+v[208]*(v[152]-v[8]/v[14]);
    rTangentTensor(0,5)=v[211]*v[83]-v[206]*v[84]+v[44]*v[86];
    rTangentTensor(1,0)=2e0*(v[209]+v[42]*v[91]+v[31]*(v[93]+v[96]));
    rTangentTensor(1,1)=2e0*(v[20]*v[212]+v[53]*v[91]+v[31]*(v[205]*v[46]*v[48]+v[96]));
    rTangentTensor(1,2)=2e0*(v[213]+v[62]*v[91]+v[31]*(v[117]+v[96]));
    rTangentTensor(1,3)=v[20]*v[214]+v[215]*v[67]+v[70]*v[91];
    rTangentTensor(1,4)=v[20]*v[216]+1e0*v[215]*v[75]+v[78]*v[91];
    rTangentTensor(1,5)=v[20]*v[217]+v[208]*(-(v[7]/v[14])+v[46]*v[83])+v[86]*v[91];
    rTangentTensor(2,0)=2e0*(v[210]+(v[114]+v[116])*v[31]+v[112]*v[42]);
    rTangentTensor(2,1)=2e0*(v[213]+(v[116]+v[117])*v[31]+v[112]*v[53]);
    rTangentTensor(2,2)=2e0*(v[218]*v[24]+v[31]*(v[116]+v[205]*v[58]*v[59])+v[112]*v[62]);
    rTangentTensor(2,3)=-(v[214]*v[24])+v[208]*(v[165]-v[4]/v[14])+v[112]*v[70];
    rTangentTensor(2,4)=v[216]*v[24]+v[219]*v[75]+v[112]*v[78];
    rTangentTensor(2,5)=-(v[217]*v[24])+v[219]*v[83]+v[112]*v[86];
    rTangentTensor(3,0)=2e0*(v[17]*v[220]+v[135]*v[221]);
    rTangentTensor(3,1)=2e0*(v[135]*v[222]+v[17]*v[223]);
    rTangentTensor(3,2)=2e0*((v[176]-v[2]/v[14])*v[208]+v[17]*(v[142]-v[218]));
    rTangentTensor(3,3)=v[146]*v[208]+v[17]*v[224];
    rTangentTensor(3,4)=v[17]*(4e0*v[148]+v[216]);
    rTangentTensor(3,5)=v[152]*v[208]+v[17]*v[225];
    rTangentTensor(4,0)=2e0*(v[21]*(v[134]-v[226])+v[208]*(v[174]-v[6]/v[14]));
    rTangentTensor(4,1)=2e0*(v[154]*v[222]+v[21]*v[223]);
    rTangentTensor(4,2)=2e0*(v[21]*v[227]+v[154]*v[228]);
    rTangentTensor(4,3)=v[21]*(4e0*v[145]+v[214]);
    rTangentTensor(4,4)=v[163]*v[208]+v[21]*v[229];
    rTangentTensor(4,5)=v[165]*v[208]+v[21]*v[225];
    rTangentTensor(5,0)=2e0*(v[18]*v[220]+v[167]*v[221]);
    rTangentTensor(5,1)=2e0*(v[18]*(v[139]-v[212])+v[208]*(-(v[3]/v[14])+v[167]*v[48]));
    rTangentTensor(5,2)=2e0*(v[18]*v[227]+v[167]*v[228]);
    rTangentTensor(5,3)=v[174]*v[208]+v[18]*v[224];
    rTangentTensor(5,4)=v[176]*v[208]+v[18]*v[229];
    rTangentTensor(5,5)=v[178]*v[208]+v[18]*v[225];

}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangianQ1P0MixedElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, TotalLagrangian );
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangianQ1P0MixedElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, TotalLagrangian );
}

} // Namespace Kratos


