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

    // Here we compute the elemental pressure
    Vector displ, displ_old;
    GetValuesVector(displ, 0);
    GetValuesVector(displ_old, 1);
    if (mKpp != 0.0)
        mPressure -= (mFp + inner_prod(mKup, displ - displ_old)) / mKpp;
    mKpp = 0.0;
    mFp = 0.0;
    noalias(mKup) = ZeroVector(mat_size);

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
        noalias(mKup) -= int_to_reference_weight * this_kinematic_variables.detF * prod(trans(this_kinematic_variables.B), inv_c_voigt);
        mKpp          -= int_to_reference_weight / bulk_modulus;
        mFp           -= int_to_reference_weight * ((this_kinematic_variables.detF - 1.0) + (mPressure / bulk_modulus));

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
        noalias(rLeftHandSideMatrix) -= outer_prod(mKup, mKup) / mKpp;
    if (CalculateResidualVectorFlag)
        noalias(rRightHandSideVector) += mKup * mFp / mKpp;

    // mPressure -= (mFp + inner_prod(mKup, displ - displ_old)) / mKpp;
    // if (Id() == 15) {
    //     KRATOS_WATCH(mPressure)
    // }
    // KRATOS_WATCH(mPressure)

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangianQ1P0MixedElement::InitializeMaterial()
{
    KRATOS_TRY

    BaseType::InitializeMaterial();

    // now we initialize the Kup vector
    const auto &r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;
    mKup.resize(mat_size, false);
    noalias(mKup) = ZeroVector(mat_size);

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
    Matrix C_modif = rC;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3;j++)
            if (i != j)
                C_modif(i, j) *= 2.0;
    double v[308];
    v[1]=C_modif(0,0);
    v[2]=C_modif(0,1);
    v[3]=C_modif(0,2);
    v[4]=C_modif(1,0);
    v[5]=C_modif(1,1);
    v[40]=-(v[2]*v[4])+v[1]*v[5];
    v[6] = C_modif(1, 2);
    v[104]=v[3]*v[4]-v[1]*v[6];
    v[91]=-(v[3]*v[5])+v[2]*v[6];
    v[7] = C_modif(2, 0);
    v[8] = C_modif(2, 1);
    v[43]=-(v[5]*v[7])+v[4]*v[8];
    v[42]=v[2]*v[7]-v[1]*v[8];
    v[9]=C_modif(2,2);
    v[84]=v[3]*v[8]-v[2]*v[9];
    v[41]=v[6]*v[7]-v[4]*v[9];
    v[39]=-(v[3]*v[7])+v[1]*v[9];
    v[179]=1e0*v[39];
    v[38]=-(v[6]*v[8])+v[5]*v[9];
    v[14]=v[1]*v[38]+v[2]*v[41]+v[3]*v[43];
    v[168]=0.5e0/std::pow(v[14],0.5e0);
    v[169]=1e0*v[168];
    v[57]=v[168]*v[43];
    v[56]=v[168]*v[42];
    v[55]=v[169]*v[41];
    v[54]=v[169]*v[40];
    v[53]=v[169]*v[179];
    v[51]=1e0*v[169]*v[38];
    v[45]=1e0/(v[14]*v[14]);
    v[50]=-(v[43]*v[45]);
    v[111]=v[104]*v[50];
    v[49]=-(v[42]*v[45]);
    v[90]=v[38]*v[49];
    v[48]=-(v[41]*v[45]);
    v[103]=v[48]*v[91];
    v[47]=-(v[40]*v[45]);
    v[108]=v[39]*v[47];
    v[97]=v[38]*v[47];
    v[96]=v[47]*v[84];
    v[46]=-(v[39]*v[45]);
    v[88]=v[38]*v[46];
    v[44]=-(v[38]*v[45]);
    v[99]=v[108]+v[1]/v[14];
    v[80]=v[5]/v[14]+v[97];
    v[79]=v[88]+v[9]/v[14];
    v[10]=C1;
    v[11]=Pressure;
    v[195]=v[11]*v[51];
    v[187]=v[11]*v[54];
    v[186]=-(v[11]*v[57]);
    v[185]=v[11]*v[56];
    v[184]=v[11]*v[55];
    v[182]=v[11]*v[53];
    v[12]=std::pow(v[14],0.5e0);
    v[170]=-(v[11]*v[12]);
    v[171]=((-4e0/3e0)*v[10])/std::pow(v[12],0.16666666666666669e1);
    v[70]=v[171]*v[57];
    v[69]=v[171]*v[56];
    v[68]=v[171]*v[55];
    v[67]=v[171]*v[54];
    v[66]=v[171]*v[53];
    v[65]=v[171]*v[51];
    v[31]=(2e0*v[10])/std::pow(v[12],0.6666666666666666e0);
    v[72]=-0.3333333333333333e0*v[31];
    v[13]=v[1]+v[5]+v[9];
    v[175]=1e0*v[13];
    v[172]=-0.3333333333333333e0*v[13];
    v[173]=1e0*v[172];
    v[160]=v[172]*v[99];
    v[155]=v[172]*v[80];
    v[138]=v[173]*v[79];
    v[120]=v[175]*v[72];
    v[188]=v[120]*v[49];
    v[181]=1e0*v[120]*v[50];
    v[180]=v[120]*v[48];
    v[194]=v[186]+v[173]*v[70];
    v[76]=1e0*v[173]*v[69];
    v[198]=-v[185]+v[76];
    v[75]=1e0*v[173]*v[68];
    v[193]=-v[184]+v[75];
    v[74]=1e0*v[173]*v[67]+v[72];
    v[196]=v[187]+v[74];
    v[73]=1e0*v[173]*v[66]+v[72];
    v[191]=v[182]+v[73];
    v[71]=1e0*v[173]*v[65]+v[72];
    v[189]=v[195]+v[71];
    v[192]=v[120]*v[46];
    v[197]=v[120]*v[47];
    v[190]=v[120]*v[44];
    v[176]=v[120]+v[170];
    v[15]=v[38]/v[14];
    v[174]=-(v[11]*v[15]);
    v[178]=v[174]*v[54]+v[170]*v[80];
    v[177]=v[174]*v[53]+v[170]*v[79];
    v[116]=-0.3333333333333333e0*v[15];
    v[114]=1e0+v[116]*v[175];
    v[17]=v[84]/v[14];
    v[18]=v[91]/v[14];
    v[20]=v[39]/v[14];
    v[183]=-(v[187]*v[20])+v[170]*v[99];
    v[140]=-0.3333333333333333e0*v[20];
    v[136]=1e0+1e0*v[140]*v[175];
    v[21]=v[104]/v[14];
    v[24]=v[40]/v[14];
    v[159]=-0.3333333333333333e0*v[24];
    v[153]=1e0+1e0*v[159]*v[175];
    rStress[0]=v[12]*v[174]+v[114]*v[31];
    rStress[1]=v[170]*v[20]+v[136]*v[31];
    rStress[2]=v[170]*v[24]+v[153]*v[31];
    rStress[3]=v[17]*v[176];
    rStress[4]=v[176]*v[21];
    rStress[5]=v[176]*v[18];

    rTangentTensor(0,0)=v[31]*(v[116]+1e0*v[173]*v[38]*v[44])-v[174]*v[51]+v[114]*v[65];
    rTangentTensor(0,1)=v[177]+(v[116]+v[138])*v[31]+v[114]*v[66];
    rTangentTensor(0,2)=v[178]+(v[116]+v[155])*v[31]+v[114]*v[67];
    rTangentTensor(0,3)=v[180]*v[38]-v[174]*v[55]+v[114]*v[68];
    rTangentTensor(0,4)=v[174]*v[56]+v[114]*v[69]+v[176]*(-(v[8]/v[14])+v[90]);
    rTangentTensor(0,5)=v[181]*v[38]-v[174]*v[57]+v[114]*v[70];
    rTangentTensor(1,0)=v[177]+(v[138]+v[140])*v[31]+v[136]*v[65];
    rTangentTensor(1,1)=v[182]*v[20]+v[31]*(v[140]+1e0*v[172]*v[179]*v[46])+v[136]*v[66];
    rTangentTensor(1,2)=v[183]+(v[140]+v[160])*v[31]+v[136]*v[67];
    rTangentTensor(1,3)=v[179]*v[180]+v[184]*v[20]+v[136]*v[68];
    rTangentTensor(1,4)=1e0*v[179]*v[188]+v[185]*v[20]+v[136]*v[69];
    rTangentTensor(1,5)=v[186]*v[20]+v[176]*(v[39]*v[50]-v[7]/v[14])+v[136]*v[70];
    rTangentTensor(2,0)=v[178]+(v[155]+v[159])*v[31]+v[153]*v[65];
    rTangentTensor(2,1)=v[183]+(v[159]+v[160])*v[31]+v[153]*v[66];
    rTangentTensor(2,2)=v[187]*v[24]+v[31]*(v[159]+v[173]*v[40]*v[47])+v[153]*v[67];
    rTangentTensor(2,3)=v[111]*v[176]-v[184]*v[24]+v[153]*v[68];
    rTangentTensor(2,4)=v[185]*v[24]+v[188]*v[40]+v[153]*v[69];
    rTangentTensor(2,5)=-(v[186]*v[24])+v[181]*v[40]+v[153]*v[70];
    rTangentTensor(3,0)=v[17]*v[189]+v[190]*v[84];
    rTangentTensor(3,1)=v[17]*v[191]+v[192]*v[84];
    rTangentTensor(3,2)=v[17]*(-v[187]+v[74])+v[176]*(-(v[2]/v[14])+v[96]);
    rTangentTensor(3,3)=v[17]*v[193]+v[176]*v[88];
    rTangentTensor(3,4)=v[17]*(v[185]+4e0*v[76]);
    rTangentTensor(3,5)=v[17]*v[194]+v[176]*v[90];
    rTangentTensor(4,0)=v[103]*v[176]+v[21]*(-v[195]+v[71]);
    rTangentTensor(4,1)=v[104]*v[192]+v[191]*v[21];
    rTangentTensor(4,2)=v[104]*v[197]+v[196]*v[21];
    rTangentTensor(4,3)=v[21]*(v[184]+4e0*v[75]);
    rTangentTensor(4,4)=v[108]*v[176]+v[198]*v[21];
    rTangentTensor(4,5)=v[194]*v[21]+v[176]*(v[111]+v[4]/v[14]);
    rTangentTensor(5,0)=v[18]*v[189]+v[190]*v[91];
    rTangentTensor(5,1)=v[18]*(-v[182]+v[73])+v[176]*(-(v[3]/v[14])+v[46]*v[91]);
    rTangentTensor(5,2)=v[18]*v[196]+v[197]*v[91];
    rTangentTensor(5,3)=v[18]*v[193]+v[176]*(v[103]+v[6]/v[14]);
    rTangentTensor(5,4)=v[18]*v[198]+v[176]*v[96];
    rTangentTensor(5,5)=v[18]*v[194]+v[176]*v[97];
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


