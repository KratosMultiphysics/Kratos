//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//


// System includes
#include <omp.h>

// External includes

// Project includes
#include "custom_elements/mpm_updated_lagrangian.hpp"
#include "includes/define.h"
#include "custom_elements/mpm_updated_lagrangian_axisym.hpp"
#include "includes/constitutive_law.h"
#include "includes/checks.h"
#include "custom_utilities/mpm_math_utilities.h"

namespace Kratos
{


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

MPMUpdatedLagrangianAxisym::MPMUpdatedLagrangianAxisym( )
    : MPMUpdatedLagrangian()
{
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************
MPMUpdatedLagrangianAxisym::MPMUpdatedLagrangianAxisym( IndexType NewId, GeometryType::Pointer pGeometry )
    : MPMUpdatedLagrangian( NewId, pGeometry )
{
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************
MPMUpdatedLagrangianAxisym::MPMUpdatedLagrangianAxisym( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : MPMUpdatedLagrangian( NewId, pGeometry, pProperties )
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************
MPMUpdatedLagrangianAxisym::MPMUpdatedLagrangianAxisym( MPMUpdatedLagrangianAxisym const& rOther)
    :MPMUpdatedLagrangian(rOther)
{
}

//******************************ASSIGNMENT OPERATOR***********************************
//************************************************************************************
MPMUpdatedLagrangianAxisym&  MPMUpdatedLagrangianAxisym::operator=(MPMUpdatedLagrangianAxisym const& rOther)
{
    Element::operator=(rOther);

    mMP = rOther.mMP;

    mDeformationGradientF0.clear();
    mDeformationGradientF0 = rOther.mDeformationGradientF0;

    mDeterminantF0 = rOther.mDeterminantF0;
    mConstitutiveLawVector = rOther.mConstitutiveLawVector;

    return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************
Element::Pointer MPMUpdatedLagrangianAxisym::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
) const
{
    return Element::Pointer( new MPMUpdatedLagrangianAxisym( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

Element::Pointer MPMUpdatedLagrangianAxisym::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
) const
{
    return Kratos::make_intrusive< MPMUpdatedLagrangianAxisym >(NewId, pGeom, pProperties);
}

//************************************CLONE*******************************************
//************************************************************************************
Element::Pointer MPMUpdatedLagrangianAxisym::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
) const
{
    MPMUpdatedLagrangianAxisym NewElement (NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewElement.mMP = mMP;

    NewElement.mConstitutiveLawVector = mConstitutiveLawVector->Clone();

    NewElement.mDeformationGradientF0 = mDeformationGradientF0;

    NewElement.mDeterminantF0 = mDeterminantF0;

    return Element::Pointer( new MPMUpdatedLagrangianAxisym(NewElement) );
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************
MPMUpdatedLagrangianAxisym::~MPMUpdatedLagrangianAxisym()
{
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianAxisym::InitializeGeneralVariables(
    GeneralVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo
)
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    const SizeType def_grad_dim = 3;

    rVariables.detF  = 1;

    rVariables.detF0 = 1;

    rVariables.detFT = 1;

    rVariables.B.resize(strain_size, number_of_nodes * dimension, false );

    rVariables.F.resize(def_grad_dim, def_grad_dim, false );

    rVariables.F0.resize(def_grad_dim, def_grad_dim, false );

    rVariables.FT.resize(def_grad_dim, def_grad_dim, false );

    rVariables.ConstitutiveMatrix.resize(strain_size, strain_size, false );

    rVariables.StrainVector.resize(strain_size, false );

    rVariables.StressVector.resize(strain_size, false );

    rVariables.DN_DX.resize( number_of_nodes, dimension, false );

    // CurrentDisp is the unknown variable. It represents the nodal delta displacement. When it is predicted is equal to zero.
    rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);
}

//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************
// AISYMM
void MPMUpdatedLagrangianAxisym::CalculateKinematics(
    GeneralVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY

    MPMUpdatedLagrangian::CalculateKinematics(rVariables, rCurrentProcessInfo);

    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

    rVariables.CurrentRadius = MPMMathUtilities<double>::CalculateRadius(r_N, GetGeometry());
    rVariables.ReferenceRadius = MPMMathUtilities<double>::CalculateRadius(r_N, GetGeometry(), Initial);

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianAxisym::CalculateDeformationMatrix(
    Matrix& rB,
    const Matrix& rDN_DX
)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    const Matrix& rN = GetGeometry().ShapeFunctionsValues();

    rB.clear(); // Set all components to zero

    const double radius = MPMMathUtilities<double>::CalculateRadius(rN, GetGeometry());

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const unsigned int index = dimension * i;

        rB(0, index + 0) = rDN_DX(i, 0);
        rB(1, index + 1) = rDN_DX(i, 1);
        rB(2, index + 0) = rN(0, i) / radius;
        rB(3, index + 0) = rDN_DX(i, 1);
        rB(3, index + 1) = rDN_DX(i, 0);
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianAxisym::CalculateAndAddKuug(
    MatrixType& rLeftHandSideMatrix,
    GeneralVariables& rVariables,
    const double& rIntegrationWeight
)
{
    KRATOS_TRY

    // Axisymmetric geometric matrix
    double alpha_1;
    double alpha_2;
    double alpha_3;

    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

    const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int index_i = 0;
    const double radius = MPMMathUtilities<double>::CalculateRadius(r_N, GetGeometry());

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        unsigned int index_j = 0;
        for (unsigned int j = 0; j < number_of_nodes; j++)
        {
            alpha_1 = rVariables.DN_DX(j, 0) * (rVariables.DN_DX(i, 0) * rVariables.StressVector[0] + rVariables.DN_DX(i, 1) * rVariables.StressVector[3]);
            alpha_2 = rVariables.DN_DX(j, 1) * (rVariables.DN_DX(i, 0) * rVariables.StressVector[3] + rVariables.DN_DX(i, 1) * rVariables.StressVector[1]);
            alpha_3 = r_N(0, i) * r_N(0, j) * rVariables.StressVector[2] * (1.0 / radius * radius);

            rLeftHandSideMatrix(index_i, index_j) += (alpha_1 + alpha_2 + alpha_3) * rIntegrationWeight;
            rLeftHandSideMatrix(index_i + 1, index_j + 1) += (alpha_1 + alpha_2) * rIntegrationWeight;

            index_j += 2;
        }
        index_i += 2;
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianAxisym::CalculateDeformationGradient(
    const Matrix& rDN_DX,
    Matrix& rF,
    Matrix& rDisplacement
)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(GetGeometry().WorkingSpaceDimension() != 2) << "Dimension given is wrong!" << std::endl;

    // Compute radius
    const double current_radius = MPMMathUtilities<double>::CalculateRadius(GetGeometry().ShapeFunctionsValues(), GetGeometry());
    const double initial_radius = MPMMathUtilities<double>::CalculateRadius(GetGeometry().ShapeFunctionsValues(), GetGeometry(), Initial);

    rF = IdentityMatrix(3);

    for (IndexType i = 0; i < GetGeometry().PointsNumber(); ++i)
    {
        rF(0, 0) += rDisplacement(i, 0) * rDN_DX(i, 0);
        rF(0, 1) += rDisplacement(i, 0) * rDN_DX(i, 1);
        rF(1, 0) += rDisplacement(i, 1) * rDN_DX(i, 0);
        rF(1, 1) += rDisplacement(i, 1) * rDN_DX(i, 1);
    }

    rF(2, 2) = current_radius / initial_radius;

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianAxisym::Initialize(
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        // Initialize parameters
        mDeterminantF0 = 1;
        mDeformationGradientF0 = IdentityMatrix(3);

        // Initialize constitutive law and materials
        InitializeMaterial(rCurrentProcessInfo);
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianAxisym::save(
    Serializer& rSerializer
) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MPMUpdatedLagrangianAxisym )
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianAxisym::load(
    Serializer& rSerializer
)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MPMUpdatedLagrangianAxisym )
}

} // Namespace Kratos
