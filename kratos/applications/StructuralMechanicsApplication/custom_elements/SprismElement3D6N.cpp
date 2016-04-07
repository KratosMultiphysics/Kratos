// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

// System includes

// External includes
#include <cmath> // or #include <math.h>

// Project includes
#include "includes/define.h"
#include "custom_elements/SprismElement3D6N.hpp"
#include "custom_utilities/structural_mechanics_math_utilities.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "structural_mechanics_application.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
    // ------------------------------------------------------------------------- //
    // ----------------------------- UTILITIES --------------------------------- //
    // ------------------------------------------------------------------------- //
    namespace Utilities
    {

    } // Namespace Utilities.

/**
 * Flags related to the element computation
 */
KRATOS_CREATE_LOCAL_FLAG( SprismElement3D6N, COMPUTE_RHS_VECTOR,                 0 );
KRATOS_CREATE_LOCAL_FLAG( SprismElement3D6N, COMPUTE_LHS_MATRIX,                 1 );
KRATOS_CREATE_LOCAL_FLAG( SprismElement3D6N, COMPUTE_RHS_VECTOR_WITH_COMPONENTS, 2 );
KRATOS_CREATE_LOCAL_FLAG( SprismElement3D6N, COMPUTE_LHS_MATRIX_WITH_COMPONENTS, 3 );

// ------------------------------------------------------------------------- //
// ------------------------------ PUBLIC ----------------------------------- //
// ------------------------------------------------------------------------- //

/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

SprismElement3D6N::SprismElement3D6N( )
        : Element( )
{
    //DO NOT CALL IT: only needed for Register and Serialization!!!
}

/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

SprismElement3D6N::SprismElement3D6N(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

SprismElement3D6N::SprismElement3D6N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
    mFinalizedStep = true; // the creation is out of the time step, it must be true

    if( GetProperties().Has(NINT_TRANS) )
    {
        if (GetProperties()[NINT_TRANS] == 2)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
        }
        else if (GetProperties()[NINT_TRANS] == 3)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
        }
        else if (GetProperties()[NINT_TRANS] == 5)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
        }
        else if (GetProperties()[NINT_TRANS] == 7)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_4;
        }
        else if (GetProperties()[NINT_TRANS] == 11)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_5;
        }
        else
        {
            std::cout << "The number of integration points is not defined.  NINT_TRANS: "<< GetProperties()[NINT_TRANS] << std::endl;
            std::cout << "Options are: 2, 3, 5, 7, 11  " << std::endl;
            std::cout << "Taking default number of integration points (NINT_TRANS = 2)  " << std::endl;
            mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
        }
    }
    else
    {
        mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
    }

    //DO NOT ADD DOFS HERE!!!
}

/*********************************** COPY CONSTRUCTOR ******************************/
/***********************************************************************************/

SprismElement3D6N::SprismElement3D6N( SprismElement3D6N const& rOther)
    :Element(rOther)
    ,mThisIntegrationMethod(rOther.mThisIntegrationMethod)
    ,mConstitutiveLawVector(rOther.mConstitutiveLawVector)
    ,mFinalizedStep(rOther.mFinalizedStep)
    ,mid_vec(rOther.mid_vec)
    ,mTotalDomainInitialSize(rOther.mTotalDomainInitialSize)
    ,mInvJ0(rOther.mInvJ0)
    ,mDetJ0(rOther.mDetJ0)
{
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

SprismElement3D6N::~SprismElement3D6N()
{
}

/********************************** ASSIGMENT OPERATOR *****************************/
/***********************************************************************************/

SprismElement3D6N&  SprismElement3D6N::operator=(SprismElement3D6N const& rOther)
{
    SprismElement3D6N::operator=(rOther);

    mThisIntegrationMethod = rOther.mThisIntegrationMethod;

    mConstitutiveLawVector.clear();
    mConstitutiveLawVector.resize( rOther.mConstitutiveLawVector.size(), false );

    mInvJ0.clear();
    mInvJ0.resize( rOther.mInvJ0.size());

    for(unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
    {
        mConstitutiveLawVector[i] = rOther.mConstitutiveLawVector[i];
        mInvJ0[i]=rOther.mInvJ0[i];
    }

    mTotalDomainInitialSize = rOther.mTotalDomainInitialSize;
    mDetJ0 = rOther.mDetJ0;

    mid_vec = rOther.mid_vec;

    return *this;
}

/************************************* OPERATIONS **********************************/
/***********************************************************************************/

Element::Pointer SprismElement3D6N::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SprismElement3D6N(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

/*********************************** CLONE ******************************************/
/************************************************************************************/

Element::Pointer SprismElement3D6N::Clone(
        IndexType NewId,
        NodesArrayType const& rThisNodes) const
{
    SprismElement3D6N NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewElement.mThisIntegrationMethod = mThisIntegrationMethod;

    if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size())
    {
        NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size(), false);
    }

    if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
    {
        KRATOS_THROW_ERROR( std::logic_error, "Constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size() );
    }
    
    for(unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
    {
      NewElement.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
    }

    //-----------//

    if ( NewElement.mInvJ0.size() != mInvJ0.size() )
    {
        NewElement.mInvJ0.resize(mInvJ0.size());
    }

    for(unsigned int i = 0; i < mInvJ0.size(); i++)
    {
        NewElement.mInvJ0[i] = mInvJ0[i];
    }

    NewElement.mTotalDomainInitialSize = mTotalDomainInitialSize;
    NewElement.mDetJ0 = mDetJ0;

    NewElement.mid_vec = mid_vec;

    return Element::Pointer( new SprismElement3D6N(NewElement));
}

//******************************* GETTING METHODS *********************************//
/***********************************************************************************/
/***********************************************************************************/

SprismElement3D6N::IntegrationMethod SprismElement3D6N::GetIntegrationMethod() const
{
    return mThisIntegrationMethod;
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
        )
{
    KRATOS_TRY;

    WeakPointerVector< Node < 3 > >& nodal_neigb = this->GetValue(NEIGHBOUR_NODES);

    unsigned int number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(nodal_neigb);
    unsigned int dim = number_of_nodes * 3;

    if (rResult.size() != dim)
    {
        rResult.resize(dim, false);
    }

    // Nodes of the central element
    unsigned int index = 0;
    for (unsigned int i = 0; i < 6; i++)
    {
        rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
        index += 3;
    }

    // Adding the ids of the neighbouring nodes
    for (unsigned int i = 0; i < 6; i++)
    {
        if (HasNeighbour(i, nodal_neigb[i]))
        {
            rResult[index]     = nodal_neigb[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = nodal_neigb[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = nodal_neigb[i].GetDof(DISPLACEMENT_Z).EquationId();
            index += 3;
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& CurrentProcessInfo
        )
{
    KRATOS_TRY;

    WeakPointerVector< Node < 3 > >& nodal_neigb = this->GetValue(NEIGHBOUR_NODES);
    rElementalDofList.resize(0);

    // Nodes of the central element
    for (unsigned int i = 0; i < GetGeometry().size(); i++)
    {
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
    }

    // Adding the dofs of the neighbouring nodes
    for (unsigned int i = 0; i < 6; i++)
    {
        if (HasNeighbour(i, nodal_neigb[i]))
        {
            rElementalDofList.push_back(nodal_neigb[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(nodal_neigb[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(nodal_neigb[i].pGetDof(DISPLACEMENT_Z));
        }
    }

    KRATOS_CATCH("");
}

/******************************** DISPLACEMENT **************************************/
/************************************************************************************/

void SprismElement3D6N::GetValuesVector(
        Vector& rValues,
        int Step
        )
{
    WeakPointerVector< Node < 3 > >& nodal_neigb = this->GetValue(NEIGHBOUR_NODES);
    unsigned int number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(nodal_neigb);

    const unsigned int MatSize = number_of_nodes * 3;
    if (rValues.size() != MatSize)
    {
        rValues.resize(MatSize, false);
    }

    unsigned int index = 0;

    // Nodes of the central element
    for (unsigned int i = 0; i < 6; i++)
    {
        const array_1d<double, 3 > & disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        rValues[index]     = disp[0];
        rValues[index + 1] = disp[1];
        rValues[index + 2] = disp[2];
        index += 3;
    }

    // Neighbour nodes
    for (int i = 0; i < 6; i++)
    {
        if (HasNeighbour(i, nodal_neigb[i]))
        {
            const array_1d<double, 3 > & disp = nodal_neigb[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            rValues[index]     = disp[0];
            rValues[index + 1] = disp[1];
            rValues[index + 2] = disp[2];
            index += 3;
        }
    }
}

/********************************** VELOCITY ****************************************/
/************************************************************************************/

void SprismElement3D6N::GetFirstDerivativesVector(
        Vector& rValues,
        int Step
        )
{
    WeakPointerVector< Node < 3 > >& nodal_neigb = this->GetValue(NEIGHBOUR_NODES);
    unsigned int number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(nodal_neigb);

    const unsigned int MatSize = number_of_nodes * 3;
    if (rValues.size() != MatSize)
    {
        rValues.resize(MatSize, false);
    }

    unsigned int index = 0;

    // Nodes of the central element
    for (unsigned int i = 0; i < 6; i++)
    {
        const array_1d<double, 3 > & vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
        rValues[index]     = vel[0];
        rValues[index + 1] = vel[1];
        rValues[index + 2] = vel[2];
        index += 3;
    }

    // Neighbour nodes
    for (int i = 0; i < 6; i++)
    {
        if (HasNeighbour(i, nodal_neigb[i]))
        {
            const array_1d<double, 3 > & vel = nodal_neigb[i].FastGetSolutionStepValue(VELOCITY, Step);
            rValues[index]     = vel[0];
            rValues[index + 1] = vel[1];
            rValues[index + 2] = vel[2];
            index += 3;
        }
    }
}

/******************************** ACCELERATION **************************************/
/************************************************************************************/

void SprismElement3D6N::GetSecondDerivativesVector(
        Vector& rValues,
        int Step
        )
{
    WeakPointerVector< Node < 3 > >& nodal_neigb = this->GetValue(NEIGHBOUR_NODES);
    unsigned int number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(nodal_neigb);

    const unsigned int MatSize = number_of_nodes * 3;
    if (rValues.size() != MatSize)
    {
        rValues.resize(MatSize, false);
    }

    unsigned int index = 0;

    // Nodes of the central element
    for (unsigned int i = 0; i < 6; i++)
    {
        const array_1d<double, 3 > & acc = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
        rValues[index]     = acc[0];
        rValues[index + 1] = acc[1];
        rValues[index + 2] = acc[2];
        index += 3;
    }

    // Neighbour nodes
    for (int i = 0; i < 6; i++)
    {
        if (HasNeighbour(i, nodal_neigb[i]))
        {
            const array_1d<double, 3 > & acc = nodal_neigb[i].FastGetSolutionStepValue(ACCELERATION, Step);
            rValues[index]     = acc[0];
            rValues[index + 1] = acc[1];
            rValues[index + 2] = acc[2];
            index += 3;
        }
    }
}

//****************************** COMPUTING METHODS ********************************//
/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo
        )
{
    KRATOS_TRY;

    /* Create local system components */
    LocalSystemComponents LocalSystem;

    /* Calculation flags */
    LocalSystem.CalculationFlags.Set(SprismElement3D6N::COMPUTE_RHS_VECTOR);

    MatrixType LeftHandSideMatrix = Matrix();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(rRightHandSideVector);

    /* Calculate elemental system */
    CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateRightHandSide(
        std::vector< VectorType >& rRightHandSideVectors,
        const std::vector< Variable< VectorType > >& rRHSVariables,
        ProcessInfo& rCurrentProcessInfo
        )
{
    KRATOS_TRY;

    /* Create local system components */
    LocalSystemComponents LocalSystem;

    /* Calculation flags */
    LocalSystem.CalculationFlags.Set(SprismElement3D6N::COMPUTE_RHS_VECTOR);
    LocalSystem.CalculationFlags.Set(SprismElement3D6N::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);

    MatrixType LeftHandSideMatrix = Matrix();

    /* Initialize sizes for the system components: */
    if( rRHSVariables.size() != rRightHandSideVectors.size() )
    {
        rRightHandSideVectors.resize(rRHSVariables.size());
    }

    for( unsigned int i = 0; i < rRightHandSideVectors.size(); i++ )
    {
        this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVectors[i], LocalSystem.CalculationFlags );
    }

    /* Set Variables to Local system components */
    LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
    LocalSystem.SetRightHandSideVectors(rRightHandSideVectors);

    LocalSystem.SetRightHandSideVariables(rRHSVariables);

    /* Calculate elemental system */
    CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo
        )
{
    KRATOS_TRY;

    /* Create local system components */
    LocalSystemComponents LocalSystem;

    /* Calculation flags */
    LocalSystem.CalculationFlags.Set(SprismElement3D6N::COMPUTE_LHS_MATRIX, true);
    LocalSystem.CalculationFlags.Set(SprismElement3D6N::COMPUTE_RHS_VECTOR, false);

    VectorType RightHandSideVector = Vector();

    /* Initialize sizes for the system components: */
    this->InitializeSystemMatrices( rLeftHandSideMatrix, RightHandSideVector, LocalSystem.CalculationFlags );

    /* Set Variables to Local system components */
    LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(RightHandSideVector);

    /* Calculate elemental system */
    CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo
        )
{
    KRATOS_TRY;

    /* Create local system components */
    LocalSystemComponents LocalSystem;

    /* Calculation flags */
    LocalSystem.CalculationFlags.Set(SprismElement3D6N::COMPUTE_LHS_MATRIX, true);
    LocalSystem.CalculationFlags.Set(SprismElement3D6N::COMPUTE_RHS_VECTOR, true);

    /* Initialize sizes for the system components: */
    this->InitializeSystemMatrices( rLeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

    /* Set Variables to Local system components */
    LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(rRightHandSideVector);

    /* Calculate elemental system */
    CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateLocalSystem(
        std::vector< MatrixType >& rLeftHandSideMatrices,
        const std::vector< Variable< MatrixType > >& rLHSVariables,
        std::vector< VectorType >& rRightHandSideVectors,
        const std::vector< Variable< VectorType > >& rRHSVariables,
        ProcessInfo& rCurrentProcessInfo
        )
{
    KRATOS_TRY;

    /* Create local system components */
    LocalSystemComponents LocalSystem;

    /* Calculation flags*/
    LocalSystem.CalculationFlags.Set(SprismElement3D6N::COMPUTE_LHS_MATRIX_WITH_COMPONENTS);
    LocalSystem.CalculationFlags.Set(SprismElement3D6N::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);

    /* Initialize sizes for the system components: */
    if( rLHSVariables.size() != rLeftHandSideMatrices.size() )
    {
        rLeftHandSideMatrices.resize(rLHSVariables.size());
    }

    if( rRHSVariables.size() != rRightHandSideVectors.size() )
    {
        rRightHandSideVectors.resize(rRHSVariables.size());
    }

    LocalSystem.CalculationFlags.Set(SprismElement3D6N::COMPUTE_LHS_MATRIX);
    for( unsigned int i = 0; i < rLeftHandSideMatrices.size(); i++ )
    {
        this->InitializeSystemMatrices( rLeftHandSideMatrices[i], rRightHandSideVectors[0], LocalSystem.CalculationFlags );
    }

    LocalSystem.CalculationFlags.Set(SprismElement3D6N::COMPUTE_RHS_VECTOR, true);
    LocalSystem.CalculationFlags.Set(SprismElement3D6N::COMPUTE_LHS_MATRIX, false);

    for( unsigned int i = 0; i < rRightHandSideVectors.size(); i++ )
    {
        this->InitializeSystemMatrices( rLeftHandSideMatrices[0], rRightHandSideVectors[i], LocalSystem.CalculationFlags );
    }

    LocalSystem.CalculationFlags.Set(SprismElement3D6N::COMPUTE_LHS_MATRIX, true);

    /* Set Variables to Local system components */
    LocalSystem.SetLeftHandSideMatrices(rLeftHandSideMatrices);
    LocalSystem.SetRightHandSideVectors(rRightHandSideVectors);

    LocalSystem.SetLeftHandSideVariables(rLHSVariables);
    LocalSystem.SetRightHandSideVariables(rRHSVariables);

    /* Calculate elemental system */
    CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo
        )
{
    KRATOS_TRY;

    WeakPointerVector< Node < 3 > >& nodal_neigb = this->GetValue(NEIGHBOUR_NODES);
        
    double Density = GetProperties()[DENSITY];

    unsigned int number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(nodal_neigb);
    unsigned int MatSize = number_of_nodes * 3;

    if (rMassMatrix.size1() != MatSize)
    {
        rMassMatrix.resize(MatSize, MatSize, false);
    }
    
    noalias(rMassMatrix) = ZeroMatrix(MatSize, MatSize);
    
    double Volume = GetGeometry().Volume();
    double TotalMass = Volume * Density;

    bool ComputeLumpedMassMatrix = false;
    if( rCurrentProcessInfo.Has(COMPUTE_LUMPED_MASS_MATRIX) )
    {
        if(rCurrentProcessInfo[COMPUTE_LUMPED_MASS_MATRIX] == true)
        {
            ComputeLumpedMassMatrix = true;
        }
    }

    // LUMPED MASS MATRIX, this one is easy because each node receives the same proportion of mass
    if (ComputeLumpedMassMatrix == true)
    {
        Vector LumpFact;
        GetGeometry().LumpingFactors(LumpFact);
        for (unsigned int i = 0; i < 3; i++)
        {
            double temp = LumpFact[i] * TotalMass;
            for (unsigned int j = 0; j < 6; j++)
            {
                unsigned int index = i * 6 + j;
                rMassMatrix(index, index) = temp;
            }
        }
    }
    // CONSISTENT MASS MATRIX
    else
    {
        // Create local system components
        LocalSystemComponents LocalSystem;

        // Calculation flags
        LocalSystem.CalculationFlags.Set(SprismElement3D6N::COMPUTE_LHS_MATRIX);

        VectorType RightHandSideVector = Vector();

        // Initialize sizes for the system components:
        this->InitializeSystemMatrices( rMassMatrix, RightHandSideVector,  LocalSystem.CalculationFlags );

        // Set Variables to Local system components
        LocalSystem.SetLeftHandSideMatrix(rMassMatrix);
        LocalSystem.SetRightHandSideVector(RightHandSideVector);

        // Calculate elemental system
        CalculateDynamicSystem( LocalSystem, rCurrentProcessInfo );
    }
    
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        ProcessInfo& rCurrentProcessInfo
        )
{
    KRATOS_TRY;

    WeakPointerVector< Node < 3 > >& nodal_neigb = this->GetValue(NEIGHBOUR_NODES);

    // 0.-Initialize the DampingMatrix:
    const unsigned int number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(nodal_neigb);
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * dimension;

    if ( rDampingMatrix.size1() != MatSize )
    {
        rDampingMatrix.resize( MatSize, MatSize, false );
    }
    
    noalias( rDampingMatrix ) = ZeroMatrix( MatSize, MatSize );

    // 1.-Calculate StiffnessMatrix:

    MatrixType StiffnessMatrix  = Matrix();

    this->CalculateLeftHandSide( StiffnessMatrix, rCurrentProcessInfo );

    // 2.-Calculate MassMatrix:

    MatrixType MassMatrix  = Matrix();

    this->CalculateMassMatrix ( MassMatrix, rCurrentProcessInfo );

    // 3.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
    double alpha = 0.0;
    if( GetProperties().Has(RAYLEIGH_ALPHA) )
    {
        alpha = GetProperties()[RAYLEIGH_ALPHA];
    }
    else if( rCurrentProcessInfo.Has(RAYLEIGH_ALPHA) )
    {
        alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];
    }

    double beta  = 0.0;
    if( GetProperties().Has(RAYLEIGH_BETA) )
    {
        beta = GetProperties()[RAYLEIGH_BETA];
    }
    else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) )
    {
        beta = rCurrentProcessInfo[RAYLEIGH_BETA];
    }

    // 4.-Compose the Damping Matrix:

    // Rayleigh Damping Matrix: alpha*M + beta*K
    noalias( rDampingMatrix ) += alpha * MassMatrix;
    noalias( rDampingMatrix ) += beta  * StiffnessMatrix;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const MatrixType& rStiffnessMatrix,
        const MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo
        )
{
    KRATOS_TRY;

    WeakPointerVector< Node < 3 > >& nodal_neigb = this->GetValue(NEIGHBOUR_NODES);

    // 0.-Initialize the DampingMatrix:
    const unsigned int number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(nodal_neigb);
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * dimension;

    if ( rDampingMatrix.size1() != MatSize )
    {
        rDampingMatrix.resize( MatSize, MatSize, false );
    }
    
    noalias( rDampingMatrix ) = ZeroMatrix( MatSize, MatSize );

    // 1.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
    double alpha = 0.0;
    if( GetProperties().Has(RAYLEIGH_ALPHA) )
    {
        alpha = GetProperties()[RAYLEIGH_ALPHA];
    }
    else if( rCurrentProcessInfo.Has(RAYLEIGH_ALPHA) )
    {
        alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];
    }

    double beta  = 0.0;
    if( GetProperties().Has(RAYLEIGH_BETA) )
    {
        beta = GetProperties()[RAYLEIGH_BETA];
    }
    else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) )
    {
        beta = rCurrentProcessInfo[RAYLEIGH_BETA];
    }

    // 2.-Compose the Damping Matrix:

    // Rayleigh Damping Matrix: alpha*M + beta*K
    noalias( rDampingMatrix ) += alpha * rMassMatrix;
    noalias( rDampingMatrix ) += beta  * rStiffnessMatrix;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::AddExplicitContribution(
        const VectorType& rRHSVector,
        const Variable<VectorType>& rRHSVariable,
        Variable<array_1d<double,3> >& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo
        )
{
    KRATOS_TRY;

    WeakPointerVector< Node < 3 > >& nodal_neigb = this->GetValue(NEIGHBOUR_NODES);

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();// + NumberOfActiveNeighbours(nodal_neigb);
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    if( rRHSVariable == EXTERNAL_FORCES_VECTOR && rDestinationVariable == EXTERNAL_FORCE )
    {
        for(unsigned int i = 0; i < number_of_nodes; i++)
        {
            int index = dimension * i;

            if (i < 6)
            {
                GetGeometry()[i].SetLock();
                array_1d<double, 3 > &ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_FORCE);

                for(unsigned int j = 0; j < dimension; j++)
                {
                    ExternalForce[j] += rRHSVector[index + j];
                }

                GetGeometry()[i].UnSetLock();
            }
            else
            {
                nodal_neigb[i].SetLock();
                array_1d<double, 3 > &ExternalForce = nodal_neigb[i].FastGetSolutionStepValue(EXTERNAL_FORCE);

                for(unsigned int j = 0; j < dimension; j++)
                {
                    ExternalForce[j] += rRHSVector[index + j];
                }

                nodal_neigb[i].UnSetLock();
            }
        }
    }

    if( rRHSVariable == INTERNAL_FORCES_VECTOR && rDestinationVariable == INTERNAL_FORCE )
    {
        for(unsigned int i = 0; i < number_of_nodes; i++)
        {
            int index = dimension * i;

            if (i < 6)
            {
                GetGeometry()[i].SetLock();
                array_1d<double, 3 > &InternalForce = GetGeometry()[i].FastGetSolutionStepValue(INTERNAL_FORCE);

                for(unsigned int j = 0; j < dimension; j++)
                {
                    InternalForce[j] += rRHSVector[index + j];
                }

                GetGeometry()[i].UnSetLock();
            }
            else
            {
                nodal_neigb[i].SetLock();
                array_1d<double, 3 > &InternalForce = nodal_neigb[i].FastGetSolutionStepValue(INTERNAL_FORCE);

                for(unsigned int j = 0; j < dimension; j++)
                {
                    InternalForce[j] += rRHSVector[index + j];
                }

                nodal_neigb[i].UnSetLock();
            }
        }
    }

    if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL )
    {
        for(unsigned int i = 0; i < number_of_nodes; i++)
        {
            int index = dimension * i;

            if (i < 6)
            {
                GetGeometry()[i].SetLock();
                array_1d<double, 3 > &ForceResidual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

                for(unsigned int j = 0; j < dimension; j++)
                {
                    ForceResidual[j] += rRHSVector[index + j];
                }

                GetGeometry()[i].UnSetLock();
            }
            else
            {
                nodal_neigb[i].SetLock();
                array_1d<double, 3 > &ForceResidual = nodal_neigb[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

                for(unsigned int j = 0; j < dimension; j++)
                {
                    ForceResidual[j] += rRHSVector[index + j];
                }

                nodal_neigb[i].UnSetLock();
            }
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        )
{
    KRATOS_TRY;

    const unsigned int& integration_point_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != integration_point_number )
    {
        rOutput.resize( integration_point_number, false );
    }

    if ( rVariable == VON_MISES_STRESS )
    {
        /* Create and initialize element variables: */
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables);

        /* Create constitutive law parameters: */
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        /* Set constitutive law flags: */
        Flags &ConstitutiveLawOptions = Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        double& alpha_eas = Element::GetValue(ALPHA_EAS);

        /* Calculate common components (B, C) */
        this->CalculateCommonComponents();

        // Reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            double zeta_gauss = 2.0 * integration_points[PointNumber].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(Variables,PointNumber, alpha_eas, zeta_gauss);
            this->CbartoFbar(Variables, PointNumber);

            // To take in account previous step writing
            if( mFinalizedStep )
            {
                this->GetHistoricalVariables(Variables,PointNumber);
            }

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(Variables,Values,PointNumber);

            // Call the constitutive law to update material variables
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy (Values);

            ComparisonUtilities EquivalentStress;
            rOutput[PointNumber] =  EquivalentStress.CalculateVonMises(Variables.StressVector);
        }
    }
    else if ( rVariable == NORM_ISOCHORIC_STRESS )
    {
        // Create and initialize element variables:
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY, true);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        double& alpha_eas = Element::GetValue(ALPHA_EAS);

        /* Calculate common components (B, C) */
        this->CalculateCommonComponents();

        // Reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            double zeta_gauss = 2.0 * integration_points[PointNumber].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(Variables,PointNumber, alpha_eas, zeta_gauss);
            this->CbartoFbar(Variables, PointNumber);

            // To take in account previous step writing
            if( mFinalizedStep )
            {
                this->GetHistoricalVariables(Variables,PointNumber);
            }

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(Variables,Values,PointNumber);

            // Call the constitutive law to update material variables
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy (Values);

            ComparisonUtilities EquivalentStress;
            rOutput[PointNumber] =  EquivalentStress.CalculateStressNorm(Variables.StressVector);
        }
    }
    else
    {
        for ( unsigned int ii = 0; ii < integration_point_number; ii++ )
        {
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rOutput[ii] );
        }
    }

    if ( rOutput.size() != 6 )
    {
        std::vector<double> rOutput_aux;
        rOutput_aux = rOutput;

        rOutput.resize( 6, false );
        Matrix interpol = StructuralMechanicsMathUtilities::interpol_PrismGiD(integration_point_number);

        for (unsigned int iii = 0; iii < 6; iii++)
        {
            rOutput[iii] = 0.0;

            for (unsigned int Gauss_Point = 0; Gauss_Point < integration_point_number; Gauss_Point++)
            {
                rOutput[iii] += interpol(Gauss_Point, iii) * rOutput_aux[Gauss_Point];
            }
        }
    }

//    if ( rOutput.size() != 3 * integration_point_number )
//    {
//        std::vector<double> rOutput_aux;
//        rOutput_aux = rOutput;

//        rOutput.resize( 3 * integration_point_number, false );

//        for (unsigned int Gauss_Point = 0; Gauss_Point < integration_point_number; Gauss_Point++)
//        {
//            for (unsigned int iii = 0; iii < 3; iii++)
//            {
//                rOutput[Gauss_Point * 3 + iii] = rOutput_aux[Gauss_Point];
//            }
//        }
//    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        )
{
    KRATOS_TRY;

    const unsigned int& integration_point_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != integration_point_number )
    {
        rOutput.resize( integration_point_number );
    }

    if ( rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR )
    {
        // Create and initialize element variables:
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

//        /* TEST */
//        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
//        /* TEST */

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        double& alpha_eas = Element::GetValue(ALPHA_EAS);

        /* Calculate common components (B, C) */
        this->CalculateCommonComponents();

        // Reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            double zeta_gauss = 2.0 * integration_points[PointNumber].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(Variables,PointNumber, alpha_eas, zeta_gauss);
            this->CbartoFbar(Variables, PointNumber);

            // To take in account previous step writing
            if( mFinalizedStep )
            {
                this->GetHistoricalVariables(Variables,PointNumber);
            }

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(Variables, Values, PointNumber);

            // Call the constitutive law to update material variables
            if( rVariable == CAUCHY_STRESS_VECTOR)
            {
                Variables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;
            }
            else
            {
                Variables.StressMeasure = ConstitutiveLaw::StressMeasure_PK2;
            }

            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(Values, Variables.StressMeasure);

//            /* TEST */
////            this->CalculateGreenLagrangeStrain(Variables.C,Variables.StrainVector);
////            this->CalculateHenckyStrain(Variables.C,Variables.StrainVector);
////            this->CalculateLinearStress(Variables);
////            this->CalculateLinearIsotropicStress(Variables);
////            this->LinearConstitutiveMatrix(Variables);
////            this->CalculateLinearStress(Variables);
//            this->CalculateGreenLagrangeStrain(Variables.C,Variables.StrainVector);
//            this->CalculateHyperelasticNeoHookeanStress(Variables);
////            this->CalculateHenckyStrain(Variables.C,Variables.StrainVector);
////            this->CalculateLogStress(Variables);
//            /* TEST */

            if (rOutput[PointNumber].size() != Variables.StressVector.size())
            {
                rOutput[PointNumber].resize( Variables.StressVector.size(), false);
            }
            rOutput[PointNumber] = Variables.StressVector;
        }
    }
    else if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR  || rVariable == ALMANSI_STRAIN_VECTOR || rVariable == HENCKY_STRAIN_VECTOR)
    {
        // Create and initialize element variables:
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        double& alpha_eas = Element::GetValue(ALPHA_EAS);

        /* Calculate common components (B, C) */
        this->CalculateCommonComponents();

        // Reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            double zeta_gauss = 2.0 * integration_points[PointNumber].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(Variables,PointNumber, alpha_eas, zeta_gauss);

            // To take in account previous step writing
            if( mFinalizedStep )
            {
                this->GetHistoricalVariables(Variables,PointNumber);
            }

            // Compute Green-Lagrange Strain
            if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR )
            {
                this->CalculateGreenLagrangeStrain( Variables.C, Variables.StrainVector );
            }
            else if( rVariable == ALMANSI_STRAIN_VECTOR )
            {
                this->CbartoFbar(Variables, PointNumber);
                this->CalculateAlmansiStrain( Variables.F, Variables.StrainVector );
            }
            else if( rVariable == HENCKY_STRAIN_VECTOR )
            {
                this->CalculateHenckyStrain( Variables.C, Variables.StrainVector );
            }

            if (rOutput[PointNumber].size() != Variables.StrainVector.size())
            {
                rOutput[PointNumber].resize( Variables.StrainVector.size(), false );
            }

            rOutput[PointNumber] = Variables.StrainVector;
        }
    }
    else
    {
        for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
        {
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue(rVariable , rOutput[ii]);
        }
    }

    if ( rOutput.size() != 6 )
    {
        std::vector<Vector> rOutput_aux;
        rOutput_aux = rOutput;

        rOutput.resize( 6 );
        Matrix interpol = StructuralMechanicsMathUtilities::interpol_PrismGiD(integration_point_number);

        for (unsigned int iii = 0; iii < 6; iii++)
        {
            rOutput[iii] = ZeroVector(rOutput[0].size());

            for (unsigned int Gauss_Point = 0; Gauss_Point < integration_point_number; Gauss_Point++)
            {
                rOutput[iii] += interpol(Gauss_Point, iii) * rOutput_aux[Gauss_Point];
            }
        }
    }

//    if ( rOutput.size() != 3  *integration_point_number )
//    {
//        std::vector<Vector> rOutput_aux;
//        rOutput_aux = rOutput;

//        rOutput.resize( 3 * integration_point_number );

//        for (unsigned int Gauss_Point = 0; Gauss_Point < integration_point_number; Gauss_Point++)
//        {
//            for (unsigned int iii = 0; iii < 3; iii++)
//            {
//                rOutput[Gauss_Point * 3 + iii] =  rOutput_aux[Gauss_Point];
//            }
//        }
//    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateOnIntegrationPoints(
        const Variable<Matrix >& rVariable,
        std::vector< Matrix >& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        )
{
    KRATOS_TRY;

    const unsigned int& integration_point_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != integration_point_number )
    {
        rOutput.resize( integration_point_number );
    }

    if ( rVariable == CAUCHY_STRESS_TENSOR || rVariable == PK2_STRESS_TENSOR )
    {
        std::vector<Vector> StressVector;
        if( rVariable == CAUCHY_STRESS_TENSOR )
        {
            this->CalculateOnIntegrationPoints( CAUCHY_STRESS_VECTOR, StressVector, rCurrentProcessInfo );
        }
        else
        {
            this->CalculateOnIntegrationPoints( PK2_STRESS_VECTOR, StressVector, rCurrentProcessInfo );
        }

        // Loop integration points
        if ( rOutput.size() != StressVector.size() )
        {
            rOutput.resize( StressVector.size() );
        }

        for ( unsigned int PointNumber = 0; PointNumber < rOutput.size(); PointNumber++ )
        {
            if (rOutput[PointNumber].size2() != 3)
            {
                rOutput[PointNumber].resize(3, 3, false);
            }
            rOutput[PointNumber] = MathUtils<double>::StressVectorToTensor(StressVector[PointNumber]);
        }
    }
    else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR  || rVariable == ALMANSI_STRAIN_TENSOR || rVariable == HENCKY_STRAIN_TENSOR)
    {
        std::vector<Vector> StrainVector;
        if( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
        {
            CalculateOnIntegrationPoints( GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );
        }
        else if ( rVariable == ALMANSI_STRAIN_TENSOR )
        {
            CalculateOnIntegrationPoints( ALMANSI_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );
        }
        else if ( rVariable == HENCKY_STRAIN_TENSOR )
        {
            CalculateOnIntegrationPoints( HENCKY_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );
        }

        // Loop integration points
        if ( rOutput.size() != StrainVector.size() )
        {
            rOutput.resize( StrainVector.size() );
        }

        for ( unsigned int PointNumber = 0; PointNumber < rOutput.size(); PointNumber++ )
        {
            if (rOutput[PointNumber].size2() != 3)
            {
                rOutput[PointNumber].resize(3, 3, false);
            }
            rOutput[PointNumber] = MathUtils<double>::StrainVectorToTensor(StrainVector[PointNumber]);
        }
    }
    else if ( rVariable == CONSTITUTIVE_MATRIX )
    {
        // Create and initialize element variables:
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        double& alpha_eas = Element::GetValue(ALPHA_EAS);

        /* Calculate common components (B, C) */
        this->CalculateCommonComponents();

        // Reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            double zeta_gauss = 2.0 * integration_points[PointNumber].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(Variables,PointNumber, alpha_eas, zeta_gauss);
            this->CbartoFbar(Variables, PointNumber);

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(Variables,Values,PointNumber);

            // Call the constitutive law to update material variables
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);

            if( rOutput[PointNumber].size2() != Variables.ConstitutiveMatrix.size2() )
            {
                rOutput[PointNumber].resize( Variables.ConstitutiveMatrix.size1() , Variables.ConstitutiveMatrix.size2() , false );
            }
            rOutput[PointNumber] = Variables.ConstitutiveMatrix;
        }
    }
    else if ( rVariable == DEFORMATION_GRADIENT )
    {
        // Create and initialize element variables:
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        double& alpha_eas = Element::GetValue(ALPHA_EAS);

        /* Calculate common components (B, C) */
        this->CalculateCommonComponents();

        // Reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            double zeta_gauss = 2.0 * integration_points[PointNumber].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(Variables,PointNumber, alpha_eas, zeta_gauss);
            this->CbartoFbar(Variables, PointNumber);

            if( rOutput[PointNumber].size2() != Variables.F.size2() )
            {
                rOutput[PointNumber].resize( Variables.F.size1() , Variables.F.size2() , false );
            }
            rOutput[PointNumber] = Variables.F;
        }
    }
    else
    {
        for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
        {
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable , rOutput[ii] );
        }
    }

    if ( rOutput.size() != 6 )
    {
        std::vector<Matrix> rOutput_aux;
        rOutput_aux = rOutput;

        rOutput.resize( 6 );
        Matrix interpol = StructuralMechanicsMathUtilities::interpol_PrismGiD(integration_point_number);

        for (unsigned int iii = 0; iii < 6; iii++)
        {
            rOutput[iii] = ZeroMatrix(rOutput[0].size1(), rOutput[0].size2());

            for (unsigned int Gauss_Point = 0; Gauss_Point < integration_point_number; Gauss_Point++)
            {
                rOutput[iii] += interpol(Gauss_Point, iii) * rOutput_aux[Gauss_Point];
            }
        }
    }

//    if ( rOutput.size() != 3 * integration_point_number )
//    {
//        std::vector<Matrix> rOutput_aux;
//        rOutput_aux = rOutput;

//        rOutput.resize( 3 * integration_point_number );

//        for (unsigned int Gauss_Point = 0; Gauss_Point < integration_point_number; Gauss_Point++)
//        {
//            for (unsigned int iii = 0; iii < 3; iii++)
//            {
//                rOutput[Gauss_Point * 3 + iii] = rOutput_aux[Gauss_Point];
//            }
//        }
//    }

    KRATOS_CATCH( "" );
}

//**************************** ON INTEGRATION POINTS ******************************//
/******************************** SET DOUBLE VALUE *********************************/
/***********************************************************************************/

void SprismElement3D6N::SetValueOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        )
{
    for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
    {
        mConstitutiveLawVector[PointNumber]->SetValue( rVariable, rValues[PointNumber], rCurrentProcessInfo );
    }
}

/******************************** SET VECTOR VALUE *********************************/
/***********************************************************************************/

void SprismElement3D6N::SetValueOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        )
{
    for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
    {
        mConstitutiveLawVector[PointNumber]->SetValue( rVariable, rValues[PointNumber], rCurrentProcessInfo );
    }

}

/******************************** SET MATRIX VALUE *********************************/
/***********************************************************************************/

void SprismElement3D6N::SetValueOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        )
{
    for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
    {
        mConstitutiveLawVector[PointNumber]->SetValue( rVariable, rValues[PointNumber], rCurrentProcessInfo );
    }
}

/****************************** SET CONSTITUTIVE VALUE *****************************/
/***********************************************************************************/

void SprismElement3D6N::SetValueOnIntegrationPoints(
        const Variable<ConstitutiveLaw::Pointer>& rVariable,
        std::vector<ConstitutiveLaw::Pointer>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        )
{
    if(rVariable == CONSTITUTIVE_LAW)
    {
        if ( mConstitutiveLawVector.size() != rValues.size() )
        {
            mConstitutiveLawVector.resize(rValues.size(), false);
            if( mConstitutiveLawVector.size() != GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod ) )
            {
                KRATOS_THROW_ERROR( std::logic_error, "Constitutive law not has the correct size ", mConstitutiveLawVector.size() );
            }
        }

        for(unsigned int i = 0; i < rValues.size(); i++)
        {
            mConstitutiveLawVector[i] = rValues[i]->Clone();
        }
    }
    if(rVariable == CONSTITUTIVE_LAW_POINTER)
    {
        if ( mConstitutiveLawVector.size() != rValues.size() )
        {
            mConstitutiveLawVector.resize(rValues.size(), false);
            if( mConstitutiveLawVector.size() != GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod ) )
            {
                KRATOS_THROW_ERROR( std::logic_error, "Constitutive law not has the correct size ", mConstitutiveLawVector.size() );
            }
        }
        for(unsigned int i = 0; i < rValues.size(); i++)
        {
            mConstitutiveLawVector[i] = rValues[i];
        }
    }
}

/******************************** GET DOUBLE VALUE *********************************/
/***********************************************************************************/

void SprismElement3D6N::GetValueOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        )
{
    if ( rVariable == VON_MISES_STRESS || rVariable == NORM_ISOCHORIC_STRESS )
    {
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    }
    else
    {
        const unsigned int& integration_point_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );
        if ( rValues.size() != integration_point_number )
        {
            rValues.resize( integration_point_number, false );
        }
        for ( unsigned int ii = 0; ii < integration_point_number; ii++ )
        {
          rValues[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rValues[ii] );
        }
    }
}

/********************************** GET VECTOR VALUE *******************************/
/***********************************************************************************/

void SprismElement3D6N::GetValueOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        )
{
    const unsigned int& integration_point_number = mConstitutiveLawVector.size();

    if ( rValues.size() != integration_point_number )
    {
        rValues.resize( integration_point_number );
    }

    if ( rVariable == PK2_STRESS_TENSOR ||  rVariable == CAUCHY_STRESS_TENSOR )
    {
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    }
    else if ( rVariable == PK2_STRESS_VECTOR ||  rVariable == CAUCHY_STRESS_VECTOR )
    {
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    }
    else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR ||  rVariable == ALMANSI_STRAIN_TENSOR ||  rVariable == HENCKY_STRAIN_TENSOR)
    {
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    }
    else
    {
        for ( unsigned int PointNumber = 0;  PointNumber < integration_point_number; PointNumber++ )
        {
            rValues[PointNumber] = mConstitutiveLawVector[PointNumber]->GetValue( rVariable, rValues[PointNumber] );
        }
    }
}

/*********************************** GET MATRIX VALUE ******************************/
/***********************************************************************************/

void SprismElement3D6N::GetValueOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        )
{
    const unsigned int& integration_point_number = mConstitutiveLawVector.size();

    if ( rValues.size() != integration_point_number )
    {
        rValues.resize( integration_point_number );
    }

    if ( rVariable == PK2_STRESS_TENSOR ||  rVariable == CAUCHY_STRESS_TENSOR )
    {
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    }
    else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR ||  rVariable == ALMANSI_STRAIN_TENSOR ||  rVariable == HENCKY_STRAIN_TENSOR)
    {
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    }
    else
    {
        for ( unsigned int PointNumber = 0;  PointNumber < integration_point_number; PointNumber++ )
        {
            rValues[PointNumber] = mConstitutiveLawVector[PointNumber]->GetValue( rVariable, rValues[PointNumber] );
        }
    }
}

/******************************** GET CONSTITUTIVE VALUE ***************************/
/***********************************************************************************/

void SprismElement3D6N::GetValueOnIntegrationPoints(
        const Variable<ConstitutiveLaw::Pointer>& rVariable,
        std::vector<ConstitutiveLaw::Pointer>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        )
{
    if(rVariable == CONSTITUTIVE_LAW || rVariable == CONSTITUTIVE_LAW_POINTER)
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
        {
            rValues.resize(mConstitutiveLawVector.size(), false);
        }
        for(unsigned int i = 0; i < rValues.size(); i++)
        {
            rValues[i] = mConstitutiveLawVector[i];
        }
    }
}

//********************************* CHECK VALUES **********************************//
/***********************************************************************************/
/***********************************************************************************/

int  SprismElement3D6N::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    /* Check the neighbours have been calculated */
    // Neighbour elements
    WeakPointerVector< Element >& elem_neigb = this->GetValue(NEIGHBOUR_ELEMENTS);
    if (elem_neigb.size() == 0)
    {
        KRATOS_THROW_ERROR(std::logic_error, "The neighbour elements are not calculated", "");
    }

    // Neighbour nodes
    WeakPointerVector< Node < 3 > >& nodal_neigb = this->GetValue(NEIGHBOUR_NODES);
    if (nodal_neigb.size() == 0)
    {
        KRATOS_THROW_ERROR(std::logic_error, "The neighbour nodes are not calculated", "");
    }

    /* Verify compatibility with the constitutive law */
    ConstitutiveLaw::Features LawFeatures;
    this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetLawFeatures(LawFeatures);

    bool correct_strain_measure = false;
    for(unsigned int i = 0; i < LawFeatures.mStrainMeasures.size(); i++)
    {
        if(LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Deformation_Gradient)
        {
            correct_strain_measure = true;
        }
    }

    if(correct_strain_measure == false)
    {
        KRATOS_THROW_ERROR( std::logic_error, "Constitutive law is not compatible with the element type ", " Large Displacements " );
    }

    // Verify that nodal variables are correctly initialized
    if (DISPLACEMENT.Key() == 0)
    {
        KRATOS_THROW_ERROR( std::invalid_argument, "DISPLACEMENT has Key zero! (check if the application is correctly registered", "" );
    }

    if (VELOCITY.Key() == 0)
    {
        KRATOS_THROW_ERROR( std::invalid_argument, "VELOCITY has Key zero! (check if the application is correctly registered", "" );
    }

    if (ACCELERATION.Key() == 0)
    {
        KRATOS_THROW_ERROR( std::invalid_argument, "ACCELERATION has Key zero! (check if the application is correctly registered", "" );
    }

    if (DENSITY.Key() == 0)
    {
        KRATOS_THROW_ERROR( std::invalid_argument, "DENSITY has Key zero! (check if the application is correctly registered", "" );
    }

    if (VOLUME_ACCELERATION.Key() == 0)
    {
        KRATOS_THROW_ERROR( std::invalid_argument,"VOLUME_ACCELERATION has Key zero! (check if the application is correctly registered", "" );
    }

    // Verify that elemental variables are correctly initialized
    if ( VON_MISES_STRESS.Key() == 0 )
    {
        KRATOS_THROW_ERROR( std::invalid_argument,"VON_MISES_STRESS has Key zero! (check if the application is correctly registered", "" );
    }

    if ( NORM_ISOCHORIC_STRESS.Key() == 0 )
    {
        KRATOS_THROW_ERROR( std::invalid_argument,"NORM_ISOCHORIC_STRESS has Key zero! (check if the application is correctly registered", "" );
    }

    if ( CAUCHY_STRESS_TENSOR.Key() == 0 )
    {
        KRATOS_THROW_ERROR( std::invalid_argument,"CAUCHY_STRESS_TENSOR has Key zero! (check if the application is correctly registered", "" );
    }

    if ( CAUCHY_STRESS_VECTOR.Key() == 0 )
    {
        KRATOS_THROW_ERROR( std::invalid_argument,"CAUCHY_STRESS_VECTOR has Key zero! (check if the application is correctly registered", "" );
    }

    if ( PK2_STRESS_TENSOR.Key() == 0 )
    {
        KRATOS_THROW_ERROR( std::invalid_argument,"PK2_STRESS_TENSOR has Key zero! (check if the application is correctly registered", "" );
    }

    if ( PK2_STRESS_VECTOR.Key() == 0 )
    {
        KRATOS_THROW_ERROR( std::invalid_argument,"PK2_STRESS_VECTOR has Key zero! (check if the application is correctly registered", "" );
    }

    if ( GREEN_LAGRANGE_STRAIN_TENSOR.Key() == 0 )
    {
        KRATOS_THROW_ERROR( std::invalid_argument,"GREEN_LAGRANGE_STRAIN_TENSOR has Key zero! (check if the application is correctly registered", "" );
    }

    if ( GREEN_LAGRANGE_STRAIN_VECTOR.Key() == 0 )
    {
        KRATOS_THROW_ERROR( std::invalid_argument,"GREEN_LAGRANGE_STRAIN_VECTOR has Key zero! (check if the application is correctly registered", "" );
    }

    if ( ALMANSI_STRAIN_TENSOR.Key() == 0 )
    {
        KRATOS_THROW_ERROR( std::invalid_argument,"ALMANSI_STRAIN_TENSOR has Key zero! (check if the application is correctly registered", "" );
    }

    if ( ALMANSI_STRAIN_VECTOR.Key() == 0 )
    {
        KRATOS_THROW_ERROR( std::invalid_argument,"ALMANSI_STRAIN_VECTOR has Key zero! (check if the application is correctly registered", "" );
    }

    if ( HENCKY_STRAIN_TENSOR.Key() == 0 )
    {
        KRATOS_THROW_ERROR( std::invalid_argument,"HENCKY_STRAIN_TENSOR has Key zero! (check if the application is correctly registered", "" );
    }

    if ( HENCKY_STRAIN_VECTOR.Key() == 0 )
    {
        KRATOS_THROW_ERROR( std::invalid_argument,"HENCKY_STRAIN_VECTOR has Key zero! (check if the application is correctly registered", "" );
    }

    if ( CONSTITUTIVE_MATRIX.Key() == 0 )
    {
        KRATOS_THROW_ERROR( std::invalid_argument,"CONSTITUTIVE_MATRIX has Key zero! (check if the application is correctly registered", "" );
    }

    if ( DEFORMATION_GRADIENT.Key() == 0 )
    {
        KRATOS_THROW_ERROR( std::invalid_argument,"DEFORMATION_GRADIENT has Key zero! (check if the application is correctly registered", "" );
    }

    // Verify that the dofs exist
    for ( unsigned int i = 0; i < this->GetGeometry().size(); i++)
    {
        if ( this->GetGeometry()[i].SolutionStepsDataHas( DISPLACEMENT ) == false)
        {
            KRATOS_THROW_ERROR( std::invalid_argument, "Missing variable DISPLACEMENT on node ", this->GetGeometry()[i].Id());
        }

        if ( this->GetGeometry()[i].HasDofFor( DISPLACEMENT_X ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Y ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Z ) == false)
        {
            KRATOS_THROW_ERROR( std::invalid_argument, "Missing one of the dofs for the variable DISPLACEMENT on node ", GetGeometry()[i].Id());
        }
    }

    // Verify that the constitutive law exists
    if (this->GetProperties().Has( CONSTITUTIVE_LAW ) == false)
    {
        KRATOS_THROW_ERROR( std::logic_error, "Constitutive law not provided for property ", this->GetProperties().Id());
    }

    if (this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() != 6)
    {
        KRATOS_THROW_ERROR( std::logic_error, "Wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) ", this->Id());
    }

    // Check constitutive law
    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
    {
        return mConstitutiveLawVector[i]->Check( GetProperties(), GetGeometry(), rCurrentProcessInfo);
    }

    return 0;

    KRATOS_CATCH( "" );
}

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::InitializeSolutionStep(
        ProcessInfo& rCurrentProcessInfo
        )
{
    ClearNodalForces();

    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        mConstitutiveLawVector[i]->InitializeSolutionStep( GetProperties(),
                GetGeometry(),
                row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ),
                rCurrentProcessInfo );

    mFinalizedStep = false;
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Create and initialize element variables:
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables);

    // StressMeasure_PK1             //stress related to reference configuration non-symmetric
    // StressMeasure_PK2             //stress related to reference configuration
    // StressMeasure_Kirchhoff       //stress related to current   configuration
    // StressMeasure_Cauchy          //stress related to current   configuration

    Variables.StressMeasure = ConstitutiveLaw::StressMeasure_PK2;

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    // Get constitutive law flags:
    Flags &ConstitutiveLawOptions=Values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    double& alpha_eas = Element::GetValue(ALPHA_EAS);

    /* Calculate common components (B, C) */
    this->CalculateCommonComponents();

    // Reading integration points
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        double zeta_gauss = 2.0 * integration_points[PointNumber].Z() - 1.0;

        // Compute element kinematics C, F ...
        this->CalculateKinematics(Variables,PointNumber, alpha_eas, zeta_gauss);
        this->CbartoFbar(Variables, PointNumber);

        // Set general variables to constitutivelaw parameters
        this->SetGeneralVariables(Variables,Values,PointNumber);

        // Call the constitutive law to update material variables
        mConstitutiveLawVector[PointNumber]->FinalizeMaterialResponse(Values, Variables.StressMeasure);

        // Call the constitutive law to finalize the solution step
        mConstitutiveLawVector[PointNumber]->FinalizeSolutionStep( GetProperties(),
                GetGeometry(),
                row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), PointNumber ),
                rCurrentProcessInfo );

        // Call the element internal variables update
        this->FinalizeStepVariables(Variables, PointNumber);
    }

    mFinalizedStep = true;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    ClearNodalForces();
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{

}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::SwitchFlagImplicitExplicit(bool& flag)
{
    KRATOS_TRY;

    bool& eas_imp = GetProperties()[EAS_IMP];
    eas_imp = flag;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::Initialize()
{
    KRATOS_TRY;

    boost::numeric::ublas::bounded_matrix<double, 12, 3 > nodes_coord = ZeroMatrix(12, 3);

    //**********************************************************************************
    /* Fill the aux matrix of coordinates */
    for (unsigned int i = 0; i < 6; i++)
    {
        const array_1d<double, 3> &CurrentPosition  = GetGeometry()[i].Coordinates();
        nodes_coord(i, 0) = CurrentPosition[0];
        nodes_coord(i, 1) = CurrentPosition[1];
        nodes_coord(i, 2) = CurrentPosition[2];
    }

    WeakPointerVector< Node < 3 > >& nodal_neigb = this->GetValue(NEIGHBOUR_NODES);
    int number_neigb = NumberOfActiveNeighbours(nodal_neigb);

    if (number_neigb == 6) // All the possible neighours
    {
        for (unsigned int i = 0; i < 6; i++)
        {
            const array_1d<double, 3> &CurrentPosition  = nodal_neigb[i].Coordinates();
            nodes_coord(i + 6, 0) = CurrentPosition[0];
            nodes_coord(i + 6, 1) = CurrentPosition[1];
            nodes_coord(i + 6, 2) = CurrentPosition[2];
        }
    }
    else
    {
        for (unsigned int i = 0; i < 6; i++)
        {
            if (HasNeighbour(i, nodal_neigb[i]))
            {
                const array_1d<double, 3> &CurrentPosition  = nodal_neigb[i].Coordinates();
                nodes_coord(i + 6, 0) = CurrentPosition[0];
                nodes_coord(i + 6, 1) = CurrentPosition[1];
                nodes_coord(i + 6, 2) = CurrentPosition[2];
            }
            else
            {
                nodes_coord(i + 6, 0) = 0.0;
                nodes_coord(i + 6, 1) = 0.0;
                nodes_coord(i + 6, 2) = 0.0;
            }
        }
    }

    //**********************************************************************************

    /* Calculate local system of coordinates of the element */
    // 0- If X is the prefered normal vector
    // 1- If Y is the prefered normal vector
    // 2- If Z is the prefered normal vector

    double ang_rot = 0.0;
    if( GetProperties().Has(ANG_ROT) )
    {
        ang_rot = GetProperties()[ANG_ROT];
    }

    this->CalculateLocalCoordinateSystem(2, ang_rot);

    //******************************** CENTRAL POINT ******************************
    // Calculate cartesian derivatives
    boost::numeric::ublas::bounded_matrix<double, 2, 4 > CartesianDerivativesCenterLower;
    boost::numeric::ublas::bounded_matrix<double, 2, 4 > CartesianDerivativesCenterUpper;

    // Lower face
    CalculateCartesianDerOnCenter_plane(0, nodes_coord, CartesianDerivativesCenterLower);
    // Upperr face
    CalculateCartesianDerOnCenter_plane(3, nodes_coord, CartesianDerivativesCenterUpper );

    /* Transversal derivative */
    CalculateCartesianDerOnCenter_trans(nodes_coord, 0); // Center
    CalculateCartesianDerOnCenter_trans(nodes_coord, 1); // Lower part
    CalculateCartesianDerOnCenter_trans(nodes_coord, 2); // Upper part

    //******************************** GAUSS POINTS *******************************

    /* Transversal derivative */
    CalculateCartesianDerOnGauss_trans(nodes_coord, mCC.mTransversalCartesianDerivativesGauss1, 0.5, 0.5, - 1.0);
    CalculateCartesianDerOnGauss_trans(nodes_coord, mCC.mTransversalCartesianDerivativesGauss4, 0.5, 0.5,   1.0);
    /* In-plane derivative */
    if (HasNeighbour(0, nodal_neigb[0])) // Assuming that if the upper element has neighbours the lower has too
    {
        CalculateCartesianDerOnGauss_plane(0, 0, nodes_coord, mCC.mInPlaneCartesianDerivativesGauss1);
        CalculateCartesianDerOnGauss_plane(0, 3, nodes_coord, mCC.mInPlaneCartesianDerivativesGauss4);
    }
    else
    {
        noalias(mCC.mInPlaneCartesianDerivativesGauss1) = CartesianDerivativesCenterLower;
        noalias(mCC.mInPlaneCartesianDerivativesGauss4) = CartesianDerivativesCenterUpper;
    }

    /* Transversal derivative */
    CalculateCartesianDerOnGauss_trans(nodes_coord, mCC.mTransversalCartesianDerivativesGauss2, 0.0, 0.5, - 1.0);
    CalculateCartesianDerOnGauss_trans(nodes_coord, mCC.mTransversalCartesianDerivativesGauss5, 0.0, 0.5,   1.0);
    /* In-plane derivative */
    if (HasNeighbour(1, nodal_neigb[1])) //Idem
    {
        CalculateCartesianDerOnGauss_plane(1, 0, nodes_coord, mCC.mInPlaneCartesianDerivativesGauss2);
        CalculateCartesianDerOnGauss_plane(1, 3, nodes_coord, mCC.mInPlaneCartesianDerivativesGauss5);
    }
    else
    {
        noalias(mCC.mInPlaneCartesianDerivativesGauss2) = CartesianDerivativesCenterLower;
        noalias(mCC.mInPlaneCartesianDerivativesGauss5) = CartesianDerivativesCenterUpper;
    }

    /* Transversal derivative */
    CalculateCartesianDerOnGauss_trans(nodes_coord, mCC.mTransversalCartesianDerivativesGauss3, 0.5, 0.0, - 1.0);
    CalculateCartesianDerOnGauss_trans(nodes_coord, mCC.mTransversalCartesianDerivativesGauss6, 0.5, 0.0,   1.0);

    /* In-plane derivative */
    if (HasNeighbour(2, nodal_neigb[2])) // Idem
    {
        CalculateCartesianDerOnGauss_plane(2, 0, nodes_coord, mCC.mInPlaneCartesianDerivativesGauss3);
        CalculateCartesianDerOnGauss_plane(2, 3, nodes_coord, mCC.mInPlaneCartesianDerivativesGauss6);
    }
    else
    {
        noalias(mCC.mInPlaneCartesianDerivativesGauss3) = CartesianDerivativesCenterLower;
        noalias(mCC.mInPlaneCartesianDerivativesGauss6) = CartesianDerivativesCenterUpper;
    }

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    /* Constitutive Law initialisation */
    if ( mConstitutiveLawVector.size() != integration_points.size() )
    {
        mConstitutiveLawVector.resize( integration_points.size(), false );
    }

    // Resizing jacobian inverses container
    mInvJ0.resize( integration_points.size() );
    mDetJ0.resize( integration_points.size(), false );

    // Compute jacobian inverses and set the domain initial size:
    GeometryType::JacobiansType J0;
    J0 = GetGeometry().Jacobian(J0, mThisIntegrationMethod);
    mTotalDomainInitialSize = 0.0;

    /* Calculating the inverse J0 */
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        // Calculating and storing inverse of the jacobian and the parameters needed
        MathUtils<double>::InvertMatrix( J0[PointNumber], mInvJ0[PointNumber], mDetJ0[PointNumber] );

        // Getting informations for integration
        double IntegrationWeight = integration_points[PointNumber].Weight();

        // Calculating the total volume
        mTotalDomainInitialSize += mDetJ0[PointNumber] * IntegrationWeight;
    }

    /* Initilize alpha_eas */
    double& alpha_eas = Element::GetValue(ALPHA_EAS);
    alpha_eas = 0.0;

    /* Reset the PK2 integrated components */
    mPK2.clear();

    /* Reset the EAS integrated components */
    mEAS.clear();

    /* Material initialisation */
    InitializeMaterial();

    /* Calculate ID vector*/
    this->CalculateIdVect();

    /* Calculate common components (B, C) */
    this->CalculateCommonComponents();

    /* Initialize previous coordinates */
    mPreviousCoor = GetVectorCurrentPosition();

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

// ------------------------------------------------------------------------- //
// ----------------------------- PROTECTED --------------------------------- //
// ------------------------------------------------------------------------- //

void SprismElement3D6N::CalculateElementalSystem(
        LocalSystemComponents& rLocalSystem,
        ProcessInfo& rCurrentProcessInfo
        )
{
    KRATOS_TRY;

    /* Create and initialize element variables: */
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables);

    // StressMeasure_PK1             //stress related to reference configuration non-symmetric
    // StressMeasure_PK2             //stress related to reference configuration
    // StressMeasure_Kirchhoff       //stress related to current   configuration
    // StressMeasure_Cauchy          //stress related to current   configuration

    Variables.StressMeasure = ConstitutiveLaw::StressMeasure_PK2;

    /* Create constitutive law parameters: */
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    /* Set constitutive law flags: */
    Flags &ConstitutiveLawOptions=Values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    /* Reading integration points */
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    /* Getting the alpha parameter of the EAS improvement */
    double& alpha_eas = Element::GetValue(ALPHA_EAS);

    /* Calculate the RHS */
    if ( rLocalSystem.CalculationFlags.Is(SprismElement3D6N::COMPUTE_RHS_VECTOR) == true ) // Update just if RHS is calculated
    {
        /* Getting the increase of displacements */
        boost::numeric::ublas::bounded_matrix<double, 36, 1 > delta_disp;

        delta_disp = GetVectorCurrentPosition() - mPreviousCoor; // Calculates the increase of displacements
        mPreviousCoor = GetVectorCurrentPosition(); // Update previous coordinates // Note: Save coordinates in an auxiliar variable

        /* Update alpha EAS */
        if (mEAS.mstiff_alpha > 1.0e-12) // Avoid division by zero
        {
            alpha_eas -= prod(mEAS.mH_EAS, delta_disp)(0, 0) / mEAS.mstiff_alpha;
        }
    }

    /* Calculate common components (B, C) */
    this->CalculateCommonComponents();

    /* Reset the PK2 integrated components */
    mPK2.clear();

    /* Reset the EAS integrated components */
    mEAS.clear();

    // Reading integration points
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        double zeta_gauss = 2.0 * integration_points[PointNumber].Z() - 1.0;

        /* Assemble B */
        this->CalculateDeformationMatrix(Variables.B, zeta_gauss, alpha_eas);

        // Compute element kinematics C, F ...
        this->CalculateKinematics(Variables, PointNumber, alpha_eas, zeta_gauss);
        this->CbartoFbar(Variables, PointNumber);

        // Set general variables to constitutivelaw parameters
        this->SetGeneralVariables(Variables, Values, PointNumber);

        // Compute stresses and constitutive parameters
        mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(Values, Variables.StressMeasure);

//        /* TEST */
//        this->CalculateGreenLagrangeStrain(Variables.C,Variables.StrainVector);
//        this->CalculateHyperelasticNeoHookeanStress(Variables);
////        this->CalculateHenckyStrain(Variables.C,Variables.StrainVector);
////        this->CalculateLogStress(Variables);
//        /* TEST */

        // Calculating weights for integration on the "reference configuration"
        double IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;

        /* Integrate in Zeta */
        IntegrateInZeta(Variables, alpha_eas, zeta_gauss, IntegrationWeight);
    }

    /* Auxiliary terms: Allocating the VolumeForce*/
    Vector VolumeForce;

    /* Calculate the RHS */
    if ( rLocalSystem.CalculationFlags.Is(SprismElement3D6N::COMPUTE_RHS_VECTOR) == true ) // Calculation of the vector is required
    {
        /* Volume forces */
        VolumeForce  = this->CalculateVolumeForce( VolumeForce, Variables );

        /* Contribution to external and internal forces */
        this->CalculateAndAddRHS ( rLocalSystem, Variables, VolumeForce, alpha_eas );
    }

    if ( rLocalSystem.CalculationFlags.Is(SprismElement3D6N::COMPUTE_LHS_MATRIX) == true ) // Calculation of the matrix is required
    {
        /* Contribution to the tangen stiffness matrix */
        this->CalculateAndAddLHS ( rLocalSystem, Variables, Values, alpha_eas );
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateDynamicSystem(
        LocalSystemComponents& rLocalSystem,
        ProcessInfo& rCurrentProcessInfo
        )
{
    KRATOS_TRY;

    // Create and initialize element variables:
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables);

    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    MatrixType  LocalLeftHandSideMatrix;
    VectorType  LocalRightHandSideVector;

    // Initialize sizes for the system components:
    this->InitializeSystemMatrices( LocalLeftHandSideMatrix, LocalRightHandSideVector, rLocalSystem.CalculationFlags );

    /* Getting the alpha parameter of the EAS improvement */
    double& alpha_eas = Element::GetValue(ALPHA_EAS);

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        double zeta_gauss = 2.0 * integration_points[PointNumber].Z() - 1.0;

        /* Assemble B */
        this->CalculateDeformationMatrix(Variables.B, zeta_gauss, alpha_eas);

        // Compute element kinematics C, F ...
        this->CalculateKinematics(Variables, PointNumber, alpha_eas, zeta_gauss);
        this->CbartoFbar(Variables, PointNumber);

        // Calculating weights for integration on the "reference configuration"
        double IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;

        if ( rLocalSystem.CalculationFlags.Is(SprismElement3D6N::COMPUTE_LHS_MATRIX) ) // Calculation of the matrix is required
        {
            LocalLeftHandSideMatrix.clear();

            this->CalculateAndAddDynamicLHS ( LocalLeftHandSideMatrix );

            MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();
            rLeftHandSideMatrix += LocalLeftHandSideMatrix;
        }
        if ( rLocalSystem.CalculationFlags.Is(SprismElement3D6N::COMPUTE_RHS_VECTOR) ) // Calculation of the vector is required
        {
            LocalRightHandSideVector.clear();

            this->CalculateAndAddDynamicRHS ( LocalRightHandSideVector, Variables, rCurrentProcessInfo, IntegrationWeight );

            VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();
            rRightHandSideVector += LocalRightHandSideVector;
        }

    // For debugging purposes
    // this->PrintElementCalculation(Variables);

    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::PrintElementCalculation(
        LocalSystemComponents& rLocalSystem,
        GeneralVariables& rVariables
        )
{
    KRATOS_TRY;

    std::cout << " Element: " << this->Id() << std::endl;

    WeakPointerVector< Node < 3 > >& nodal_neigb = this->GetValue(NEIGHBOUR_NODES);
    unsigned int number_neigb = NumberOfActiveNeighbours(nodal_neigb);

    for ( unsigned int i = 0; i < 6; i++ )
    {
        array_1d<double, 3> &CurrentPosition  = GetGeometry()[i].Coordinates();
        array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
        array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
        array_1d<double, 3> PreviousPosition  = CurrentPosition - (CurrentDisplacement-PreviousDisplacement);
        std::cout << " Previous  Position  node[" << GetGeometry()[i].Id() << "]: "<<PreviousPosition << std::endl;
    }

    for ( unsigned int i = 0; i < number_neigb; i++ )
    {
        array_1d<double, 3> &CurrentPosition  = nodal_neigb[i].Coordinates();
        array_1d<double, 3 > & CurrentDisplacement  = nodal_neigb[i].FastGetSolutionStepValue(DISPLACEMENT);
        array_1d<double, 3 > & PreviousDisplacement = nodal_neigb[i].FastGetSolutionStepValue(DISPLACEMENT,1);
        array_1d<double, 3> PreviousPosition  = CurrentPosition - (CurrentDisplacement-PreviousDisplacement);
        std::cout << " Previous  Position  neighbour node[" << nodal_neigb[i].Id() << "]: "<<PreviousPosition << std::endl;
    }

    for ( unsigned int i = 0; i < 6; i++ )
    {
        array_1d<double, 3> & CurrentPosition  = GetGeometry()[i].Coordinates();
        std::cout << " Current  Position  node[" << GetGeometry()[i].Id()<<"]: " << CurrentPosition << std::endl;
    }

    for ( unsigned int i = 0; i < number_neigb; i++ )
    {
        array_1d<double, 3> & CurrentPosition  = nodal_neigb[i].Coordinates();
        std::cout << " Current  Position neighbour node[" << nodal_neigb[i].Id()<<"]: " << CurrentPosition << std::endl;
    }

    for ( unsigned int i = 0; i < 6; i++ )
    {
        array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
        std::cout << " Previous Displacement node[" << GetGeometry()[i].Id() << "]: " << PreviousDisplacement << std::endl;
    }

    for ( unsigned int i = 0; i < number_neigb; i++ )
    {
        array_1d<double, 3 > & PreviousDisplacement = nodal_neigb[i].FastGetSolutionStepValue(DISPLACEMENT,1);
        std::cout << " Previous Displacement neighbour node[" << nodal_neigb[i].Id() << "]: " << PreviousDisplacement << std::endl;
    }

    for ( unsigned int i = 0; i < 6; i++ )
    {
        array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
        std::cout << " Current  Displacement  node[" << GetGeometry()[i].Id() << "]: " << CurrentDisplacement << std::endl;
    }

    for ( unsigned int i = 0; i < number_neigb; i++ )
    {
        array_1d<double, 3 > & CurrentDisplacement  = nodal_neigb[i].FastGetSolutionStepValue(DISPLACEMENT);
        std::cout << " Current  Displacement  node[" << nodal_neigb[i].Id() << "]: " << CurrentDisplacement << std::endl;
    }

    std::cout << " Stress " << rVariables.StressVector << std::endl;
    std::cout << " Strain " << rVariables.StrainVector << std::endl;
    std::cout << " F  " << rVariables.F<<std::endl;
    std::cout << " ConstitutiveMatrix " <<rVariables.ConstitutiveMatrix << std::endl;
    std::cout << " K " << rLocalSystem.GetLeftHandSideMatrix() << std::endl;
    std::cout << " f " << rLocalSystem.GetRightHandSideVector() << std::endl;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

bool SprismElement3D6N::HasNeighbour(unsigned int index, const Node < 3 > & neighb)
{
    bool quad_on = true;
    if( GetProperties().Has(QUAD_ON) )
    {
        quad_on = GetProperties()[QUAD_ON];
    }

    if (neighb.Id() == GetGeometry()[index].Id())
    {
        return false;
    }
    else
    {
        if (quad_on == true)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

unsigned int SprismElement3D6N::NumberOfActiveNeighbours(WeakPointerVector< Node < 3 > >& neighbs)
{
    unsigned int active_neighbours = 0;
    for (unsigned int i = 0; i < neighbs.size(); i++)
    {
        if (HasNeighbour(i, neighbs[i]))
        {
            active_neighbours++;
        }
    }
    return active_neighbours;
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateCommonComponents()
{
    KRATOS_TRY;

    WeakPointerVector< Node < 3 > >& nodal_neigb = this->GetValue(NEIGHBOUR_NODES);
    boost::numeric::ublas::bounded_matrix<double, 12, 3 > nodes_coord = ZeroMatrix(12, 3); // Coordinates of the nodes

    /* Fill the aux matrix of coordinates */
    for (unsigned int i = 0; i < 6; i++)
    {
        const array_1d<double, 3> &CurrentPosition  = GetGeometry()[i].Coordinates();
        nodes_coord(i, 0) = CurrentPosition[0];
        nodes_coord(i, 1) = CurrentPosition[1];
        nodes_coord(i, 2) = CurrentPosition[2];
    }

    int number_neigb = NumberOfActiveNeighbours(nodal_neigb);

    if (number_neigb == 6) // All the possible neighours
    {
        for (unsigned int i = 0; i < 6; i++)
        {
            const array_1d<double, 3> &CurrentPosition  = nodal_neigb[i].Coordinates();
            nodes_coord(i + 6, 0) = CurrentPosition[0];
            nodes_coord(i + 6, 1) = CurrentPosition[1];
            nodes_coord(i + 6, 2) = CurrentPosition[2];
        }
    }
    else
    {
        for (unsigned int i = 0; i < 6; i++)
        {
            if (HasNeighbour(i, nodal_neigb[i]))
            {
                const array_1d<double, 3> &CurrentPosition  = nodal_neigb[i].Coordinates();
                nodes_coord(i + 6, 0) = CurrentPosition[0];
                nodes_coord(i + 6, 1) = CurrentPosition[1];
                nodes_coord(i + 6, 2) = CurrentPosition[2];
            }
            else
            {
                nodes_coord(i + 6, 0) = 0.0;
                nodes_coord(i + 6, 1) = 0.0;
                nodes_coord(i + 6, 2) = 0.0;
            }
        }
    }

    /* Declare deformation Gradient F components */
    // In plane components
    boost::numeric::ublas::bounded_matrix<double, 3, 2 > InPlaneGradientFGauss;
    // Transversal components
    array_1d<double, 3 > TransverseGradientF0, TransverseGradientF1, TransverseGradientF2;

    //*****************************************************************************

    /* COMPUTATION OF B TANGENTS MATRIX */

    /* MEMBRANE CONTRIBUTION */
    /* Calculating the membrane strain-displacement matrix */
    // Lower face
    mCC.mB_membrane_lower = ZeroMatrix(3, 18);
    mCC.mC_membrane_lower = ZeroMatrix(3, 1);

    // Gauss point 1
    CalculateInPlaneGradientFGauss(InPlaneGradientFGauss, mCC.mInPlaneCartesianDerivativesGauss1, nodes_coord, 0, 0);
    CalculateAndAdd_B_Membrane(mCC.mB_membrane_lower, mCC.mC_membrane_lower, mCC.mInPlaneCartesianDerivativesGauss1, InPlaneGradientFGauss, 0);

    // Gauss point 2
    CalculateInPlaneGradientFGauss(InPlaneGradientFGauss, mCC.mInPlaneCartesianDerivativesGauss2, nodes_coord, 1, 0);
    CalculateAndAdd_B_Membrane(mCC.mB_membrane_lower, mCC.mC_membrane_lower, mCC.mInPlaneCartesianDerivativesGauss2, InPlaneGradientFGauss, 1);

    // Gauss point 3
    CalculateInPlaneGradientFGauss(InPlaneGradientFGauss, mCC.mInPlaneCartesianDerivativesGauss3, nodes_coord, 2, 0);
    CalculateAndAdd_B_Membrane(mCC.mB_membrane_lower, mCC.mC_membrane_lower, mCC.mInPlaneCartesianDerivativesGauss3, InPlaneGradientFGauss, 2);

    mCC.mB_membrane_lower *= 0.333333333333333333333333333333333;
    mCC.mC_membrane_lower *= 0.333333333333333333333333333333333;

    // Upper face
    mCC.mB_membrane_upper = ZeroMatrix(3, 18);
    mCC.mC_membrane_upper = ZeroMatrix(3, 1);

    // Gauss point 4
    CalculateInPlaneGradientFGauss(InPlaneGradientFGauss, mCC.mInPlaneCartesianDerivativesGauss4, nodes_coord, 0, 3);
    CalculateAndAdd_B_Membrane(mCC.mB_membrane_upper, mCC.mC_membrane_upper, mCC.mInPlaneCartesianDerivativesGauss4, InPlaneGradientFGauss, 0);

    // Gauss point 5
    CalculateInPlaneGradientFGauss(InPlaneGradientFGauss, mCC.mInPlaneCartesianDerivativesGauss5, nodes_coord, 1, 3);
    CalculateAndAdd_B_Membrane(mCC.mB_membrane_upper, mCC.mC_membrane_upper, mCC.mInPlaneCartesianDerivativesGauss5, InPlaneGradientFGauss, 1);

    // Gauss point 6
    CalculateInPlaneGradientFGauss(InPlaneGradientFGauss, mCC.mInPlaneCartesianDerivativesGauss6, nodes_coord, 2, 3);
    CalculateAndAdd_B_Membrane(mCC.mB_membrane_upper, mCC.mC_membrane_upper, mCC.mInPlaneCartesianDerivativesGauss6, InPlaneGradientFGauss, 2);

    mCC.mB_membrane_upper *= 0.333333333333333333333333333333333;
    mCC.mC_membrane_upper *= 0.333333333333333333333333333333333;

    /* SHEAR CONTRIBUTION */
    /* Calculating the shear strain-displacement matrix */

    // Declaring variables
    array_1d<double, 3 > TransverseGradientFt, TransverseGradientFxi, TransverseGradientFeta;

    // Lower face

    /* Calculate f components in the face */
    CalculateTransverseGradientFinP(TransverseGradientFt, TransverseGradientFxi, TransverseGradientFeta, nodes_coord, 0);

    /* Calculate f transverse components */
    CalculateTransverseGradientF(TransverseGradientF0, mCC.mTransversalCartesianDerivativesGauss1, nodes_coord);
    CalculateTransverseGradientF(TransverseGradientF1, mCC.mTransversalCartesianDerivativesGauss2, nodes_coord);
    CalculateTransverseGradientF(TransverseGradientF2, mCC.mTransversalCartesianDerivativesGauss3, nodes_coord);

    /* Shear contribution to the deformation matrix */
    CalculateAndAdd_B_Shear(mCC.mB_shear_lower, mCC.mC_shear_lower, mCC.mTransversalCartesianDerivativesGauss1, mCC.mTransversalCartesianDerivativesGauss2,\
                            mCC.mTransversalCartesianDerivativesGauss3, TransverseGradientF0, TransverseGradientF1, TransverseGradientF2,\
                            TransverseGradientFt, TransverseGradientFxi, TransverseGradientFeta, mCC.mJinv_plane_lower, 0);

    // Upper face

    /* Calculate f components in the face */
    CalculateTransverseGradientFinP(TransverseGradientFt, TransverseGradientFxi, TransverseGradientFeta, nodes_coord, 3);

    /* Calculate f transverse components */
    CalculateTransverseGradientF(TransverseGradientF0, mCC.mTransversalCartesianDerivativesGauss4, nodes_coord);
    CalculateTransverseGradientF(TransverseGradientF1, mCC.mTransversalCartesianDerivativesGauss5, nodes_coord);
    CalculateTransverseGradientF(TransverseGradientF2, mCC.mTransversalCartesianDerivativesGauss6, nodes_coord);

    /* Shear contribution to the deformation matrix */
    CalculateAndAdd_B_Shear(mCC.mB_shear_upper, mCC.mC_shear_upper, mCC.mTransversalCartesianDerivativesGauss4, mCC.mTransversalCartesianDerivativesGauss5,\
                            mCC.mTransversalCartesianDerivativesGauss6, TransverseGradientF0, TransverseGradientF1, TransverseGradientF2,\
                            TransverseGradientFt, TransverseGradientFxi, TransverseGradientFeta, mCC.mJinv_plane_upper, 9);

    /* NORMAL TRANSVERSE */
    /* Calculate f normal components */
    CalculateTransverseGradientF(TransverseGradientF0, mCC.mTransversalCartesianDerivativesCenter, nodes_coord);

    /* Calculating the normal transverse strain-displacement matrix */
    CalculateAndAdd_B_Normal(mCC.mB_normal, mCC.mC_normal, mCC.mTransversalCartesianDerivativesCenter, TransverseGradientF0);

//    if (this->Id() == 1)
//    {
//        std::cout << "\n "<< std::endl;
////        KRATOS_WATCH(mCC.mB_membrane_lower);
////        KRATOS_WATCH(mCC.mB_membrane_upper);
////        KRATOS_WATCH(mCC.mB_shear_lower);
////        KRATOS_WATCH(mCC.mB_shear_upper);
////        KRATOS_WATCH(mCC.mB_normal);
//        KRATOS_WATCH(mCC.mC_membrane_lower);
//        KRATOS_WATCH(mCC.mC_membrane_upper);
//        KRATOS_WATCH(mCC.mC_shear_lower);
//        KRATOS_WATCH(mCC.mC_shear_upper);
//        KRATOS_WATCH(mCC.mC_normal);
//    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateLocalCoordinateSystem(
        const int& choose,
        const double& ang
        )
{
    KRATOS_TRY;

    /* Mid-surface vectors */
    double norm;
    array_1d<double, 3 > vxe, vye;
    vxe[0] = 0.5 * ((GetGeometry()[2].X0() + GetGeometry()[5].X0()) - (GetGeometry()[1].X0() + GetGeometry()[4].X0()));
    vxe[1] = 0.5 * ((GetGeometry()[2].Y0() + GetGeometry()[5].Y0()) - (GetGeometry()[1].Y0() + GetGeometry()[4].Y0()));
    vxe[2] = 0.5 * ((GetGeometry()[2].Z0() + GetGeometry()[5].Z0()) - (GetGeometry()[1].Z0() + GetGeometry()[4].Z0()));

    vye[0] = 0.5 * ((GetGeometry()[0].X0() + GetGeometry()[3].X0()) - (GetGeometry()[2].X0() + GetGeometry()[5].X0()));
    vye[1] = 0.5 * ((GetGeometry()[0].Y0() + GetGeometry()[3].Y0()) - (GetGeometry()[2].Y0() + GetGeometry()[5].Y0()));
    vye[2] = 0.5 * ((GetGeometry()[0].Z0() + GetGeometry()[3].Z0()) - (GetGeometry()[2].Z0() + GetGeometry()[5].Z0()));

    MathUtils<double>::CrossProduct(mvze, vxe, vye);
    norm = norm_2(mvze);
    mvze /= norm;

    double threshold = 1e-5;
    double ortho_comp;

    /* Performing the calculation */
    // 0- If X is the prefered normal vector
    // 1- If Y is the prefered normal vector
    // 2- If Z is the prefered normal vector

    if (choose == 0)
    {
        ortho_comp = mvze[1] * mvze[1] + mvze[2] * mvze[2]; // Component in th Y-Z plane
        if (ortho_comp < threshold) // If mvze is almost orthogonal to  Y-Z plane
        {
            mvye[0] = - mvze[2]; // Choose mvxe orthogonal to global Y direction
            mvye[1] =       0.0;
            mvye[2] =   mvze[0];

            norm = norm_2(mvxe);
            mvxe /= norm;
            MathUtils<double>::CrossProduct(mvxe, mvye, mvze);
        }
        else // SELECT local y=mvxe in the global YZ plane
        {
            mvxe[0] =       0.0;
            mvxe[1] =   mvze[2];
            mvxe[2] = - mvze[1];

            norm = norm_2(mvxe);
            mvxe /= norm;

            mvye[0] =          ortho_comp; // Choose mvxe orthogonal to global X direction
            mvye[1] = - mvze[0] * mvze[1];
            mvye[2] = - mvze[0] * mvze[2];

            norm = norm_2(mvye);
            mvye /= norm;
        }
    }
    else if (choose == 1)
    {
        ortho_comp = mvze[0] * mvze[0] + mvze[2] * mvze[2]; // Component in th Z-X plane
        if (ortho_comp < threshold) // If vze is almost orthogonal to  Z-X plane
        {
            mvye[0] =       0.0; // Choose mvxe orthogonal to global X direction
            mvye[1] =   mvze[2];
            mvye[2] = - mvze[1];

            norm = norm_2(mvye);
            mvye /= norm;
            MathUtils<double>::CrossProduct(mvxe, mvye, mvze);
        }
        else // SELECT local z=mvxe in the global ZX plane
        {
            mvxe[0] = - mvze[2]; // Choose mvxe orthogonal to global Y direction
            mvxe[1] =       0.0;
            mvxe[2] = - mvze[0];

            norm = norm_2(mvxe);
            mvxe /= norm;

            mvye[0] = - mvze[0] * mvze[1];
            mvye[1] =          ortho_comp;
            mvye[2] = - mvze[2] * mvze[1];

            norm = norm_2(mvye);
            mvye /= norm;
        }
    }
    else if (choose == 2)
    {
        ortho_comp = mvze[0] * mvze[0] + mvze[1] * mvze[1]; // Component in th X-Y plane
        if (ortho_comp < threshold) // If vze is almost orthogonal to  X-Y plane
        {
            mvye[0] =       0.0; // Choose mvxe orthogonal to global X direction
            mvye[1] =   mvze[2];
            mvye[2] = - mvze[1];

            norm = norm_2(mvye);
            mvye /= norm;
            MathUtils<double>::CrossProduct(mvxe, mvye, mvze);
        }
        else // SELECT local x=mvxe in the global XY plane
        {
            mvxe[0] = - mvze[1];
            mvxe[1] =   mvze[0];
            mvxe[2] =       0.0;

            norm = norm_2(mvxe);
            mvxe /= norm;

            mvye[0] = - mvze[0] * mvze[2]; // Choose mvxe orthogonal to global Z direction
            mvye[1] = - mvze[1] * mvze[2];
            mvye[2] =          ortho_comp;

            norm = norm_2(mvye);
            mvye /= norm;
        }
    }
    else
    {
        mvxe[0] = 1.0;
        mvxe[1] = 0.0;
        mvxe[2] = 0.0;

        mvye[0] = 0.0;
        mvye[1] = 1.0;
        mvye[2] = 0.0;
    }

    if (ang != 0.0)
    {
        // Compute angle between local system mvxe-mvye and L1
        double cosa = cos(ang);
        double sina = sin(ang);
        // Rotate local system mvxe-mvye to best fit L1-L2
        mvze = mvxe; // Reusing as auxiliar value
        mvxe =   cosa * mvxe + sina * mvye;
        mvye = - sina * mvze  + cosa * mvye;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateIdVect()
{
    KRATOS_TRY;

    WeakPointerVector< Node < 3 > >& nodal_neigb = this->GetValue(NEIGHBOUR_NODES);

    /* Compute mid_vec */
    for (unsigned int i = 0; i < 18; i++)
    {
        mid_vec[i] = i;
    }
    unsigned int index = 18;
    for (unsigned int i = 0; i < 6; i++)
    {
        if (HasNeighbour(i, nodal_neigb[i]))
        {
            for (unsigned int j = 0; j < 3; j++)
            {
                mid_vec[18 + i * 3 + j] = index + j;
            }
            index += 3;
        }
        else
        {
            for (unsigned int j = 0; j < 3; j++)
            {
                mid_vec[18 + i * 3 + j] = 1000;
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::ComputeLocalDerivatives(
        boost::numeric::ublas::bounded_matrix<double, 6, 3 > & local_der_patch,
        const double xi,
        const double eta,
        const double zeta
        )
{
    double L_1 = 0.5 * (1.0 - zeta);
    double L_2 = 0.5 * (1.0 + zeta);

    /* Derivative in direction nu and xi */
    // Lower face
    local_der_patch(0, 0) = - L_1;
    local_der_patch(1, 0) =   L_1;
    local_der_patch(2, 0) =   0.0;

    local_der_patch(0, 1) = - L_1;
    local_der_patch(1, 1) =   0.0;
    local_der_patch(2, 1) =   L_1;

    // Upper face
    local_der_patch(3, 0) = - L_2;
    local_der_patch(4, 0) =   L_2;
    local_der_patch(5, 0) =   0.0;

    local_der_patch(3, 1) = - L_2;
    local_der_patch(4, 1) =   0.0;
    local_der_patch(5, 1) =   L_2;

    /* Derivative in direction zeta */
//    local_der_patch(0, 2) = - 0.5 * (1.0 - eta - xi);
//    local_der_patch(1, 2) = - 0.5 * xi;
//    local_der_patch(2, 2) = - 0.5 * eta;
//    local_der_patch(3, 2) =   0.5 * (1.0 - eta - xi);
//    local_der_patch(4, 2) =   0.5 * xi;
//    local_der_patch(5, 2) =   0.5 * eta;

    local_der_patch(0, 2) = - 1.0 + eta + xi;
    local_der_patch(1, 2) = - xi;
    local_der_patch(2, 2) = - eta;
    local_der_patch(3, 2) =   1.0 - eta - xi;
    local_der_patch(4, 2) =   xi;
    local_der_patch(5, 2) =   eta;
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::ComputeLocalDerivativesQuadratic(
        boost::numeric::ublas::bounded_matrix<double, 4, 2 > & local_der_patch,
        const int node_gauss
        )
{
    /* Local coordinates */
    double xi  = 0.0;
    double eta = 0.0;

    if (node_gauss == 0)
    {
        xi  = 0.50000000000000000000;
        eta = 0.50000000000000000000;
    }
    else if (node_gauss == 1)
    {
        xi  = 0.00000000000000000000;
        eta = 0.50000000000000000000;
    }
    else if (node_gauss == 2)
    {
        xi  = 0.50000000000000000000;
        eta = 0.00000000000000000000;
    }

    local_der_patch = ZeroMatrix(4, 2);

    /* Derivative in main nodes */
    local_der_patch(0, 0) = - 1.0 + eta;
    local_der_patch(0, 1) = - 1.0 + xi;
    local_der_patch(1, 0) =   1.0 - eta;
    local_der_patch(1, 1) =   1.0 - xi - 2.0 * eta;
    local_der_patch(2, 0) =   1.0 - 2.0 * xi - eta;
    local_der_patch(2, 1) =   1.0 - xi;

    /* Derivative in neighbour nodes */
    if (node_gauss == 0)
    {
        local_der_patch(3, 0) = xi + eta - 0.5;
        local_der_patch(3, 1) = xi + eta - 0.5;
    }
    else if (node_gauss == 1)
    {
        local_der_patch(3, 0) = xi - 0.5;
        local_der_patch(3, 1) = 0.0;
    }
    else if (node_gauss == 2)
    {
        local_der_patch(3, 0) = 0.0;
        local_der_patch(3, 1) = eta - 0.5;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateJacobianCenterGauss(
        GeometryType::JacobiansType& J,
        std::vector< Matrix >& Jinv,
        Vector& detJ,
        const int& rPointNumber,
        const double& zeta
        )
{
    /* Fill the aux matrix of coordinates */
    boost::numeric::ublas::bounded_matrix<double, 3, 6 > nodes_coord;
    for (unsigned int i = 0; i < 6; i++)
    {
        const array_1d<double, 3> &CurrentPosition  = GetGeometry()[i].Coordinates();
        nodes_coord(0, i) = CurrentPosition[0];
        nodes_coord(1, i) = CurrentPosition[1];
        nodes_coord(2, i) = CurrentPosition[2];
    }

    double xi  = 0.333333333333333333333333333333333;
    double eta = 0.333333333333333333333333333333333;

    /* Local derivatives patch */
    boost::numeric::ublas::bounded_matrix<double, 6, 3 > local_der_patch;
    ComputeLocalDerivatives(local_der_patch, xi, eta, zeta);

    /* Compute Jacobian */
    noalias(J[rPointNumber]) = prod(nodes_coord, local_der_patch);

    /* Compute inverse of the Jaccobian */
    MathUtils<double>::InvertMatrix( J[rPointNumber], Jinv[rPointNumber], detJ[rPointNumber] );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateJacobian(
        double & detJ,
        boost::numeric::ublas::bounded_matrix<double, 3, 3 > & J,
        boost::numeric::ublas::bounded_matrix<double, 6, 3 > & local_der_patch,
        const boost::numeric::ublas::bounded_matrix<double, 12, 3 > & nodes_coord,
        const double xi,
        const double eta,
        const double zeta
        )
{
    /* Auxiliar coordinates of the nodes */
    boost::numeric::ublas::bounded_matrix<double, 3, 6 > nodes_coord_aux;

    for (unsigned int i = 0; i < 6; i++)
    {
        nodes_coord_aux(0, i) = nodes_coord(i, 0);
        nodes_coord_aux(1, i) = nodes_coord(i, 1);
        nodes_coord_aux(2, i) = nodes_coord(i, 2);
    }

    /* Local derivatives patch */
    ComputeLocalDerivatives(local_der_patch, xi, eta, zeta);

    /* Compute Jacobian */
    noalias(J) = prod(nodes_coord_aux, local_der_patch);

    /* Compute determinant */
    StructuralMechanicsMathUtilities::DetMat3x3(J, detJ);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateJacobianAndInv(
        boost::numeric::ublas::bounded_matrix<double, 3, 3 > & J,
        boost::numeric::ublas::bounded_matrix<double, 3, 3 > & Jinv,
        boost::numeric::ublas::bounded_matrix<double, 6, 3 > & local_der_patch,
        const boost::numeric::ublas::bounded_matrix<double, 3, 6 > & nodes_coord,
        const double xi,
        const double eta,
        const double zeta
        )
{
    /* Local derivatives patch */
    ComputeLocalDerivatives(local_der_patch, xi, eta, zeta);

    /* Compute Jacobian */
    noalias(J) = prod(nodes_coord, local_der_patch);

    /* Compute inverse of the Jaccobian */
    StructuralMechanicsMathUtilities::InvMat3x3(J, Jinv);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateJacobianAndInv(
        boost::numeric::ublas::bounded_matrix<double, 3, 3 > & J,
        boost::numeric::ublas::bounded_matrix<double, 3, 3 > & Jinv,
        const boost::numeric::ublas::bounded_matrix<double, 3, 6 > & nodes_coord,
        const double xi,
        const double eta,
        const double zeta
        )
{
    /* Local derivatives patch */
    boost::numeric::ublas::bounded_matrix<double, 6, 3 > local_der_patch;
    ComputeLocalDerivatives(local_der_patch, xi, eta, zeta);

    /* Compute Jacobian */
    noalias(J) = prod(nodes_coord, local_der_patch);

    /* Compute inverse of the Jaccobian */
    StructuralMechanicsMathUtilities::InvMat3x3(J, Jinv);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateCartesianDerOnCenter_plane(
        const int index,
        const boost::numeric::ublas::bounded_matrix<double, 12, 3 > & nodes_coord,
        boost::numeric::ublas::bounded_matrix<double, 2, 4 > & CartesianDerivativesCenter
        )
{
    double norm0, norm;
    array_1d<double, 3 > vxe, vye;
    vxe[0] = GetGeometry()[2 + index].X0() - GetGeometry()[1 + index].X0();
    vxe[1] = GetGeometry()[2 + index].Y0() - GetGeometry()[1 + index].Y0();
    vxe[2] = GetGeometry()[2 + index].Z0() - GetGeometry()[1 + index].Z0();

    vye[0] = GetGeometry()[0 + index].X0() - GetGeometry()[2 + index].X0();
    vye[1] = GetGeometry()[0 + index].Y0() - GetGeometry()[2 + index].Y0();
    vye[2] = GetGeometry()[0 + index].Z0() - GetGeometry()[2 + index].Z0();

    array_1d<double, 3 > t1g, t2g, t3g;
    MathUtils<double>::CrossProduct(t3g, vxe, vye);
    norm0 = norm_2(t3g);
    t3g /= norm0;

    MathUtils<double>::CrossProduct(t2g, t3g, mvxe);
    norm = norm_2(t2g);
    t2g /= norm;

    MathUtils<double>::CrossProduct(t1g, t2g, t3g);
    norm = norm_2(t1g);
    t1g /= norm;

    array_1d<double, 3 > a, b;

    a[0] = inner_prod(vxe, t1g)/norm0;
    a[1] = inner_prod(vye, t1g)/norm0;
    a[2] = -(a[0] + a[1]);
    b[0] = inner_prod(vxe, t2g)/norm0;
    b[1] = inner_prod(vye, t2g)/norm0;
    b[2] = -(b[0] + b[1]);

    CartesianDerivativesCenter = ZeroMatrix(2, 4);
    for (unsigned int i = 0; i < 3; i++)
    {
       CartesianDerivativesCenter(0, i) = - b[i];
       CartesianDerivativesCenter(1, i) =   a[i];
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateCartesianDerOnGauss_plane(
        const int node_gauss,
        const int index,
        const boost::numeric::ublas::bounded_matrix<double, 12, 3 > & nodes_coord,
        boost::numeric::ublas::bounded_matrix<double, 2, 4 > & InPlaneCartesianDerivativesGauss
        )
{
    /* Local derivatives patch */
    boost::numeric::ublas::bounded_matrix<double, 4, 2 > local_der_patch;
    ComputeLocalDerivativesQuadratic(local_der_patch,node_gauss);

    /* Auxiliar coordinates of the nodes */
    boost::numeric::ublas::bounded_matrix<double, 3, 4 > nodes_coord_aux;

    for (unsigned int i = 0; i < 3; i++)
    {
        nodes_coord_aux(0, i) = nodes_coord(i + index, 0);
        nodes_coord_aux(1, i) = nodes_coord(i + index, 1);
        nodes_coord_aux(2, i) = nodes_coord(i + index, 2);
    }

    nodes_coord_aux(0, 3) = nodes_coord(node_gauss + 6 + index, 0);
    nodes_coord_aux(1, 3) = nodes_coord(node_gauss + 6 + index, 1);
    nodes_coord_aux(2, 3) = nodes_coord(node_gauss + 6 + index, 2);

    /* Compute local derivatives */
    boost::numeric::ublas::bounded_matrix<double, 3, 2 > Xd;
    noalias(Xd) = prod(nodes_coord_aux, local_der_patch);

    /* Split local derivatives */
    array_1d<double, 3 > Xdxi, Xdeta;
    Xdxi[0]  = Xd(0, 0);
    Xdxi[1]  = Xd(1, 0);
    Xdxi[2]  = Xd(2, 0);
    Xdeta[0] = Xd(0, 1);
    Xdeta[1] = Xd(1, 1);
    Xdeta[2] = Xd(2, 1);

    /* Compute orthonormal vectors */
    array_1d<double, 3 > t1g, t2g, t3g;
    StructuralMechanicsMathUtilities::Comp_Orthonor_Base(t1g, t2g, t3g, mvxe, Xdxi, Xdeta);

    /* Compute Jacobian */
    boost::numeric::ublas::bounded_matrix<double, 2, 2 > jac;
    jac(0, 0) = inner_prod(Xdxi,  t1g);
    jac(0, 1) = inner_prod(Xdxi,  t2g);
    jac(1, 0) = inner_prod(Xdeta, t1g);
    jac(1, 1) = inner_prod(Xdeta, t2g);

    /* Compute the inverse of the Jacobian */
    boost::numeric::ublas::bounded_matrix<double, 2, 2 > Jinv_plane;
    StructuralMechanicsMathUtilities::InvMat2x2(jac, Jinv_plane);

    /* Compute the Cartesian derivatives */
    noalias(InPlaneCartesianDerivativesGauss) = prod(Jinv_plane, trans(local_der_patch));
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateCartesianDerOnGauss_trans(
        const boost::numeric::ublas::bounded_matrix<double, 12, 3 > & nodes_coord,
        boost::numeric::ublas::bounded_matrix<double, 6, 1 > & TransversalCartesianDerivativesGauss,
        const double xi,
        const double eta,
        const double zeta
        )
{
    /* Compute local derivatives */
    double det;
    boost::numeric::ublas::bounded_matrix<double, 3, 3 > Xd;
    boost::numeric::ublas::bounded_matrix<double, 6, 3 > local_der_patch;
    CalculateJacobian(det, Xd, local_der_patch, nodes_coord, xi, eta, zeta);

    /* Split local derivatives */
    array_1d<double, 3 > Xdxi, Xdeta;
    Xdxi[0]  = Xd(0, 0);
    Xdxi[1]  = Xd(1, 0);
    Xdxi[2]  = Xd(2, 0);
    Xdeta[0] = Xd(0, 1);
    Xdeta[1] = Xd(1, 1);
    Xdeta[2] = Xd(2, 1);

    /* Compute orthonormal vectors */
    array_1d<double, 3 > t1g, t2g, t3g;
    boost::numeric::ublas::bounded_matrix<double, 3, 3 > t = ZeroMatrix(3, 3);
    StructuralMechanicsMathUtilities::Comp_Orthonor_Base(t, t1g, t2g, t3g, mvxe, Xdxi, Xdeta);

    /* Compute Jacobian */
    boost::numeric::ublas::bounded_matrix<double, 3, 3 > jac;
    noalias(jac) = prod(t, Xd);

    /* Compute inverse of the Jaccobian (just third column) */
    boost::numeric::ublas::bounded_matrix<double, 3, 1 > Jinv_trans;
    StructuralMechanicsMathUtilities::InvMat3x3Col(jac, 3, det,  Jinv_trans);

    /* Compute Cartesian derivatives */
    noalias(TransversalCartesianDerivativesGauss) = prod(local_der_patch, Jinv_trans);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateCartesianDerOnCenter_trans(
        const boost::numeric::ublas::bounded_matrix<double, 12, 3 > & nodes_coord,
        const int part
        )
{
    double xi  = 0.333333333333333333333333333333333;
    double eta = 0.333333333333333333333333333333333;
    double zeta;

    if (part == 0)
    {
        zeta =   0.0;
    }
    else if (part == 1)
    {
        zeta = - 1.0;
    }
    else if (part == 2)
    {
        zeta =   1.0;
    }
    else
    {
        KRATOS_THROW_ERROR( std::invalid_argument," This part id is not possible, just 0, 1 or 2  ", part );
    }

    /* Auxiliar coordinates of the nodes */
    boost::numeric::ublas::bounded_matrix<double, 3, 6 > nodes_coord_aux;
    for (unsigned int i = 0; i < 6; i++)
    {
        nodes_coord_aux(0, i) = nodes_coord(i, 0);
        nodes_coord_aux(1, i) = nodes_coord(i, 1);
        nodes_coord_aux(2, i) = nodes_coord(i, 2);
    }

    /* Auxiliar components to calculate the Jacobian and his inverse */
    boost::numeric::ublas::bounded_matrix<double, 3, 3 > J;
    boost::numeric::ublas::bounded_matrix<double, 3, 3 > Jinv;

    if (part == 0)
    {
        /* Calculate the Jacobian and his inverse */
        boost::numeric::ublas::bounded_matrix<double, 6, 3 > local_der_patch;
        CalculateJacobianAndInv(J, Jinv, local_der_patch, nodes_coord_aux, xi, eta, zeta);

        // Compute cartesian (y3) derivatives of the shape functions necessary to compute f_3
        /* Compute Cartesian derivatives */
        boost::numeric::ublas::bounded_matrix<double, 6, 3 > TransversalCartesianDerivativesGauss_aux;
        noalias(TransversalCartesianDerivativesGauss_aux) = prod(local_der_patch, Jinv);

        for (unsigned int i = 0; i < 6 ; i++)
        {
            mCC.mTransversalCartesianDerivativesCenter(i, 0) = \
                      mvze[0] * TransversalCartesianDerivativesGauss_aux(i, 0) + \
                      mvze[1] * TransversalCartesianDerivativesGauss_aux(i, 1) + \
                      mvze[2] * TransversalCartesianDerivativesGauss_aux(i, 2);
        }
     }
     else
     {
        /* Calculate the Jacobian and his inverse */
        CalculateJacobianAndInv(J, Jinv, nodes_coord_aux, xi, eta, zeta);

         /* Split local derivatives */
         array_1d<double, 3 > Xdxi, Xdeta;
         Xdxi[0]   = Jinv(0, 0);
         Xdxi[1]   = Jinv(0, 1);
         Xdxi[2]   = Jinv(0, 2);
         Xdeta[0]  = Jinv(1, 0);
         Xdeta[1]  = Jinv(1, 1);
         Xdeta[2]  = Jinv(1, 2);

         /* Compute inverse of the Jaccobian (just in plane components)*/
         if (part == 1)
         {
             mCC.mJinv_plane_lower(0, 0) = inner_prod(Xdxi,  mvxe);
             mCC.mJinv_plane_lower(0, 1) = inner_prod(Xdeta, mvxe);
             mCC.mJinv_plane_lower(1, 0) = inner_prod(Xdxi,  mvye);
             mCC.mJinv_plane_lower(1, 1) = inner_prod(Xdeta, mvye);
         }
         else if (part == 2)
         {
             mCC.mJinv_plane_upper(0, 0) = inner_prod(Xdxi,  mvxe);
             mCC.mJinv_plane_upper(0, 1) = inner_prod(Xdeta, mvxe);
             mCC.mJinv_plane_upper(1, 0) = inner_prod(Xdxi,  mvye);
             mCC.mJinv_plane_upper(1, 1) = inner_prod(Xdeta, mvye);
         }
     }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateInPlaneGradientFGauss(
        boost::numeric::ublas::bounded_matrix<double, 3, 2 > & InPlaneGradientFGauss,
        const boost::numeric::ublas::bounded_matrix<double, 2, 4 > & InPlaneCartesianDerivativesGauss,
        const boost::numeric::ublas::bounded_matrix<double, 12, 3 > & nodes_coord,
        const int node_gauss,
        const int index
        )
{
    /* Auxiliar operators */
    boost::numeric::ublas::bounded_matrix<double, 3, 3 > nodes_coord_aux;
    boost::numeric::ublas::bounded_matrix<double, 3, 2 > InPlaneCartesianDerivativesGauss_aux;

    for (unsigned int i = 0; i < 3; i++)
    {
        nodes_coord_aux(0, i) = nodes_coord(i + index, 0);
        nodes_coord_aux(1, i) = nodes_coord(i + index, 1);
        nodes_coord_aux(2, i) = nodes_coord(i + index, 2);

        InPlaneCartesianDerivativesGauss_aux(i, 0) = InPlaneCartesianDerivativesGauss(0, i);
        InPlaneCartesianDerivativesGauss_aux(i, 1) = InPlaneCartesianDerivativesGauss(1, i);
    }

    noalias(InPlaneGradientFGauss) = prod(nodes_coord_aux, InPlaneCartesianDerivativesGauss_aux);

    WeakPointerVector< Node < 3 > >& nodal_neigb = this->GetValue(NEIGHBOUR_NODES);
    if (HasNeighbour(node_gauss, nodal_neigb[node_gauss]))
    {
        for (unsigned int j = 0; j < 3 ; j++)
        {
            InPlaneGradientFGauss(j, 0) += nodes_coord(node_gauss + 6 + index, j) * InPlaneCartesianDerivativesGauss(0, 3);
            InPlaneGradientFGauss(j, 1) += nodes_coord(node_gauss + 6 + index, j) * InPlaneCartesianDerivativesGauss(1, 3);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateTransverseGradientF(
        array_1d<double, 3 > & TransverseGradientF,
        const boost::numeric::ublas::bounded_matrix<double, 1, 6 > & TransversalCartesianDerivativesGauss,
        const boost::numeric::ublas::bounded_matrix<double, 12, 3 > & nodes_coord
        )
{
    noalias(TransverseGradientF) = ZeroVector(3);

    for (unsigned int i = 0; i < 6; i++)
    {
        for (unsigned int j = 0; j < 3; j++)
        {
            TransverseGradientF[j] += TransversalCartesianDerivativesGauss(0, i) * nodes_coord(i, j);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateTransverseGradientFinP(
        array_1d<double, 3 > & TransverseGradientFt,
        array_1d<double, 3 > & TransverseGradientFxi,
        array_1d<double, 3 > & TransverseGradientFeta,
        const boost::numeric::ublas::bounded_matrix<double, 12, 3 > & nodes_coord,
        const int index
        )
{
    TransverseGradientFt[0]   = nodes_coord(2 + index, 0) - nodes_coord(1 + index, 0);
    TransverseGradientFt[1]   = nodes_coord(2 + index, 1) - nodes_coord(1 + index, 1);
    TransverseGradientFt[2]   = nodes_coord(2 + index, 2) - nodes_coord(1 + index, 2);

    TransverseGradientFxi[0]  = nodes_coord(0 + index, 0) - nodes_coord(2 + index, 0);
    TransverseGradientFxi[1]  = nodes_coord(0 + index, 1) - nodes_coord(2 + index, 1);
    TransverseGradientFxi[2]  = nodes_coord(0 + index, 2) - nodes_coord(2 + index, 2);

    TransverseGradientFeta[0] = nodes_coord(1 + index, 0) - nodes_coord(0 + index, 0);
    TransverseGradientFeta[1] = nodes_coord(1 + index, 1) - nodes_coord(0 + index, 1);
    TransverseGradientFeta[2] = nodes_coord(1 + index, 2) - nodes_coord(0 + index, 2);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAdd_B_Membrane(
        boost::numeric::ublas::bounded_matrix<double, 3, 18 > & mB_membrane,
        boost::numeric::ublas::bounded_matrix<double, 3, 1  > & mC_membrane,
        const boost::numeric::ublas::bounded_matrix<double, 2, 4 > & InPlaneCartesianDerivativesGauss,
        const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & InPlaneGradientFGauss,
        const int node_gauss
        )
{
    for (unsigned int i = 0; i < 4; i++)
    {
        unsigned int base = i * 3;
        if (i == 3)
        {
            base += node_gauss * 3;
        }
        for (unsigned int j = 0; j < 3; j++)
        {
            mB_membrane(0, base + j) += InPlaneCartesianDerivativesGauss(0, i) * InPlaneGradientFGauss(j, 0);
            mB_membrane(1, base + j) += InPlaneCartesianDerivativesGauss(1, i) * InPlaneGradientFGauss(j, 1);
            mB_membrane(2, base + j) += InPlaneCartesianDerivativesGauss(1, i) * InPlaneGradientFGauss(j, 0)
                                     +  InPlaneCartesianDerivativesGauss(0, i) * InPlaneGradientFGauss(j, 1);
        }
    }

    /* Calculate de componets of Cauchy tensor */
    // In plane auxiliar components
    array_1d<double, 3 > auxDeformationGradientF1;
    array_1d<double, 3 > auxDeformationGradientF2;

    auxDeformationGradientF1[0] = InPlaneGradientFGauss(0, 0);
    auxDeformationGradientF1[1] = InPlaneGradientFGauss(1, 0);
    auxDeformationGradientF1[2] = InPlaneGradientFGauss(2, 0);
    auxDeformationGradientF2[0] = InPlaneGradientFGauss(0, 1);
    auxDeformationGradientF2[1] = InPlaneGradientFGauss(1, 1);
    auxDeformationGradientF2[2] = InPlaneGradientFGauss(2, 1);

    mC_membrane(0, 0) += inner_prod(auxDeformationGradientF1, auxDeformationGradientF1);
    mC_membrane(1, 0) += inner_prod(auxDeformationGradientF2, auxDeformationGradientF2);
    mC_membrane(2, 0) += inner_prod(auxDeformationGradientF1, auxDeformationGradientF2);

}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAdd_Membrane_Kgeometric(
        boost::numeric::ublas::bounded_matrix<double, 36, 36 > & Kgeometricmembrane,
        const boost::numeric::ublas::bounded_matrix<double, 2, 4 > & mInPlaneCartesianDerivativesGauss1,
        const boost::numeric::ublas::bounded_matrix<double, 2, 4 > & mInPlaneCartesianDerivativesGauss2,
        const boost::numeric::ublas::bounded_matrix<double, 2, 4 > & mInPlaneCartesianDerivativesGauss3,
        const array_1d<double, 3 > & S_membrane,
        const int index
        )
{
    boost::numeric::ublas::bounded_matrix<double, 6, 6 > H = ZeroMatrix(6, 6);

    unsigned int ii;
    unsigned int jj;
    for (unsigned int i = 0; i < 4; i++)
    {
        for (unsigned int j = 0; j < 4; j++)
        {
            // Gauss 1
            ii = i;
            jj = j;
            H(ii, jj) += S_membrane[0] *  mInPlaneCartesianDerivativesGauss1(0, i) * mInPlaneCartesianDerivativesGauss1(0, j)
                       + S_membrane[1] *  mInPlaneCartesianDerivativesGauss1(1, i) * mInPlaneCartesianDerivativesGauss1(1, j)
                       + S_membrane[2] * (mInPlaneCartesianDerivativesGauss1(0, i) * mInPlaneCartesianDerivativesGauss1(1, j)
                                        + mInPlaneCartesianDerivativesGauss1(1, i) * mInPlaneCartesianDerivativesGauss1(0, j));

            // Gauss 2
            if (i ==  3)
            {
                ii = 4;
            }
            else
            {
                ii = i;
            }
            if (j ==  3)
            {
                jj = 4;
            }
            else
            {
                jj = j;
            }

            H(ii, jj) += S_membrane[0] *  mInPlaneCartesianDerivativesGauss2(0, i) * mInPlaneCartesianDerivativesGauss2(0, j)
                       + S_membrane[1] *  mInPlaneCartesianDerivativesGauss2(1, i) * mInPlaneCartesianDerivativesGauss2(1, j)
                       + S_membrane[2] * (mInPlaneCartesianDerivativesGauss2(0, i) * mInPlaneCartesianDerivativesGauss2(1, j)
                                        + mInPlaneCartesianDerivativesGauss2(1, i) * mInPlaneCartesianDerivativesGauss2(0, j));

            // Gauss 3
            if (i ==  3)
            {
                ii = 5;
            }
            else
            {
                ii = i;
            }
            if (j ==  3)
            {
                jj = 5;
            }
            else
            {
                jj = j;
            }

            H(ii, jj) += S_membrane[0] *  mInPlaneCartesianDerivativesGauss3(0, i) * mInPlaneCartesianDerivativesGauss3(0, j)
                       + S_membrane[1] *  mInPlaneCartesianDerivativesGauss3(1, i) * mInPlaneCartesianDerivativesGauss3(1, j)
                       + S_membrane[2] * (mInPlaneCartesianDerivativesGauss3(0, i) * mInPlaneCartesianDerivativesGauss3(1, j)
                                        + mInPlaneCartesianDerivativesGauss3(1, i) * mInPlaneCartesianDerivativesGauss3(0, j));
        }
    }

    H *= 0.333333333333333333333333333;

    // Assembling in Kgeometricmembrane
    unsigned int rowindex;
    unsigned int colindex;
    for (unsigned int i = 0; i < 6; i++)
    {
        if (i < 3)
        {
            rowindex = i * 3 + index;
        }
        else
        {
            rowindex = i * 3 + index + 9;
        }
        for (unsigned int j = i; j < 6; j++)
        {
            if (j < 3)
            {
                colindex = j * 3 + index;
            }
            else
            {
                colindex = j * 3 + index + 9;
            }
            for(unsigned int ii = 0; ii < 3; ii++)
            {
                Kgeometricmembrane(rowindex + ii,colindex + ii) += H (i, j);
                if (rowindex != colindex) // Skip diagonal
                {
                    Kgeometricmembrane(colindex + ii, rowindex + ii) += H (i, j); // Symmetric part
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAdd_B_Shear(
        boost::numeric::ublas::bounded_matrix<double, 2, 18 > & mB_shear,
        boost::numeric::ublas::bounded_matrix<double, 2, 1 > & C_shear,
        const boost::numeric::ublas::bounded_matrix<double, 6, 1 > & mTransversalCartesianDerivativesGauss1,
        const boost::numeric::ublas::bounded_matrix<double, 6, 1 > & mTransversalCartesianDerivativesGauss2,
        const boost::numeric::ublas::bounded_matrix<double, 6, 1 > & mTransversalCartesianDerivativesGauss3,
        const array_1d<double, 3 > & TransverseGradientFGauss1,
        const array_1d<double, 3 > & TransverseGradientFGauss2,
        const array_1d<double, 3 > & TransverseGradientFGauss3,
        const array_1d<double, 3 > & TransverseGradientFt,
        const array_1d<double, 3 > & TransverseGradientFxi,
        const array_1d<double, 3 > & TransverseGradientFeta,
        const boost::numeric::ublas::bounded_matrix<double, 2, 2 > & Jinv_plane,
        const int index
        )
{
    // Considering the Gauss point in the middle of the element
    double eta_p = 0.33333333333333333333;
    double xi_p  = 0.33333333333333333333;
    boost::numeric::ublas::bounded_matrix<double, 2, 3 > Pa;
    Pa(0, 0) = - xi_p;
    Pa(0, 1) = - xi_p;
    Pa(0, 2) = 1.0 - xi_p;
    Pa(1, 0) = eta_p;
    Pa(1, 1) = eta_p - 1.0;
    Pa(1, 2) = eta_p;

    boost::numeric::ublas::bounded_matrix<double, 3, 18 > aux_B_shear = ZeroMatrix(3, 18);

    /* First contribution*/
    for (unsigned int i = 0; i < 6; i++)
    {
        unsigned int base = i * 3;
        for (unsigned int j = 0; j < 3; j++)
        {
            aux_B_shear(0, base + j) += mTransversalCartesianDerivativesGauss1(i, 0) * TransverseGradientFt[j];
            aux_B_shear(1, base + j) += mTransversalCartesianDerivativesGauss2(i, 0) * TransverseGradientFxi[j];
            aux_B_shear(2, base + j) += mTransversalCartesianDerivativesGauss3(i, 0) * TransverseGradientFeta[j];
        }
    }

    /* Second contibution */
    for (unsigned int i = 0; i < 3; i++)
    {
        /* First row */
        aux_B_shear(0, i + index + 3) -= TransverseGradientFGauss1[i];
        aux_B_shear(0, i + index + 6) += TransverseGradientFGauss1[i];

        /* Second row */
        aux_B_shear(1, i + index)     += TransverseGradientFGauss2[i];
        aux_B_shear(1, i + index + 6) -= TransverseGradientFGauss2[i];

        /* Third row */
        aux_B_shear(2, i + index)     -= TransverseGradientFGauss3[i];
        aux_B_shear(2, i + index + 3) += TransverseGradientFGauss3[i];
    }

    boost::numeric::ublas::bounded_matrix<double, 2, 3 > aux_prod;
    noalias(aux_prod) = prod(Jinv_plane, Pa);
    noalias(mB_shear) = prod(aux_prod, aux_B_shear);

    // Calculating the components of C
    boost::numeric::ublas::bounded_matrix<double, 3, 1 > aux_C_shear;
    aux_C_shear(0, 0) =   inner_prod(TransverseGradientFt  , TransverseGradientFGauss1);
    aux_C_shear(1, 0) =   inner_prod(TransverseGradientFxi , TransverseGradientFGauss2);
    aux_C_shear(2, 0) =   inner_prod(TransverseGradientFeta, TransverseGradientFGauss3);

    noalias(C_shear) = prod(aux_prod, aux_C_shear);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAdd_Shear_Kgeometric(
    boost::numeric::ublas::bounded_matrix<double, 18, 18 > & Kgeometricshear,
    const boost::numeric::ublas::bounded_matrix<double, 6, 1 > & mTransversalCartesianDerivativesGauss1,
    const boost::numeric::ublas::bounded_matrix<double, 6, 1 > & mTransversalCartesianDerivativesGauss2,
    const boost::numeric::ublas::bounded_matrix<double, 6, 1 > & mTransversalCartesianDerivativesGauss3,
    const boost::numeric::ublas::bounded_matrix<double, 2, 2 > & Jinv_plane,
    const array_1d<double, 2 > & S_shear,
    const int index
    )
{
    double Q1 = 0.333333333333333333333333333 * (S_shear[0] * Jinv_plane(0, 0) + S_shear[1] * Jinv_plane(0, 1));
    double Q2 = 0.333333333333333333333333333 * (S_shear[0] * Jinv_plane(1, 0) + S_shear[1] * Jinv_plane(1, 1));

//    array_1d<double, 3 > q;
//    q[0] = -Q1 + Q2;
//    q[1] = -(Q1 + 2.0 * Q2);
//    q[2] = (2.0 * Q1 + Q2);

//    int delta;
//    if (index == 9)
//    {
//        delta = 3;
//    }
//    else
//    {
//        delta = 0;
//    }
//    for (unsigned int i = 0; i < 3; i++) // For each DOF
//    {
//        /* First assembling */
//        Kgeometricshear(i + index + 3, i + index + 3) -= q[0] * mTransversalCartesianDerivativesGauss1(1 + delta, 0);
//        Kgeometricshear(i + index + 6, i + index + 6) += q[0] * mTransversalCartesianDerivativesGauss1(2 + delta, 0);

//        /* Second assembling */
//        Kgeometricshear(i + index, i + index)         += q[1] * mTransversalCartesianDerivativesGauss2(0 + delta, 0);
//        Kgeometricshear(i + index + 6, i + index + 6) -= q[1] * mTransversalCartesianDerivativesGauss2(2 + delta, 0);
//        /* Third assembling */
//        Kgeometricshear(i + index, i + index)         -= q[2] * mTransversalCartesianDerivativesGauss3(0 + delta, 0);
//        Kgeometricshear(i + index + 3, i + index + 3) += q[2] * mTransversalCartesianDerivativesGauss3(1 + delta, 0);
//    }

    array_1d<double, 3 > n1; // Side node with + contribution (previous DOF position)
    array_1d<double, 3 > n2; // Side node with - contribution (previous DOF position)

//    if (index == 0)
//    {
//        n1[0] = 6;
//        n1[1] = 3;
//        n1[2] = 0;

//        n2[0] = 6;
//        n2[1] = 3;
//        n2[2] = 0;
//    }
//    else
//    {
//        n1[0] = 15;
//        n1[1] = 12;
//        n1[2] = 9;

//        n2[0] = 15;
//        n2[1] = 12;
//        n2[2] = 9;
//    }

    // Note: Technically this is the correct one
    if (index == 0)
    {
        n1[0] = 6;
        n1[1] = 0;
        n1[2] = 3;

        n2[0] = 3;
        n2[1] = 6;
        n2[2] = 0;
    }
    else
    {
        n1[0] = 15;
        n1[1] = 9;
        n1[2] = 12;

        n2[0] = 12;
        n2[1] = 15;
        n2[2] = 9;
    }

    double value = 0.0;
    for (unsigned int k = 0; k < 3; k++)
    {
        unsigned int l = 0; // Initializes DOF associated to N_3
        for (unsigned int i = 0; i < 6; i++) //  For each node
        {
            if (k == 0)
            {
                value = (-Q1 + Q2) *  mTransversalCartesianDerivativesGauss1(i, 0);
            }
            else if (k == 1)
            {
                value = -(Q1 + 2.0 * Q2) * mTransversalCartesianDerivativesGauss2(i, 0);
            }
            else if (k == 2)
            {
                value = (2.0 * Q1 + Q2) * mTransversalCartesianDerivativesGauss3(i, 0);
            }

            for (unsigned j = 0; j < 3; j++) // For each DOF (diagonal only)
            {
                Kgeometricshear(n1[k] + j, l + j) += value;
                Kgeometricshear(l + j, n1[k] + j) += value;
            }

            for (unsigned j = 0; j < 3; j++)  // For each DOF (diagonal only)
            {
                Kgeometricshear(n2[k] + j, l + j) -= value;
                Kgeometricshear(l + j, n2[k] + j) -= value;
            }

            l += 3; // Increment DOF position I
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAdd_B_Normal(
        boost::numeric::ublas::bounded_matrix<double, 1, 18 > & mB_normal,
        double & mC_normal,
        const boost::numeric::ublas::bounded_matrix<double, 6, 1 > & mTransversalCartesianDerivativesCenter,
        const array_1d<double, 3 > & TransversalDeformationGradientF
        )
{
        for (unsigned int i = 0; i < 6; i++)
        {
            unsigned int base = i * 3;
            mB_normal(0, base)     = mTransversalCartesianDerivativesCenter(i, 0) * TransversalDeformationGradientF[0];
            mB_normal(0, base + 1) = mTransversalCartesianDerivativesCenter(i, 0) * TransversalDeformationGradientF[1];
            mB_normal(0, base + 2) = mTransversalCartesianDerivativesCenter(i, 0) * TransversalDeformationGradientF[2];
        }

        mC_normal = inner_prod(TransversalDeformationGradientF, TransversalDeformationGradientF);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAdd_Normal_Kgeometric(
        boost::numeric::ublas::bounded_matrix<double, 18, 18 > & Kgeometricnormal,
        const boost::numeric::ublas::bounded_matrix<double, 6, 1 > & mTransversalCartesianDerivativesCenter,
        const double mS_normal
        )
{
    boost::numeric::ublas::bounded_matrix<double, 6, 6 > H = ZeroMatrix(6, 6);
    for (unsigned int i = 0; i < 6; i++)
    {
        double aux = mS_normal * mTransversalCartesianDerivativesCenter(i, 0);
        for (unsigned int j = 0; j < 6; j++)
        {
            H(i, j) =  aux * mTransversalCartesianDerivativesCenter(j, 0);
        }
    }

    noalias(H) = mS_normal * prod(mTransversalCartesianDerivativesCenter, trans(mTransversalCartesianDerivativesCenter));

    unsigned int rowindex;
    unsigned int colindex;
    for (unsigned int i = 0; i < 6; i++)
    {
        rowindex = i * 3;
        for (unsigned int j = 0; j < 6; j++)
        {
            colindex = j * 3;
            for(unsigned int ii = 0; ii < 3; ii++)
            {
                Kgeometricnormal(rowindex + ii,colindex + ii) += H(i, j);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

boost::numeric::ublas::bounded_matrix<double, 36, 1 > SprismElement3D6N::CalculateDisp(const int& step)
{
    KRATOS_TRY;

    boost::numeric::ublas::bounded_matrix<double, 36, 1 > disp_vec;

    WeakPointerVector< Node < 3 > >& nodal_neigb = this->GetValue(NEIGHBOUR_NODES);

    /* Element nodes */
    for (unsigned int index = 0; index < 6; index++)
    {
        array_1d<double,3> disp     = GetGeometry()[index].FastGetSolutionStepValue(DISPLACEMENT, step);
        disp_vec(index * 3,     0) = disp[0];
        disp_vec(index * 3 + 1, 0) = disp[1];
        disp_vec(index * 3 + 2, 0) = disp[2];
    }

    /* Neighbour nodes */
    int number_neigb = NumberOfActiveNeighbours(nodal_neigb);

    if (number_neigb == 6) // All the possible neighours
    {
        for (unsigned int index = 0; index < 6; index++)
        {
            array_1d<double,3> disp     = nodal_neigb[index].FastGetSolutionStepValue(DISPLACEMENT, step);
            disp_vec(18 + index * 3    , 0) = disp[0];
            disp_vec(18 + index * 3 + 1, 0) = disp[1];
            disp_vec(18 + index * 3 + 2, 0) = disp[2];
        }
    }
    else
    {
        for (unsigned int index = 0; index < 6; index++)
        {
            if (HasNeighbour(index, nodal_neigb[index]))
            {
                array_1d<double,3> disp     = nodal_neigb[index].FastGetSolutionStepValue(DISPLACEMENT, step);
                disp_vec(18 + index * 3    , 0) = disp[0];
                disp_vec(18 + index * 3 + 1, 0) = disp[1];
                disp_vec(18 + index * 3 + 2, 0) = disp[2];
            }
            else
            {
                disp_vec(18 + index * 3    , 0) = 0.0;
                disp_vec(18 + index * 3 + 1, 0) = 0.0;
                disp_vec(18 + index * 3 + 2, 0) = 0.0;
            }
        }
    }


    return disp_vec;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

boost::numeric::ublas::bounded_matrix<double, 36, 1 > SprismElement3D6N::GetVectorCurrentPosition()
{
    KRATOS_TRY;

    boost::numeric::ublas::bounded_matrix<double, 36, 1 > VectorCurrentPosition;

    WeakPointerVector< Node < 3 > >& nodal_neigb = this->GetValue(NEIGHBOUR_NODES);

    /* Element nodes */
    for (unsigned int index = 0; index < 6; index++)
    {
        array_1d<double,3> CurrentPosition = GetGeometry()[index].Coordinates();
        VectorCurrentPosition(index * 3,     0) = CurrentPosition[0];
        VectorCurrentPosition(index * 3 + 1, 0) = CurrentPosition[1];
        VectorCurrentPosition(index * 3 + 2, 0) = CurrentPosition[2];
    }

    /* Neighbour nodes */
    int number_neigb = NumberOfActiveNeighbours(nodal_neigb);

    if (number_neigb == 6) // All the possible neighours
    {
        for (unsigned int index = 0; index < 6; index++)
        {
            array_1d<double,3> CurrentPosition = nodal_neigb[index].Coordinates();
            VectorCurrentPosition(18 + index * 3    , 0) = CurrentPosition[0];
            VectorCurrentPosition(18 + index * 3 + 1, 0) = CurrentPosition[1];
            VectorCurrentPosition(18 + index * 3 + 2, 0) = CurrentPosition[2];
        }
    }
    else
    {
        for (unsigned int index = 0; index < 6; index++)
        {
            if (HasNeighbour(index, nodal_neigb[index]))
            {
                array_1d<double,3> CurrentPosition = nodal_neigb[index].Coordinates();
                VectorCurrentPosition(18 + index * 3    , 0) = CurrentPosition[0];
                VectorCurrentPosition(18 + index * 3 + 1, 0) = CurrentPosition[1];
                VectorCurrentPosition(18 + index * 3 + 2, 0) = CurrentPosition[2];
            }
            else
            {
                VectorCurrentPosition(18 + index * 3    , 0) = 0.0;
                VectorCurrentPosition(18 + index * 3 + 1, 0) = 0.0;
                VectorCurrentPosition(18 + index * 3 + 2, 0) = 0.0;
            }
        }
    }

    return VectorCurrentPosition;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::IntegrateInZeta(
        GeneralVariables& rVariables,
        const double& alpha_eas,
        const double& zeta_gauss,
        const double& rIntegrationWeight
        )
{
    KRATOS_TRY;
    double L_1 = 0.5 * (1.0 - zeta_gauss);
    double L_2 = 0.5 * (1.0 + zeta_gauss);

    double fact_eas = exp(2.0 * alpha_eas * zeta_gauss);

    /* INTEGRATE PK2 IN ZETA */
    // Integrate stresses in the reference configuration
    /* In plane stresses */
    // Lower
    mPK2.mS_membrane_lower(0) +=  L_1 * rIntegrationWeight * rVariables.StressVector[0]; // xx
    mPK2.mS_membrane_lower(1) +=  L_1 * rIntegrationWeight * rVariables.StressVector[1]; // yy
    mPK2.mS_membrane_lower(2) +=  L_1 * rIntegrationWeight * rVariables.StressVector[3]; // xy
    // Upper
    mPK2.mS_membrane_upper(0) +=  L_2 * rIntegrationWeight * rVariables.StressVector[0]; // xx
    mPK2.mS_membrane_upper(1) +=  L_2 * rIntegrationWeight * rVariables.StressVector[1]; // yy
    mPK2.mS_membrane_upper(2) +=  L_2 * rIntegrationWeight * rVariables.StressVector[3]; // xy

    /* Transversal stresses */ // Note: Order according to the Voigt Notation in the Wiki
    // Lower face
    mPK2.mS_shear_lower(0)    +=  L_1 * rIntegrationWeight * rVariables.StressVector[5]; // xz
    mPK2.mS_shear_lower(1)    +=  L_1 * rIntegrationWeight * rVariables.StressVector[4]; // yz
    // Upper face
    mPK2.mS_shear_upper(0)    +=  L_2 * rIntegrationWeight * rVariables.StressVector[5]; // xz
    mPK2.mS_shear_upper(1)    +=  L_2 * rIntegrationWeight * rVariables.StressVector[4]; // yz

    /* Normal stress */
    mPK2.mS_normal            +=  fact_eas * rIntegrationWeight * rVariables.StressVector[2]; // zz

    /* INTEGRATE EAS IN ZETA */
    // Calculate EAS residual
    mEAS.mrhs_alpha += rIntegrationWeight * zeta_gauss * rVariables.StressVector[2] * rVariables.C[2];

    // Calculate EAS stiffness
    mEAS.mstiff_alpha += rIntegrationWeight * zeta_gauss * zeta_gauss * rVariables.C[2]
            * (rVariables.ConstitutiveMatrix(2, 2) * rVariables.C[2] + 2.0 * rVariables.StressVector[2]);

    boost::numeric::ublas::bounded_matrix<double, 1, 36 > B3;
    boost::numeric::ublas::bounded_matrix<double, 1,  6 > D3;

    for (unsigned int i = 0; i < 6; i++)
    {
        D3(0, i) = rVariables.ConstitutiveMatrix(2, i);
    }
    for (unsigned int i = 0; i < 36; i++)
    {
        B3(0, i) = rVariables.B(2, i);
    }

    // Calculate H operator
    noalias(mEAS.mH_EAS) += rIntegrationWeight * zeta_gauss
            * (rVariables.C[2] * prod(D3, rVariables.B) + 2.0 * rVariables.StressVector[2] * B3);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAddLHS(
        LocalSystemComponents& rLocalSystem,
        GeneralVariables& rVariables,
        ConstitutiveLaw::Parameters& rValues,
        double& alpha_eas
        )
{
    /* Contributions of the stiffness matrix calculated on the reference configuration */
    if( rLocalSystem.CalculationFlags.Is( SprismElement3D6N::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
    {
        std::vector<MatrixType>& rLeftHandSideMatrices = rLocalSystem.GetLeftHandSideMatrices();
        const std::vector< Variable< MatrixType > >& rLeftHandSideVariables = rLocalSystem.GetLeftHandSideVariables();

        for( unsigned int i = 0; i < rLeftHandSideVariables.size(); i++ )
        {
            bool calculated = false;
            /* Calculate the Material Stiffness Matrix */
            if( rLeftHandSideVariables[i] == MATERIAL_STIFFNESS_MATRIX )
            {
                /* Reading integration points */
                const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

                for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
                {
                    double zeta_gauss = 2.0 * integration_points[PointNumber].Z() - 1.0;

                    /* Assemble B */
                    this->CalculateDeformationMatrix(rVariables.B, zeta_gauss, alpha_eas);

                    // Compute element kinematics C, F ...
                    this->CalculateKinematics(rVariables, PointNumber, alpha_eas, zeta_gauss);
                    this->CbartoFbar(rVariables, PointNumber);

                    // Set general variables to constitutivelaw parameters
                    this->SetGeneralVariables(rVariables, rValues, PointNumber);

                    // Compute stresses and constitutive parameters
                    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, rVariables.StressMeasure);

                    // Calculating weights for integration on the "reference configuration"
                    double IntegrationWeight = integration_points[PointNumber].Weight() * rVariables.detJ;

                    /* Operation performed: add Km to the LefsHandSideMatrix */
                    this->CalculateAndAddKuum( rLeftHandSideMatrices[i], rVariables, IntegrationWeight);
                }
                calculated = true;
            }

            /* Calculate the Geometric Stiffness Matrix */
            if( rLeftHandSideVariables[i] == GEOMETRIC_STIFFNESS_MATRIX )
            {
                /* Operation performed: add Kg to the LefsHandSideMatrix */
                this->CalculateAndAddKuug( rLeftHandSideMatrices[i] );
                calculated = true;
            }

            /* Implicit or explicit EAS update*/
            bool eas_imp = true;
            if( GetProperties().Has(EAS_IMP) )
            {
                eas_imp = GetProperties()[EAS_IMP];
            }

            if (eas_imp == true)
            {
                /* Apply EAS stabilization */
                ApplyEASLHS(rLeftHandSideMatrices[i]);
            }

            if(calculated == false)
            {
                KRATOS_THROW_ERROR(std::logic_error, " ELEMENT can not supply the required local system variable: ",rLeftHandSideVariables[i]);
            }
        }
    }
    else
    {
        MatrixType& LeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        /* Calculate the Material Stiffness Matrix */
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            double zeta_gauss = 2.0 * integration_points[PointNumber].Z() - 1.0;

            /* Assemble B */
            this->CalculateDeformationMatrix(rVariables.B, zeta_gauss, alpha_eas);

            // Compute element kinematics C, F ...
            this->CalculateKinematics(rVariables, PointNumber, alpha_eas, zeta_gauss);
            this->CbartoFbar(rVariables, PointNumber);

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(rVariables, rValues, PointNumber);

            // Compute stresses and constitutive parameters
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, rVariables.StressMeasure);

            // Calculating weights for integration on the "reference configuration"
            double IntegrationWeight = integration_points[PointNumber].Weight() * rVariables.detJ;

            /* Operation performed: add Km to the LefsHandSideMatrix */
            this->CalculateAndAddKuum( LeftHandSideMatrix, rVariables, IntegrationWeight);
        }

        /* Calculate the Geometric Stiffness Matrix */
        /* Operation performed: add Kg to the LefsHandSideMatrix */
        this->CalculateAndAddKuug( LeftHandSideMatrix );

        /* Implicit or explicit EAS update*/
        bool eas_imp = true;
        if( GetProperties().Has(EAS_IMP) )
        {
            eas_imp = GetProperties()[EAS_IMP];
        }

        if (eas_imp == true)
        {
            /* Apply EAS stabilization */
            ApplyEASLHS(LeftHandSideMatrix);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAddDynamicLHS(MatrixType& rLeftHandSideMatrix)
{

  // Mass matrix
  WeakPointerVector< Node < 3 > >& nodal_neigb = this->GetValue(NEIGHBOUR_NODES);

  double Density = GetProperties()[DENSITY];

  unsigned int number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(nodal_neigb);
  unsigned int MatSize = number_of_nodes * 3;

  if (rLeftHandSideMatrix.size1() != MatSize)
  {
      rLeftHandSideMatrix.resize(MatSize, MatSize, false);
  }

  noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize, MatSize);

  double Volume = GetGeometry().Volume();
  double TotalMass = Volume * Density;

  TotalMass /= 72.0; // Dividing for the coefficient
  for (unsigned int i = 0; i < 6; i++) // Main nodes
  {
      for (unsigned int j = 0; j < 3; j++) // DOF (X, Y, Z)
      {
          unsigned int index = i * 3 + j;
          if (i == 0)
          {
              // Superior band
              rLeftHandSideMatrix(index, index +  3) = 2.0 * TotalMass;
              rLeftHandSideMatrix(index, index +  6) = 2.0 * TotalMass;
              rLeftHandSideMatrix(index, index +  9) = 2.0 * TotalMass;
              rLeftHandSideMatrix(index, index + 12) = TotalMass;
              rLeftHandSideMatrix(index, index + 15) = TotalMass;
              // Symmetric part
              rLeftHandSideMatrix(index +  3, index) = 2.0 * TotalMass;
              rLeftHandSideMatrix(index +  6, index) = 2.0 * TotalMass;
              rLeftHandSideMatrix(index +  9, index) = 2.0 * TotalMass;
              rLeftHandSideMatrix(index + 12, index) = TotalMass;
              rLeftHandSideMatrix(index + 15, index) = TotalMass;
          }
          else if (i == 1)
          {
              // Superior band
              rLeftHandSideMatrix(index, index +  3) = 2.0 * TotalMass;
              rLeftHandSideMatrix(index, index +  6) = TotalMass;
              rLeftHandSideMatrix(index, index +  9) = 2.0 * TotalMass;
              rLeftHandSideMatrix(index, index + 12) = TotalMass;
              // Symmetric part
              rLeftHandSideMatrix(index +  3, index) = 2.0 * TotalMass;
              rLeftHandSideMatrix(index +  6, index) = TotalMass;
              rLeftHandSideMatrix(index +  9, index) = 2.0 * TotalMass;
              rLeftHandSideMatrix(index + 12, index) = TotalMass;
          }
          else if (i == 2)
          {
              // Superior band
              rLeftHandSideMatrix(index, index + 3) = TotalMass;
              rLeftHandSideMatrix(index, index + 6) = TotalMass;
              rLeftHandSideMatrix(index, index + 9) = 2.0 * TotalMass;
              // Symmetric part
              rLeftHandSideMatrix(index + 3, index) = TotalMass;
              rLeftHandSideMatrix(index + 6, index) = TotalMass;
              rLeftHandSideMatrix(index + 9, index) = 2.0 * TotalMass;
          }
          else if (i == 3)
          {
              // Superior band
              rLeftHandSideMatrix(index, index + 3) = 2.0 * TotalMass;
              rLeftHandSideMatrix(index, index + 6) = 2.0 * TotalMass;
              // Symmetric part
              rLeftHandSideMatrix(index + 3, index) = 2.0 * TotalMass;
              rLeftHandSideMatrix(index + 6, index) = 2.0 * TotalMass;
          }
          else if (i == 4)
          {
              // Superior band
              rLeftHandSideMatrix(index, index + 3) = 2.0 * TotalMass;
              // Symmetric part
              rLeftHandSideMatrix(index + 3, index) = 2.0 * TotalMass;
          }

          // Diagonal part
          rLeftHandSideMatrix(index, index) = 4.0 * TotalMass;
      }
  }

//  KRATOS_WATCH( rLeftHandSideMatrix );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAddRHS(
        LocalSystemComponents& rLocalSystem,
        GeneralVariables& rVariables,
        Vector& rVolumeForce,
        double& alpha_eas
        )
{
    /* Contribution of the internal and external forces */
    if( rLocalSystem.CalculationFlags.Is( SprismElement3D6N::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
    {
        std::vector<VectorType>& RightHandSideVectors = rLocalSystem.GetRightHandSideVectors();
        const std::vector< Variable< VectorType > >& rRightHandSideVariables = rLocalSystem.GetRightHandSideVariables();
        for( unsigned int i = 0; i < rRightHandSideVariables.size(); i++ )
        {
            bool calculated = false;
            if( rRightHandSideVariables[i] == EXTERNAL_FORCES_VECTOR )
            {
                /* Operation performed: RightHandSideVector += ExtForce */
                this->CalculateAndAddExternalForces( RightHandSideVectors[i], rVariables, rVolumeForce );
                calculated = true;
            }

            if( rRightHandSideVariables[i] == INTERNAL_FORCES_VECTOR )
            {
                /* Operation performed: RightHandSideVector -= IntForce */
                this->CalculateAndAddInternalForces( RightHandSideVectors[i], alpha_eas );
                calculated = true;
            }

            if(calculated == false)
            {
                KRATOS_THROW_ERROR( std::logic_error, " ELEMENT can not supply the required local system variable: ", rRightHandSideVariables[i] );
            }
        }
    }
    else
    {
        VectorType& RightHandSideVector = rLocalSystem.GetRightHandSideVector();

        /* Operation performed: RightHandSideVector += ExtForce */
        this->CalculateAndAddExternalForces( RightHandSideVector, rVariables, rVolumeForce );

        /* Operation performed: RightHandSideVector -= IntForce */
        this->CalculateAndAddInternalForces( RightHandSideVector, alpha_eas );
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAddDynamicRHS(
        VectorType& rRightHandSideVector,
        GeneralVariables& rVariables,
        ProcessInfo& rCurrentProcessInfo,
        double& rIntegrationWeight
        )
{
    // DO NOTHING (right now)
    //KRATOS_WATCH( rRightHandSideVector )
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAddKuum(
        MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight
        )
{
    KRATOS_TRY;
    
    boost::numeric::ublas::bounded_matrix<double, 36, 36 > mK; // Local stiffness matrix
    boost::numeric::ublas::bounded_matrix<double, 6, 36 > aux_CB;

    /* Calculate K */
    noalias(aux_CB) = prod(rVariables.ConstitutiveMatrix, rVariables.B);
    noalias(mK) = rIntegrationWeight * prod(trans(rVariables.B), aux_CB);

    for (unsigned int i = 0; i < 36; i++)
    {
        if (mid_vec[i] < 1000)
        {
            for (unsigned int j = 0; j < 36; j++)
            {
                if (mid_vec[j] < 1000)
                {
                    rLeftHandSideMatrix(mid_vec[i], mid_vec[j]) += mK(i, j);
                }
            }
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix)
{
    KRATOS_TRY;

    /* The stress is already integrated, we just calculate it once */

    /* Auxiliar stiffness matrix */
    boost::numeric::ublas::bounded_matrix<double, 18, 18 > aux_mK = ZeroMatrix(18, 18); // Auxiliar stiffness matrix
    boost::numeric::ublas::bounded_matrix<double, 36, 36 > mK     = ZeroMatrix(36, 36); // Stiffness matrix

    /* COMPUTATION OF GEOMETRIC STIFFNESS MATRIX */

    /* MEMBRANE CONTRIBUTION */
    /* Adding the geometric membrane stiffness */
    // Lower face
    CalculateAndAdd_Membrane_Kgeometric(mK, mCC.mInPlaneCartesianDerivativesGauss1, mCC.mInPlaneCartesianDerivativesGauss2,\
                                                mCC.mInPlaneCartesianDerivativesGauss3, mPK2.mS_membrane_lower, 0);
    // Upper face
    CalculateAndAdd_Membrane_Kgeometric(mK, mCC.mInPlaneCartesianDerivativesGauss4, mCC.mInPlaneCartesianDerivativesGauss5,\
                                                mCC.mInPlaneCartesianDerivativesGauss6, mPK2.mS_membrane_upper, 9);

//    /* SHEAR CONTRIBUTION */
//    /* Adding the geometric shear stiffness */
//    // Lower face
//    CalculateAndAdd_Shear_Kgeometric(aux_mK, mCC.mTransversalCartesianDerivativesGauss1,\
//                                          mCC.mTransversalCartesianDerivativesGauss2, mCC.mTransversalCartesianDerivativesGauss3,\
//                                          mCC.mJinv_plane_lower, mPK2.mS_shear_lower, 0);
//    // Upper face
//    CalculateAndAdd_Shear_Kgeometric(aux_mK, mCC.mTransversalCartesianDerivativesGauss4,\
//                                          mCC.mTransversalCartesianDerivativesGauss5, mCC.mTransversalCartesianDerivativesGauss6,\
//                                          mCC.mJinv_plane_upper, mPK2.mS_shear_upper, 9);

    /* NORMAL TRANSVERSE */
    /* Adding the geometric normal stiffness */
    CalculateAndAdd_Normal_Kgeometric(aux_mK, mCC.mTransversalCartesianDerivativesCenter, mPK2.mS_normal);

    // Transfering to the complete stiffness matrix
    for (unsigned int i = 0; i < 18; i++)
    {
        for (unsigned int j = 0; j < 18; j++)
        {
            mK(i, j) += aux_mK(i, j);
        }
    }

    for (unsigned int i = 0; i < 36; i++)
    {
        if (mid_vec[i] < 1000)
        {
            for (unsigned int j = 0; j < 36; j++)
            {
                if (mid_vec[j] < 1000)
                {
                    rLeftHandSideMatrix(mid_vec[i], mid_vec[j]) += mK(i, j);
                }
            }
        }
    }

//    std::cout<<std::endl;
//    std::cout<<" Kmat + Kgeo "<<rLeftHandSideMatrix<<std::endl;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::ApplyEASLHS(MatrixType& rLeftHandSideMatrix )
{
    KRATOS_TRY;

    // Allocate auxiliar vector and matrix
    boost::numeric::ublas::bounded_matrix<double, 36, 36 > lhs_aux = ZeroMatrix(36, 36);

    noalias(lhs_aux) -= prod(trans(mEAS.mH_EAS), mEAS.mH_EAS) / mEAS.mstiff_alpha;

    // Note: rIntegrationWeight already considered in the integration
    for (unsigned int i = 0; i < 36; i++)
    {
        if (mid_vec[i] < 1000)
        {
            for (unsigned int j = 0; j < 36; j++)
            {
                if (mid_vec[j] < 1000)
                {
                    rLeftHandSideMatrix(mid_vec[i], mid_vec[j]) += lhs_aux(i, j);
                }
            }
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::ApplyEASRHS(
        boost::numeric::ublas::bounded_matrix<double, 36, 1 > & rhs_full,
        double& alpha_eas
        )
{
    KRATOS_TRY;

    /* Calculate the RHS */
    noalias(rhs_full) -= trans(mEAS.mH_EAS) * mEAS.mrhs_alpha / mEAS.mstiff_alpha;

//////    std::cout << "\n Element: " << this->Id() << std::endl;
//    if (this->Id() == 1)
//    {
//        std::cout << "Before update" << std::endl;
//        KRATOS_WATCH(alpha_eas);
//    }

    // Update ALPHA_EAS
    alpha_eas -= mEAS.mrhs_alpha / mEAS.mstiff_alpha;

//    if (this->Id() == 1)
//    {
//        std::cout << "After update" << std::endl;
//        KRATOS_WATCH(mEAS.mrhs_alpha);
//        KRATOS_WATCH(mEAS.mstiff_alpha);
//        KRATOS_WATCH(mEAS.mH_EAS);
//        KRATOS_WATCH(alpha_eas);
//    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAddExternalForces(
        VectorType& rRightHandSideVector,
        GeneralVariables& rVariables,
        Vector& rVolumeForce
        )
{
    KRATOS_TRY;

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    double aux_div = number_of_nodes;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = 3 * i;
        for ( unsigned int j = 0; j < 3; j++ )
        {
            rRightHandSideVector[index + j] += rVolumeForce[j]/aux_div;
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAddInternalForces(
        VectorType& rRightHandSideVector,
        double& alpha_eas
        )
{
    KRATOS_TRY;

    boost::numeric::ublas::bounded_matrix<double, 36, 1 > rhs_full = ZeroMatrix(36, 1);

    int aux_index = 0;
    for (unsigned int i = 0; i < 18; i++)
    {
        if (i == 9)
        {
            aux_index += 9;
        }

        /* Calculate residual forces */
        /* Apply membrane stress, adding the in-plane nodal force contribution */
        /* Nodes 1-3  and 7-9 */
        rhs_full(aux_index + i, 0)     += mPK2.mS_membrane_lower[0] * mCC.mB_membrane_lower(0, i); // xx
        rhs_full(aux_index + i, 0)     += mPK2.mS_membrane_lower[1] * mCC.mB_membrane_lower(1, i); // yy
        rhs_full(aux_index + i, 0)     += mPK2.mS_membrane_lower[2] * mCC.mB_membrane_lower(2, i); // xy

        /* Nodes 4-6  and 10-12 */
        rhs_full(aux_index + i + 9, 0) += mPK2.mS_membrane_upper[0] * mCC.mB_membrane_upper(0, i); // xx
        rhs_full(aux_index + i + 9, 0) += mPK2.mS_membrane_upper[1] * mCC.mB_membrane_upper(1, i); // yy
        rhs_full(aux_index + i + 9, 0) += mPK2.mS_membrane_upper[2] * mCC.mB_membrane_upper(2, i); // xy

        /* Apply transversal forces */
        /* Apply shear stress, adding the transverse nodal force contribution */
        rhs_full(i, 0) += mPK2.mS_shear_lower[0] * mCC.mB_shear_lower(0, i); // xz
        rhs_full(i, 0) += mPK2.mS_shear_lower[1] * mCC.mB_shear_lower(1, i); // yz
        rhs_full(i, 0) += mPK2.mS_shear_upper[0] * mCC.mB_shear_upper(0, i); // xz
        rhs_full(i, 0) += mPK2.mS_shear_upper[1] * mCC.mB_shear_upper(1, i); // yz

        /* Apply normal transverse stress */
        rhs_full(i, 0) += mPK2.mS_normal * mCC.mB_normal(0, i); // zz
    }

    /* Apply EAS stabilization */
    ApplyEASRHS(rhs_full, alpha_eas);

    for (unsigned int i = 0; i < 36; i++)
    {
        if (mid_vec[i] < 1000)
        {
            rRightHandSideVector[mid_vec[i]] -= rhs_full(i, 0);
        }
    }

//    if (this->Id() == 1)
//    {
//        KRATOS_WATCH(rhs_full);
//    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::SetGeneralVariables(
        GeneralVariables& rVariables,
        ConstitutiveLaw::Parameters& rValues,
        const int & rPointNumber
        )
{
    if(rVariables.detF < 0)
    {
        std::cout<<" Element: "<<this->Id()<<std::endl;
        unsigned int number_of_nodes = GetGeometry().PointsNumber();

        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            array_1d<double, 3> &CurrentPosition  = GetGeometry()[i].Coordinates();
            array_1d<double, 3> & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3> & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
            array_1d<double, 3> PreviousPosition  = CurrentPosition - (CurrentDisplacement-PreviousDisplacement);
            std::cout<<" NODE ["<<GetGeometry()[i].Id()<<"]: "<<PreviousPosition<<" (Cur: "<<CurrentPosition<<") "<<std::endl;
            std::cout<<" ---Disp: "<<CurrentDisplacement<<" (Pre: "<<PreviousDisplacement<<")"<<std::endl;
        }

        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            if(GetGeometry()[i].SolutionStepsDataHas(CONTACT_FORCE))
            {
                array_1d<double, 3 > & PreContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE,1);
                array_1d<double, 3 > & ContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE);
                std::cout<<" ---Contact_Force: (Pre:"<<PreContactForce<<", Cur:"<<ContactForce<<") "<<std::endl;
            }
            else
            {
                std::cout<<" ---Contact_Force: NULL "<<std::endl;
            }
        }
        KRATOS_WATCH(rVariables.F);
        KRATOS_THROW_ERROR( std::invalid_argument," SPRISM ELEMENT INVERTED: |F| < 0  detF = ", rVariables.detF );
    }

    rValues.SetDeterminantF(rVariables.detF);
    rValues.SetDeformationGradientF(rVariables.F);
    rValues.SetStrainVector(rVariables.StrainVector);
    rValues.SetStressVector(rVariables.StressVector);
    rValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);

    // Adding the standard prism shape functions
    rValues.SetShapeFunctionsValues(rVariables.N);
    rValues.SetShapeFunctionsDerivatives(rVariables.DN_DX);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::InitializeSystemMatrices(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        Flags& rCalculationFlags
        )
{
    // Resizing as needed the LHS
    WeakPointerVector< Node < 3 > >& nodal_neigb = this->GetValue(NEIGHBOUR_NODES);
    const unsigned int number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(nodal_neigb);
    unsigned int MatSize = number_of_nodes * 3;

    if ( rCalculationFlags.Is(SprismElement3D6N::COMPUTE_LHS_MATRIX) ) // Calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
        {
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );
        }

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); // Resetting LHS
    }

    // Resizing as needed the RHS
    if ( rCalculationFlags.Is(SprismElement3D6N::COMPUTE_RHS_VECTOR) ) // Calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != MatSize )
        {
            rRightHandSideVector.resize( MatSize, false );
        }

    rRightHandSideVector = ZeroVector( MatSize ); // Resetting RHS
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::InitializeMaterial()
{
    KRATOS_TRY;

    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(),
                    row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
        }
    }
    else
    {
        KRATOS_THROW_ERROR( std::logic_error, "A constitutive law needs to be specified for the element with ID ", this->Id() );
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::ResetConstitutiveLaw()
{
    KRATOS_TRY;

    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            mConstitutiveLawVector[i]->ResetMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
        }
    }
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::ClearNodalForces()
{
    KRATOS_TRY;

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
      if( GetGeometry()[i].SolutionStepsDataHas(EXTERNAL_FORCE) && GetGeometry()[i].SolutionStepsDataHas(INTERNAL_FORCE) ){

        array_1d<double, 3 > & ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_FORCE);
        array_1d<double, 3 > & InternalForce = GetGeometry()[i].FastGetSolutionStepValue(INTERNAL_FORCE);

        GetGeometry()[i].SetLock();
        ExternalForce.clear();
        InternalForce.clear();
        GetGeometry()[i].UnSetLock();
      }
    }

    KRATOS_CATCH( "" );
}

/******************************* COMPUTE KINEMATICS ********************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateKinematics(
        GeneralVariables& rVariables,
        const int& rPointNumber,
        const double& alpha_eas,
        const double& zeta_gauss
        )
{
    KRATOS_TRY;

    // Jacobian Determinant for the isoparametric and numerical integration
    rVariables.detJ = mDetJ0[rPointNumber];

    double L_1 = 0.5 * (1.0 - zeta_gauss);
    double L_2 = 0.5 * (1.0 + zeta_gauss);

    double fact_eas = exp(2.0 * alpha_eas * zeta_gauss);  // EAS factor

    /* Assemble C */
    rVariables.C[0] = L_1 * mCC.mC_membrane_lower(0, 0) + L_2 * mCC.mC_membrane_upper(0, 0); // xx
    rVariables.C[1] = L_1 * mCC.mC_membrane_lower(1, 0) + L_2 * mCC.mC_membrane_upper(1, 0); // yy
    rVariables.C[2] = fact_eas * mCC.mC_normal;                                              // zz
    rVariables.C[3] = L_1 * mCC.mC_membrane_lower(2, 0) + L_2 * mCC.mC_membrane_upper(2, 0); // xy
    rVariables.C[4] = L_1 * mCC.mC_shear_lower(1, 0)    + L_2 * mCC.mC_shear_upper(1, 0);    // yz
    rVariables.C[5] = L_1 * mCC.mC_shear_lower(0, 0)    + L_2 * mCC.mC_shear_upper(0, 0);    // xz

    rVariables.detF = rVariables.C[0] * rVariables.C[1] * rVariables.C[2] + 2 * rVariables.C[3] * rVariables.C[4] * rVariables.C[5]\
                    - rVariables.C[5] * rVariables.C[5] * rVariables.C[1] -     rVariables.C[4] * rVariables.C[4] * rVariables.C[0]\
                    - rVariables.C[3] * rVariables.C[3] * rVariables.C[2];

    if (rVariables.detF < 1.0e-8)
    {
        KRATOS_WATCH(rVariables.C);
        KRATOS_THROW_ERROR( std::logic_error, "The determinant of C is zero or negative.  det(C): ", rVariables.detF );
    }

    rVariables.detF = sqrt(rVariables.detF);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CbartoFbar(
        GeneralVariables& rVariables,
        const int& rPointNumber
        )
{
    KRATOS_TRY;

    /* We perform a polar decomposition of the C_bar and F(regular) to obtain F_bar */

    /* Decompose C_bar */
    Matrix EigenValuesMatrix;
    EigenValuesMatrix.resize(3, 3, false);
    Matrix EigenVectorsMatrix;
    EigenVectorsMatrix.resize(3, 3, false);

    // Assemble matrix C_bar
    Matrix C_bar;
    C_bar.resize(3, 3, false);
    StructuralMechanicsMathUtilities::VectorToTensor(rVariables.C, C_bar);

//    if (this->Id() == 1)
//    {
//        KRATOS_WATCH(rVariables.C);
//        KRATOS_WATCH(C_bar - IdentityMatrix(3));
//    }

    // Decompose matrix C_bar
    StructuralMechanicsMathUtilities::EigenVectors(C_bar, EigenVectorsMatrix, EigenValuesMatrix, 1e-24, 100);

    for (unsigned int i = 0; i < 3; i++)
    {
        EigenValuesMatrix(i, i) = sqrt(EigenValuesMatrix(i, i));
    }

    boost::numeric::ublas::bounded_matrix<double, 3, 3 > U_bar;
    noalias(U_bar) = prod( EigenValuesMatrix, EigenVectorsMatrix );

//    if (this->Id() == 1)
//    {
//        KRATOS_WATCH(U_bar);
//        KRATOS_WATCH(prod(trans(U_bar), U_bar)-IdentityMatrix(3));
//    }

    /* Decompose F */
    boost::numeric::ublas::bounded_matrix<double, 3, 3 > F;
    noalias(F) = prod( rVariables.J[rPointNumber], mInvJ0[rPointNumber] );

//    if (this->Id() == 1)
//    {
//        std::cout << "F of Kratos" << std::endl;
//        KRATOS_WATCH(F - IdentityMatrix(3));
//    }

    boost::numeric::ublas::bounded_matrix<double, 3, 3 > C;
    noalias(C) = prod( trans(F), F );

    // Decompose matrix C
    StructuralMechanicsMathUtilities::EigenVectors(C, EigenVectorsMatrix, EigenValuesMatrix, 1e-24, 100);

    for (unsigned int i = 0; i < 3; i++)
    {
        EigenValuesMatrix(i, i) = sqrt(EigenValuesMatrix(i, i));
    }

    boost::numeric::ublas::bounded_matrix<double, 3, 3 > U;
    boost::numeric::ublas::bounded_matrix<double, 3, 3 > invU;
    boost::numeric::ublas::bounded_matrix<double, 3, 3 > R;
    noalias(U) = prod( EigenValuesMatrix, EigenVectorsMatrix );

    StructuralMechanicsMathUtilities::InvMat3x3(U, invU);
    noalias(R) = prod( F, invU );

    StructuralMechanicsMathUtilities::InvMat3x3(U_bar, invU);

//    if (this->Id() == 1)
//    {
//        KRATOS_WATCH(R);
//        KRATOS_WATCH(prod(rVariables.F, invU));

//        KRATOS_WATCH(prod(trans(R),R)-IdentityMatrix(3));
//    }

    /* Calculate F_bar */
    noalias(rVariables.F) = prod(R, U_bar);

//    if (this->Id() == 1)
//    {
//        KRATOS_WATCH(rVariables.F -IdentityMatrix(3) );
//        KRATOS_WATCH(prod(trans(rVariables.F),rVariables.F) - IdentityMatrix(3));
//    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateDeformationMatrix(
        Matrix& rB,
        const double& zeta_gauss,
        const double& alpha_eas
        )
{
    KRATOS_TRY;

    rB.clear(); // Set all components to zero

    double L_1 = 0.5 * (1.0 - zeta_gauss);
    double L_2 = 0.5 * (1.0 + zeta_gauss);

    double fact_eas = exp(2.0 * alpha_eas * zeta_gauss); // EAS factor

    for (unsigned int index = 0; index < 9; index++)
    {
        /* Element nodes */ // Note: It's important to consider the Voigt notation order considered in Kratos
        // Lower face
        rB(0, index)      = L_1 * mCC.mB_membrane_lower(0, index);  // xx
        rB(1, index)      = L_1 * mCC.mB_membrane_lower(1, index);  // yy
        rB(2, index)      = fact_eas * mCC.mB_normal(0, index);     // zz
        rB(3, index)      = L_1 * mCC.mB_membrane_lower(2, index);  // xy
        rB(4, index)      = L_1 * mCC.mB_shear_lower(1, index) + L_2 * mCC.mB_shear_upper(1, index); // yz
        rB(5, index)      = L_1 * mCC.mB_shear_lower(0, index) + L_2 * mCC.mB_shear_upper(0, index); // xz
        // Upper face
        rB(0, index + 9)  = L_2 * mCC.mB_membrane_upper(0, index);  // xx
        rB(1, index + 9)  = L_2 * mCC.mB_membrane_upper(1, index);  // yy
        rB(2, index + 9)  = fact_eas * mCC.mB_normal(0, index + 9); // zz
        rB(3, index + 9)  = L_2 * mCC.mB_membrane_upper(2, index);  // xy
        rB(4, index + 9)  = L_1 * mCC.mB_shear_lower(1, index + 9) + L_2 * mCC.mB_shear_upper(1, index + 9); // yz
        rB(5, index + 9)  = L_1 * mCC.mB_shear_lower(0, index + 9) + L_2 * mCC.mB_shear_upper(0, index + 9); // xz

        /* Neighbour nodes */
        // Lower face
        rB(0, index + 18) = L_1 * mCC.mB_membrane_lower(0, index + 9); // xx
        rB(1, index + 18) = L_1 * mCC.mB_membrane_lower(1, index + 9); // yy
        rB(3, index + 18) = L_1 * mCC.mB_membrane_lower(2, index + 9); // xy
        // Upper face
        rB(0, index + 27) = L_2 * mCC.mB_membrane_upper(0, index + 9); // xx
        rB(1, index + 27) = L_2 * mCC.mB_membrane_upper(1, index + 9); // yy
        rB(3, index + 27) = L_2 * mCC.mB_membrane_upper(2, index + 9); // xy
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::InitializeGeneralVariables(GeneralVariables& rVariables)
{
    rVariables.detF = 1.0;
    rVariables.detJ = 1.0;
    rVariables.F.resize(3, 3, false);
    rVariables.B.resize(6, 36, false);
    rVariables.F = IdentityMatrix(3);
    rVariables.C.resize(6, false);
    rVariables.ConstitutiveMatrix.resize(6, 6, false);
    rVariables.StrainVector.resize(6, false);
    rVariables.StressVector.resize(6, false);

    rVariables.DN_DX.resize( 6, 3, false);

    // Reading shape functions
    rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

    // Calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
    rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::FinalizeStepVariables(
        GeneralVariables & rVariables,
        const int& rPointNumber
        )
{
   // Note: For the future, to include internal variables
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::GetHistoricalVariables(
        GeneralVariables& rVariables,
        const int& rPointNumber
        )
{
//    /* Deformation Gradient F ( set to identity ) */
//    unsigned int size =  rVariables.F.size1();

//    rVariables.detF  = 1;
//    rVariables.F     = IdentityMatrix(size);

}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateLinearStress(GeneralVariables& rVariables)
{
    KRATOS_TRY;
    
    noalias(rVariables.StressVector) = prod(rVariables.ConstitutiveMatrix, rVariables.StrainVector);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateLinearIsotropicStress(GeneralVariables& rVariables)
{
    KRATOS_TRY;

    double E = GetProperties()[YOUNG_MODULUS];
    double NU = GetProperties()[POISSON_RATIO];

    double K = E /(3.0 * (1.0 - 2.0 * NU));
    double G = E /(2.0 * (1.0 +       NU));

    double delta;
    delta = ( rVariables.C[0] + rVariables.C[1] + rVariables.C[2] - 3.0) / 6.0;

    rVariables.StrainVector = rVariables.C;
    for (unsigned int i = 0; i < 3 ; i++)
    {
        rVariables.StrainVector[i] -= (delta * 2.0 + 1.0);
    }

    rVariables.StressVector = G * rVariables.StrainVector;
    for (unsigned int i = 0; i < 3 ; i++)
    {
        rVariables.StressVector[i] += delta * 3.0 * K;
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateHyperelasticNeoHookeanStress(GeneralVariables& rVariables)
{
    KRATOS_TRY;

    double E = GetProperties()[YOUNG_MODULUS];
    double NU = GetProperties()[POISSON_RATIO];

    double LameLambda = E * NU/((1.0 + NU) * (1.0 - 2.0 * NU));
    double LameMu = E /(2.0 * (1.0 +       NU));

    // Assemble matrix C_bar
    Matrix C_bar;
    C_bar.resize(3, 3, false);
    StructuralMechanicsMathUtilities::VectorToTensor(rVariables.C, C_bar);

    Matrix invC_bar;
    invC_bar.resize(3, 3, false);
    double detC_bar;
    MathUtils<double>::InvertMatrix( C_bar, invC_bar, detC_bar);

    Matrix S_mat;
    S_mat.resize(3, 3, false);

    double factor = log(rVariables.detF);

    noalias(S_mat) = LameLambda * factor * invC_bar;
    noalias(S_mat) += LameMu * ( IdentityMatrix(3) - invC_bar);

    StructuralMechanicsMathUtilities::TensorToVector(S_mat, rVariables.StressVector);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateLogStress(GeneralVariables& rVariables)
{
    KRATOS_TRY;

    double E = GetProperties()[YOUNG_MODULUS];
    double NU = GetProperties()[POISSON_RATIO];

    double LameLambda = E * NU/((1.0 + NU) * (1.0 - 2.0 * NU));
    double LameMu = E /(2.0 * (1.0 +       NU));

    Matrix S_mat;
    S_mat.resize(3, 3, false);
    S_mat = StructuralMechanicsMathUtilities::StrainVectorToTensor(rVariables.StrainVector);

    double trace = S_mat(0, 0) + S_mat(1, 1) + S_mat(2, 2);
    noalias(S_mat) = 2.0 * LameMu * S_mat;
    noalias(S_mat) += LameLambda * trace * IdentityMatrix(3);

//    noalias(rVariables.StressVector) = prod(rVariables.ConstitutiveMatrix, rVariables.StrainVector);
//    StructuralMechanicsMathUtilities::VectorToTensor(rVariables.StressVector, S_mat);

//    Matrix invF;
//    invF.resize(3, 3, false);
//    double detF;
//    MathUtils<double>::InvertMatrix( rVariables.F, invF, detF);

//    noalias(mat) = prod(invF, mat);
//    noalias(mat) = prod(mat, trans(invF));

//    StructuralMechanicsMathUtilities::TensorToVector(mat, rVariables.StressVector);

//    noalias(rVariables.StressVector) = prod(rVariables.ConstitutiveMatrix, rVariables.StrainVector);
////    rVariables.StressVector *= rVariables.detF;

    // Declare the different matrix
    Matrix CMatrix;
    CMatrix.resize(3, 3, false);
    Matrix EigenValuesMatrix;
    EigenValuesMatrix.resize(3, 3, false);
    Matrix EigenVectorsMatrix;
    EigenVectorsMatrix.resize(3, 3, false);

    // Assemble matrix C
    StructuralMechanicsMathUtilities::VectorToTensor(rVariables.C, CMatrix);

    // Decompose matrix
    StructuralMechanicsMathUtilities::EigenVectors(CMatrix, EigenVectorsMatrix, EigenValuesMatrix, 1.0e-24, 10);

////    Matrix R;
////    R.resize(3, 3, false);
////    for (unsigned int i = 0; i < 3; i++)
////    {
////        R(i,0) = mvxe[i];
////        R(i,1) = mvye[i];
////        R(i,2) = mvze[i];
////    }
////    noalias(EigenVectorsMatrix) = prod(EigenVectorsMatrix, R);

//    // Define the rotated tensor of Kirchoff
//    Matrix S_mat;
//    S_mat.resize(3, 3,false);
//    StructuralMechanicsMathUtilities::VectorToTensor(rVariables.StressVector, S_mat);

    noalias(S_mat) = prod(trans(EigenVectorsMatrix), S_mat);
    noalias(S_mat) = prod(S_mat, EigenVectorsMatrix);

    for (unsigned int i = 0; i < 3; i++)
    {
        for (unsigned int j = i; j < 3 ; j++)
        {
            if (j == i)
            {
                 S_mat(i, j) /= EigenValuesMatrix(i, i);
            }
            else
            {
                if (std::abs(EigenValuesMatrix(i, i) - EigenValuesMatrix(j, j)) > 1.0e-30) // Avoid division by 0
                {
                    double fact = 2.0 * log(sqrt(EigenValuesMatrix(i, i))/sqrt(EigenValuesMatrix(j, j)));
                    fact /= (EigenValuesMatrix(i, i) - EigenValuesMatrix(j, j));
                    S_mat(i, j) *= fact;
                    S_mat(j, i) *= fact;
                }
            }
        }
    }

//    for (unsigned int i = 0; i < 3; i++)
//    {
//        for (unsigned int j = i; j < 3 ; j++)
//        {
//            if (j == i)
//            {
//                 S_mat(i, j) /= pow(EigenValuesMatrix(i, i), 2.0);
//            }
//            else
//            {
//                if (std::abs(EigenValuesMatrix(i, i) - EigenValuesMatrix(j, j)) > 1.0e-30) // Avoid division by 0
//                {
//                    double fact = 2.0 * log(EigenValuesMatrix(i, i)/EigenValuesMatrix(j, j));
//                    fact /= (pow(EigenValuesMatrix(i, i), 2.0) - pow(EigenValuesMatrix(j, j), 2.0));
//                    S_mat(i, j) *= fact;
//                    S_mat(j, i) *= fact;
//                }
//            }
//        }
//    }

    noalias(S_mat) = prod(EigenVectorsMatrix, S_mat);
    noalias(S_mat) = prod(S_mat, trans(EigenVectorsMatrix));

    StructuralMechanicsMathUtilities::TensorToVector(S_mat, rVariables.StressVector);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::LinearConstitutiveMatrix(GeneralVariables& rVariables)
{
    KRATOS_TRY;

    double E = GetProperties()[YOUNG_MODULUS];
    double NU = GetProperties()[POISSON_RATIO];

    double K = E /(3.0 * (1.0 - 2.0 * NU));
    double G = E /(2.0 * (1.0 +       NU));

    rVariables.ConstitutiveMatrix = ZeroMatrix(6, 6);

    double ddiag = K + 4.0 * G /3.0;
    double dndia = K - 2.0 * G /3.0;

    rVariables.ConstitutiveMatrix(0, 0) = ddiag;
    rVariables.ConstitutiveMatrix(0, 1) = dndia;
    rVariables.ConstitutiveMatrix(0, 2) = dndia;
    rVariables.ConstitutiveMatrix(1, 1) = ddiag;
    rVariables.ConstitutiveMatrix(1, 2) = dndia;
    rVariables.ConstitutiveMatrix(2, 2) = ddiag;
    rVariables.ConstitutiveMatrix(3, 3) = G;
    rVariables.ConstitutiveMatrix(4, 4) = G;
    rVariables.ConstitutiveMatrix(5, 5) = G;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateGreenLagrangeStrain(
        const Vector& rC,
        Vector& rStrainVector
        )
{
    KRATOS_TRY;

    //Green Lagrange Strain Calculation
    if (rStrainVector.size() != 6)
    {
        rStrainVector.resize(6, false);
    }

    rStrainVector[0] = 0.5 * (rC[0] - 1.00); // xx
    rStrainVector[1] = 0.5 * (rC[1] - 1.00); // yy
    rStrainVector[2] = 0.5 * (rC[2] - 1.00); // zz
    rStrainVector[3] = rC[3]; // xy
    rStrainVector[4] = rC[4]; // yz
    rStrainVector[5] = rC[5]; // xz

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateGreenLagrangeStrain(
        const Matrix& rF,
        Vector& rStrainVector
        )
{
    KRATOS_TRY;

    // Right Cauchy-Green Calculation
    Matrix C ( 3, 3 );

    noalias( C ) = prod( trans( rF ), rF );

    // Green Lagrange Strain Calculation
    if ( rStrainVector.size() != 6 ) rStrainVector.resize( 6, false );

    rStrainVector[0] = 0.5 * (C(0, 0) - 1.0); // xx
    rStrainVector[1] = 0.5 * (C(1, 1) - 1.0); // yy
    rStrainVector[2] = 0.5 * (C(2, 2) - 1.0); // zz
    rStrainVector[3] = C(0, 1); // xy
    rStrainVector[4] = C(1, 2); // yz
    rStrainVector[5] = C(0, 2); // xz

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateHenckyStrain(
        const Vector& rC,
        Vector& rStrainVector
        )
{
    KRATOS_TRY;

    // Declare the different matrix
    Matrix CMatrix;
    CMatrix.resize(3, 3, false);
    Matrix EigenValuesMatrix;
    EigenValuesMatrix.resize(3, 3, false);
    Matrix EigenVectorsMatrix;
    EigenVectorsMatrix.resize(3, 3, false);

    // Assemble matrix C
    StructuralMechanicsMathUtilities::VectorToTensor(rC, CMatrix);

    // Decompose matrix
    StructuralMechanicsMathUtilities::EigenVectors(CMatrix, EigenVectorsMatrix, EigenValuesMatrix, 1e-24, 10);

    // Calculate the eigenvalues of the E matrix
    EigenValuesMatrix(0, 0) = 0.5 * log(EigenValuesMatrix(0, 0));
    EigenValuesMatrix(1, 1) = 0.5 * log(EigenValuesMatrix(1, 1));
    EigenValuesMatrix(2, 2) = 0.5 * log(EigenValuesMatrix(2, 2));

    // Calculate E matrix
    boost::numeric::ublas::bounded_matrix<double, 3, 3 > EMatrix;
    noalias(EMatrix) = prod(trans(EigenVectorsMatrix), EigenValuesMatrix);
    noalias(EMatrix) = prod(EMatrix, EigenVectorsMatrix);

    // Hencky Strain Calculation
    if (rStrainVector.size() != 6)
    {
        rStrainVector.resize(6, false);
    }

    rStrainVector[0] = EMatrix(0, 0); // xx
    rStrainVector[1] = EMatrix(1, 1); // yy
    rStrainVector[2] = EMatrix(2, 2); // zz
    rStrainVector[3] = 2.0 * EMatrix(0, 1); // xy
    rStrainVector[4] = 2.0 * EMatrix(1, 2); // yz
    rStrainVector[5] = 2.0 * EMatrix(0, 2); // xz

    KRATOS_CATCH( "" );
}

//************************************************************************************
//************************************************************************************

void SprismElement3D6N::CalculateAlmansiStrain(
        const Matrix& rF,
        Vector& rStrainVector
        )
{
    KRATOS_TRY;

    // Tensor Cauchy-Green Calculation
    Matrix TensorCauchyGreen = prod(rF, trans(rF));

    // Calculating the inverse of the jacobian
    Matrix InverseTensorCauchyGreen (3, 3);
    double det_b = 0.0;
    MathUtils<double>::InvertMatrix(TensorCauchyGreen, InverseTensorCauchyGreen, det_b);

    // Almansi Strain Calculation
    if ( rStrainVector.size() != 6)
    {
        rStrainVector.resize(6, false);
    }

    rStrainVector[0] = 0.5 * (1.00 - InverseTensorCauchyGreen(0, 0));
    rStrainVector[1] = 0.5 * (1.00 - InverseTensorCauchyGreen(1, 1));
    rStrainVector[2] = 0.5 * (1.00 - InverseTensorCauchyGreen(2, 2));
    rStrainVector[3] = - InverseTensorCauchyGreen(0, 1); // xy
    rStrainVector[4] = - InverseTensorCauchyGreen(1, 2); // yz
    rStrainVector[5] = - InverseTensorCauchyGreen(0, 2); // xz

    KRATOS_CATCH( "" );
}

/************************* CALCULATE VOLUME ACCELERATION ***************************/
/***********************************************************************************/

Vector& SprismElement3D6N::CalculateVolumeForce(
        Vector& rVolumeForce,
        GeneralVariables& rVariables
        )
{
    KRATOS_TRY;

    /* Reading integration points */
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    rVolumeForce = ZeroVector(3);

    array_1d<double,3> accel = ZeroVector(3);

    if( GetGeometry()[0].SolutionStepsDataHas(VOLUME_ACCELERATION) )
    {
        accel = GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    }

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        double IntegrationWeight = integration_points[PointNumber].Weight() * rVariables.detJ;
        rVolumeForce += IntegrationWeight * accel;
    }

    rVolumeForce *= GetProperties()[DENSITY];

    return rVolumeForce;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    int IntMethod = int(mThisIntegrationMethod);
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.save("FinalizedStep",mFinalizedStep);
    rSerializer.save("id_vec",mid_vec);
    rSerializer.save("mTotalDomainInitialSize",mTotalDomainInitialSize);
    rSerializer.save("InvJ0",mInvJ0);
    rSerializer.save("DetJ0",mDetJ0);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.load("FinalizedStep",mFinalizedStep);
    rSerializer.load("id_vec",mid_vec);
    rSerializer.load("mTotalDomainInitialSize",mTotalDomainInitialSize);
    rSerializer.load("InvJ0",mInvJ0);
    rSerializer.load("DetJ0",mDetJ0);
}

} // Namespace Kratos.
