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

// Project includes
#include "custom_elements/SprismElement3D6N.hpp"

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
KRATOS_CREATE_LOCAL_FLAG( SprismElement3D6N, EAS_IMPLICIT_EXPLICIT,              4 ); // True means implicit // TODO: change this using templates!!!
KRATOS_CREATE_LOCAL_FLAG( SprismElement3D6N, TOTAL_UPDATED_LAGRANGIAN,           5 ); // True means total lagrangian // TODO: change this using templates!!!
KRATOS_CREATE_LOCAL_FLAG( SprismElement3D6N, QUADRATIC_ELEMENT,                  6 ); // True means quadratic in-plane behaviour // TODO: Idem

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
            mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_1;
        }
        else if (GetProperties()[NINT_TRANS] == 3)
        {
            mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_2;
        }
        else if (GetProperties()[NINT_TRANS] == 5)
        {
            mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_3;
        }
        else if (GetProperties()[NINT_TRANS] == 7)
        {
            mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_4;
        }
        else if (GetProperties()[NINT_TRANS] == 11)
        {
            mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_5;
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
    ,mAuxMatCont(rOther.mAuxMatCont)
    ,mAuxCont(rOther.mAuxCont)
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
    Element::operator=(rOther);

    mThisIntegrationMethod = rOther.mThisIntegrationMethod;

    mConstitutiveLawVector.clear();
    mConstitutiveLawVector.resize( rOther.mConstitutiveLawVector.size() );

    mAuxMatCont.clear();
    mAuxMatCont.resize( rOther.mAuxMatCont.size());

    for(unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
    {
        mConstitutiveLawVector[i] = rOther.mConstitutiveLawVector[i];
        mAuxMatCont[i]=rOther.mAuxMatCont[i];
    }

    mTotalDomainInitialSize = rOther.mTotalDomainInitialSize;
    mAuxCont = rOther.mAuxCont;

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
    return Kratos::make_shared<SprismElement3D6N>(NewId, GetGeometry().Create(ThisNodes), pProperties);
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
        NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());
    }

    if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
    {
        KRATOS_ERROR << "Constitutive law not has the correct size " << NewElement.mConstitutiveLawVector.size() << std::endl;
    }
    
    for(unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
    {
      NewElement.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
    }

    //-----------//

    if ( NewElement.mAuxMatCont.size() != mAuxMatCont.size() )
    {
        NewElement.mAuxMatCont.resize(mAuxMatCont.size());
    }

    for(unsigned int i = 0; i < mAuxMatCont.size(); i++)
    {
        NewElement.mAuxMatCont[i] = mAuxMatCont[i];
    }

    NewElement.mTotalDomainInitialSize = mTotalDomainInitialSize;
    NewElement.mAuxCont = mAuxCont;

    NewElement.mid_vec = mid_vec;

    return Kratos::make_shared<SprismElement3D6N>(NewElement);
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

    WeakPointerVector< Node < 3 > >& NeighbourNodes = this->GetValue(NEIGHBOUR_NODES);

    const unsigned int NumberOfNodes = GetGeometry().size() + NumberOfActiveNeighbours(NeighbourNodes);
    const unsigned int dim = NumberOfNodes * 3;

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
        if (HasNeighbour(i, NeighbourNodes[i]))
        {
            rResult[index]     = NeighbourNodes[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = NeighbourNodes[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = NeighbourNodes[i].GetDof(DISPLACEMENT_Z).EquationId();
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

    WeakPointerVector< Node < 3 > >& NeighbourNodes = this->GetValue(NEIGHBOUR_NODES);
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
        if (HasNeighbour(i, NeighbourNodes[i]))
        {
            rElementalDofList.push_back(NeighbourNodes[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(NeighbourNodes[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(NeighbourNodes[i].pGetDof(DISPLACEMENT_Z));
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
    WeakPointerVector< Node < 3 > >& NeighbourNodes = this->GetValue(NEIGHBOUR_NODES);
    const unsigned int NumberOfNodes = GetGeometry().size() + NumberOfActiveNeighbours(NeighbourNodes);

    const unsigned int MatSize = NumberOfNodes * 3;
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
        if (HasNeighbour(i, NeighbourNodes[i]))
        {
            const array_1d<double, 3 > & disp = NeighbourNodes[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
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
    WeakPointerVector< Node < 3 > >& NeighbourNodes = this->GetValue(NEIGHBOUR_NODES);
    const unsigned int NumberOfNodes = GetGeometry().size() + NumberOfActiveNeighbours(NeighbourNodes);

    const unsigned int MatSize = NumberOfNodes * 3;
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
        if (HasNeighbour(i, NeighbourNodes[i]))
        {
            const array_1d<double, 3 > & vel = NeighbourNodes[i].FastGetSolutionStepValue(VELOCITY, Step);
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
    WeakPointerVector< Node < 3 > >& NeighbourNodes = this->GetValue(NEIGHBOUR_NODES);
    const unsigned int NumberOfNodes = GetGeometry().size() + NumberOfActiveNeighbours(NeighbourNodes);

    const unsigned int MatSize = NumberOfNodes * 3;
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
        if (HasNeighbour(i, NeighbourNodes[i]))
        {
            const array_1d<double, 3 > & acc = NeighbourNodes[i].FastGetSolutionStepValue(ACCELERATION, Step);
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

    WeakPointerVector< Node < 3 > >& NeighbourNodes = this->GetValue(NEIGHBOUR_NODES);
        
    const double Density = GetProperties()[DENSITY];

    const unsigned int NumberOfNodes = GetGeometry().size() + NumberOfActiveNeighbours(NeighbourNodes);
    const unsigned int MatSize = NumberOfNodes * 3;

    if (rMassMatrix.size1() != MatSize)
    {
        rMassMatrix.resize(MatSize, MatSize, false);
    }
    
    noalias(rMassMatrix) = ZeroMatrix(MatSize, MatSize);
    
    const double Volume = GetGeometry().Volume();
    const double TotalMass = Volume * Density;

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

    WeakPointerVector< Node < 3 > >& NeighbourNodes = this->GetValue(NEIGHBOUR_NODES);

    // 0.-Initialize the DampingMatrix:
    const unsigned int NumberOfNodes = GetGeometry().size() + NumberOfActiveNeighbours(NeighbourNodes);
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the LHS
    const unsigned int MatSize = NumberOfNodes * dimension;

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

    WeakPointerVector< Node < 3 > >& NeighbourNodes = this->GetValue(NEIGHBOUR_NODES);

    // 0.-Initialize the DampingMatrix:
    const unsigned int NumberOfNodes = GetGeometry().size() + NumberOfActiveNeighbours(NeighbourNodes);
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the LHS
    const unsigned int MatSize = NumberOfNodes * dimension;

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

void SprismElement3D6N::CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        )
{
    KRATOS_TRY;

    const unsigned int& IntegrationPointNumber = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != IntegrationPointNumber )
    {
        rOutput.resize( IntegrationPointNumber, false );
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

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        double& AlphaEAS = Element::GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives CartDeriv;
        this->CalculateCartesianDerivatives(CartDeriv);

        /* Calculate common components (B, C) */
        CommonComponents CC;
        CC.clear();
        this->CalculateCommonComponents(CC, CartDeriv);

        // Reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < IntegrationPoints.size(); PointNumber++ )
        {
            const double ZetaGauss = 2.0 * IntegrationPoints[PointNumber].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(Variables, CC, PointNumber, AlphaEAS, ZetaGauss);

            // To take in account previous step writing
            if( mFinalizedStep )
            {
                this->GetHistoricalVariables(Variables,PointNumber);
            }

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(Variables,Values,PointNumber);

            // Call the constitutive law to update material variables
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy (Values);
            
            Matrix StressTensor  = MathUtils<double>::StressVectorToTensor(Variables.StressVector); //reduced dimension stress tensor


            // In general coordinates:
            double SigmaEquivalent =  (0.5)*((StressTensor(0,0)-StressTensor(1,1))*((StressTensor(0,0)-StressTensor(1,1)))+
                                            (StressTensor(1,1)-StressTensor(2,2))*((StressTensor(1,1)-StressTensor(2,2)))+
                                            (StressTensor(2,2)-StressTensor(0,0))*((StressTensor(2,2)-StressTensor(0,0)))+
                                            6*(StressTensor(0,1)*StressTensor(1,0)+StressTensor(1,2)*StressTensor(2,1)+StressTensor(2,0)*StressTensor(0,2)));

            if( SigmaEquivalent < 0 )
            {
                SigmaEquivalent = 0;
            }

            SigmaEquivalent = std::sqrt(SigmaEquivalent);

            rOutput[PointNumber] =  SigmaEquivalent;
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

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY, true);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        double& AlphaEAS = Element::GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives CartDeriv;
        this->CalculateCartesianDerivatives(CartDeriv);

        /* Calculate common components (B, C) */
        CommonComponents CC;
        CC.clear();
        this->CalculateCommonComponents(CC, CartDeriv);

        // Reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < IntegrationPoints.size(); PointNumber++ )
        {
            const double ZetaGauss = 2.0 * IntegrationPoints[PointNumber].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(Variables, CC,PointNumber, AlphaEAS, ZetaGauss);

            // To take in account previous step writing
            if( mFinalizedStep )
            {
                this->GetHistoricalVariables(Variables,PointNumber);
            }

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(Variables,Values,PointNumber);

            // Call the constitutive law to update material variables
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy (Values);
            
            Matrix StressTensor  = MathUtils<double>::StressVectorToTensor(Variables.StressVector); //reduced dimension stress tensor

            double StressNorm =  ((StressTensor(0,0)*StressTensor(0,0))+(StressTensor(1,1)*StressTensor(1,1))+(StressTensor(2,2)*StressTensor(2,2))+
                                (StressTensor(0,1)*StressTensor(0,1))+(StressTensor(0,2)*StressTensor(0,2))+(StressTensor(1,2)*StressTensor(1,2))+
                                (StressTensor(1,0)*StressTensor(1,0))+(StressTensor(2,0)*StressTensor(2,0))+(StressTensor(2,1)*StressTensor(2,1)));

            StressNorm = sqrt(StressNorm);

            rOutput[PointNumber] =  StressNorm;
        }
    }
    else if ( rVariable == STRAIN_ENERGY )
    {
        // Create and initialize element variables:
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY, true);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        double& AlphaEAS = Element::GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives CartDeriv;
        this->CalculateCartesianDerivatives(CartDeriv);

        /* Calculate common components (B, C) */
        CommonComponents CC;
        CC.clear();
        this->CalculateCommonComponents(CC, CartDeriv);

        // Reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < IntegrationPoints.size(); PointNumber++ )
        {
            const double ZetaGauss = 2.0 * IntegrationPoints[PointNumber].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(Variables, CC,PointNumber, AlphaEAS, ZetaGauss);

            // To take in account previous step writing
            if( mFinalizedStep )
            {
                this->GetHistoricalVariables(Variables,PointNumber);
            }

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(Variables,Values,PointNumber);

            double StrainEnergy = 0.0;

            // Compute stresses and constitutive parameters
            if ( mELementalFlags.Is(SprismElement3D6N::TOTAL_UPDATED_LAGRANGIAN) == true )
            {
                mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseKirchhoff(Values);
            }
            else
            {
                mConstitutiveLawVector[PointNumber]->CalculateMaterialResponsePK2(Values);
            }
            mConstitutiveLawVector[PointNumber]->GetValue(STRAIN_ENERGY, StrainEnergy);

            rOutput[PointNumber] = Variables.detJ * IntegrationPoints[PointNumber].Weight() * StrainEnergy;  // 1/2 * sigma * epsilon
        }
    }
    else
    {
        for ( unsigned int ii = 0; ii < IntegrationPointNumber; ii++ )
        {
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rOutput[ii] );
        }
    }

    if ( rOutput.size() != 6 )
    {
        std::vector<double> rOutput_aux;
        rOutput_aux = rOutput;

        rOutput.resize( 6, false );
        Matrix interpol = StructuralMechanicsMathUtilities::InterpolPrismGiD(IntegrationPointNumber);

        for (unsigned int iii = 0; iii < 6; iii++)
        {
            rOutput[iii] = 0.0;

            for (unsigned int Gauss_Point = 0; Gauss_Point < IntegrationPointNumber; Gauss_Point++)
            {
                rOutput[iii] += interpol(Gauss_Point, iii) * rOutput_aux[Gauss_Point];
            }
        }
    }

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

    const unsigned int& IntegrationPointNumber = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != IntegrationPointNumber )
    {
        rOutput.resize( IntegrationPointNumber );
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

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        double& AlphaEAS = Element::GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives CartDeriv;
        this->CalculateCartesianDerivatives(CartDeriv);

        /* Calculate common components (B, C) */
        CommonComponents CC;
        CC.clear();
        this->CalculateCommonComponents(CC, CartDeriv);

        // Reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < IntegrationPoints.size(); PointNumber++ )
        {
            const double ZetaGauss = 2.0 * IntegrationPoints[PointNumber].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(Variables, CC, PointNumber, AlphaEAS, ZetaGauss);

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
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        double& AlphaEAS = Element::GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives CartDeriv;
        this->CalculateCartesianDerivatives(CartDeriv);

        /* Calculate common components (B, C) */
        CommonComponents CC;
        CC.clear();
        this->CalculateCommonComponents(CC, CartDeriv);

        // Reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < IntegrationPoints.size(); PointNumber++ )
        {
            const double ZetaGauss = 2.0 * IntegrationPoints[PointNumber].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(Variables, CC, PointNumber, AlphaEAS, ZetaGauss);

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
        Matrix interpol = StructuralMechanicsMathUtilities::InterpolPrismGiD(IntegrationPointNumber);

        for (unsigned int iii = 0; iii < 6; iii++)
        {
            rOutput[iii] = ZeroVector(rOutput[0].size());

            for (unsigned int Gauss_Point = 0; Gauss_Point < IntegrationPointNumber; Gauss_Point++)
            {
                rOutput[iii] += interpol(Gauss_Point, iii) * rOutput_aux[Gauss_Point];
            }
        }
    }

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

    const unsigned int& IntegrationPointNumber = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != IntegrationPointNumber )
    {
        rOutput.resize( IntegrationPointNumber );
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
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        double& AlphaEAS = Element::GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives CartDeriv;
        this->CalculateCartesianDerivatives(CartDeriv);

        /* Calculate common components (B, C) */
        CommonComponents CC;
        CC.clear();
        this->CalculateCommonComponents(CC, CartDeriv);

        // Reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < IntegrationPoints.size(); PointNumber++ )
        {
            const double& ZetaGauss = 2.0 * IntegrationPoints[PointNumber].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(Variables, CC, PointNumber, AlphaEAS, ZetaGauss);

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
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        double& AlphaEAS = Element::GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives CartDeriv;
        this->CalculateCartesianDerivatives(CartDeriv);

        /* Calculate common components (B, C) */
        CommonComponents CC;
        CC.clear();
        this->CalculateCommonComponents(CC, CartDeriv);

        // Reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < IntegrationPoints.size(); PointNumber++ )
        {
            const double & ZetaGauss = 2.0 * IntegrationPoints[PointNumber].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(Variables, CC, PointNumber, AlphaEAS, ZetaGauss);

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
        Matrix interpol = StructuralMechanicsMathUtilities::InterpolPrismGiD(IntegrationPointNumber);

        for (unsigned int iii = 0; iii < 6; iii++)
        {
            rOutput[iii] = ZeroMatrix(rOutput[0].size1(), rOutput[0].size2());

            for (unsigned int Gauss_Point = 0; Gauss_Point < IntegrationPointNumber; Gauss_Point++)
            {
                rOutput[iii] += interpol(Gauss_Point, iii) * rOutput_aux[Gauss_Point];
            }
        }
    }

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
        if (rVariable == DETERMINANT_F)
        {
            mAuxCont[PointNumber] = rValues[PointNumber];
        }

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
            mConstitutiveLawVector.resize(rValues.size());
            if( mConstitutiveLawVector.size() != GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod ) )
            {
                KRATOS_ERROR << "Constitutive law not has the correct size " << mConstitutiveLawVector.size() << std::endl;
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
        const unsigned int& IntegrationPointNumber = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );
        if ( rValues.size() != IntegrationPointNumber )
        {
            rValues.resize( IntegrationPointNumber, false );
        }
        for ( unsigned int ii = 0; ii < IntegrationPointNumber; ii++ )
        {
            if (rVariable == DETERMINANT_F)
            {
                rValues[ii] = mAuxCont[ii];
            }
            else
            {
                rValues[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rValues[ii] );
            }
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
    const unsigned int& IntegrationPointNumber = mConstitutiveLawVector.size();

    if ( rValues.size() != IntegrationPointNumber )
    {
        rValues.resize( IntegrationPointNumber );
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
        for ( unsigned int PointNumber = 0;  PointNumber < IntegrationPointNumber; PointNumber++ )
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
    const unsigned int& IntegrationPointNumber = mConstitutiveLawVector.size();

    if ( rValues.size() != IntegrationPointNumber )
    {
        rValues.resize( IntegrationPointNumber );
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
        for ( unsigned int PointNumber = 0;  PointNumber < IntegrationPointNumber; PointNumber++ )
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
    if(rVariable == CONSTITUTIVE_LAW)
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
        {
            rValues.resize(mConstitutiveLawVector.size());
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
    WeakPointerVector< Element >& NeighbourElements = this->GetValue(NEIGHBOUR_ELEMENTS);
    if (NeighbourElements.size() == 0)
    {
        KRATOS_ERROR << "The neighbour elements are not calculated" << std::endl;
    }

    // Neighbour nodes
    WeakPointerVector< Node < 3 > >& NeighbourNodes = this->GetValue(NEIGHBOUR_NODES);
    if (NeighbourNodes.size() == 0)
    {
        KRATOS_ERROR << "The neighbour nodes are not calculated" << std::endl;
    }

    /* Verify compatibility with the constitutive law */
    ConstitutiveLaw::Features LawFeatures;
    this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetLawFeatures(LawFeatures);

    bool CorrectStrainMeasure = false;
    for(unsigned int i = 0; i < LawFeatures.mStrainMeasures.size(); i++)
    {
        if(LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Deformation_Gradient)
        {
            CorrectStrainMeasure = true;
        }
    }

    if(CorrectStrainMeasure == false)
    {
        KRATOS_ERROR << "Constitutive law is not compatible with the element type SprismElement3D6N" << std::endl;
    }

    // Verify that nodal variables are correctly initialized
    if (DISPLACEMENT.Key() == 0)
    {
        KRATOS_ERROR << "DISPLACEMENT has Key zero! (check if the application is correctly registered"<< std::endl;
    }

    if (VELOCITY.Key() == 0)
    {
        KRATOS_ERROR << "VELOCITY has Key zero! (check if the application is correctly registered"<< std::endl;
    }

    if (ACCELERATION.Key() == 0)
    {
        KRATOS_ERROR << "ACCELERATION has Key zero! (check if the application is correctly registered"<< std::endl;
    }

    if (DENSITY.Key() == 0)
    {
        KRATOS_ERROR << "DENSITY has Key zero! (check if the application is correctly registered"<< std::endl;
    }

    if (VOLUME_ACCELERATION.Key() == 0)
    {
        KRATOS_ERROR << "VOLUME_ACCELERATION has Key zero! (check if the application is correctly registered"<< std::endl;
    }

    // Verify that elemental variables are correctly initialized
    if ( VON_MISES_STRESS.Key() == 0 )
    {
        KRATOS_ERROR << "VON_MISES_STRESS has Key zero! (check if the application is correctly registered"<< std::endl;
    }

    if ( NORM_ISOCHORIC_STRESS.Key() == 0 )
    {
        KRATOS_ERROR << "NORM_ISOCHORIC_STRESS has Key zero! (check if the application is correctly registered"<< std::endl;
    }

    if ( CAUCHY_STRESS_TENSOR.Key() == 0 )
    {
        KRATOS_ERROR << "CAUCHY_STRESS_TENSOR has Key zero! (check if the application is correctly registered"<< std::endl;
    }

    if ( CAUCHY_STRESS_VECTOR.Key() == 0 )
    {
        KRATOS_ERROR << "CAUCHY_STRESS_VECTOR has Key zero! (check if the application is correctly registered"<< std::endl;
    }

    if ( PK2_STRESS_TENSOR.Key() == 0 )
    {
        KRATOS_ERROR << "PK2_STRESS_TENSOR has Key zero! (check if the application is correctly registered"<< std::endl;
    }

    if ( PK2_STRESS_VECTOR.Key() == 0 )
    {
        KRATOS_ERROR << "PK2_STRESS_VECTOR has Key zero! (check if the application is correctly registered" << std::endl;
    }

    if ( GREEN_LAGRANGE_STRAIN_TENSOR.Key() == 0 )
    {
        KRATOS_ERROR << "GREEN_LAGRANGE_STRAIN_TENSOR has Key zero! (check if the application is correctly registered" << std::endl;
    }

    if ( GREEN_LAGRANGE_STRAIN_VECTOR.Key() == 0 )
    {
        KRATOS_ERROR << "GREEN_LAGRANGE_STRAIN_VECTOR has Key zero! (check if the application is correctly registered" << std::endl;
    }

    if ( ALMANSI_STRAIN_TENSOR.Key() == 0 )
    {
        KRATOS_ERROR << "ALMANSI_STRAIN_TENSOR has Key zero! (check if the application is correctly registered" << std::endl;
    }

    if ( ALMANSI_STRAIN_VECTOR.Key() == 0 )
    {
        KRATOS_ERROR << "ALMANSI_STRAIN_VECTOR has Key zero! (check if the application is correctly registered" << std::endl;
    }

    if ( HENCKY_STRAIN_TENSOR.Key() == 0 )
    {
        KRATOS_ERROR << "HENCKY_STRAIN_TENSOR has Key zero! (check if the application is correctly registered" << std::endl;
    }

    if ( HENCKY_STRAIN_VECTOR.Key() == 0 )
    {
        KRATOS_ERROR << "HENCKY_STRAIN_VECTOR has Key zero! (check if the application is correctly registered" << std::endl;
    }

    if ( CONSTITUTIVE_MATRIX.Key() == 0 )
    {
        KRATOS_ERROR << "CONSTITUTIVE_MATRIX has Key zero! (check if the application is correctly registered" << std::endl;
    }

    if ( DEFORMATION_GRADIENT.Key() == 0 )
    {
        KRATOS_ERROR << "DEFORMATION_GRADIENT has Key zero! (check if the application is correctly registered" << std::endl;
    }

    // Verify that the dofs exist
    for ( unsigned int i = 0; i < this->GetGeometry().size(); i++)
    {
        if ( this->GetGeometry()[i].SolutionStepsDataHas( DISPLACEMENT ) == false)
        {
            KRATOS_ERROR << "Missing variable DISPLACEMENT on node " << this->GetGeometry()[i].Id() << std::endl;
        }

        if ( this->GetGeometry()[i].HasDofFor( DISPLACEMENT_X ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Y ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Z ) == false)
        {
            KRATOS_ERROR << "Missing one of the dofs for the variable DISPLACEMENT on node " << GetGeometry()[i].Id() << std::endl;
        }
    }

    // Verify that the constitutive law exists
    if (this->GetProperties().Has( CONSTITUTIVE_LAW ) == false)
    {
        KRATOS_ERROR << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;
    }

    if (this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() != 6)
    {
        KRATOS_ERROR << "Wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) " << this->Id() << std::endl;
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

void SprismElement3D6N::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
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

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    // Get constitutive law flags:
    Flags &ConstitutiveLawOptions=Values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    double& AlphaEAS = Element::GetValue(ALPHA_EAS);

    /* Calculate the cartesian derivatives */
    CartesianDerivatives CartDeriv;
    this->CalculateCartesianDerivatives(CartDeriv);

    /* Calculate common components (B, C) */
    CommonComponents CC;
    CC.clear();
    this->CalculateCommonComponents(CC, CartDeriv);

    // Reading integration points
    for ( unsigned int PointNumber = 0; PointNumber < IntegrationPoints.size(); PointNumber++ )
    {
        const double& ZetaGauss = 2.0 * IntegrationPoints[PointNumber].Z() - 1.0;

        // Compute element kinematics C, F ...
        this->CalculateKinematics(Variables, CC, PointNumber, AlphaEAS, ZetaGauss);

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

void SprismElement3D6N::Initialize()
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    /* Constitutive Law initialisation */
    if ( mConstitutiveLawVector.size() != IntegrationPoints.size() )
    {
        mConstitutiveLawVector.resize( IntegrationPoints.size() );
    }

    /* Implicit or explicit EAS update */
    if( GetProperties().Has(EAS_IMP) )
    {
        mELementalFlags.Set(SprismElement3D6N::EAS_IMPLICIT_EXPLICIT, GetProperties()[EAS_IMP]);
    }
    else
    {
        mELementalFlags.Set(SprismElement3D6N::EAS_IMPLICIT_EXPLICIT, true);
    }

    /* Total or updated lagrangian */
    if( GetProperties().Has(SPRISM_TL_UL) )
    {
        mELementalFlags.Set(SprismElement3D6N::TOTAL_UPDATED_LAGRANGIAN, GetProperties()[SPRISM_TL_UL]);
    }
    else
    {
        mELementalFlags.Set(SprismElement3D6N::TOTAL_UPDATED_LAGRANGIAN, true);
    }

    /* Quadratic or linear element */
    if( GetProperties().Has(QUAD_ON) )
    {
        mELementalFlags.Set(SprismElement3D6N::QUADRATIC_ELEMENT, GetProperties()[QUAD_ON]);
    }
    else
    {
        mELementalFlags.Set(SprismElement3D6N::QUADRATIC_ELEMENT, true);
    }

    // Resizing the containers
    mAuxMatCont.resize( IntegrationPoints.size() );
    mAuxCont.resize( IntegrationPoints.size(), false );

    if ( mELementalFlags.Is(SprismElement3D6N::TOTAL_UPDATED_LAGRANGIAN) == true ) // Jacobian inverses
    {
        // Compute jacobian inverses and set the domain initial size:
        GeometryType::JacobiansType J0;
        J0 = GetGeometry().Jacobian(J0, mThisIntegrationMethod);
        mTotalDomainInitialSize = 0.0;

        /* Calculating the inverse J0 */
        for ( unsigned int PointNumber = 0; PointNumber < IntegrationPoints.size(); PointNumber++ )
        {
            // Calculating and storing inverse of the jacobian and the parameters needed
            MathUtils<double>::InvertMatrix( J0[PointNumber], mAuxMatCont[PointNumber], mAuxCont[PointNumber] );

            // Getting informations for integration
            double IntegrationWeight = IntegrationPoints[PointNumber].Weight();

            // Calculating the total volume
            mTotalDomainInitialSize += mAuxCont[PointNumber] * IntegrationWeight;
        }
    }
    else // Historic deformation gradient
    {
        mTotalDomainInitialSize = 0.0; // Just initialize, not used in UL

        for ( unsigned int PointNumber = 0; PointNumber < IntegrationPoints.size(); PointNumber++ )
        {
            mAuxCont[PointNumber] = 1.00;
            mAuxMatCont[PointNumber] = IdentityMatrix(3);
        }
    }

    /* Initialize AlphaEAS */
    double& AlphaEAS = Element::GetValue(ALPHA_EAS);
    AlphaEAS = 0.0;

    /* Initialize EAS parameters*/
    mEAS.clear();

    /* Material initialisation */
    InitializeMaterial();

    /* Calculate ID vector*/
    this->CalculateIdVect();

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

    /* Create constitutive law parameters: */
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    /* Set constitutive law flags: */
    Flags &ConstitutiveLawOptions=Values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    /* Reading integration points */
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    /* Getting the alpha parameter of the EAS improvement */
    double& AlphaEAS = Element::GetValue(ALPHA_EAS);

    /* Calculate the RHS */
    if ( rLocalSystem.CalculationFlags.Is(SprismElement3D6N::COMPUTE_RHS_VECTOR) == true ) // Update just if RHS is calculated
    {
        /* Getting the increase of displacements */
        bounded_matrix<double, 36, 1 > delta_disp;

        delta_disp = GetVectorCurrentPosition() - mPreviousCoor; // Calculates the increase of displacements
        mPreviousCoor = GetVectorCurrentPosition(); // Update previous coordinates // Note: Save coordinates in an auxiliar variable

        /* Update alpha EAS */
        if (mEAS.stiff_alpha > 1.0e-12) // Avoid division by zero
        {
            AlphaEAS -= prod(mEAS.H_EAS, delta_disp)(0, 0) / mEAS.stiff_alpha;
        }
    }

    /* Calculate the cartesian derivatives */
    CartesianDerivatives CartDeriv;
    this->CalculateCartesianDerivatives(CartDeriv);

    /* Calculate common components (B, C) */
    CommonComponents CC;
    CC.clear();
    this->CalculateCommonComponents(CC, CartDeriv);

    /* Reset the integrated stress components */
    StressIntegratedComponents IntStress;
    IntStress.clear();

    /* Reset the EAS integrated components */
    mEAS.clear();

    // Reading integration points
    for ( unsigned int PointNumber = 0; PointNumber < IntegrationPoints.size(); PointNumber++ )
    {
        const double ZetaGauss = 2.0 * IntegrationPoints[PointNumber].Z() - 1.0;

        /* Assemble B */
        this->CalculateDeformationMatrix(Variables.B, CC, ZetaGauss, AlphaEAS);

        // Compute element kinematics C, F ...
        this->CalculateKinematics(Variables, CC, PointNumber, AlphaEAS, ZetaGauss);

        // Set general variables to constitutivelaw parameters
        this->SetGeneralVariables(Variables, Values, PointNumber);

        // Compute stresses and constitutive parameters
        mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(Values, Variables.StressMeasure);

        // Calculating weights for integration on the "reference configuration"
        const double IntegrationWeight = IntegrationPoints[PointNumber].Weight() * Variables.detJ;

        /* Integrate in Zeta */
        IntegrateInZeta(Variables, IntStress, AlphaEAS, ZetaGauss, IntegrationWeight);
    }

    /* Auxiliary terms: Allocating the VolumeForce*/
    Vector VolumeForce;

    /* Calculate the RHS */
    if ( rLocalSystem.CalculationFlags.Is(SprismElement3D6N::COMPUTE_RHS_VECTOR) == true ) // Calculation of the vector is required
    {
        /* Volume forces */
        this->CalculateVolumeForce( VolumeForce, Variables );

        /* Contribution to external and internal forces */
        this->CalculateAndAddRHS ( rLocalSystem, Variables, VolumeForce, IntStress, CC, AlphaEAS );
    }

    if ( rLocalSystem.CalculationFlags.Is(SprismElement3D6N::COMPUTE_LHS_MATRIX) == true ) // Calculation of the matrix is required
    {
        /* Contribution to the tangent stiffness matrix */
        this->CalculateAndAddLHS( rLocalSystem, Variables, Values, IntStress, CC, CartDeriv, AlphaEAS );
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
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    MatrixType  LocalLeftHandSideMatrix;
    VectorType  LocalRightHandSideVector;

    // Initialize sizes for the system components:
    this->InitializeSystemMatrices( LocalLeftHandSideMatrix, LocalRightHandSideVector, rLocalSystem.CalculationFlags );

    /* Getting the alpha parameter of the EAS improvement */
    double& AlphaEAS = Element::GetValue(ALPHA_EAS);

    /* Calculate the cartesian derivatives */
    CartesianDerivatives CartDeriv;
    this->CalculateCartesianDerivatives(CartDeriv);

    /* Calculate common components (B, C) */
    CommonComponents CC;
    CC.clear();
    this->CalculateCommonComponents(CC, CartDeriv);

    for ( unsigned int PointNumber = 0; PointNumber < IntegrationPoints.size(); PointNumber++ )
    {
        const double ZetaGauss = 2.0 * IntegrationPoints[PointNumber].Z() - 1.0;

        /* Assemble B */
        this->CalculateDeformationMatrix(Variables.B, CC, ZetaGauss, AlphaEAS);

        // Compute element kinematics C, F ...
        this->CalculateKinematics(Variables, CC, PointNumber, AlphaEAS, ZetaGauss);

        // Calculating weights for integration on the "reference configuration"
        const double IntegrationWeight = IntegrationPoints[PointNumber].Weight() * Variables.detJ;

        if ( rLocalSystem.CalculationFlags.Is(SprismElement3D6N::COMPUTE_LHS_MATRIX) ) // Calculation of the matrix is required
        {
            LocalLeftHandSideMatrix.clear();

            this->CalculateAndAddDynamicLHS ( LocalLeftHandSideMatrix, Variables, IntegrationWeight);

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

    WeakPointerVector< Node < 3 > >& NeighbourNodes = this->GetValue(NEIGHBOUR_NODES);
    const unsigned int NumberNeighbours = NumberOfActiveNeighbours(NeighbourNodes);

    for ( unsigned int i = 0; i < 6; i++ )
    {
        array_1d<double, 3> &CurrentPosition  = GetGeometry()[i].Coordinates();
        array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
        array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
        array_1d<double, 3> PreviousPosition  = CurrentPosition - (CurrentDisplacement-PreviousDisplacement);
        std::cout << " Previous  Position  node[" << GetGeometry()[i].Id() << "]: "<<PreviousPosition << std::endl;
    }

    for ( unsigned int i = 0; i < NumberNeighbours; i++ )
    {
        array_1d<double, 3> &CurrentPosition  = NeighbourNodes[i].Coordinates();
        array_1d<double, 3 > & CurrentDisplacement  = NeighbourNodes[i].FastGetSolutionStepValue(DISPLACEMENT);
        array_1d<double, 3 > & PreviousDisplacement = NeighbourNodes[i].FastGetSolutionStepValue(DISPLACEMENT,1);
        array_1d<double, 3> PreviousPosition  = CurrentPosition - (CurrentDisplacement-PreviousDisplacement);
        std::cout << " Previous  Position  neighbour node[" << NeighbourNodes[i].Id() << "]: "<<PreviousPosition << std::endl;
    }

    for ( unsigned int i = 0; i < 6; i++ )
    {
        array_1d<double, 3> & CurrentPosition  = GetGeometry()[i].Coordinates();
        std::cout << " Current  Position  node[" << GetGeometry()[i].Id()<<"]: " << CurrentPosition << std::endl;
    }

    for ( unsigned int i = 0; i < NumberNeighbours; i++ )
    {
        array_1d<double, 3> & CurrentPosition  = NeighbourNodes[i].Coordinates();
        std::cout << " Current  Position neighbour node[" << NeighbourNodes[i].Id()<<"]: " << CurrentPosition << std::endl;
    }

    for ( unsigned int i = 0; i < 6; i++ )
    {
        array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
        std::cout << " Previous Displacement node[" << GetGeometry()[i].Id() << "]: " << PreviousDisplacement << std::endl;
    }

    for ( unsigned int i = 0; i < NumberNeighbours; i++ )
    {
        array_1d<double, 3 > & PreviousDisplacement = NeighbourNodes[i].FastGetSolutionStepValue(DISPLACEMENT,1);
        std::cout << " Previous Displacement neighbour node[" << NeighbourNodes[i].Id() << "]: " << PreviousDisplacement << std::endl;
    }

    for ( unsigned int i = 0; i < 6; i++ )
    {
        array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
        std::cout << " Current  Displacement  node[" << GetGeometry()[i].Id() << "]: " << CurrentDisplacement << std::endl;
    }

    for ( unsigned int i = 0; i < NumberNeighbours; i++ )
    {
        array_1d<double, 3 > & CurrentDisplacement  = NeighbourNodes[i].FastGetSolutionStepValue(DISPLACEMENT);
        std::cout << " Current  Displacement  node[" << NeighbourNodes[i].Id() << "]: " << CurrentDisplacement << std::endl;
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
    if (neighb.Id() == GetGeometry()[index].Id())
    {
        return false;
    }
    else
    {
        if ( mELementalFlags.Is(SprismElement3D6N::QUADRATIC_ELEMENT) == true )
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

void SprismElement3D6N::GetNodalCoordinates(
        bounded_matrix<double, 12, 3 > & NodesCoord,
        WeakPointerVector< Node < 3 > >& NeighbourNodes,
        const Configuration ThisConfiguration
        )
{
     NodesCoord = ZeroMatrix(12, 3);
     const unsigned int NumberNeighbours = NumberOfActiveNeighbours(NeighbourNodes);

     if (ThisConfiguration == Initial)
     {
         /* Fill the aux matrix of coordinates */
         for (unsigned int i = 0; i < 6; i++)
         {
             NodesCoord(i, 0) = GetGeometry()[i].X0();
             NodesCoord(i, 1) = GetGeometry()[i].Y0();
             NodesCoord(i, 2) = GetGeometry()[i].Z0();
         }

         if (NumberNeighbours == 6) // All the possible neighours
         {
             for (unsigned int i = 0; i < 6; i++)
             {
                 NodesCoord(i + 6, 0) = NeighbourNodes[i].X0();
                 NodesCoord(i + 6, 1) = NeighbourNodes[i].Y0();
                 NodesCoord(i + 6, 2) = NeighbourNodes[i].Z0();
             }
         }
         else
         {
             for (unsigned int i = 0; i < 6; i++)
             {
                 if (HasNeighbour(i, NeighbourNodes[i]))
                 {
                     NodesCoord(i + 6, 0) = NeighbourNodes[i].X0();
                     NodesCoord(i + 6, 1) = NeighbourNodes[i].Y0();
                     NodesCoord(i + 6, 2) = NeighbourNodes[i].Z0();
                 }
                 else
                 {
                     NodesCoord(i + 6, 0) = 0.0;
                     NodesCoord(i + 6, 1) = 0.0;
                     NodesCoord(i + 6, 2) = 0.0;
                 }
             }
         }
     }
     else if (ThisConfiguration == Current)
     {
         /* Fill the aux matrix of coordinates */
         for (unsigned int i = 0; i < 6; i++)
         {
             const array_1d<double, 3> &CurrentPosition  = GetGeometry()[i].Coordinates();
             NodesCoord(i, 0) = CurrentPosition[0];
             NodesCoord(i, 1) = CurrentPosition[1];
             NodesCoord(i, 2) = CurrentPosition[2];
         }

         if (NumberNeighbours == 6) // All the possible neighours
         {
             for (unsigned int i = 0; i < 6; i++)
             {
                 const array_1d<double, 3> &CurrentPosition  = NeighbourNodes[i].Coordinates();
                 NodesCoord(i + 6, 0) = CurrentPosition[0];
                 NodesCoord(i + 6, 1) = CurrentPosition[1];
                 NodesCoord(i + 6, 2) = CurrentPosition[2];
             }
         }
         else
         {
             for (unsigned int i = 0; i < 6; i++)
             {
                 if (HasNeighbour(i, NeighbourNodes[i]))
                 {
                     const array_1d<double, 3> &CurrentPosition  = NeighbourNodes[i].Coordinates();
                     NodesCoord(i + 6, 0) = CurrentPosition[0];
                     NodesCoord(i + 6, 1) = CurrentPosition[1];
                     NodesCoord(i + 6, 2) = CurrentPosition[2];
                 }
                 else
                 {
                     NodesCoord(i + 6, 0) = 0.0;
                     NodesCoord(i + 6, 1) = 0.0;
                     NodesCoord(i + 6, 2) = 0.0;
                 }
             }
         }
     }
     else
     {
         std::string Config = (ThisConfiguration == Initial) ? "Initial" : "Current";
         KRATOS_ERROR << " The configuration is not possible, the posibilities are Current and Initial: " << Config << std::endl;
     }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateCartesianDerivatives(CartesianDerivatives& CartDeriv)
{
    bounded_matrix<double, 12, 3 > NodesCoord; // Coordinates of the nodes
    WeakPointerVector< Node < 3 > >& NeighbourNodes = this->GetValue(NEIGHBOUR_NODES);
    if ( mELementalFlags.Is(SprismElement3D6N::TOTAL_UPDATED_LAGRANGIAN) == true )
    {
        this->GetNodalCoordinates(NodesCoord, NeighbourNodes, Initial);
    }
    else
    {
        this->GetNodalCoordinates(NodesCoord, NeighbourNodes, Current);
    }

    /* Calculate local system of coordinates of the element */
    // 0- If X is the prefered normal vector
    // 1- If Y is the prefered normal vector
    // 2- If Z is the prefered normal vector

    double AngRot = 0.0; // TODO: Change to consider multiple plies
    if( GetProperties().Has(ANG_ROT) )
    {
        AngRot = GetProperties()[ANG_ROT];
    }

    this->CalculateLocalCoordinateSystem(2, AngRot);

    //******************************** CENTRAL POINT ******************************
    // Calculate cartesian derivatives
    bounded_matrix<double, 2, 4 > CartesianDerivativesCenterLower;
    bounded_matrix<double, 2, 4 > CartesianDerivativesCenterUpper;

    // Lower face
    CalculateCartesianDerOnCenter_plane(0, NodesCoord, CartesianDerivativesCenterLower);
    // Upperr face
    CalculateCartesianDerOnCenter_plane(3, NodesCoord, CartesianDerivativesCenterUpper );

    /* Transversal derivative */
    CalculateCartesianDerOnCenter_trans(CartDeriv, NodesCoord, 0); // Center
    CalculateCartesianDerOnCenter_trans(CartDeriv, NodesCoord, 1); // Lower part
    CalculateCartesianDerOnCenter_trans(CartDeriv, NodesCoord, 2); // Upper part

    //******************************** GAUSS POINTS *******************************

    /* Transversal derivative */
    CalculateCartesianDerOnGauss_trans(NodesCoord, CartDeriv.TransversalCartesianDerivativesGauss1, 0.5, 0.5, - 1.0);
    CalculateCartesianDerOnGauss_trans(NodesCoord, CartDeriv.TransversalCartesianDerivativesGauss4, 0.5, 0.5,   1.0);

    /* In-plane derivative */
    if (HasNeighbour(0, NeighbourNodes[0])) // Assuming that if the upper element has neighbours the lower has too
    {
        CalculateCartesianDerOnGauss_plane(0, 0, NodesCoord, CartDeriv.InPlaneCartesianDerivativesGauss1);
        CalculateCartesianDerOnGauss_plane(0, 3, NodesCoord, CartDeriv.InPlaneCartesianDerivativesGauss4);
    }
    else
    {
        noalias(CartDeriv.InPlaneCartesianDerivativesGauss1) = CartesianDerivativesCenterLower;
        noalias(CartDeriv.InPlaneCartesianDerivativesGauss4) = CartesianDerivativesCenterUpper;
    }

    /* Transversal derivative */
    CalculateCartesianDerOnGauss_trans(NodesCoord, CartDeriv.TransversalCartesianDerivativesGauss2, 0.0, 0.5, - 1.0);
    CalculateCartesianDerOnGauss_trans(NodesCoord, CartDeriv.TransversalCartesianDerivativesGauss5, 0.0, 0.5,   1.0);

    /* In-plane derivative */
    if (HasNeighbour(1, NeighbourNodes[1])) //Idem
    {
        CalculateCartesianDerOnGauss_plane(1, 0, NodesCoord, CartDeriv.InPlaneCartesianDerivativesGauss2);
        CalculateCartesianDerOnGauss_plane(1, 3, NodesCoord, CartDeriv.InPlaneCartesianDerivativesGauss5);
    }
    else
    {
        noalias(CartDeriv.InPlaneCartesianDerivativesGauss2) = CartesianDerivativesCenterLower;
        noalias(CartDeriv.InPlaneCartesianDerivativesGauss5) = CartesianDerivativesCenterUpper;
    }

    /* Transversal derivative */
    CalculateCartesianDerOnGauss_trans(NodesCoord, CartDeriv.TransversalCartesianDerivativesGauss3, 0.5, 0.0, - 1.0);
    CalculateCartesianDerOnGauss_trans(NodesCoord, CartDeriv.TransversalCartesianDerivativesGauss6, 0.5, 0.0,   1.0);

    /* In-plane derivative */
    if (HasNeighbour(2, NeighbourNodes[2])) // Idem
    {
        CalculateCartesianDerOnGauss_plane(2, 0, NodesCoord, CartDeriv.InPlaneCartesianDerivativesGauss3);
        CalculateCartesianDerOnGauss_plane(2, 3, NodesCoord, CartDeriv.InPlaneCartesianDerivativesGauss6);
    }
    else
    {
        noalias(CartDeriv.InPlaneCartesianDerivativesGauss3) = CartesianDerivativesCenterLower;
        noalias(CartDeriv.InPlaneCartesianDerivativesGauss6) = CartesianDerivativesCenterUpper;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateCommonComponents(
        CommonComponents& CC,
        const CartesianDerivatives& CartDeriv
        )
{
    KRATOS_TRY;

    bounded_matrix<double, 12, 3 > NodesCoord; // Coordinates of the nodes
    WeakPointerVector< Node < 3 > >& NeighbourNodes = this->GetValue(NEIGHBOUR_NODES);
    this->GetNodalCoordinates(NodesCoord, NeighbourNodes, Current);

    /* Declare deformation Gradient F components */
    // In plane components
    bounded_matrix<double, 3, 2 > InPlaneGradientFGauss;
    // Transversal components
    array_1d<double, 3 > TransverseGradientF0, TransverseGradientF1, TransverseGradientF2;

    //*****************************************************************************

    /* COMPUTATION OF B TANGENTS MATRIX */

    /* MEMBRANE CONTRIBUTION */
    /* Calculating the membrane strain-displacement matrix */
    // Lower face

    // Gauss point 1
    CalculateInPlaneGradientFGauss(InPlaneGradientFGauss, CartDeriv.InPlaneCartesianDerivativesGauss1, NodesCoord, 0, 0);
    CalculateAndAdd_B_Membrane(CC.B_membrane_lower, CC.C_membrane_lower, CartDeriv.InPlaneCartesianDerivativesGauss1, InPlaneGradientFGauss, 0);

    // Gauss point 2
    CalculateInPlaneGradientFGauss(InPlaneGradientFGauss, CartDeriv.InPlaneCartesianDerivativesGauss2, NodesCoord, 1, 0);
    CalculateAndAdd_B_Membrane(CC.B_membrane_lower, CC.C_membrane_lower, CartDeriv.InPlaneCartesianDerivativesGauss2, InPlaneGradientFGauss, 1);

    // Gauss point 3
    CalculateInPlaneGradientFGauss(InPlaneGradientFGauss, CartDeriv.InPlaneCartesianDerivativesGauss3, NodesCoord, 2, 0);
    CalculateAndAdd_B_Membrane(CC.B_membrane_lower, CC.C_membrane_lower, CartDeriv.InPlaneCartesianDerivativesGauss3, InPlaneGradientFGauss, 2);

    CC.B_membrane_lower *= 0.333333333333333333333333333333333;
    CC.C_membrane_lower *= 0.333333333333333333333333333333333;

    // Upper face

    // Gauss point 4
    CalculateInPlaneGradientFGauss(InPlaneGradientFGauss, CartDeriv.InPlaneCartesianDerivativesGauss4, NodesCoord, 0, 3);
    CalculateAndAdd_B_Membrane(CC.B_membrane_upper, CC.C_membrane_upper, CartDeriv.InPlaneCartesianDerivativesGauss4, InPlaneGradientFGauss, 0);

    // Gauss point 5
    CalculateInPlaneGradientFGauss(InPlaneGradientFGauss, CartDeriv.InPlaneCartesianDerivativesGauss5, NodesCoord, 1, 3);
    CalculateAndAdd_B_Membrane(CC.B_membrane_upper, CC.C_membrane_upper, CartDeriv.InPlaneCartesianDerivativesGauss5, InPlaneGradientFGauss, 1);

    // Gauss point 6
    CalculateInPlaneGradientFGauss(InPlaneGradientFGauss, CartDeriv.InPlaneCartesianDerivativesGauss6, NodesCoord, 2, 3);
    CalculateAndAdd_B_Membrane(CC.B_membrane_upper, CC.C_membrane_upper, CartDeriv.InPlaneCartesianDerivativesGauss6, InPlaneGradientFGauss, 2);

    CC.B_membrane_upper *= 0.333333333333333333333333333333333;
    CC.C_membrane_upper *= 0.333333333333333333333333333333333;

    /* SHEAR CONTRIBUTION */
    /* Calculating the shear strain-displacement matrix */

    // Declaring variables
    array_1d<double, 3 > TransverseGradientFt, TransverseGradientFxi, TransverseGradientFeta;

    // Lower face

    /* Calculate f components in the face */
    CalculateTransverseGradientFinP(TransverseGradientFt, TransverseGradientFxi, TransverseGradientFeta, NodesCoord, 0);

    /* Calculate f transverse components */
    CalculateTransverseGradientF(TransverseGradientF0, CartDeriv.TransversalCartesianDerivativesGauss1, NodesCoord);
    CalculateTransverseGradientF(TransverseGradientF1, CartDeriv.TransversalCartesianDerivativesGauss2, NodesCoord);
    CalculateTransverseGradientF(TransverseGradientF2, CartDeriv.TransversalCartesianDerivativesGauss3, NodesCoord);

    /* Shear contribution to the deformation matrix */
    CalculateAndAdd_B_Shear(CC.B_shear_lower, CC.C_shear_lower, CartDeriv.TransversalCartesianDerivativesGauss1,
                            CartDeriv.TransversalCartesianDerivativesGauss2, CartDeriv.TransversalCartesianDerivativesGauss3,
                            TransverseGradientF0, TransverseGradientF1, TransverseGradientF2, TransverseGradientFt,
                            TransverseGradientFxi, TransverseGradientFeta, CartDeriv.Jinv_plane_lower, 0);

    // Upper face

    /* Calculate f components in the face */
    CalculateTransverseGradientFinP(TransverseGradientFt, TransverseGradientFxi, TransverseGradientFeta, NodesCoord, 3);

    /* Calculate f transverse components */
    CalculateTransverseGradientF(TransverseGradientF0, CartDeriv.TransversalCartesianDerivativesGauss4, NodesCoord);
    CalculateTransverseGradientF(TransverseGradientF1, CartDeriv.TransversalCartesianDerivativesGauss5, NodesCoord);
    CalculateTransverseGradientF(TransverseGradientF2, CartDeriv.TransversalCartesianDerivativesGauss6, NodesCoord);

    /* Shear contribution to the deformation matrix */
    CalculateAndAdd_B_Shear(CC.B_shear_upper, CC.C_shear_upper, CartDeriv.TransversalCartesianDerivativesGauss4,
                            CartDeriv.TransversalCartesianDerivativesGauss5, CartDeriv.TransversalCartesianDerivativesGauss6,
                            TransverseGradientF0, TransverseGradientF1, TransverseGradientF2, TransverseGradientFt,
                            TransverseGradientFxi, TransverseGradientFeta, CartDeriv.Jinv_plane_upper, 9);

    /* NORMAL TRANSVERSE */
    /* Calculate f normal components */
    CalculateTransverseGradientF(TransverseGradientF0, CartDeriv.TransversalCartesianDerivativesCenter, NodesCoord);

    /* Calculating the normal transverse strain-displacement matrix */
    CalculateAndAdd_B_Normal(CC.B_normal, CC.C_normal, CartDeriv.TransversalCartesianDerivativesCenter, TransverseGradientF0);

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
    double norm; // TODO: Use the geometry normal when avalaible
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
    double OrthoComp;

    /* Performing the calculation */
    // 0- If X is the prefered normal vector
    // 1- If Y is the prefered normal vector
    // 2- If Z is the prefered normal vector

    if (choose == 0)
    {
        OrthoComp = mvze[1] * mvze[1] + mvze[2] * mvze[2]; // Component in th Y-Z plane
        if (OrthoComp < threshold) // If mvze is almost orthogonal to  Y-Z plane
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

            mvye[0] =          OrthoComp; // Choose mvxe orthogonal to global X direction
            mvye[1] = - mvze[0] * mvze[1];
            mvye[2] = - mvze[0] * mvze[2];

            norm = norm_2(mvye);
            mvye /= norm;
        }
    }
    else if (choose == 1)
    {
        OrthoComp = mvze[0] * mvze[0] + mvze[2] * mvze[2]; // Component in th Z-X plane
        if (OrthoComp < threshold) // If vze is almost orthogonal to  Z-X plane
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
            mvye[1] =          OrthoComp;
            mvye[2] = - mvze[2] * mvze[1];

            norm = norm_2(mvye);
            mvye /= norm;
        }
    }
    else if (choose == 2)
    {
        OrthoComp = mvze[0] * mvze[0] + mvze[1] * mvze[1]; // Component in th X-Y plane
        if (OrthoComp < threshold) // If vze is almost orthogonal to  X-Y plane
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
            mvye[2] =          OrthoComp;

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

    WeakPointerVector< Node < 3 > >& NeighbourNodes = this->GetValue(NEIGHBOUR_NODES);

    /* Compute mid_vec */ // TODO: Optimze this
    for (unsigned int i = 0; i < 18; i++)
    {
        mid_vec[i] = i;
    }
    unsigned int index = 18;
    for (unsigned int i = 0; i < 6; i++)
    {
        if (HasNeighbour(i, NeighbourNodes[i]))
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
        bounded_matrix<double, 6, 3 > & LocalDerivativePatch,
        const double xi,
        const double eta,
        const double zeta
        )
{
    double L_1 = 0.5 * (1.0 - zeta);
    double L_2 = 0.5 * (1.0 + zeta);

    /* Derivative in direction nu and xi */
    // Lower face
    LocalDerivativePatch(0, 0) = - L_1;
    LocalDerivativePatch(1, 0) =   L_1;
    LocalDerivativePatch(2, 0) =   0.0;

    LocalDerivativePatch(0, 1) = - L_1;
    LocalDerivativePatch(1, 1) =   0.0;
    LocalDerivativePatch(2, 1) =   L_1;

    // Upper face
    LocalDerivativePatch(3, 0) = - L_2;
    LocalDerivativePatch(4, 0) =   L_2;
    LocalDerivativePatch(5, 0) =   0.0;

    LocalDerivativePatch(3, 1) = - L_2;
    LocalDerivativePatch(4, 1) =   0.0;
    LocalDerivativePatch(5, 1) =   L_2;

    /* Derivative in direction zeta */
    LocalDerivativePatch(0, 2) = - 1.0 + eta + xi;
    LocalDerivativePatch(1, 2) = - xi;
    LocalDerivativePatch(2, 2) = - eta;
    LocalDerivativePatch(3, 2) =   1.0 - eta - xi;
    LocalDerivativePatch(4, 2) =   xi;
    LocalDerivativePatch(5, 2) =   eta;
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::ComputeLocalDerivativesQuadratic(
        bounded_matrix<double, 4, 2 > & LocalDerivativePatch,
        const int NodeGauss
        )
{
    /* Local coordinates */
    double xi  = 0.0;
    double eta = 0.0;

    if (NodeGauss == 0)
    {
        xi  = 0.5;
        eta = 0.5;
    }
    else if (NodeGauss == 1)
    {
        xi  = 0.0;
        eta = 0.5;
    }
    else if (NodeGauss == 2)
    {
        xi  = 0.5;
        eta = 0.0;
    }

    LocalDerivativePatch = ZeroMatrix(4, 2);

    /* Derivative in main nodes */
    LocalDerivativePatch(0, 0) = - 1.0 + eta;
    LocalDerivativePatch(0, 1) = - 1.0 + xi;
    LocalDerivativePatch(1, 0) =   1.0 - eta;
    LocalDerivativePatch(1, 1) =   1.0 - xi - 2.0 * eta;
    LocalDerivativePatch(2, 0) =   1.0 - 2.0 * xi - eta;
    LocalDerivativePatch(2, 1) =   1.0 - xi;

    /* Derivative in neighbour nodes */
    if (NodeGauss == 0)
    {
        LocalDerivativePatch(3, 0) = xi + eta - 0.5;
        LocalDerivativePatch(3, 1) = xi + eta - 0.5;
    }
    else if (NodeGauss == 1)
    {
        LocalDerivativePatch(3, 0) = xi - 0.5;
        LocalDerivativePatch(3, 1) = 0.0;
    }
    else if (NodeGauss == 2)
    {
        LocalDerivativePatch(3, 0) = 0.0;
        LocalDerivativePatch(3, 1) = eta - 0.5;
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
    bounded_matrix<double, 3, 6 > NodesCoord;
    for (unsigned int i = 0; i < 6; i++)
    {
        const array_1d<double, 3> &CurrentPosition  = GetGeometry()[i].Coordinates();
        NodesCoord(0, i) = CurrentPosition[0];
        NodesCoord(1, i) = CurrentPosition[1];
        NodesCoord(2, i) = CurrentPosition[2];
    }

    double xi  = 0.333333333333333333333333333333333;
    double eta = 0.333333333333333333333333333333333;

    /* Local derivatives patch */
    bounded_matrix<double, 6, 3 > LocalDerivativePatch;
    ComputeLocalDerivatives(LocalDerivativePatch, xi, eta, zeta);

    /* Compute Jacobian */
    noalias(J[rPointNumber]) = prod(NodesCoord, LocalDerivativePatch);

    /* Compute inverse of the Jaccobian */
    MathUtils<double>::InvertMatrix( J[rPointNumber], Jinv[rPointNumber], detJ[rPointNumber] );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateJacobian(
        double & detJ,
        bounded_matrix<double, 3, 3 > & J,
        bounded_matrix<double, 6, 3 > & LocalDerivativePatch,
        const bounded_matrix<double, 12, 3 > & NodesCoord,
        const double xi,
        const double eta,
        const double zeta
        )
{
    /* Auxiliar coordinates of the nodes */
    bounded_matrix<double, 3, 6 > NodesCoordAux;

    for (unsigned int i = 0; i < 6; i++)
    {
        NodesCoordAux(0, i) = NodesCoord(i, 0);
        NodesCoordAux(1, i) = NodesCoord(i, 1);
        NodesCoordAux(2, i) = NodesCoord(i, 2);
    }

    /* Local derivatives patch */
    ComputeLocalDerivatives(LocalDerivativePatch, xi, eta, zeta);

    /* Compute Jacobian */
    noalias(J) = prod(NodesCoordAux, LocalDerivativePatch);

    /* Compute determinant */
    detJ = MathUtils<double>::Det3(J);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateJacobianAndInv(
        bounded_matrix<double, 3, 3 > & J,
        bounded_matrix<double, 3, 3 > & Jinv,
        bounded_matrix<double, 6, 3 > & LocalDerivativePatch,
        const bounded_matrix<double, 3, 6 > & NodesCoord,
        const double xi,
        const double eta,
        const double zeta
        )
{
    /* Local derivatives patch */
    ComputeLocalDerivatives(LocalDerivativePatch, xi, eta, zeta);

    /* Compute Jacobian */
    noalias(J) = prod(NodesCoord, LocalDerivativePatch);

    /* Compute inverse of the Jaccobian */
    double detJ;
    Jinv = MathUtils<double>::InvertMatrix<3>(J, detJ);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateJacobianAndInv(
        bounded_matrix<double, 3, 3 > & J,
        bounded_matrix<double, 3, 3 > & Jinv,
        const bounded_matrix<double, 3, 6 > & NodesCoord,
        const double xi,
        const double eta,
        const double zeta
        )
{
    /* Local derivatives patch */
    bounded_matrix<double, 6, 3 > LocalDerivativePatch;
    ComputeLocalDerivatives(LocalDerivativePatch, xi, eta, zeta);

    /* Compute Jacobian */
    noalias(J) = prod(NodesCoord, LocalDerivativePatch);

    /* Compute inverse of the Jaccobian */
    double detJ;
    Jinv = MathUtils<double>::InvertMatrix<3>(J, detJ);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateCartesianDerOnCenter_plane(
        const int index,
        const bounded_matrix<double, 12, 3 > & NodesCoord,
        bounded_matrix<double, 2, 4 > & CartesianDerivativesCenter
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
        const int NodeGauss,
        const int index,
        const bounded_matrix<double, 12, 3 > & NodesCoord,
        bounded_matrix<double, 2, 4 > & InPlaneCartesianDerivativesGauss
        )
{
    /* Local derivatives patch */
    bounded_matrix<double, 4, 2 > LocalDerivativePatch;
    ComputeLocalDerivativesQuadratic(LocalDerivativePatch,NodeGauss);

    /* Auxiliar coordinates of the nodes */
    bounded_matrix<double, 3, 4 > NodesCoordAux;

    for (unsigned int i = 0; i < 3; i++)
    {
        NodesCoordAux(0, i) = NodesCoord(i + index, 0);
        NodesCoordAux(1, i) = NodesCoord(i + index, 1);
        NodesCoordAux(2, i) = NodesCoord(i + index, 2);
    }

    NodesCoordAux(0, 3) = NodesCoord(NodeGauss + 6 + index, 0);
    NodesCoordAux(1, 3) = NodesCoord(NodeGauss + 6 + index, 1);
    NodesCoordAux(2, 3) = NodesCoord(NodeGauss + 6 + index, 2);

    /* Compute local derivatives */
    bounded_matrix<double, 3, 2 > Xd;
    noalias(Xd) = prod(NodesCoordAux, LocalDerivativePatch);

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
    bounded_matrix<double, 2, 2 > jac;
    jac(0, 0) = inner_prod(Xdxi,  t1g);
    jac(0, 1) = inner_prod(Xdxi,  t2g);
    jac(1, 0) = inner_prod(Xdeta, t1g);
    jac(1, 1) = inner_prod(Xdeta, t2g);

    /* Compute the inverse of the Jacobian */
    double AuxDet;
    const bounded_matrix<double, 2, 2 > JinvPlane = MathUtils<double>::InvertMatrix<2>(jac, AuxDet);

    /* Compute the Cartesian derivatives */
    noalias(InPlaneCartesianDerivativesGauss) = prod(JinvPlane, trans(LocalDerivativePatch));
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateCartesianDerOnGauss_trans(
        const bounded_matrix<double, 12, 3 > & NodesCoord,
        bounded_matrix<double, 6, 1 > & TransversalCartesianDerivativesGauss,
        const double xi,
        const double eta,
        const double zeta
        )
{
    /* Compute local derivatives */
    double det;
    bounded_matrix<double, 3, 3 > Xd;
    bounded_matrix<double, 6, 3 > LocalDerivativePatch;
    CalculateJacobian(det, Xd, LocalDerivativePatch, NodesCoord, xi, eta, zeta);

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
    bounded_matrix<double, 3, 3 > t = ZeroMatrix(3, 3);
    StructuralMechanicsMathUtilities::Comp_Orthonor_Base(t, t1g, t2g, t3g, mvxe, Xdxi, Xdeta);

    /* Compute Jacobian */
    bounded_matrix<double, 3, 3 > jac;
    noalias(jac) = prod(t, Xd);

    /* Compute inverse of the Jaccobian (just third column) */
    bounded_matrix<double, 3 ,1> JinvTrans;
    JinvTrans(0, 0) =   (jac(0, 1) * jac(1, 2) - jac(0, 2) * jac(1, 1)) / det;
    JinvTrans(1, 0) = - (jac(0, 0) * jac(1, 2) - jac(0, 2) * jac(1, 0)) / det;
    JinvTrans(2, 0) =   (jac(0, 0) * jac(1, 1) - jac(1, 0) * jac(0, 1)) / det;

    /* Compute Cartesian derivatives */
    noalias(TransversalCartesianDerivativesGauss) = prod(LocalDerivativePatch, JinvTrans);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateCartesianDerOnCenter_trans(
        CartesianDerivatives& CartDeriv,
        const bounded_matrix<double, 12, 3 > & NodesCoord,
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
        KRATOS_ERROR << " This part id is not possible, just 0, 1 or 2  " << part << std::endl;
    }

    /* Auxiliar coordinates of the nodes */
    bounded_matrix<double, 3, 6 > NodesCoordAux;
    for (unsigned int i = 0; i < 6; i++)
    {
        NodesCoordAux(0, i) = NodesCoord(i, 0);
        NodesCoordAux(1, i) = NodesCoord(i, 1);
        NodesCoordAux(2, i) = NodesCoord(i, 2);
    }

    /* Auxiliar components to calculate the Jacobian and his inverse */
    bounded_matrix<double, 3, 3 > J;
    bounded_matrix<double, 3, 3 > Jinv;

    if (part == 0)
    {
        /* Calculate the Jacobian and his inverse */
        bounded_matrix<double, 6, 3 > LocalDerivativePatch;
        CalculateJacobianAndInv(J, Jinv, LocalDerivativePatch, NodesCoordAux, xi, eta, zeta);

        // Compute cartesian (y3) derivatives of the shape functions necessary to compute f_3
        /* Compute Cartesian derivatives */
        bounded_matrix<double, 6, 3 > TransversalCartesianDerivativesGauss_aux;
        noalias(TransversalCartesianDerivativesGauss_aux) = prod(LocalDerivativePatch, Jinv);

        for (unsigned int i = 0; i < 6 ; i++)
        {
            CartDeriv.TransversalCartesianDerivativesCenter(i, 0) =
                     mvze[0] * TransversalCartesianDerivativesGauss_aux(i, 0) +
                     mvze[1] * TransversalCartesianDerivativesGauss_aux(i, 1) +
                     mvze[2] * TransversalCartesianDerivativesGauss_aux(i, 2);
        }
     }
     else
     {
        /* Calculate the Jacobian and his inverse */
        CalculateJacobianAndInv(J, Jinv, NodesCoordAux, xi, eta, zeta);

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
             CartDeriv.Jinv_plane_lower(0, 0) = inner_prod(Xdxi,  mvxe);
             CartDeriv.Jinv_plane_lower(0, 1) = inner_prod(Xdeta, mvxe);
             CartDeriv.Jinv_plane_lower(1, 0) = inner_prod(Xdxi,  mvye);
             CartDeriv.Jinv_plane_lower(1, 1) = inner_prod(Xdeta, mvye);
         }
         else if (part == 2)
         {
             CartDeriv.Jinv_plane_upper(0, 0) = inner_prod(Xdxi,  mvxe);
             CartDeriv.Jinv_plane_upper(0, 1) = inner_prod(Xdeta, mvxe);
             CartDeriv.Jinv_plane_upper(1, 0) = inner_prod(Xdxi,  mvye);
             CartDeriv.Jinv_plane_upper(1, 1) = inner_prod(Xdeta, mvye);
         }
     }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateInPlaneGradientFGauss(
        bounded_matrix<double, 3, 2 > & InPlaneGradientFGauss,
        const bounded_matrix<double, 2, 4 > & InPlaneCartesianDerivativesGauss,
        const bounded_matrix<double, 12, 3 > & NodesCoord,
        const int NodeGauss,
        const int index
        )
{
    /* Auxiliar operators */
    bounded_matrix<double, 3, 3 > NodesCoordAux;
    bounded_matrix<double, 3, 2 > InPlaneCartesianDerivativesGauss_aux;

    for (unsigned int i = 0; i < 3; i++)
    {
        NodesCoordAux(0, i) = NodesCoord(i + index, 0);
        NodesCoordAux(1, i) = NodesCoord(i + index, 1);
        NodesCoordAux(2, i) = NodesCoord(i + index, 2);

        InPlaneCartesianDerivativesGauss_aux(i, 0) = InPlaneCartesianDerivativesGauss(0, i);
        InPlaneCartesianDerivativesGauss_aux(i, 1) = InPlaneCartesianDerivativesGauss(1, i);
    }

    noalias(InPlaneGradientFGauss) = prod(NodesCoordAux, InPlaneCartesianDerivativesGauss_aux);

    WeakPointerVector< Node < 3 > >& NeighbourNodes = this->GetValue(NEIGHBOUR_NODES);
    if (HasNeighbour(NodeGauss, NeighbourNodes[NodeGauss]))
    {
        for (unsigned int j = 0; j < 3 ; j++)
        {
            InPlaneGradientFGauss(j, 0) += NodesCoord(NodeGauss + 6 + index, j) * InPlaneCartesianDerivativesGauss(0, 3);
            InPlaneGradientFGauss(j, 1) += NodesCoord(NodeGauss + 6 + index, j) * InPlaneCartesianDerivativesGauss(1, 3);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateTransverseGradientF(
        array_1d<double, 3 > & TransverseGradientF,
        const bounded_matrix<double, 6, 1 > & TransversalCartesianDerivativesGauss,
        const bounded_matrix<double, 12, 3 > & NodesCoord
        )
{
    noalias(TransverseGradientF) = ZeroVector(3);

    for (unsigned int i = 0; i < 6; i++)
    {
        for (unsigned int j = 0; j < 3; j++)
        {
            TransverseGradientF[j] += TransversalCartesianDerivativesGauss(i, 0) * NodesCoord(i, j);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateTransverseGradientFinP(
        array_1d<double, 3 > & TransverseGradientFt,
        array_1d<double, 3 > & TransverseGradientFxi,
        array_1d<double, 3 > & TransverseGradientFeta,
        const bounded_matrix<double, 12, 3 > & NodesCoord,
        const int index
        )
{
    TransverseGradientFt[0]   = NodesCoord(2 + index, 0) - NodesCoord(1 + index, 0);
    TransverseGradientFt[1]   = NodesCoord(2 + index, 1) - NodesCoord(1 + index, 1);
    TransverseGradientFt[2]   = NodesCoord(2 + index, 2) - NodesCoord(1 + index, 2);

    TransverseGradientFxi[0]  = NodesCoord(0 + index, 0) - NodesCoord(2 + index, 0);
    TransverseGradientFxi[1]  = NodesCoord(0 + index, 1) - NodesCoord(2 + index, 1);
    TransverseGradientFxi[2]  = NodesCoord(0 + index, 2) - NodesCoord(2 + index, 2);

    TransverseGradientFeta[0] = NodesCoord(1 + index, 0) - NodesCoord(0 + index, 0);
    TransverseGradientFeta[1] = NodesCoord(1 + index, 1) - NodesCoord(0 + index, 1);
    TransverseGradientFeta[2] = NodesCoord(1 + index, 2) - NodesCoord(0 + index, 2);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAdd_B_Membrane(
        bounded_matrix<double, 3, 18 > & B_membrane,
        bounded_matrix<double, 3, 1  > & C_membrane,
        const bounded_matrix<double, 2, 4 > & InPlaneCartesianDerivativesGauss,
        const bounded_matrix<double, 3, 2 > & InPlaneGradientFGauss,
        const int NodeGauss
        )
{
    for (unsigned int i = 0; i < 4; i++)
    {
        unsigned int base = i * 3;
        if (i == 3)
        {
            base += NodeGauss * 3;
        }
        for (unsigned int j = 0; j < 3; j++)
        {
            B_membrane(0, base + j) += InPlaneCartesianDerivativesGauss(0, i) * InPlaneGradientFGauss(j, 0);
            B_membrane(1, base + j) += InPlaneCartesianDerivativesGauss(1, i) * InPlaneGradientFGauss(j, 1);
            B_membrane(2, base + j) += InPlaneCartesianDerivativesGauss(1, i) * InPlaneGradientFGauss(j, 0)
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

    C_membrane(0, 0) += inner_prod(auxDeformationGradientF1, auxDeformationGradientF1);
    C_membrane(1, 0) += inner_prod(auxDeformationGradientF2, auxDeformationGradientF2);
    C_membrane(2, 0) += inner_prod(auxDeformationGradientF1, auxDeformationGradientF2);

}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAdd_Membrane_Kgeometric(
        bounded_matrix<double, 36, 36 > & Kgeometricmembrane,
        const bounded_matrix<double, 2, 4 > & InPlaneCartesianDerivativesGauss1,
        const bounded_matrix<double, 2, 4 > & InPlaneCartesianDerivativesGauss2,
        const bounded_matrix<double, 2, 4 > & InPlaneCartesianDerivativesGauss3,
        const array_1d<double, 3 > & S_membrane,
        const int index
        )
{
    bounded_matrix<double, 6, 6 > H = ZeroMatrix(6, 6);

    unsigned int ii;
    unsigned int jj;
    for (unsigned int i = 0; i < 4; i++)
    {
        for (unsigned int j = 0; j < 4; j++)
        {
            // Gauss 1
            ii = i;
            jj = j;
            H(ii, jj) += S_membrane[0] *  InPlaneCartesianDerivativesGauss1(0, i) * InPlaneCartesianDerivativesGauss1(0, j)
                       + S_membrane[1] *  InPlaneCartesianDerivativesGauss1(1, i) * InPlaneCartesianDerivativesGauss1(1, j)
                       + S_membrane[2] * (InPlaneCartesianDerivativesGauss1(0, i) * InPlaneCartesianDerivativesGauss1(1, j)
                                        + InPlaneCartesianDerivativesGauss1(1, i) * InPlaneCartesianDerivativesGauss1(0, j));

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

            H(ii, jj) += S_membrane[0] *  InPlaneCartesianDerivativesGauss2(0, i) * InPlaneCartesianDerivativesGauss2(0, j)
                       + S_membrane[1] *  InPlaneCartesianDerivativesGauss2(1, i) * InPlaneCartesianDerivativesGauss2(1, j)
                       + S_membrane[2] * (InPlaneCartesianDerivativesGauss2(0, i) * InPlaneCartesianDerivativesGauss2(1, j)
                                        + InPlaneCartesianDerivativesGauss2(1, i) * InPlaneCartesianDerivativesGauss2(0, j));

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

            H(ii, jj) += S_membrane[0] *  InPlaneCartesianDerivativesGauss3(0, i) * InPlaneCartesianDerivativesGauss3(0, j)
                       + S_membrane[1] *  InPlaneCartesianDerivativesGauss3(1, i) * InPlaneCartesianDerivativesGauss3(1, j)
                       + S_membrane[2] * (InPlaneCartesianDerivativesGauss3(0, i) * InPlaneCartesianDerivativesGauss3(1, j)
                                        + InPlaneCartesianDerivativesGauss3(1, i) * InPlaneCartesianDerivativesGauss3(0, j));
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
        bounded_matrix<double, 2, 18 > & B_shear,
        bounded_matrix<double, 2, 1 > & C_shear,
        const bounded_matrix<double, 6, 1 > & TransversalCartesianDerivativesGauss1,
        const bounded_matrix<double, 6, 1 > & TransversalCartesianDerivativesGauss2,
        const bounded_matrix<double, 6, 1 > & TransversalCartesianDerivativesGauss3,
        const array_1d<double, 3 > & TransverseGradientFGauss1,
        const array_1d<double, 3 > & TransverseGradientFGauss2,
        const array_1d<double, 3 > & TransverseGradientFGauss3,
        const array_1d<double, 3 > & TransverseGradientFt,
        const array_1d<double, 3 > & TransverseGradientFxi,
        const array_1d<double, 3 > & TransverseGradientFeta,
        const bounded_matrix<double, 2, 2 > & Jinv_plane,
        const int index
        )
{
    // Considering the Gauss point in the middle of the element
    double eta_p = 0.33333333333333333333;
    double xi_p  = 0.33333333333333333333;
    bounded_matrix<double, 2, 3 > Pa;
    Pa(0, 0) = - xi_p;
    Pa(0, 1) = - xi_p;
    Pa(0, 2) = 1.0 - xi_p;
    Pa(1, 0) = eta_p;
    Pa(1, 1) = eta_p - 1.0;
    Pa(1, 2) = eta_p;

    bounded_matrix<double, 3, 18 > aux_B_shear = ZeroMatrix(3, 18);

    /* First contribution*/
    for (unsigned int i = 0; i < 6; i++)
    {
        unsigned int base = i * 3;
        for (unsigned int j = 0; j < 3; j++)
        {
            aux_B_shear(0, base + j) += TransversalCartesianDerivativesGauss1(i, 0) * TransverseGradientFt[j];
            aux_B_shear(1, base + j) += TransversalCartesianDerivativesGauss2(i, 0) * TransverseGradientFxi[j];
            aux_B_shear(2, base + j) += TransversalCartesianDerivativesGauss3(i, 0) * TransverseGradientFeta[j];
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

    bounded_matrix<double, 2, 3 > aux_prod;
    noalias(aux_prod) = prod(Jinv_plane, Pa);
    noalias(B_shear) = prod(aux_prod, aux_B_shear);

    // Calculating the components of C
    bounded_matrix<double, 3, 1 > aux_C_shear;
    aux_C_shear(0, 0) =   inner_prod(TransverseGradientFt  , TransverseGradientFGauss1);
    aux_C_shear(1, 0) =   inner_prod(TransverseGradientFxi , TransverseGradientFGauss2);
    aux_C_shear(2, 0) =   inner_prod(TransverseGradientFeta, TransverseGradientFGauss3);

    noalias(C_shear) = prod(aux_prod, aux_C_shear);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAdd_Shear_Kgeometric(
    bounded_matrix<double, 18, 18 > & Kgeometricshear,
    const bounded_matrix<double, 6, 1 > & TransversalCartesianDerivativesGauss1,
    const bounded_matrix<double, 6, 1 > & TransversalCartesianDerivativesGauss2,
    const bounded_matrix<double, 6, 1 > & TransversalCartesianDerivativesGauss3,
    const bounded_matrix<double, 2, 2 > & Jinv_plane,
    const array_1d<double, 2 > & S_shear,
    const int index
    )
{
    const double Q1 = 0.333333333333333333333333333 * (S_shear[0] * Jinv_plane(0, 0) + S_shear[1] * Jinv_plane(0, 1));
    const double Q2 = 0.333333333333333333333333333 * (S_shear[0] * Jinv_plane(1, 0) + S_shear[1] * Jinv_plane(1, 1));

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
//        Kgeometricshear(i + index + 3, i + index + 3) -= q[0] * TransversalCartesianDerivativesGauss1(1 + delta, 0);
//        Kgeometricshear(i + index + 6, i + index + 6) += q[0] * TransversalCartesianDerivativesGauss1(2 + delta, 0);

//        /* Second assembling */
//        Kgeometricshear(i + index, i + index)         += q[1] * TransversalCartesianDerivativesGauss2(0 + delta, 0);
//        Kgeometricshear(i + index + 6, i + index + 6) -= q[1] * TransversalCartesianDerivativesGauss2(2 + delta, 0);
//        /* Third assembling */
//        Kgeometricshear(i + index, i + index)         -= q[2] * TransversalCartesianDerivativesGauss3(0 + delta, 0);
//        Kgeometricshear(i + index + 3, i + index + 3) += q[2] * TransversalCartesianDerivativesGauss3(1 + delta, 0);
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
                value = (-Q1 + Q2) *  TransversalCartesianDerivativesGauss1(i, 0);
            }
            else if (k == 1)
            {
                value = -(Q1 + 2.0 * Q2) * TransversalCartesianDerivativesGauss2(i, 0);
            }
            else if (k == 2)
            {
                value = (2.0 * Q1 + Q2) * TransversalCartesianDerivativesGauss3(i, 0);
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
        bounded_matrix<double, 1, 18 > & B_normal,
        double & C_normal,
        const bounded_matrix<double, 6, 1 > & TransversalCartesianDerivativesCenter,
        const array_1d<double, 3 > & TransversalDeformationGradientF
        )
{
        for (unsigned int i = 0; i < 6; i++)
        {
            unsigned int base = i * 3;
            B_normal(0, base)     = TransversalCartesianDerivativesCenter(i, 0) * TransversalDeformationGradientF[0];
            B_normal(0, base + 1) = TransversalCartesianDerivativesCenter(i, 0) * TransversalDeformationGradientF[1];
            B_normal(0, base + 2) = TransversalCartesianDerivativesCenter(i, 0) * TransversalDeformationGradientF[2];
        }

        C_normal = inner_prod(TransversalDeformationGradientF, TransversalDeformationGradientF);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAdd_Normal_Kgeometric(
        bounded_matrix<double, 18, 18 > & Kgeometricnormal,
        const bounded_matrix<double, 6, 1 > & TransversalCartesianDerivativesCenter,
        const double S_normal
        )
{
    bounded_matrix<double, 6, 6 > H = ZeroMatrix(6, 6);
    for (unsigned int i = 0; i < 6; i++)
    {
        double aux = S_normal * TransversalCartesianDerivativesCenter(i, 0);
        for (unsigned int j = 0; j < 6; j++)
        {
            H(i, j) =  aux * TransversalCartesianDerivativesCenter(j, 0);
        }
    }

    noalias(H) = S_normal * prod(TransversalCartesianDerivativesCenter, trans(TransversalCartesianDerivativesCenter));

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

bounded_matrix<double, 36, 1 > SprismElement3D6N::CalculateDisp(const int& step)
{
    KRATOS_TRY;

    bounded_matrix<double, 36, 1 > DisplacementsArray;

    WeakPointerVector< Node < 3 > >& NeighbourNodes = this->GetValue(NEIGHBOUR_NODES);

    /* Element nodes */
    for (unsigned int index = 0; index < 6; index++)
    {
        array_1d<double,3> disp = GetGeometry()[index].FastGetSolutionStepValue(DISPLACEMENT, step);
        DisplacementsArray(index * 3,     0) = disp[0];
        DisplacementsArray(index * 3 + 1, 0) = disp[1];
        DisplacementsArray(index * 3 + 2, 0) = disp[2];
    }

    /* Neighbour nodes */
    int NumberNeighbours = NumberOfActiveNeighbours(NeighbourNodes);

    if (NumberNeighbours == 6) // All the possible neighours
    {
        for (unsigned int index = 0; index < 6; index++)
        {
            array_1d<double,3> disp = NeighbourNodes[index].FastGetSolutionStepValue(DISPLACEMENT, step);
            DisplacementsArray(18 + index * 3    , 0) = disp[0];
            DisplacementsArray(18 + index * 3 + 1, 0) = disp[1];
            DisplacementsArray(18 + index * 3 + 2, 0) = disp[2];
        }
    }
    else
    {
        for (unsigned int index = 0; index < 6; index++)
        {
            if (HasNeighbour(index, NeighbourNodes[index]))
            {
                array_1d<double,3> disp = NeighbourNodes[index].FastGetSolutionStepValue(DISPLACEMENT, step);
                DisplacementsArray(18 + index * 3    , 0) = disp[0];
                DisplacementsArray(18 + index * 3 + 1, 0) = disp[1];
                DisplacementsArray(18 + index * 3 + 2, 0) = disp[2];
            }
            else
            {
                DisplacementsArray(18 + index * 3    , 0) = 0.0;
                DisplacementsArray(18 + index * 3 + 1, 0) = 0.0;
                DisplacementsArray(18 + index * 3 + 2, 0) = 0.0;
            }
        }
    }

    return DisplacementsArray;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

bounded_matrix<double, 36, 1 > SprismElement3D6N::GetVectorCurrentPosition()
{
    KRATOS_TRY;

    bounded_matrix<double, 36, 1 > VectorCurrentPosition;

    WeakPointerVector< Node < 3 > >& NeighbourNodes = this->GetValue(NEIGHBOUR_NODES);

    /* Element nodes */
    for (unsigned int index = 0; index < 6; index++)
    {
        array_1d<double,3> CurrentPosition = GetGeometry()[index].Coordinates();
        VectorCurrentPosition(index * 3,     0) = CurrentPosition[0];
        VectorCurrentPosition(index * 3 + 1, 0) = CurrentPosition[1];
        VectorCurrentPosition(index * 3 + 2, 0) = CurrentPosition[2];
    }

    /* Neighbour nodes */
    int NumberNeighbours = NumberOfActiveNeighbours(NeighbourNodes);

    if (NumberNeighbours == 6) // All the possible neighours
    {
        for (unsigned int index = 0; index < 6; index++)
        {
            array_1d<double,3> CurrentPosition = NeighbourNodes[index].Coordinates();
            VectorCurrentPosition(18 + index * 3    , 0) = CurrentPosition[0];
            VectorCurrentPosition(18 + index * 3 + 1, 0) = CurrentPosition[1];
            VectorCurrentPosition(18 + index * 3 + 2, 0) = CurrentPosition[2];
        }
    }
    else
    {
        for (unsigned int index = 0; index < 6; index++)
        {
            if (HasNeighbour(index, NeighbourNodes[index]))
            {
                array_1d<double,3> CurrentPosition = NeighbourNodes[index].Coordinates();
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
        StressIntegratedComponents& IntStress,
        const double& AlphaEAS,
        const double& ZetaGauss,
        const double& rIntegrationWeight
        )
{
    KRATOS_TRY;
    
    const double L1 = 0.5 * (1.0 - ZetaGauss);
    const double L2 = 0.5 * (1.0 + ZetaGauss);

    const double FactorEAS = std::exp(2.0 * AlphaEAS * ZetaGauss);

    /* INTEGRATE PK2 IN ZETA */
    // Integrate stresses in the reference configuration
    /* In plane stresses */
    // Lower
    IntStress.S_membrane_lower(0) +=  L1 * rIntegrationWeight * rVariables.StressVector[0]; // xx
    IntStress.S_membrane_lower(1) +=  L1 * rIntegrationWeight * rVariables.StressVector[1]; // yy
    IntStress.S_membrane_lower(2) +=  L1 * rIntegrationWeight * rVariables.StressVector[3]; // xy
    // Upper
    IntStress.S_membrane_upper(0) +=  L2 * rIntegrationWeight * rVariables.StressVector[0]; // xx
    IntStress.S_membrane_upper(1) +=  L2 * rIntegrationWeight * rVariables.StressVector[1]; // yy
    IntStress.S_membrane_upper(2) +=  L2 * rIntegrationWeight * rVariables.StressVector[3]; // xy

    /* Transversal stresses */ // Note: Order according to the Voigt Notation in the Wiki
    // Lower face
    IntStress.S_shear_lower(0)    +=  L1 * rIntegrationWeight * rVariables.StressVector[5]; // xz
    IntStress.S_shear_lower(1)    +=  L1 * rIntegrationWeight * rVariables.StressVector[4]; // yz
    // Upper face
    IntStress.S_shear_upper(0)    +=  L2 * rIntegrationWeight * rVariables.StressVector[5]; // xz
    IntStress.S_shear_upper(1)    +=  L2 * rIntegrationWeight * rVariables.StressVector[4]; // yz

    /* Normal stress */
    IntStress.S_normal            +=  FactorEAS * rIntegrationWeight * rVariables.StressVector[2]; // zz

    /* INTEGRATE EAS IN ZETA */
    // Calculate EAS residual
    mEAS.rhs_alpha += rIntegrationWeight * ZetaGauss * rVariables.StressVector[2] * rVariables.C[2];

    // Calculate EAS stiffness
    mEAS.stiff_alpha += rIntegrationWeight * ZetaGauss * ZetaGauss * rVariables.C[2]
            * (rVariables.ConstitutiveMatrix(2, 2) * rVariables.C[2] + 2.0 * rVariables.StressVector[2]);

    bounded_matrix<double, 1, 36 > B3;
    bounded_matrix<double, 1,  6 > D3;

    for (unsigned int i = 0; i < 6; i++)
    {
        D3(0, i) = rVariables.ConstitutiveMatrix(2, i);
    }
    for (unsigned int i = 0; i < 36; i++)
    {
        B3(0, i) = rVariables.B(2, i);
    }

    // Calculate H operator
    noalias(mEAS.H_EAS) += rIntegrationWeight * ZetaGauss
            * (rVariables.C[2] * prod(D3, rVariables.B) + 2.0 * rVariables.StressVector[2] * B3);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAddLHS(
        LocalSystemComponents& rLocalSystem,
        GeneralVariables& rVariables,
        ConstitutiveLaw::Parameters& rValues,
        const StressIntegratedComponents& IntStress,
        const CommonComponents& CC,
        const CartesianDerivatives& CartDeriv,
        double& AlphaEAS
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
                const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

                for ( unsigned int PointNumber = 0; PointNumber < IntegrationPoints.size(); PointNumber++ )
                {
                    const double ZetaGauss = 2.0 * IntegrationPoints[PointNumber].Z() - 1.0;

                    /* Assemble B */
                    this->CalculateDeformationMatrix(rVariables.B, CC, ZetaGauss, AlphaEAS);

                    // Compute element kinematics C, F ...
                    this->CalculateKinematics(rVariables, CC, PointNumber, AlphaEAS, ZetaGauss);

                    // Set general variables to constitutivelaw parameters
                    this->SetGeneralVariables(rVariables, rValues, PointNumber);

                    // Compute stresses and constitutive parameters
                    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, rVariables.StressMeasure);

                    // Calculating weights for integration on the "reference configuration"
                    const double IntegrationWeight = IntegrationPoints[PointNumber].Weight() * rVariables.detJ;

                    /* Operation performed: add Km to the LefsHandSideMatrix */
                    this->CalculateAndAddKuum( rLeftHandSideMatrices[i], rVariables, IntegrationWeight);
                }
                calculated = true;
            }

            /* Calculate the Geometric Stiffness Matrix */
            if( rLeftHandSideVariables[i] == GEOMETRIC_STIFFNESS_MATRIX )
            {
                /* Operation performed: add Kg to the LefsHandSideMatrix */
                this->CalculateAndAddKuug( rLeftHandSideMatrices[i], IntStress, CartDeriv );
                calculated = true;
            }

            /* Implicit or explicit EAS update*/
            if ( mELementalFlags.Is(SprismElement3D6N::EAS_IMPLICIT_EXPLICIT) == true )
            {
                /* Apply EAS stabilization */
                ApplyEASLHS(rLeftHandSideMatrices[i]);
            }

            if(calculated == false)
            {
                KRATOS_ERROR << " ELEMENT can not supply the required local system variable: " << rLeftHandSideVariables[i] << std::endl;
            }
        }
    }
    else
    {
        MatrixType& LeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        /* Calculate the Material Stiffness Matrix */
        for ( unsigned int PointNumber = 0; PointNumber < IntegrationPoints.size(); PointNumber++ )
        {
            const double ZetaGauss = 2.0 * IntegrationPoints[PointNumber].Z() - 1.0;

            /* Assemble B */
            this->CalculateDeformationMatrix(rVariables.B, CC, ZetaGauss, AlphaEAS);

            // Compute element kinematics C, F ...
            this->CalculateKinematics(rVariables, CC, PointNumber, AlphaEAS, ZetaGauss);

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(rVariables, rValues, PointNumber);

            // Compute stresses and constitutive parameters
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, rVariables.StressMeasure);

            // Calculating weights for integration on the "reference configuration"
            const double IntegrationWeight = IntegrationPoints[PointNumber].Weight() * rVariables.detJ;

            /* Operation performed: add Km to the LefsHandSideMatrix */
            this->CalculateAndAddKuum( LeftHandSideMatrix, rVariables, IntegrationWeight);
        }

        /* Calculate the Geometric Stiffness Matrix */
        /* Operation performed: add Kg to the LefsHandSideMatrix */
        this->CalculateAndAddKuug( LeftHandSideMatrix, IntStress, CartDeriv );

        /* Implicit or explicit EAS update*/
        if ( mELementalFlags.Is(SprismElement3D6N::EAS_IMPLICIT_EXPLICIT) == true )
        {
            /* Apply EAS stabilization */
            ApplyEASLHS(LeftHandSideMatrix);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAddDynamicLHS(
        MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight
        )
{
  // Mass matrix
  WeakPointerVector< Node < 3 > >& NeighbourNodes = this->GetValue(NEIGHBOUR_NODES);

  // Compute volume change
  double VolumeChange;
  this->CalculateVolumeChange( VolumeChange, rVariables );

  // Compute density
  double Density = VolumeChange * GetProperties()[DENSITY];

  unsigned int NumberOfNodes = GetGeometry().size() + NumberOfActiveNeighbours(NeighbourNodes);
  unsigned int MatSize = NumberOfNodes * 3;

  if (rLeftHandSideMatrix.size1() != MatSize)
  {
      rLeftHandSideMatrix.resize(MatSize, MatSize, false);
  }

  noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize, MatSize);

  double TotalMass = rIntegrationWeight * Density;

  // Manually
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
              rLeftHandSideMatrix(index, index + 12) =       TotalMass;
              rLeftHandSideMatrix(index, index + 15) =       TotalMass;
              // Symmetric part
              rLeftHandSideMatrix(index +  3, index) = 2.0 * TotalMass;
              rLeftHandSideMatrix(index +  6, index) = 2.0 * TotalMass;
              rLeftHandSideMatrix(index +  9, index) = 2.0 * TotalMass;
              rLeftHandSideMatrix(index + 12, index) =       TotalMass;
              rLeftHandSideMatrix(index + 15, index) =       TotalMass;
          }
          else if (i == 1)
          {
              // Superior band
              rLeftHandSideMatrix(index, index +  3) = 2.0 * TotalMass;
              rLeftHandSideMatrix(index, index +  6) =       TotalMass;
              rLeftHandSideMatrix(index, index +  9) = 2.0 * TotalMass;
              rLeftHandSideMatrix(index, index + 12) =       TotalMass;
              // Symmetric part
              rLeftHandSideMatrix(index +  3, index) = 2.0 * TotalMass;
              rLeftHandSideMatrix(index +  6, index) =       TotalMass;
              rLeftHandSideMatrix(index +  9, index) = 2.0 * TotalMass;
              rLeftHandSideMatrix(index + 12, index) =       TotalMass;
          }
          else if (i == 2)
          {
              // Superior band
              rLeftHandSideMatrix(index, index + 3) =       TotalMass;
              rLeftHandSideMatrix(index, index + 6) =       TotalMass;
              rLeftHandSideMatrix(index, index + 9) = 2.0 * TotalMass;
              // Symmetric part
              rLeftHandSideMatrix(index + 3, index) =       TotalMass;
              rLeftHandSideMatrix(index + 6, index) =       TotalMass;
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
          rLeftHandSideMatrix(index, index)         = 4.0 * TotalMass;
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
        const StressIntegratedComponents& IntStress,
        const CommonComponents& CC,
        double& AlphaEAS
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
                this->CalculateAndAddInternalForces( RightHandSideVectors[i], IntStress, CC, AlphaEAS );
                calculated = true;
            }

            if(calculated == false)
            {
                KRATOS_ERROR << " ELEMENT can not supply the required local system variable: " << rRightHandSideVariables[i] << std::endl;
            }
        }
    }
    else
    {
        VectorType& RightHandSideVector = rLocalSystem.GetRightHandSideVector();

        /* Operation performed: RightHandSideVector += ExtForce */
        this->CalculateAndAddExternalForces( RightHandSideVector, rVariables, rVolumeForce );

        /* Operation performed: RightHandSideVector -= IntForce */
        this->CalculateAndAddInternalForces( RightHandSideVector, IntStress, CC, AlphaEAS );
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAddDynamicRHS(
        VectorType& rRightHandSideVector,
        GeneralVariables& rVariables,
        ProcessInfo& rCurrentProcessInfo,
        const double& rIntegrationWeight
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
    
    bounded_matrix<double, 36, 36 >      K; // Local stiffness matrix
    bounded_matrix<double,  6, 36 > aux_CB;

    /* Calculate K */
    noalias(aux_CB) = prod(rVariables.ConstitutiveMatrix, rVariables.B);
    noalias(K)      = rIntegrationWeight * prod(trans(rVariables.B), aux_CB);

    for (unsigned int i = 0; i < 36; i++)
    {
        if (mid_vec[i] < 1000)
        {
            for (unsigned int j = 0; j < 36; j++)
            {
                if (mid_vec[j] < 1000)
                {
                    rLeftHandSideMatrix(mid_vec[i], mid_vec[j]) += K(i, j);
                }
            }
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAddKuug(
        MatrixType& rLeftHandSideMatrix,
        const StressIntegratedComponents& IntStress,
        const CartesianDerivatives& CartDeriv
        )
{
    KRATOS_TRY;

    /* The stress is already integrated, we just calculate it once */

    /* Auxiliar stiffness matrix */
    bounded_matrix<double, 18, 18 > aux_K = ZeroMatrix(18, 18); // Auxiliar stiffness matrix
    bounded_matrix<double, 36, 36 >     K = ZeroMatrix(36, 36); // Stiffness matrix

    /* COMPUTATION OF GEOMETRIC STIFFNESS MATRIX */

    /* MEMBRANE CONTRIBUTION */
    /* Adding the geometric membrane stiffness */
    // Lower face
    CalculateAndAdd_Membrane_Kgeometric(K, CartDeriv.InPlaneCartesianDerivativesGauss1, CartDeriv.InPlaneCartesianDerivativesGauss2,
                                           CartDeriv.InPlaneCartesianDerivativesGauss3, IntStress.S_membrane_lower, 0);
    // Upper face
    CalculateAndAdd_Membrane_Kgeometric(K, CartDeriv.InPlaneCartesianDerivativesGauss4, CartDeriv.InPlaneCartesianDerivativesGauss5,
                                           CartDeriv.InPlaneCartesianDerivativesGauss6, IntStress.S_membrane_upper, 9);

//    /* SHEAR CONTRIBUTION */
//    /* Adding the geometric shear stiffness */
//    // Lower face
//    CalculateAndAdd_Shear_Kgeometric(aux_K, CartDeriv.TransversalCartesianDerivativesGauss1,
//                                            CartDeriv.TransversalCartesianDerivativesGauss2, CartDeriv.TransversalCartesianDerivativesGauss3,
//                                            CartDeriv.Jinv_plane_lower, IntStress.S_shear_lower, 0);
//    // Upper face
//    CalculateAndAdd_Shear_Kgeometric(aux_K, CartDeriv.TransversalCartesianDerivativesGauss4,
//                                            CartDeriv.TransversalCartesianDerivativesGauss5, CartDeriv.TransversalCartesianDerivativesGauss6,
//                                            CartDeriv.Jinv_plane_upper, IntStress.S_shear_upper, 9);

    /* NORMAL TRANSVERSE */
    /* Adding the geometric normal stiffness */
    CalculateAndAdd_Normal_Kgeometric(aux_K, CartDeriv.TransversalCartesianDerivativesCenter, IntStress.S_normal);

    // Transfering to the complete stiffness matrix
    for (unsigned int i = 0; i < 18; i++)
    {
        for (unsigned int j = 0; j < 18; j++)
        {
            K(i, j) += aux_K(i, j);
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
                    rLeftHandSideMatrix(mid_vec[i], mid_vec[j]) += K(i, j);
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

void SprismElement3D6N::ApplyEASLHS(MatrixType& rLeftHandSideMatrix)
{
    KRATOS_TRY;

    // Allocate auxiliar vector and matrix
    bounded_matrix<double, 36, 36 > lhs_aux = ZeroMatrix(36, 36);

    noalias(lhs_aux) -= prod(trans(mEAS.H_EAS), mEAS.H_EAS) / mEAS.stiff_alpha;

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
        bounded_matrix<double, 36, 1 > & rhs_full,
        double& AlphaEAS
        )
{
    KRATOS_TRY;

    /* Calculate the RHS */
    noalias(rhs_full) -= trans(mEAS.H_EAS) * mEAS.rhs_alpha / mEAS.stiff_alpha;

    /* Update ALPHA_EAS */
    AlphaEAS -= mEAS.rhs_alpha / mEAS.stiff_alpha;

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

    unsigned int NumberOfNodes = GetGeometry().PointsNumber();
    double aux_div = NumberOfNodes;

    for ( unsigned int i = 0; i < NumberOfNodes; i++ )
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
        const StressIntegratedComponents& IntStress,
        const CommonComponents& CC,
        double& AlphaEAS
        )
{
    KRATOS_TRY;

    bounded_matrix<double, 36, 1 > rhs_full = ZeroMatrix(36, 1);

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
        rhs_full(aux_index + i, 0)     += IntStress.S_membrane_lower[0] * CC.B_membrane_lower(0, i); // xx
        rhs_full(aux_index + i, 0)     += IntStress.S_membrane_lower[1] * CC.B_membrane_lower(1, i); // yy
        rhs_full(aux_index + i, 0)     += IntStress.S_membrane_lower[2] * CC.B_membrane_lower(2, i); // xy

        /* Nodes 4-6  and 10-12 */
        rhs_full(aux_index + i + 9, 0) += IntStress.S_membrane_upper[0] * CC.B_membrane_upper(0, i); // xx
        rhs_full(aux_index + i + 9, 0) += IntStress.S_membrane_upper[1] * CC.B_membrane_upper(1, i); // yy
        rhs_full(aux_index + i + 9, 0) += IntStress.S_membrane_upper[2] * CC.B_membrane_upper(2, i); // xy

        /* Apply transversal forces */
        /* Apply shear stress, adding the transverse nodal force contribution */
        rhs_full(i, 0) += IntStress.S_shear_lower[0] * CC.B_shear_lower(0, i); // xz
        rhs_full(i, 0) += IntStress.S_shear_lower[1] * CC.B_shear_lower(1, i); // yz
        rhs_full(i, 0) += IntStress.S_shear_upper[0] * CC.B_shear_upper(0, i); // xz
        rhs_full(i, 0) += IntStress.S_shear_upper[1] * CC.B_shear_upper(1, i); // yz

        /* Apply normal transverse stress */
        rhs_full(i, 0) += IntStress.S_normal * CC.B_normal(0, i); // zz
    }

    /* Apply EAS stabilization */
    ApplyEASRHS(rhs_full, AlphaEAS);

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
        unsigned int NumberOfNodes = GetGeometry().PointsNumber();

        for (unsigned int i = 0; i < NumberOfNodes; i++)
        {
            array_1d<double, 3> &CurrentPosition  = GetGeometry()[i].Coordinates();
            array_1d<double, 3> & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3> & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
            array_1d<double, 3> PreviousPosition  = CurrentPosition - (CurrentDisplacement-PreviousDisplacement);
            std::cout<<" NODE ["<<GetGeometry()[i].Id()<<"]: "<<PreviousPosition<<" (Cur: "<<CurrentPosition<<") "<<std::endl;
            std::cout<<" ---Disp: "<<CurrentDisplacement<<" (Pre: "<<PreviousDisplacement<<")"<<std::endl;
        }

        KRATOS_WATCH(rVariables.F);
        KRATOS_ERROR << " SPRISM ELEMENT INVERTED: |F| < 0  detF = " << rVariables.detF << std::endl;
    }

    // Compute total F: FT
    rVariables.detFT = rVariables.detF * rVariables.detF0;
    rVariables.FT    = prod( rVariables.F, rVariables.F0 );

    rValues.SetDeterminantF(rVariables.detFT);
    rValues.SetDeformationGradientF(rVariables.FT);
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
    WeakPointerVector< Node < 3 > >& NeighbourNodes = this->GetValue(NEIGHBOUR_NODES);
    const unsigned int NumberOfNodes = GetGeometry().size() + NumberOfActiveNeighbours(NeighbourNodes);
    const unsigned int MatSize = NumberOfNodes * 3;

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
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;
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

    const unsigned int NumberOfNodes = GetGeometry().PointsNumber();
    for ( unsigned int i = 0; i < NumberOfNodes; i++ )
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
        const CommonComponents& CC,
        const int& rPointNumber,
        const double& AlphaEAS,
        const double& ZetaGauss
        )
{
    KRATOS_TRY;

    const double & L_1 = 0.5 * (1.0 - ZetaGauss);
    const double & L_2 = 0.5 * (1.0 + ZetaGauss);

    const double & FactorEAS = std::exp(2.0 * AlphaEAS * ZetaGauss);  // EAS factor

    /* Assemble C */
    rVariables.C[0] = L_1 * CC.C_membrane_lower(0, 0) + L_2 * CC.C_membrane_upper(0, 0); // xx
    rVariables.C[1] = L_1 * CC.C_membrane_lower(1, 0) + L_2 * CC.C_membrane_upper(1, 0); // yy
    rVariables.C[2] = FactorEAS * CC.C_normal;                                            // zz
    rVariables.C[3] = L_1 * CC.C_membrane_lower(2, 0) + L_2 * CC.C_membrane_upper(2, 0); // xy
    rVariables.C[4] = L_1 * CC.C_shear_lower(1, 0)    + L_2 * CC.C_shear_upper(1, 0);    // yz
    rVariables.C[5] = L_1 * CC.C_shear_lower(0, 0)    + L_2 * CC.C_shear_upper(0, 0);    // xz

    rVariables.detF = rVariables.C[0] * rVariables.C[1] * rVariables.C[2] + 2 * rVariables.C[3] * rVariables.C[4] * rVariables.C[5]
                    - rVariables.C[5] * rVariables.C[5] * rVariables.C[1] -     rVariables.C[4] * rVariables.C[4] * rVariables.C[0]
                    - rVariables.C[3] * rVariables.C[3] * rVariables.C[2];

    if (rVariables.detF < 1.0e-8)
    {
        KRATOS_WATCH(rVariables.C);
        KRATOS_ERROR << "The determinant of C is zero or negative.  det(C): " << rVariables.detF << std::endl;
    }

    rVariables.detF = std::sqrt(rVariables.detF);

    if ( mELementalFlags.Is(SprismElement3D6N::TOTAL_UPDATED_LAGRANGIAN) == true )
    {
        // PK2 stress measure
        rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_PK2;

        // Jacobian Determinant for the isoparametric and numerical integration
        rVariables.detJ = mAuxCont[rPointNumber];
    }
    else
    {
        // Cauchy stress measure
        rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

        //Determinant of the Deformation Gradient F0
        rVariables.detF0 = mAuxCont[rPointNumber];
        rVariables.F0    = mAuxMatCont[rPointNumber];
    }

    this->CbartoFbar(rVariables, rPointNumber);

    // Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    // Set Shape Functions Values for this integration point
    rVariables.N = row( Ncontainer, rPointNumber);

    KRATOS_CATCH( "" );
}

/***************************** COMPUTE DELTA POSITION ******************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateDeltaPosition(Matrix & rDeltaPosition)
{
    KRATOS_TRY;

    rDeltaPosition = ZeroMatrix( 6 , 3);

    for ( unsigned int i = 0; i < 6; i++ )
    {
        const array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
        const array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, 1);

        for ( unsigned int j = 0; j < 3; j++ )
        {
            rDeltaPosition(i,j) = CurrentDisplacement[j] - PreviousDisplacement[j];
        }
    }

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

    /* We perform a polar decomposition of the CBar and F(regular) to obtain F_bar */

    /* Decompose CBar */
    bounded_matrix<double, 3, 3> EigenVectorsMatrix;
    bounded_matrix<double, 3, 3> EigenValuesMatrix;

    // Assemble matrix CBar
    const Matrix CBar = MathUtils<double>::VectorToSymmetricTensor(rVariables.C);

    // Decompose matrix CBar
    MathUtils<double>::EigenSystem<3>(CBar, EigenVectorsMatrix, EigenValuesMatrix, 1e-24, 100);

    for (unsigned int i = 0; i < 3; i++)
    {
        EigenValuesMatrix(i, i) = std::sqrt(EigenValuesMatrix(i, i));
    }

    const Matrix UBar = prod( EigenValuesMatrix, EigenVectorsMatrix );

    /* Decompose F */
    Matrix F = ZeroMatrix(3, 3);
    if ( mELementalFlags.Is(SprismElement3D6N::TOTAL_UPDATED_LAGRANGIAN) == true )
    {
        // Deformation Gradient F [dx_n+1/dx_n]
        noalias(F) = prod( rVariables.j[rPointNumber], mAuxMatCont[rPointNumber] );
    }
    else
    {
        // Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
        Matrix InvJ(3, 3);
        MathUtils<double>::InvertMatrix( rVariables.J[rPointNumber], InvJ, rVariables.detJ);

        // Deformation Gradient F [dx_n+1/dx_n]
        noalias(F) = prod( rVariables.j[rPointNumber], InvJ );
    }

    const Matrix C = prod( trans(F), F );

    // Decompose matrix C
    MathUtils<double>::EigenSystem<3>(C, EigenVectorsMatrix, EigenValuesMatrix, 1e-24, 100);

    for (unsigned int i = 0; i < 3; i++)
    {
        EigenValuesMatrix(i, i) = std::sqrt(EigenValuesMatrix(i, i));
    }

    const Matrix U  = prod( EigenValuesMatrix, EigenVectorsMatrix );

    double AuxDet;
    Matrix invU(3, 3);
    MathUtils<double>::InvertMatrix(U, invU, AuxDet);
    const Matrix R  = prod( F, invU );

    /* Calculate F_bar */
    noalias(rVariables.F) = prod(R, UBar);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateDeformationMatrix(
        Matrix& rB,
        const CommonComponents& CC,
        const double& ZetaGauss,
        const double& AlphaEAS
        )
{
    KRATOS_TRY;

    rB.clear(); // Set all components to zero

    const double L1 = 0.5 * (1.0 - ZetaGauss);
    const double L2 = 0.5 * (1.0 + ZetaGauss);

    const double FactorEAS = std::exp(2.0 * AlphaEAS * ZetaGauss); // EAS factor

    for (unsigned int index = 0; index < 9; index++)
    {
        /* Element nodes */ // Note: It's important to consider the Voigt notation order considered in Kratos
        // Lower face
        rB(0, index)      = L1 * CC.B_membrane_lower(0, index);  // xx
        rB(1, index)      = L1 * CC.B_membrane_lower(1, index);  // yy
        rB(2, index)      = FactorEAS * CC.B_normal(0, index);     // zz
        rB(3, index)      = L1 * CC.B_membrane_lower(2, index);  // xy
        rB(4, index)      = L1 * CC.B_shear_lower(1, index) + L2 * CC.B_shear_upper(1, index); // yz
        rB(5, index)      = L1 * CC.B_shear_lower(0, index) + L2 * CC.B_shear_upper(0, index); // xz
        // Upper face
        rB(0, index + 9)  = L2 * CC.B_membrane_upper(0, index);  // xx
        rB(1, index + 9)  = L2 * CC.B_membrane_upper(1, index);  // yy
        rB(2, index + 9)  = FactorEAS * CC.B_normal(0, index + 9); // zz
        rB(3, index + 9)  = L2 * CC.B_membrane_upper(2, index);  // xy
        rB(4, index + 9)  = L1 * CC.B_shear_lower(1, index + 9) + L2 * CC.B_shear_upper(1, index + 9); // yz
        rB(5, index + 9)  = L1 * CC.B_shear_lower(0, index + 9) + L2 * CC.B_shear_upper(0, index + 9); // xz

        /* Neighbour nodes */
        // Lower face
        rB(0, index + 18) = L1 * CC.B_membrane_lower(0, index + 9); // xx
        rB(1, index + 18) = L1 * CC.B_membrane_lower(1, index + 9); // yy
        rB(3, index + 18) = L1 * CC.B_membrane_lower(2, index + 9); // xy
        // Upper face
        rB(0, index + 27) = L2 * CC.B_membrane_upper(0, index + 9); // xx
        rB(1, index + 27) = L2 * CC.B_membrane_upper(1, index + 9); // yy
        rB(3, index + 27) = L2 * CC.B_membrane_upper(2, index + 9); // xy
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::InitializeGeneralVariables(GeneralVariables& rVariables)
{
    // StressMeasure_PK1             //stress related to reference configuration non-symmetric
    // StressMeasure_PK2             //stress related to reference configuration
    // StressMeasure_Kirchhoff       //stress related to current   configuration
    // StressMeasure_Cauchy          //stress related to current   configuration

    // StressMeasure
    if ( mELementalFlags.Is(SprismElement3D6N::TOTAL_UPDATED_LAGRANGIAN) == true )
    {
        rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_PK2;
    }
    else
    {
        rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;
    }

    // Doubles
    rVariables.detF  = 1.0;
    rVariables.detF0 = 1.0;
    rVariables.detFT = 1.0;
    rVariables.detJ  = 1.0;

    // Vectors
    rVariables.StrainVector = ZeroVector(6);
    rVariables.StressVector = ZeroVector(6);
    rVariables.C = ZeroVector(6);
    rVariables.N = ZeroVector(6);

    // Matrices
    rVariables.F  = IdentityMatrix(3);
    rVariables.F0 = IdentityMatrix(3);
    rVariables.FT = IdentityMatrix(3);
    rVariables.B  = ZeroMatrix(6, 36);

    rVariables.DN_DX = ZeroMatrix(6, 3);
    rVariables.DeltaPosition = ZeroMatrix(6, 3);
    rVariables.ConstitutiveMatrix = ZeroMatrix(6, 6);

    // Reading shape functions
    rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

    // Jacobians
    rVariables.J.resize(1, false);
    rVariables.j.resize(1, false);
    rVariables.J[0] = ZeroMatrix(1, 1);
    rVariables.j[0] = ZeroMatrix(1, 1);

    // Calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
    rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod );

    if ( mELementalFlags.Is(SprismElement3D6N::TOTAL_UPDATED_LAGRANGIAN) == false )
    {
        //Calculate Delta Position
        this->CalculateDeltaPosition(rVariables.DeltaPosition);
        rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::FinalizeStepVariables(
        GeneralVariables & rVariables,
        const int& rPointNumber
        )
{
    if ( mELementalFlags.Is(SprismElement3D6N::TOTAL_UPDATED_LAGRANGIAN) == false )
    {
        // Update internal (historical) variables
        mAuxCont[rPointNumber]         = rVariables.detF * rVariables.detF0;
        mAuxMatCont[rPointNumber] = prod(rVariables.F, rVariables.F0);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::GetHistoricalVariables(
        GeneralVariables& rVariables,
        const int& rPointNumber
        )
{
    /* Deformation Gradient F ( set to identity ) */
    const unsigned int size =  rVariables.F.size1();

    rVariables.detF  = 1.0;
    rVariables.F     = IdentityMatrix(size);
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
    if ( rStrainVector.size() != 6 )
    {
        rStrainVector.resize( 6, false );
    }

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
    bounded_matrix<double, 3, 3> EigenValuesMatrix;
    bounded_matrix<double, 3, 3> EigenVectorsMatrix;

    // Assemble matrix C
    const Matrix CMatrix = MathUtils<double>::VectorToSymmetricTensor(rC);

    // Decompose matrix
    MathUtils<double>::EigenSystem<3>(CMatrix, EigenVectorsMatrix, EigenValuesMatrix, 1e-24, 10);

    // Calculate the eigenvalues of the E matrix
    EigenValuesMatrix(0, 0) = 0.5 * std::log(EigenValuesMatrix(0, 0));
    EigenValuesMatrix(1, 1) = 0.5 * std::log(EigenValuesMatrix(1, 1));
    EigenValuesMatrix(2, 2) = 0.5 * std::log(EigenValuesMatrix(2, 2));

    // Calculate E matrix
    bounded_matrix<double, 3, 3 > EMatrix;
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

/**************************** CALCULATE VOLUME CHANGE ******************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateVolumeChange( double& rVolumeChange, GeneralVariables& rVariables )
{
    KRATOS_TRY;

    if ( mELementalFlags.Is(SprismElement3D6N::TOTAL_UPDATED_LAGRANGIAN) == true )
    {
        rVolumeChange = 1.0;
    }
    else
    {
        rVolumeChange = 1.0 / (rVariables.detF * rVariables.detF0);
    }

    KRATOS_CATCH( "" );
}

/************************* CALCULATE VOLUME ACCELERATION ***************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateVolumeForce(
        Vector& rVolumeForce,
        GeneralVariables& rVariables
        )
{
    KRATOS_TRY;

    /* Reading integration points */
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    rVolumeForce = ZeroVector(3);

    array_1d<double,3> accel = ZeroVector(3);

    if( GetGeometry()[0].SolutionStepsDataHas(VOLUME_ACCELERATION) )
    {
        accel = GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    }

    for ( unsigned int PointNumber = 0; PointNumber < IntegrationPoints.size(); PointNumber++ )
    {
        double IntegrationWeight = IntegrationPoints[PointNumber].Weight() * rVariables.detJ;
        rVolumeForce += IntegrationWeight * accel;
    }

    // Compute volume change
    double VolumeChange;
    this->CalculateVolumeChange( VolumeChange, rVariables );

    rVolumeForce *= VolumeChange * GetProperties()[DENSITY];

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
    rSerializer.save("AuxMatCont",mAuxMatCont);
    rSerializer.save("AuxCont",mAuxCont);
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
    rSerializer.load("AuxMatCont",mAuxMatCont);
    rSerializer.load("AuxCont",mAuxCont);
}

} // Namespace Kratos.
