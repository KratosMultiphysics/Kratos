// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Quirin Aumann,
//                   Aron Noordam
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "spring_damper_element.hpp"

#include "structural_mechanics_application_variables.h"

namespace Kratos
{
//***********************DEFAULT CONSTRUCTOR******************************************
//************************************************************************************
    template<std::size_t TDim>
    SpringDamperElement<TDim>::SpringDamperElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************
    template<std::size_t TDim>
    SpringDamperElement<TDim>::SpringDamperElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : Element( NewId, pGeometry, pProperties )
{

}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************
template<std::size_t TDim>
SpringDamperElement<TDim>::SpringDamperElement(SpringDamperElement const& rOther)
    :Element(rOther)
{

}

//******************************ASSIGNMENT OPERATOR***********************************
//************************************************************************************
template<std::size_t TDim>
SpringDamperElement<TDim>& SpringDamperElement<TDim>::operator=(SpringDamperElement const& rOther)
{
    //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

    Element::operator=(rOther);

    return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************
template<std::size_t TDim>
Element::Pointer SpringDamperElement<TDim>::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
    //NEEDED TO CREATE AN ELEMENT
    return Kratos::make_intrusive<SpringDamperElement>( NewId, GetGeometry().Create( rThisNodes ), pProperties );
}

//************************************************************************************
//************************************************************************************
template<std::size_t TDim>
Element::Pointer SpringDamperElement<TDim>::Create( IndexType NewId,  GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    //NEEDED TO CREATE AN ELEMENT
    return Kratos::make_intrusive<SpringDamperElement>( NewId, pGeom, pProperties );
}

//************************************CLONE*******************************************
//************************************************************************************
template<std::size_t TDim>
Element::Pointer SpringDamperElement<TDim>::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
    //YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
    //ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

    SpringDamperElement NewElement(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    return Kratos::make_intrusive<SpringDamperElement>(NewElement);
}


template<std::size_t TDim>
SpringDamperElement<TDim>::~SpringDamperElement()
{
}

//************* GETTING METHODS
//************************************************************************************
//************************************************************************************
template<std::size_t TDim>
void SpringDamperElement<TDim>::GetDofList( DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo ) const
{
    //NEEDED TO DEFINE THE DOFS OF THE ELEMENT

    // Resizing as needed

    if (rElementalDofList.size() != msElementSize)
    {
        rElementalDofList.resize(msElementSize);
    }
    

    for ( SizeType i = 0; i < msNumNodes; ++i ) {
        const SizeType index = i * msLocalSize;
        if constexpr (TDim == 2){
            
            rElementalDofList[index] = GetGeometry()[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[index + 1] = GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
            rElementalDofList[index + 2] = GetGeometry()[i].pGetDof(ROTATION_Z);
        }
        else if constexpr (TDim == 3){
            rElementalDofList[index] = GetGeometry()[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[index + 1] = GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
            rElementalDofList[index + 2] = GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
            rElementalDofList[index + 3] = GetGeometry()[i].pGetDof(ROTATION_X);
            rElementalDofList[index + 4] = GetGeometry()[i].pGetDof(ROTATION_Y);
            rElementalDofList[index + 5] = GetGeometry()[i].pGetDof(ROTATION_Z);
        }
    }
}

//************************************************************************************
//************************************************************************************

template<std::size_t TDim>
void SpringDamperElement<TDim>::EquationIdVector( EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo ) const
{
    //NEEDED TO DEFINE GLOBAL IDS FOR THE CORRECT ASSEMBLY
    if ( rResult.size() != msElementSize )
    {
        rResult.resize( msElementSize, false );
    }

    for ( std::size_t i = 0; i < GetGeometry().size(); ++i)
    {
        const SizeType index = i * msLocalSize;
        if constexpr (TDim == 2) {
            rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(ROTATION_Z).EquationId();
        }
        else if constexpr (TDim == 3) {
            rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof(ROTATION_X).EquationId();
            rResult[index + 4] = GetGeometry()[i].GetDof(ROTATION_Y).EquationId();
            rResult[index + 5] = GetGeometry()[i].GetDof(ROTATION_Z).EquationId();
        }
    }
}

//*********************************DISPLACEMENT***************************************
//************************************************************************************

template<std::size_t TDim>
void SpringDamperElement<TDim>::GetValuesVector( Vector& rValues, int Step ) const
{
    //GIVES THE VECTOR WITH THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT DISPLACEMENTS)
    if ( rValues.size() != msElementSize )
    {
        rValues.resize( msElementSize, false );
    }

    for ( std::size_t i = 0; i < GetGeometry().size(); ++i)
    {
        const array_1d<double, 3>& disp = GetGeometry()[i].FastGetSolutionStepValue( DISPLACEMENT, Step );
        const array_1d<double, 3>& rot  = GetGeometry()[i].FastGetSolutionStepValue( ROTATION, Step );

        const SizeType index = i * msLocalSize;
        if constexpr (TDim == 2) {
            rValues[index] = disp[0];
            rValues[index + 1] = disp[1];
            rValues[index + 2] = rot[2];
        }
        else if constexpr (TDim == 3){
            rValues[index] = disp[0];
            rValues[index + 1] = disp[1];
            rValues[index + 2] = disp[2];
            rValues[index + 3] = rot[0];
            rValues[index + 4] = rot[1];
            rValues[index + 5] = rot[2];
        }

    }
}


//************************************VELOCITY****************************************
//************************************************************************************

template<std::size_t TDim>
void SpringDamperElement<TDim>::GetFirstDerivativesVector( Vector& rValues, int Step ) const
{
    //GIVES THE VECTOR WITH THE TIME DERIVATIVE OF THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT VELOCITIES)
    if ( rValues.size() != msElementSize )
    {
        rValues.resize( msElementSize, false );
    }

    for ( std::size_t i = 0; i < GetGeometry().size(); ++i)
    {
        const array_1d<double, 3>& vel = GetGeometry()[i].FastGetSolutionStepValue( VELOCITY, Step );
        const array_1d<double, 3>& avel = GetGeometry()[i].FastGetSolutionStepValue( ANGULAR_VELOCITY, Step );
        std::size_t index = i * msLocalSize;

        if constexpr (TDim == 2) {
            rValues[index] = vel[0];
            rValues[index + 1] = vel[1];
            rValues[index + 2] = avel[2];
        }
        else if constexpr (TDim == 3) {
            rValues[index] = vel[0];
            rValues[index + 1] = vel[1];
            rValues[index + 2] = vel[2];
            rValues[index + 3] = avel[0];
            rValues[index + 4] = avel[1];
            rValues[index + 5] = avel[2];
        }

    }
}

//*********************************ACCELERATION***************************************
//************************************************************************************

template<std::size_t TDim>
void SpringDamperElement<TDim>::GetSecondDerivativesVector( Vector& rValues, int Step ) const
{
    //GIVES THE VECTOR WITH THE TIME SECOND DERIVATIVE OF THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT ACCELERATIONS)
    if ( rValues.size() != msElementSize )
    {
        rValues.resize( msElementSize, false );
    }

    for ( std::size_t i = 0; i < GetGeometry().size(); ++i)
    {
        const array_1d<double, 3>& acc = GetGeometry()[i].FastGetSolutionStepValue( ACCELERATION, Step );
        const array_1d<double, 3>& aacc = GetGeometry()[i].FastGetSolutionStepValue( ANGULAR_ACCELERATION, Step );
        std::size_t index = i * msLocalSize;

        if constexpr (TDim == 2) {
            rValues[index] = acc[0];
            rValues[index + 1] = acc[1];
            rValues[index + 2] = aacc[2];
        }
        else if constexpr (TDim == 3) {
            rValues[index] = acc[0];
            rValues[index + 1] = acc[1];
            rValues[index + 2] = acc[2];
            rValues[index + 3] = aacc[0];
            rValues[index + 4] = aacc[1];
            rValues[index + 5] = aacc[2];
        }
    }
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

template<std::size_t TDim>
void SpringDamperElement<TDim>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{

    KRATOS_TRY;

    this->ConstCalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector);

    KRATOS_CATCH( "" );
}

//***********************************************************************************
//***********************************************************************************

template<std::size_t TDim>
void SpringDamperElement<TDim>::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    this->ConstCalculateRightHandSide(rRightHandSideVector);
}

//***********************************************************************************
//***********************************************************************************

template<std::size_t TDim>
void SpringDamperElement<TDim>::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    this->ConstCalculateLeftHandSide(rLeftHandSideMatrix);
}

//************************************************************************************
//************************************************************************************

template<std::size_t TDim>
void SpringDamperElement<TDim>::CalculateMassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //this is a massless element
    std::size_t system_size = msElementSize;

    if ( rMassMatrix.size1() != system_size )
    {
        rMassMatrix.resize( system_size, system_size, false );
    }

    rMassMatrix = ZeroMatrix( system_size, system_size );

    KRATOS_CATCH( "" );
}

//************************************************************************************
//************************************************************************************

template<std::size_t TDim>
void SpringDamperElement<TDim>::CalculateDampingMatrix( MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    this->ConstCalculateDampingMatrix(rDampingMatrix);

    KRATOS_CATCH( "" );
}

template<std::size_t TDim>
void SpringDamperElement<TDim>::ConstCalculateDampingMatrix(MatrixType& rDampingMatrix) const
{
    KRATOS_TRY;

    const std::size_t system_size = msElementSize;

    rDampingMatrix = ZeroMatrix(system_size, system_size);

    if (this->Has(NODAL_DAMPING_RATIO) || this->Has(NODAL_ROTATIONAL_DAMPING_RATIO)) {
        array_1d<double, msLocalSize> elemental_damping_ratio = ZeroVector(msLocalSize);
        if (this->Has(NODAL_DAMPING_RATIO)) {
            const array_1d<double, 3>& nodal_damping = this->GetValue(NODAL_DAMPING_RATIO);

            elemental_damping_ratio[0] = nodal_damping[0];
            elemental_damping_ratio[1] = nodal_damping[1];

            if constexpr (TDim == 3) {
                elemental_damping_ratio[2] = nodal_damping[2];
            }
        }

        if (this->Has(NODAL_ROTATIONAL_DAMPING_RATIO)) {
            const array_1d<double, 3>& nodal_rotational_damping = this->GetValue(NODAL_ROTATIONAL_DAMPING_RATIO);
            if constexpr (TDim == 2) {
                elemental_damping_ratio[2] = nodal_rotational_damping[2];
            }
            else if constexpr (TDim == 3) {
                elemental_damping_ratio[3] = nodal_rotational_damping[0];
                elemental_damping_ratio[4] = nodal_rotational_damping[1];
                elemental_damping_ratio[5] = nodal_rotational_damping[2];
            }
        }

        for (std::size_t i = 0; i < msLocalSize; ++i) {
            rDampingMatrix(i, i) += elemental_damping_ratio[i];
            rDampingMatrix(i + msLocalSize, i + msLocalSize) += elemental_damping_ratio[i];
            rDampingMatrix(i, i + msLocalSize) -= elemental_damping_ratio[i];
            rDampingMatrix(i + msLocalSize, i) -= elemental_damping_ratio[i];
        }
    }

    KRATOS_CATCH("");
}


template<std::size_t TDim>
void SpringDamperElement<TDim>::ConstCalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector) const
{
    KRATOS_TRY;

    /* Calculate elemental system */
    // Compute RHS (RHS = rRightHandSideVector = Fext - Fint)
    this->ConstCalculateRightHandSide(rRightHandSideVector);

    // Compute LHS
    this->ConstCalculateLeftHandSide(rLeftHandSideMatrix);

    KRATOS_CATCH("");
}


template<std::size_t TDim>
void SpringDamperElement<TDim>::ConstCalculateLeftHandSide(MatrixType& rLeftHandSideMatrix) const
{
    KRATOS_TRY;
    // Resizing the LHS
    std::size_t system_size = msElementSize;

    if (rLeftHandSideMatrix.size1() != system_size) {
        rLeftHandSideMatrix.resize(system_size, system_size, false);
    }

    noalias(rLeftHandSideMatrix) = ZeroMatrix(system_size, system_size); //resetting LHS

    // elemental_stiffness: kx, ky, kz, cpx, cpy, cpz
    array_1d<double, msLocalSize > elemental_stiffness = ZeroVector(msLocalSize);

    const array_1d<double, 3>& nodal_stiffness = this->GetValue(NODAL_DISPLACEMENT_STIFFNESS);
    const array_1d<double, 3>& nodal_rot_stiffness = this->GetValue(NODAL_ROTATIONAL_STIFFNESS);

    if constexpr (TDim == 2) {

        elemental_stiffness[0] = nodal_stiffness[0];
        elemental_stiffness[1] = nodal_stiffness[1];
        elemental_stiffness[2] = nodal_rot_stiffness[2];

    }
    else if constexpr (TDim == 3) {
        elemental_stiffness[0] = nodal_stiffness[0];
        elemental_stiffness[1] = nodal_stiffness[1];
        elemental_stiffness[2] = nodal_stiffness[2];

        elemental_stiffness[3] = nodal_rot_stiffness[0];
        elemental_stiffness[4] = nodal_rot_stiffness[1];
        elemental_stiffness[5] = nodal_rot_stiffness[2];
    }

    for (std::size_t i = 0; i < msLocalSize; ++i) {
        rLeftHandSideMatrix(i, i) += elemental_stiffness[i];
        rLeftHandSideMatrix(i + msLocalSize, i + msLocalSize) += elemental_stiffness[i];
        rLeftHandSideMatrix(i, i + msLocalSize) -= elemental_stiffness[i];
        rLeftHandSideMatrix(i + msLocalSize, i) -= elemental_stiffness[i];
    }
    KRATOS_CATCH("");
}

template<std::size_t TDim>
void SpringDamperElement<TDim>::ConstCalculateRightHandSide(VectorType& rRightHandSideVector) const
{
    KRATOS_TRY

    if (rRightHandSideVector.size() != msElementSize) {
        rRightHandSideVector.resize(msElementSize, false);
    }

    rRightHandSideVector = ZeroVector(msElementSize); //resetting RHS

    array_1d<double, msElementSize > current_displacement = ZeroVector(msElementSize);
    array_1d<double, msLocalSize> elemental_stiffness = ZeroVector(msLocalSize);

    // Getting actual values
    const array_1d<double, 3>& nodal_stiffness = this->GetValue(NODAL_DISPLACEMENT_STIFFNESS);
    const array_1d<double, 3>& nodal_rot_stiffness = this->GetValue(NODAL_ROTATIONAL_STIFFNESS);

    if constexpr (TDim == 2) {
        elemental_stiffness[0] = nodal_stiffness[0];
        elemental_stiffness[1] = nodal_stiffness[1];

        elemental_stiffness[2] = nodal_rot_stiffness[2];
    }
    else if constexpr (TDim == 3) {
        elemental_stiffness[0] = nodal_stiffness[0];
        elemental_stiffness[1] = nodal_stiffness[1];
        elemental_stiffness[2] = nodal_stiffness[2];

        elemental_stiffness[3] = nodal_rot_stiffness[0];
        elemental_stiffness[4] = nodal_rot_stiffness[1];
        elemental_stiffness[5] = nodal_rot_stiffness[2];
    }

    // Getting geometry
    auto& r_geometry = this->GetGeometry();

    // Getting values
    const array_1d<double, 3> ddisp = r_geometry[1].FastGetSolutionStepValue(DISPLACEMENT)
        - r_geometry[0].FastGetSolutionStepValue(DISPLACEMENT);

    const array_1d<double, 3> drot = r_geometry[1].FastGetSolutionStepValue(ROTATION)
        - r_geometry[0].FastGetSolutionStepValue(ROTATION);

    if constexpr (TDim == 2) {
        current_displacement[0] = -ddisp[0];
        current_displacement[1] = -ddisp[1];
        current_displacement[2] = -drot[2];
        current_displacement[3] = ddisp[0];
        current_displacement[4] = ddisp[1];
        current_displacement[5] = drot[2];
    }
    else if constexpr (TDim == 3) {
        current_displacement[0] = -ddisp[0];
        current_displacement[1] = -ddisp[1];
        current_displacement[2] = -ddisp[2];
        current_displacement[3] = -drot[0];
        current_displacement[4] = -drot[1];
        current_displacement[5] = -drot[2];
        current_displacement[6] = ddisp[0];
        current_displacement[7] = ddisp[1];
        current_displacement[8] = ddisp[2];
        current_displacement[9] = drot[0];
        current_displacement[10] = drot[1];
        current_displacement[11] = drot[2];
    }

    for (std::size_t i = 0; i < msElementSize; ++i) {
        rRightHandSideVector[i] -= elemental_stiffness[i % msLocalSize] * current_displacement[i];
    }
    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************

template<std::size_t TDim>
int SpringDamperElement<TDim>::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    // Verify that the dofs exist
    for ( std::size_t i = 0; i < this->GetGeometry().size(); i++ ) {
        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        const NodeType& rnode = this->GetGeometry()[i];

        // The displacement terms
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X,rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y,rnode)
        if constexpr (TDim == 3) {
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode)
        }

        // The rotational terms
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ROTATION,rnode)

        if constexpr (TDim == 3) {
            KRATOS_CHECK_DOF_IN_NODE(ROTATION_X, rnode)
            KRATOS_CHECK_DOF_IN_NODE(ROTATION_Y, rnode)
        }
        KRATOS_CHECK_DOF_IN_NODE(ROTATION_Z,rnode)
    }

    return 0;

    KRATOS_CATCH( "Problem in the Check in the SpringDamperElement" )
}


//************************************************************************************
//************************************************************************************
template<std::size_t TDim>
void SpringDamperElement<TDim>::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
}

template<std::size_t TDim>
void SpringDamperElement<TDim>::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
}

template class SpringDamperElement<2>;
template class SpringDamperElement<3>;


} // Namespace Kratos


