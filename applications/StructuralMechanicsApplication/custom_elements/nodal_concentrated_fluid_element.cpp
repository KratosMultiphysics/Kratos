// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Navaneeth K Narayanan
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/checks.h"
#include "structural_mechanics_application_variables.h"
#include "custom_elements/nodal_concentrated_fluid_element.h"
#include "utilities/integration_utilities.h"
#include "geometries/triangle_3d_3.h"
#include "utilities/divide_geometry.h"
#include "utilities/divide_triangle_2d_3.h"

namespace Kratos
{
//***********************DEFAULT CONSTRUCTOR******************************************
//************************************************************************************

NodalConcentratedFluidElement::NodalConcentratedFluidElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)

{
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

NodalConcentratedFluidElement::NodalConcentratedFluidElement(
    IndexType NewId, GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)

{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

NodalConcentratedFluidElement::NodalConcentratedFluidElement(NodalConcentratedFluidElement const &rOther)
    : Element(rOther)
{
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

NodalConcentratedFluidElement &NodalConcentratedFluidElement::operator=(NodalConcentratedFluidElement const &rOther)
{
    //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

    Element::operator=(rOther);

    return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer NodalConcentratedFluidElement::Create(
    IndexType NewId,
    NodesArrayType const &rThisNodes,
    PropertiesType::Pointer pProperties) const
{
    //NEEDED TO CREATE AN ELEMENT
    return Kratos::make_shared<NodalConcentratedFluidElement>(NewId, GetGeometry().Create(rThisNodes), pProperties);
}

//************************************************************************************
//************************************************************************************

Element::Pointer NodalConcentratedFluidElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    //NEEDED TO CREATE AN ELEMENT
    return Kratos::make_shared<NodalConcentratedFluidElement>(NewId, pGeom, pProperties);
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer NodalConcentratedFluidElement::Clone(
    IndexType NewId,
    NodesArrayType const &rThisNodes) const
{
    //YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
    //ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

    NodalConcentratedFluidElement new_element(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

    return Kratos::make_shared<NodalConcentratedFluidElement>(new_element);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

NodalConcentratedFluidElement::~NodalConcentratedFluidElement()
{
}

//************* GETTING METHODS
//************************************************************************************
//************************************************************************************

void NodalConcentratedFluidElement::GetDofList(
    DofsVectorType &rElementalDofList,
    ProcessInfo &rCurrentProcessInfo)
{
    //NEEDED TO DEFINE THE DOFS OF THE ELEMENT
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rElementalDofList.resize(0);

    SizeType number_of_master_nodes = mMasterNodes.size();

    for (IndexType i_node = 0; i_node < number_of_master_nodes; ++i_node)
    {
        Node<3> &node = mMasterNodes[i_node];

        rElementalDofList.push_back(node.pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back(node.pGetDof(DISPLACEMENT_Y));

        if (dimension == 3)
            rElementalDofList.push_back(node.pGetDof(DISPLACEMENT_Z));
    }
}

//************************************************************************************
//************************************************************************************

void NodalConcentratedFluidElement::EquationIdVector(
    EquationIdVectorType &rResult,
    ProcessInfo &rCurrentProcessInfo)
{
    //NEEDED TO DEFINE GLOBAL IDS FOR THE CORRECT ASSEMBLY
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    SizeType number_of_master_nodes = mMasterNodes.size();
    SizeType sys_size = number_of_master_nodes * dimension;

    if (rResult.size() != sys_size)
        rResult.resize(sys_size, false);

    for (IndexType i_node = 0; i_node < number_of_master_nodes; ++i_node)
    {
        Node<3> &node = mMasterNodes[i_node];
        IndexType index = i_node * dimension;
        rResult[index] = node.GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = node.GetDof(DISPLACEMENT_Y).EquationId();
        if (dimension == 3)
            rResult[index + 2] = node.GetDof(DISPLACEMENT_Z).EquationId();
    }
}

//*********************************DISPLACEMENT***************************************
//************************************************************************************

void NodalConcentratedFluidElement::GetValuesVector(Vector &rValues, int Step)
{
    //GIVES THE VECTOR WITH THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT DISPLACEMENTS)
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    SizeType number_of_master_nodes = mMasterNodes.size();
    SizeType sys_size = number_of_master_nodes * dimension;

    if (rValues.size() != sys_size)
        rValues.resize(sys_size, false);

    for (IndexType i_node = 0; i_node < number_of_master_nodes; ++i_node)
    {
        Node<3> &node = mMasterNodes[i_node];
        IndexType index = i_node * dimension;
        const array_1d<double, 3> &displacement = node.FastGetSolutionStepValue(DISPLACEMENT, Step);

        for (IndexType i_dim = 0; i_node < dimension; ++i_dim)
        {
            rValues[index + i_dim] = displacement[i_dim];
        }
    }
}

//************************************VELOCITY****************************************
//************************************************************************************

void NodalConcentratedFluidElement::GetFirstDerivativesVector(Vector &rValues, int Step)
{
    //GIVES THE VECTOR WITH THE TIME DERIVATIVE OF THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT VELOCITIES)
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    SizeType number_of_master_nodes = mMasterNodes.size();
    SizeType sys_size = number_of_master_nodes * dimension;

    if (rValues.size() != sys_size)
        rValues.resize(sys_size, false);

    for (IndexType i_node = 0; i_node < number_of_master_nodes; ++i_node)
    {
        Node<3> &node = mMasterNodes[i_node];
        IndexType index = i_node * dimension;
        const array_1d<double, 3> &velocity = node.FastGetSolutionStepValue(VELOCITY, Step);

        for (IndexType i_dim = 0; i_node < dimension; ++i_dim)
        {
            rValues[index + i_dim] = velocity[i_dim];
        }
    }
}

//*********************************ACCELERATION***************************************
//************************************************************************************

void NodalConcentratedFluidElement::GetSecondDerivativesVector(Vector &rValues, int Step)
{
    //GIVES THE VECTOR WITH THE TIME SECOND DERIVATIVE OF THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT ACCELERATIONS)
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    SizeType number_of_master_nodes = mMasterNodes.size();
    SizeType sys_size = number_of_master_nodes * dimension;

    if (rValues.size() != sys_size)
        rValues.resize(sys_size, false);

    for (IndexType i_node = 0; i_node < number_of_master_nodes; ++i_node)
    {
        Node<3> &node = mMasterNodes[i_node];
        IndexType index = i_node * dimension;
        const array_1d<double, 3> &acceleration = node.FastGetSolutionStepValue(ACCELERATION, Step);

        for (IndexType i_dim = 0; i_node < dimension; ++i_dim)
        {
            rValues[index + i_dim] = acceleration[i_dim];
        }
    }
}

//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************

void NodalConcentratedFluidElement::Initialize()
{

    /*     Vector normal = pGetProperties()->GetValue(FREE_SURFACE_NORMAL);

    double norm_normal = norm_2(normal);
    if (norm_normal > std::numeric_limits<double>::epsilon())
        normal = normal / norm_normal;
    else
        normal *= 0;

    unsigned int nr_of_non_zeros = 0;

    for (IndexType i = 0; i < 3; ++i)
    {
        if (fabs(normal[i]) > std::numeric_limits<double>::epsilon())
            nr_of_non_zeros++;
    }

    if (nr_of_non_zeros > 1)
        KRATOS_ERROR << "The normal should be either in X, Y or Z direction" << std::endl;

    Vector position_vector = ZeroVector(3);

    for (IndexType i = 0; i < 3; ++i)
        position_vector[i] = i;

    mNormalDirection = MathUtils<double>::Dot(normal, position_vector);

    double target_volume = pGetProperties()->GetValue(FLUID_VOLUME);

    this->GetValue(NODAL_MASS) = target_volume * pGetProperties()->GetValue(SPECIFIC_WEIGHT) / 9.81; */

    pGetProperties()->GetValue(FREE_SURFACE_CENTRE) = GetGeometry()[0].Coordinates();
    mMasterConditions = GetValue(MASTER_CONDITIONS);
    mMasterNodes = GetValue(MASTER_NODES);
}

////************************************************************************************
////************************************************************************************

void NodalConcentratedFluidElement::InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    pGetProperties()->GetValue(FREE_SURFACE_CENTRE) = GetGeometry()[0].Coordinates();

    KRATOS_CATCH("");
}

void NodalConcentratedFluidElement::InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    Vector centre = GetGeometry()[0].Coordinates();

    pGetProperties()->GetValue(FREE_SURFACE_CENTRE) = centre;
    double radius = pGetProperties()->GetValue(FREE_SURFACE_RADIUS);
    Vector normal = pGetProperties()->GetValue(FREE_SURFACE_NORMAL);
    double target_volume = pGetProperties()->GetValue(FLUID_VOLUME);

    VolumeCalculationUnderPlaneUtility VolumeCalculator(centre, radius, normal);
    VolumeCalculator.UpdatePositionOfPlaneBasedOnTargetVolumeForConditions(mMasterConditions, mMasterNodes, target_volume);

    //pGetProperties()->GetValue(CURRENT_FLUID_VOLUME) = VolumeCalculator.GetVolume();
    pGetProperties()->GetValue(FREE_SURFACE_AREA) = VolumeCalculator.GetIntersectedArea();
    pGetProperties()->GetValue(FREE_SURFACE_CENTRE) = VolumeCalculator.GetCentre();
    GetGeometry()[0].Coordinates() = VolumeCalculator.GetCentre();

    Vector displacement = GetGeometry()[0].Coordinates() - GetGeometry()[0].GetInitialPosition();

    GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT) = displacement;

    std::cout << " Displacement " << GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT) << std::endl;

    KRATOS_CATCH("");
}

////************************************************************************************
////************************************************************************************

void NodalConcentratedFluidElement::FinalizeNonLinearIteration(ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}

////************************************************************************************
////************************************************************************************

void NodalConcentratedFluidElement::FinalizeSolutionStep(ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}

//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

void NodalConcentratedFluidElement::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo)
{

    KRATOS_TRY;

    /* Calculate elemental system */

    // Compute RHS (RHS = rRightHandSideVector = Fext - Fint)
    this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    // Compute LHS
    this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

    KRATOS_CATCH("");
}

//***********************************************************************************
//***********************************************************************************

void NodalConcentratedFluidElement::CalculateRightHandSide(VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo)
{

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the RHS
    SizeType number_of_master_nodes = mMasterNodes.size();
    SizeType sys_size = number_of_master_nodes * dimension;

    if (rRightHandSideVector.size() != sys_size)
        rRightHandSideVector.resize(sys_size, false);

    rRightHandSideVector = ZeroVector(sys_size); //resetting RHS

    /* double target_volume = pGetProperties()->GetValue(FLUID_VOLUME);
    double current_volume = pGetProperties()->GetValue(CURRENT_FLUID_VOLUME);

    // Calculate forces from all the conditions
    if (mIsDynamic)
    {

        Vector force_vector = ZeroVector(3);
        SizeType number_of_master_conditions = mMasterConditions.size();

        for (IndexType i_cond = 0; i_cond < number_of_master_conditions; ++i_cond)
        {
            Condition &cond = mMasterConditions[i_cond];
            SizeType number_of_nodes = cond.GetGeometry().size();
            Vector condition_nodal_forces = ZeroVector(number_of_nodes * 3);

            cond.CalculateRightHandSide(condition_nodal_forces, rCurrentProcessInfo);

            for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node)
            {

                IndexType i = i_node * 3;

                force_vector[0] += condition_nodal_forces[i];
                force_vector[1] += condition_nodal_forces[i + 1];
                force_vector[2] += condition_nodal_forces[i + 2];
            }
        }

        std::cout << "Total vertical force (Weight) " << -force_vector[mNormalDirection] << std::endl;

        for (IndexType dim = 0; dim < 3; ++dim)
        {
            if (dim != mNormalDirection) // adding forces other than the normal direction
                rRightHandSideVector[dim] -= force_vector[dim];
        }
    }
    else
        rRightHandSideVector[mNormalDirection] += (target_volume - current_volume) * mScalingParameter; // For volume conservation in normal direction */
}

//***********************************************************************************
//***********************************************************************************

void NodalConcentratedFluidElement::CalculateLeftHandSide(MatrixType &rLeftHandSideMatrix, ProcessInfo &rCurrentProcessInfo)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    SizeType number_of_master_nodes = mMasterNodes.size();
    SizeType sys_size = number_of_master_nodes * dimension;

    std::cout << "nodal_element sys size ::  " << sys_size << std::endl;
    if (rLeftHandSideMatrix.size1() != sys_size)
        rLeftHandSideMatrix.resize(sys_size, sys_size, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(sys_size, sys_size); //resetting LHS

    /*  // We add the nodal stiffness
    double nodal_area = pGetProperties()->GetValue(FREE_SURFACE_AREA);

    if (!mIsDynamic)
    {
        for (unsigned int j = 0; j < dimension; ++j)
            rLeftHandSideMatrix(j, j) += nodal_area * mScalingParameter; // stiffness included in the other directions to avoid singularity
    }
    else
    {
        rLeftHandSideMatrix(mNormalDirection, mNormalDirection) += nodal_area * mScalingParameter;
    } */

    // Adding volume conserving stiffness on the wetted nodes
    Matrix Kij(dimension, dimension);
    double intersected_area = pGetProperties()->GetValue(FREE_SURFACE_AREA);
    double specific_weight = pGetProperties()->GetValue(SPECIFIC_WEIGHT);
    array_1d<double, 3> nodal_normal_m;
    array_1d<double, 3> nodal_normal_n;
    double coeff;

    for (IndexType m = 0; m < number_of_master_nodes; m++)
    {

        const IndexType RowIndex = m * dimension;
        Node<3> &node_m = mMasterNodes[m];
        nodal_normal_m = node_m.FastGetSolutionStepValue(NORMAL);

        for (IndexType n = 0; n < number_of_master_nodes; n++)

        {

            const IndexType ColIndex = n * dimension;
            Node<3> &node_n = mMasterNodes[n];

            nodal_normal_n = node_n.FastGetSolutionStepValue(NORMAL);

            DyadicProduct(Kij, nodal_normal_m, nodal_normal_n);

            coeff = -specific_weight / intersected_area;

            Kij *= coeff;

            MathUtils<double>::SubtractMatrix(rLeftHandSideMatrix, Kij, RowIndex, ColIndex);
        }
    }
}

void NodalConcentratedFluidElement::DyadicProduct(Matrix &M, const array_1d<double, 3> &U, const array_1d<double, 3> &V)
{
    KRATOS_TRY;

    if (U.size() == V.size())
    {

        for (IndexType i = 0; i < U.size(); ++i)
        {

            for (IndexType j = 0; j < U.size(); ++j)
            {

                M(i, j) = U[i] * V[j];
            }
        }
    }
    else
        KRATOS_ERROR << "The size of the two input vectors don't match" << std::endl;

    KRATOS_CATCH("")
}

//*************************COMPUTE DELTA POSITION*************************************
//************************************************************************************

/* Matrix &NodalConcentratedFluidElement::CalculateDeltaPosition(Matrix &rDeltaPosition)
{
    KRATOS_TRY;

    //KRATOS NODAL CURRENT POSITION (X = X0 + DISPLACEMENT_X) IS ALWAYS COMPUTED
    GeometryType &geom = GetGeometry();
    const unsigned int dimension = geom.WorkingSpaceDimension();

    rDeltaPosition = ZeroMatrix(1, dimension);

    rDeltaPosition(0, 0) = GetGeometry()[0].X() - GetGeometry()[0].X0();
    rDeltaPosition(0, 1) = GetGeometry()[0].Y() - GetGeometry()[0].Y0();
    if (dimension == 3)
        rDeltaPosition(0, 2) = GetGeometry()[0].Z() - GetGeometry()[0].Z0();

    return rDeltaPosition;

    KRATOS_CATCH("");
} */

//************************************************************************************
//************************************************************************************

void NodalConcentratedFluidElement::CalculateMassMatrix(MatrixType &rMassMatrix, ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    //lumped
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int system_size = dimension;

    if (rMassMatrix.size1() != system_size)
        rMassMatrix.resize(system_size, system_size, false);

    rMassMatrix = ZeroMatrix(system_size, system_size);

    // We get the reference
    //const auto &rconst_this = *this;

    /* // Get the nodal mass
    const double nodal_mass = rconst_this.GetValue(NODAL_MASS);

    for (unsigned int j = 0; j < dimension; ++j)
    {
        if (j != mNormalDirection)
            rMassMatrix(j, j) += nodal_mass; // Dynamics in normal direction is absent
    } */

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************

void NodalConcentratedFluidElement::CalculateDampingMatrix(
    MatrixType &rDampingMatrix,
    ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    //0.-Initialize the DampingMatrix:
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the LHS
    const unsigned int system_size = dimension;

    rDampingMatrix = ZeroMatrix(system_size, system_size);

    /* const array_1d<double, 3> &nodal_damping_ratio = this->GetValue(NODAL_DAMPING_RATIO);
    for (unsigned int j = 0; j < dimension; ++j)
        if (j != mNormalDirection)
            rDampingMatrix(j, j) += nodal_damping_ratio[j]; // Dynamics in normal direction is absent */

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************

int NodalConcentratedFluidElement::Check(const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(NODAL_MASS)
    KRATOS_CHECK_VARIABLE_KEY(NODAL_AREA)
    KRATOS_CHECK_VARIABLE_KEY(NODAL_DAMPING_RATIO)
    KRATOS_CHECK_VARIABLE_KEY(VOLUME_ACCELERATION)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for (std::size_t i = 0; i < this->GetGeometry().size(); ++i)
    {
        Node<3> &rnode = this->GetGeometry()[i];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, rnode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VOLUME_ACCELERATION, rnode)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode)
    }

    return 0;

    KRATOS_CATCH("Problem in the Check in the NodalConcentratedFluidElement")
}

//************************************************************************************
//************************************************************************************

void NodalConcentratedFluidElement::save(Serializer &rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
}

void NodalConcentratedFluidElement::load(Serializer &rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
}

} // Namespace Kratos
