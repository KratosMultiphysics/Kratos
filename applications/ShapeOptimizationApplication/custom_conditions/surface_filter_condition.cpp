// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl
//
// ==============================================================================
// System includes

// External includes

// Project includes
#include "custom_conditions/surface_filter_condition.h"
#include "../StructuralMechanicsApplication/custom_utilities/shellt3_local_coordinate_system.hpp"
#include "../StructuralMechanicsApplication/custom_utilities/shellq4_local_coordinate_system.hpp"
#include "includes/variables.h"
#include "includes/checks.h"
#include "utilities/atomic_utilities.h"

namespace Kratos
{
SurfaceFilterCondition::SurfaceFilterCondition(SurfaceFilterCondition const& rOther)
    : BaseType(rOther)
{
}

/***********************************************************************************/
/***********************************************************************************/

SurfaceFilterCondition& SurfaceFilterCondition::operator=(SurfaceFilterCondition const& rOther)
{
    //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

    Condition::operator=(rOther);

    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer SurfaceFilterCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<SurfaceFilterCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer SurfaceFilterCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<SurfaceFilterCondition>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer SurfaceFilterCondition::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<SurfaceFilterCondition>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SurfaceFilterCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    KRATOS_DEBUG_ERROR_IF_NOT(rCurrentProcessInfo.Has(HELMHOLTZ_DIRECTION))
        << "HELMHOLTZ_DIRECTION not defined in the ProcessInfo!" << std::endl;

    const GeometryType &rgeom = this->GetGeometry();
    const SizeType num_nodes = rgeom.size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if (rResult.size() != num_nodes)
        rResult.resize(num_nodes, false);

    unsigned int pos = this->GetGeometry()[0].GetDofPosition(HELMHOLTZ_VARS_X);

    if (dimension == 2) {
        for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
        if (rCurrentProcessInfo[HELMHOLTZ_DIRECTION] == 1)
            rResult[i_node] =
                rgeom[i_node].GetDof(HELMHOLTZ_VARS_X, pos).EquationId();

        else if (rCurrentProcessInfo[HELMHOLTZ_DIRECTION] == 2)
            rResult[i_node] =
                rgeom[i_node].GetDof(HELMHOLTZ_VARS_Y, pos + 1).EquationId();
        }
    } else {
        for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
        if (rCurrentProcessInfo[HELMHOLTZ_DIRECTION] == 1)
            rResult[i_node] =
                rgeom[i_node].GetDof(HELMHOLTZ_VARS_X, pos).EquationId();
        else if (rCurrentProcessInfo[HELMHOLTZ_DIRECTION] == 2)
            rResult[i_node] =
                rgeom[i_node].GetDof(HELMHOLTZ_VARS_Y, pos + 1).EquationId();
        else if (rCurrentProcessInfo[HELMHOLTZ_DIRECTION] == 3)
            rResult[i_node] =
                rgeom[i_node].GetDof(HELMHOLTZ_VARS_Z, pos + 2).EquationId();
        }
    }
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/
void SurfaceFilterCondition::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF_NOT(rCurrentProcessInfo.Has(HELMHOLTZ_DIRECTION))
        << "HELMHOLTZ_DIRECTION not defined in the ProcessInfo!" << std::endl;


    const GeometryType &rgeom = this->GetGeometry();
    const SizeType num_nodes = rgeom.size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if (rElementalDofList.size() != num_nodes)
        rElementalDofList.resize(num_nodes);

    if (dimension == 2)
        for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
        if (rCurrentProcessInfo[HELMHOLTZ_DIRECTION] == 1)
            rElementalDofList[i_node] = rgeom[i_node].pGetDof(HELMHOLTZ_VARS_X);
        else if (rCurrentProcessInfo[HELMHOLTZ_DIRECTION] == 2)
            rElementalDofList[i_node] = rgeom[i_node].pGetDof(HELMHOLTZ_VARS_Y);
        }
    else
        for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
        if (rCurrentProcessInfo[HELMHOLTZ_DIRECTION] == 1)
            rElementalDofList[i_node] = rgeom[i_node].pGetDof(HELMHOLTZ_VARS_X);
        if (rCurrentProcessInfo[HELMHOLTZ_DIRECTION] == 2)
            rElementalDofList[i_node] = rgeom[i_node].pGetDof(HELMHOLTZ_VARS_Y);
        if (rCurrentProcessInfo[HELMHOLTZ_DIRECTION] == 3)
            rElementalDofList[i_node] = rgeom[i_node].pGetDof(HELMHOLTZ_VARS_Z);
        }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SurfaceFilterCondition::GetValuesVector(
    Vector& rValues,
    int Step
    ) const
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dim = GetGeometry().WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dim;

    if (rValues.size() != mat_size) {
        rValues.resize(mat_size, false);
    }

    for (SizeType i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3 > & r_displacement = GetGeometry()[i].FastGetSolutionStepValue(HELMHOLTZ_VARS, Step);
        SizeType index = i * dim;
        for(SizeType k = 0; k < dim; ++k) {
            rValues[index + k] = r_displacement[k];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SurfaceFilterCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/
void SurfaceFilterCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{

    const auto& r_prop = GetProperties();
    // Checking some reqs
    KRATOS_ERROR_IF_NOT(r_prop.Has(HELMHOLTZ_RADIUS)) << "HELMHOLTZ_RADIUS has to be provided for the calculations of the HelmholtzElement!" << std::endl;

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(HELMHOLTZ_DIRECTION))
        << "HELMHOLTZ_DIRECTION not defined in the ProcessInfo!" << std::endl;   

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(COMPUTE_CONTROL_POINTS))
        << "COMPUTE_CONTROL_POINTS not defined in the ProcessInfo!" << std::endl;  

    const unsigned int component_index = rCurrentProcessInfo[HELMHOLTZ_DIRECTION] - 1; 

    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_points = r_geometry.size();
    const unsigned int dim = r_geometry.WorkingSpaceDimension();


    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != number_of_points)
        rLeftHandSideMatrix.resize(number_of_points,number_of_points,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points,number_of_points); //resetting LHS


    //resizing as needed the RHS
    if(rRightHandSideVector.size() != number_of_points)
        rRightHandSideVector.resize(number_of_points,false);
    noalias(rRightHandSideVector) = ZeroVector(number_of_points); //resetting RHS   

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(r_geometry.GetDefaultIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(r_geometry.GetDefaultIntegrationMethod());


    MatrixType dNxy = ZeroMatrix(3, 2);  /*!< shape function cartesian derivatives */

    //now get the coord transformation object
    ShellT3_LocalCoordinateSystem transformer(r_geometry[0].GetInitialPosition(),
                                              r_geometry[1].GetInitialPosition(),
                                              r_geometry[2].GetInitialPosition());


    const double x12 = transformer.X1() - transformer.X2();
    const double x23 = transformer.X2() - transformer.X3();
    const double x31 = transformer.X3() - transformer.X1();
    const double x21 = -x12;
    const double x32 = -x23;
    const double x13 = -x31;

    const double y12 = transformer.Y1() - transformer.Y2();
    const double y23 = transformer.Y2() - transformer.Y3();
    const double y31 = transformer.Y3() - transformer.Y1();
    const double y21 = -y12;

    const double y13 = -y31;

    const double A = 0.5*(y21*x13 - x21*y13);
    const double A2 = 2.0*A;

    // cartesian derivatives
    dNxy(0, 0) = (y13 - y12) / A2;
    dNxy(0, 1) = (x12 - x13) / A2;
    dNxy(1, 0) = -y13 / A2;
    dNxy(1, 1) = x13 / A2;
    dNxy(2, 0) = y12 / A2;
    dNxy(2, 1) = -x12 / A2;

    const double r_helmholtz = r_prop[HELMHOLTZ_RADIUS];

    MatrixType local_K = A * r_helmholtz * r_helmholtz * prod(dNxy, trans(dNxy));

    /*compute rotation matrix aligned with HELMHOLTZ_DIRECTION */

    //get the orientation matrix 
    const MatrixType& or_mat = transformer.Orientation();
    MatrixType R(3,3);
    for(int i=0;i<3;i++){
        R(i,0) = or_mat(component_index,0); 
        R(i,1) = or_mat(component_index,1);
        R(i,2) = or_mat(component_index,2);
    }
    
    MatrixType temp(3,3);
    noalias(temp) = prod(trans(R), local_K);
    noalias(rLeftHandSideMatrix) = prod(temp, R);

    // std::cout<<"Hi Reza U R calling me"<<std::endl;

    // std::cout<<"local_K : "<<local_K<<std::endl;
    // std::cout<<"R : "<<R<<std::endl;
    // std::cout<<"temp : "<<temp<<std::endl;

    // std::exit(0);
    
    // if(!rCurrentProcessInfo[COMPUTE_CONTROL_POINTS]){
    //    noalias(rLeftHandSideMatrix) = prod(temp, R);
    // }
     

}

/***********************************************************************************/
/***********************************************************************************/

int SurfaceFilterCondition::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    // Base check
    Condition::Check(rCurrentProcessInfo);

    // Check that the condition's nodes contain all required SolutionStepData and Degrees of freedom
    for (const auto& r_node : this->GetGeometry().Points()) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HELMHOLTZ_VARS,r_node)

        KRATOS_CHECK_DOF_IN_NODE(HELMHOLTZ_VARS_X, r_node)
        KRATOS_CHECK_DOF_IN_NODE(HELMHOLTZ_VARS_Y, r_node)
        KRATOS_CHECK_DOF_IN_NODE(HELMHOLTZ_VARS_Z, r_node)
    }

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

void SurfaceFilterCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
}

/***********************************************************************************/
/***********************************************************************************/

void SurfaceFilterCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
}

} // Namespace Kratos
