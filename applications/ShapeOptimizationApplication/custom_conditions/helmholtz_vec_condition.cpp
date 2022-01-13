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
#include "custom_conditions/helmholtz_vec_condition.h"
#include "includes/variables.h"
#include "includes/checks.h"
#include "utilities/atomic_utilities.h"

namespace Kratos
{
HelmholtzVecCondition::HelmholtzVecCondition(HelmholtzVecCondition const& rOther)
    : BaseType(rOther)
{
}

/***********************************************************************************/
/***********************************************************************************/

HelmholtzVecCondition& HelmholtzVecCondition::operator=(HelmholtzVecCondition const& rOther)
{
    //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

    Condition::operator=(rOther);

    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer HelmholtzVecCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<HelmholtzVecCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer HelmholtzVecCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<HelmholtzVecCondition>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer HelmholtzVecCondition::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<HelmholtzVecCondition>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void HelmholtzVecCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    if (rResult.size() != dimension * number_of_nodes)
        rResult.resize(dimension * number_of_nodes,false);

    const SizeType pos = this->GetGeometry()[0].GetDofPosition(HELMHOLTZ_VARS_X);

    if(dimension == 2) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * 2;
            rResult[index] = GetGeometry()[i].GetDof(HELMHOLTZ_VARS_X,pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(HELMHOLTZ_VARS_Y,pos+1).EquationId();
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * 3;
            rResult[index] = GetGeometry()[i].GetDof(HELMHOLTZ_VARS_X,pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(HELMHOLTZ_VARS_Y,pos+1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(HELMHOLTZ_VARS_Z,pos+2).EquationId();
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/
void HelmholtzVecCondition::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(dimension*number_of_nodes);

    if(dimension == 2) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(HELMHOLTZ_VARS_X));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(HELMHOLTZ_VARS_Y));
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(HELMHOLTZ_VARS_X));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(HELMHOLTZ_VARS_Y));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(HELMHOLTZ_VARS_Z));
        }
    }

    KRATOS_CATCH("")
}

//******************************************************************************
//******************************************************************************
void HelmholtzVecCondition::GetValuesVector(VectorType &rValues,
                                            int Step) const {
  const GeometryType &rgeom = this->GetGeometry();
  const SizeType num_nodes = rgeom.PointsNumber();
  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
  const unsigned int local_size = num_nodes * dimension;

  if (rValues.size() != local_size)
    rValues.resize(local_size, false);

  if (dimension == 2) {
    SizeType index = 0;
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VARS_X, Step);
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VARS_Y, Step);
    }
  } else if (dimension == 3) {
    SizeType index = 0;
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VARS_X, Step);
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VARS_Y, Step);
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VARS_Z, Step);
    }
  }
}

/***********************************************************************************/
/***********************************************************************************/

void HelmholtzVecCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

void HelmholtzVecCondition::Calculate(const Variable<Matrix>& rVariable, Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == HELMHOLTZ_MASS_MATRIX)
        CalculateSurfaceMassMatrix(rOutput,rCurrentProcessInfo);

}
/***********************************************************************************/
/***********************************************************************************/
void HelmholtzVecCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * dimension;

    if ( rLeftHandSideMatrix.size1() != mat_size )
        rLeftHandSideMatrix.resize( mat_size, mat_size, false );

    noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS

    // Resizing as needed the RHS 
    if ( rRightHandSideVector.size() != mat_size )
        rRightHandSideVector.resize( mat_size, false );

    rRightHandSideVector = ZeroVector( mat_size ); //resetting RHS

    // Check that parents have been computed
    // These are required to retrieve the material properties and the viscous stress
    auto& parentElement = this->GetValue(NEIGHBOUR_ELEMENTS);
    KRATOS_ERROR_IF(parentElement.size() > 1) << "A condition was assigned more than one parent element." << std::endl;
    KRATOS_ERROR_IF(parentElement.size() == 0) << "A condition was NOT assigned a parent element." << std::endl;

    MatrixType M;
    CalculateSurfaceMassMatrix(M,rCurrentProcessInfo);  
    MatrixType A;
    CalculateSurfaceStiffnessMatrix(A,rCurrentProcessInfo); 

    MatrixType K = A + M;     


    const unsigned int number_of_points = r_geometry.size();
    Vector nodal_vals(number_of_points*dimension);
    for(unsigned int node_element = 0; node_element<number_of_points; node_element++)
    {
        const VectorType &source = r_geometry[node_element].FastGetSolutionStepValue(HELMHOLTZ_SOURCE);
        auto node_weight = r_geometry[node_element].GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
        if(rCurrentProcessInfo[COMPUTE_CONTROL_POINTS])
            node_weight = 1.0;
        nodal_vals[3 * node_element + 0] = source[0]/node_weight;
        nodal_vals[3 * node_element + 1] = source[1]/node_weight;
        nodal_vals[3 * node_element + 2] = source[2]/node_weight;
    }

    if(rCurrentProcessInfo[COMPUTE_CONTROL_POINTS]){
        noalias(rLeftHandSideMatrix) += M;
        noalias(rRightHandSideVector) += prod(K,nodal_vals);
        //apply drichlet BC
        Vector temp(number_of_points*3);
        for (SizeType iNode = 0; iNode < number_of_points; ++iNode) {
            const VectorType &vars = r_geometry[iNode].FastGetSolutionStepValue(HELMHOLTZ_VARS,0);
            temp[3*iNode] = vars[0];
            temp[3*iNode+1] = vars[1];
            temp[3*iNode+2] = vars[2];
        }
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);     
    }            
    else{
        noalias(rLeftHandSideMatrix) += K;
        // noalias(rRightHandSideVector) += nodal_vals;
    } 

}

/***********************************************************************************/
/***********************************************************************************/

int HelmholtzVecCondition::Check( const ProcessInfo& rCurrentProcessInfo ) const
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

void HelmholtzVecCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
}

/***********************************************************************************/
/***********************************************************************************/

void HelmholtzVecCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
}

/***********************************************************************************/
/***********************************************************************************/

void HelmholtzVecCondition::CalculateNormal(VectorType & r_n) const
{
    const auto& r_cond_geom = GetGeometry();

    array_1d<double,3> v1,v2;
    v1[0] = r_cond_geom[1].X() - r_cond_geom[0].X();
    v1[1] = r_cond_geom[1].Y() - r_cond_geom[0].Y();
    v1[2] = r_cond_geom[1].Z() - r_cond_geom[0].Z();

    v2[0] = r_cond_geom[2].X() - r_cond_geom[0].X();
    v2[1] = r_cond_geom[2].Y() - r_cond_geom[0].Y();
    v2[2] = r_cond_geom[2].Z() - r_cond_geom[0].Z();

    r_n.resize(3);
    MathUtils<double>::CrossProduct(r_n,v1,v2);
    double norm = MathUtils<double>::Norm3(r_n);
    r_n /= norm;
}

/***********************************************************************************/
/***********************************************************************************/
void HelmholtzVecCondition::CalculateSurfaceMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_cond_geom = GetGeometry();
    SizeType dimension = r_cond_geom.WorkingSpaceDimension();
    SizeType number_of_nodes = r_cond_geom.size();
    SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rMassMatrix.size1() != mat_size || rMassMatrix.size2() != mat_size)
        rMassMatrix.resize( mat_size, mat_size, false );
    rMassMatrix = ZeroMatrix( mat_size, mat_size );

    IntegrationMethod integration_method = GeometryData::IntegrationMethod::GI_GAUSS_4;
    const GeometryType::IntegrationPointsArrayType& integration_points = r_cond_geom.IntegrationPoints( integration_method );
    MatrixType Ncontainer;
    GetParentElementShapeFunctionsValues(Ncontainer,integration_method,rCurrentProcessInfo);

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = r_cond_geom.IntegrationPoints(integration_method);
    const unsigned int NumGauss = IntegrationPoints.size();
    Vector GaussPtsJDet = ZeroVector(NumGauss);
    r_cond_geom.DeterminantOfJacobian(GaussPtsJDet, integration_method);

    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {

        const double integration_weight = integration_points[point_number].Weight() * GaussPtsJDet[point_number];
        const Vector& rN = row(Ncontainer,point_number);

        for ( IndexType i = 0; i < number_of_nodes; ++i ) {
            const SizeType index_i = i * dimension;

            for ( IndexType j = 0; j < number_of_nodes; ++j ) {
                const SizeType index_j = j * dimension;
                const double NiNj_weight = rN[i] * rN[j] * integration_weight;

                for ( IndexType k = 0; k < dimension; ++k )
                    rMassMatrix( index_i + k, index_j + k ) += NiNj_weight;
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/
void HelmholtzVecCondition::GetParentElementShapeFunctionsValues(
    MatrixType& rNMatrix,
    const IntegrationMethod& rIntegrationMethod,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_cond_geom = GetGeometry();
    SizeType cond_number_of_nodes = r_cond_geom.size();
    const GeometryType::IntegrationPointsArrayType& cond_integration_points = r_cond_geom.IntegrationPoints( rIntegrationMethod );
    SizeType mat_size1 = cond_integration_points.size();
    SizeType mat_size2 = cond_number_of_nodes;

    rNMatrix.resize( mat_size1, mat_size2, false );
    rNMatrix = ZeroMatrix( mat_size1, mat_size2 );

    auto& parentElement = this->GetValue(NEIGHBOUR_ELEMENTS);
    const auto& r_elem_geom = parentElement[0].GetGeometry();

    for ( IndexType point_number = 0; point_number < cond_integration_points.size(); ++point_number ) {    
        Point cond_gp_local_pt = Point(cond_integration_points[point_number].Coordinates());
        Point cond_gp_global_pt;
        r_cond_geom.GlobalCoordinates(cond_gp_global_pt,cond_gp_local_pt);
        Point elem_cond__gp_local_pt;
        r_elem_geom.PointLocalCoordinates(elem_cond__gp_local_pt,cond_gp_global_pt);

        for(IndexType cond_point_number = 0; cond_point_number < r_cond_geom.size(); ++cond_point_number ){
            for ( IndexType elem_point_number = 0; elem_point_number < r_elem_geom.size(); ++elem_point_number ){
                if(r_cond_geom[cond_point_number].Id()==r_elem_geom[elem_point_number].Id()){
                    double elem_shape_func_value = r_elem_geom.ShapeFunctionValue(elem_point_number,elem_cond__gp_local_pt);
                    rNMatrix(point_number,cond_point_number) = elem_shape_func_value;
                }
            }
        }
    }  

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void HelmholtzVecCondition::CalculateSurfaceStiffnessMatrix(
    MatrixType& rStiffnessMatrix,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_cond_prop = GetProperties();

    // Checking radius
    KRATOS_ERROR_IF_NOT(r_cond_prop.Has(HELMHOLTZ_RADIUS)) << "HELMHOLTZ_RADIUS has to be provided for the calculations of the HelmholtzVecCondition!" << std::endl;

    const auto& r_cond_geom = GetGeometry();
    SizeType dimension = r_cond_geom.WorkingSpaceDimension();
    SizeType number_of_nodes = r_cond_geom.size();
    SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rStiffnessMatrix.size1() != mat_size || rStiffnessMatrix.size2() != mat_size)
        rStiffnessMatrix.resize( mat_size, mat_size, false );
    rStiffnessMatrix = ZeroMatrix( mat_size, mat_size );

    IntegrationMethod integration_method = GeometryData::IntegrationMethod::GI_GAUSS_4;
    const GeometryType::IntegrationPointsArrayType& integration_points = r_cond_geom.IntegrationPoints(integration_method);
    const unsigned int NumGauss = integration_points.size();
    Vector GaussPtsJDet = ZeroVector(NumGauss);
    r_cond_geom.DeterminantOfJacobian(GaussPtsJDet, integration_method);

    // get the normal
    VectorType n_surf;
    CalculateNormal(n_surf);
    MatrixType id_matrix = IdentityMatrix(dimension,dimension);
    MatrixType tangent_projection_matrix = id_matrix - outer_prod(n_surf, n_surf);
    MatrixType A_dirc = ZeroMatrix(number_of_nodes,number_of_nodes);
    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {

        const double integration_weight = integration_points[point_number].Weight() * GaussPtsJDet[point_number];
        MatrixType DN_DX;
        GetParentElementShapeFunctionsGlobalGradients(DN_DX,point_number,integration_method,rCurrentProcessInfo);
        MatrixType DN_DX_t = prod(DN_DX,tangent_projection_matrix);
        MatrixType D = SetAndModifyConstitutiveLaw(GaussPtsJDet[point_number]);
        MatrixType B = CalculateBMatrix(DN_DX_t);
        const double r_helmholtz = r_cond_prop[HELMHOLTZ_RADIUS];
        noalias(A_dirc) += integration_weight * r_helmholtz * r_helmholtz * prod(trans(B), Matrix(prod(D, B)));
    }

    //contruct the stifness matrix in all dims
    for(IndexType i=0;i<number_of_nodes;i++)
        for(IndexType j=0;j<dimension;j++)
            for(IndexType k=0;k<number_of_nodes;k++)
                rStiffnessMatrix(dimension*i+j,dimension*k+j) = A_dirc(i,k);


    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/
void HelmholtzVecCondition::GetParentElementShapeFunctionsGlobalGradients(
    MatrixType& rDN_DX,
    const IndexType PointNumber,
    const IntegrationMethod& rIntegrationMethod,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_cond_geom = GetGeometry();
    SizeType cond_number_of_nodes = r_cond_geom.size();
    SizeType dimension = r_cond_geom.WorkingSpaceDimension();
    const GeometryType::IntegrationPointsArrayType& cond_integration_points = r_cond_geom.IntegrationPoints( rIntegrationMethod );
    SizeType mat_size1 = cond_number_of_nodes;
    SizeType mat_size2 = dimension;

    rDN_DX.resize( mat_size1, mat_size2, false );
    rDN_DX = ZeroMatrix( mat_size1, mat_size2 );

    auto& parentElement = this->GetValue(NEIGHBOUR_ELEMENTS);
    const auto& r_elem_geom = parentElement[0].GetGeometry();
    
    Point cond_gp_local_pt = Point(cond_integration_points[PointNumber].Coordinates());
    Point cond_gp_global_pt;
    r_cond_geom.GlobalCoordinates(cond_gp_global_pt,cond_gp_local_pt);
    Point elem_cond__gp_local_pt;
    r_elem_geom.PointLocalCoordinates(elem_cond__gp_local_pt,cond_gp_global_pt);

    MatrixType DN_De;
    r_elem_geom.ShapeFunctionsLocalGradients(DN_De,elem_cond__gp_local_pt);

    MatrixType InvJ0;
    r_elem_geom.InverseOfJacobian(InvJ0,elem_cond__gp_local_pt);

    MatrixType elem_DN_DX = prod(DN_De,InvJ0);

    for(IndexType cond_point_number = 0; cond_point_number < r_cond_geom.size(); ++cond_point_number ){
        for ( IndexType elem_point_number = 0; elem_point_number < r_elem_geom.size(); ++elem_point_number ){
            if(r_cond_geom[cond_point_number].Id()==r_elem_geom[elem_point_number].Id()){
                for(IndexType k = 0; k<dimension; k++)
                    rDN_DX(cond_point_number,k) = elem_DN_DX(elem_point_number,k);
            }
        }
    }    

    KRATOS_CATCH("");
}
//******************************************************************************
//******************************************************************************
HelmholtzVecCondition::MatrixType
HelmholtzVecCondition::CalculateBMatrix(const MatrixType& rDN_DX) const {
  KRATOS_TRY;

  const GeometryType &rgeom = this->GetGeometry();
  SizeType dimension = rgeom.WorkingSpaceDimension();
  const SizeType num_nodes = rgeom.PointsNumber();

  MatrixType B;

  if (dimension == 2) {
    B = ZeroMatrix(3, num_nodes * 2);

    SizeType index = 0;
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      B(0, index + 0) = rDN_DX(i_node, 0);
      B(0, index + 1) = 0.0;
      B(1, index + 0) = 0.0;
      B(1, index + 1) = rDN_DX(i_node, 1);
      B(2, index + 0) = rDN_DX(i_node, 1);
      B(2, index + 1) = rDN_DX(i_node, 0);
      index += 2;
    }
  }

  else if (dimension == 3) {
    B = ZeroMatrix(6, num_nodes * 3);

    SizeType index = 0;
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      B(0, index + 0) = rDN_DX(i_node, 0);
      B(1, index + 1) = rDN_DX(i_node, 1);
      B(2, index + 2) = rDN_DX(i_node, 2);
      B(3, index + 0) = rDN_DX(i_node, 1);
      B(3, index + 1) = rDN_DX(i_node, 0);
      B(4, index + 1) = rDN_DX(i_node, 2);
      B(4, index + 2) = rDN_DX(i_node, 1);
      B(5, index + 0) = rDN_DX(i_node, 2);
      B(5, index + 2) = rDN_DX(i_node, 0);
      index += 3;
    }
  }

  return B;

  KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
HelmholtzVecCondition::MatrixType
HelmholtzVecCondition::SetAndModifyConstitutiveLaw(const double detJ0) const {
  KRATOS_TRY;

  const GeometryType &rgeom = this->GetGeometry();
  SizeType dimension = rgeom.WorkingSpaceDimension();

  // Stiffening of elements using Jacobian determinants and exponent between
  // 0.0 and 2.0
  const double factor =
      100;               // Factor influences how far the HELMHOLTZ_VARS spreads
                         // into the fluid mesh
  const double xi = 1.5; // 1.5 Exponent influences stiffening of smaller
                         // elements; 0 = no stiffening
  const double quotient = factor / detJ0;
  const double weighting_factor = 1.0;
  const double poisson_coefficient = this->pGetProperties()->Has(HELMHOLTZ_POISSON_RATIO)
    ? this->pGetProperties()->GetValue(HELMHOLTZ_POISSON_RATIO) : 0.3;

  // The ratio between lambda and mu affects relative stiffening against
  // volume or shape change.
  const double lambda =
      weighting_factor * poisson_coefficient /
      ((1 + poisson_coefficient) * (1 - 2 * poisson_coefficient));
  const double mu = weighting_factor / (2 * (1 + poisson_coefficient));

  MatrixType constitutive_matrix;

  // stress = lambda*tr(strain tensor)*I + 2*mu*(strain tensor).
  if (dimension == 2) {
    constitutive_matrix = ZeroMatrix(3, 3);
    constitutive_matrix(0, 0) = lambda + 2 * mu;
    constitutive_matrix(1, 1) = constitutive_matrix(0, 0);
    constitutive_matrix(2, 2) = mu;
    constitutive_matrix(0, 1) = lambda;
    constitutive_matrix(1, 0) = lambda;
  }

  else if (dimension == 3) {
    constitutive_matrix = ZeroMatrix(6, 6);
    constitutive_matrix(0, 0) = lambda + 2 * mu;
    constitutive_matrix(1, 1) = constitutive_matrix(0, 0);
    constitutive_matrix(2, 2) = constitutive_matrix(0, 0);
    constitutive_matrix(3, 3) = mu;
    constitutive_matrix(4, 4) = mu;
    constitutive_matrix(5, 5) = mu;
    constitutive_matrix(0, 1) = lambda;
    constitutive_matrix(1, 0) = lambda;
    constitutive_matrix(0, 2) = lambda;
    constitutive_matrix(2, 0) = lambda;
    constitutive_matrix(1, 2) = lambda;
    constitutive_matrix(2, 1) = lambda;
  }

  return constitutive_matrix;

  KRATOS_CATCH("");
}

} // Namespace Kratos
