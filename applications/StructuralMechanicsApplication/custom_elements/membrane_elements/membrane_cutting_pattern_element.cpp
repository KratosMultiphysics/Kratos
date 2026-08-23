// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Klaus B. Sautter
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "utilities/math_utils.h"
#include "membrane_cutting_pattern_element.hpp"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_math_utilities.hpp"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "utilities/integration_utilities.h"
#include "utilities/atomic_utilities.h"


namespace Kratos
{
  
  MembraneCuttingPatternElement::MembraneCuttingPatternElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
  {

  }


  MembraneCuttingPatternElement::MembraneCuttingPatternElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
  {

  }


  Element::Pointer MembraneCuttingPatternElement::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties) const

  {
    return Kratos::make_intrusive< MembraneCuttingPatternElement >(NewId, GetGeometry().Create(rThisNodes), pProperties);
  }


  Element::Pointer MembraneCuttingPatternElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const

  {
    return Kratos::make_intrusive< MembraneCuttingPatternElement >(NewId, pGeom, pProperties);
  }


  void MembraneCuttingPatternElement::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
  {
    KRATOS_TRY;
    
    SizeType num_nodes, local_size;
    SizeType local_index = 0;
    
    num_nodes = GetGeometry().size();
    
    local_size = num_nodes * 3;
    
    const SizeType d_pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);
    
    if (rResult.size() != local_size)
      rResult.resize(local_size, false);
    
      for (SizeType i_node = 0; i_node < num_nodes; ++i_node)
      {
        rResult[local_index++] = this->GetGeometry()[i_node].GetDof(DISPLACEMENT_X, d_pos).EquationId();
        rResult[local_index++] = this->GetGeometry()[i_node].GetDof(DISPLACEMENT_Y, d_pos + 1).EquationId();
        rResult[local_index++] = this->GetGeometry()[i_node].GetDof(DISPLACEMENT_Z, d_pos + 2).EquationId();
      }
      
      KRATOS_CATCH("")
    }



void MembraneCuttingPatternElement::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    SizeType num_nodes, local_size;
    num_nodes = GetGeometry().size();
    local_size = num_nodes * 3;

    if (rElementalDofList.size() != local_size)
        rElementalDofList.resize(local_size);

    SizeType local_index = 0;

    for (SizeType i_node = 0; i_node < num_nodes; ++i_node)
    {
        rElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(DISPLACEMENT_X);
        rElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(DISPLACEMENT_Y);
        rElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(DISPLACEMENT_Z);
    }
}


void MembraneCuttingPatternElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(GetIntegrationMethod());

        //Constitutive Law initialisation
        if ( mConstitutiveLawVector.size() != integration_points.size() )
            mConstitutiveLawVector.resize( integration_points.size() );

        if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
            const GeometryType& r_geometry = GetGeometry();
            const Properties& r_properties = GetProperties();
            const auto& N_values = r_geometry.ShapeFunctionsValues(GetIntegrationMethod());
            for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
                mConstitutiveLawVector[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
                mConstitutiveLawVector[point_number]->InitializeMaterial( r_properties, r_geometry, row(N_values , point_number ));
            }
        } else {
            KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;
        }

    }
    KRATOS_CATCH( "" )
}


  void MembraneCuttingPatternElement::Relaxation(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) {

    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType number_dofs = dimension * number_of_nodes;

    const IntegrationMethod integration_method = r_geom.GetDefaultIntegrationMethod();

    Vector internal_forces = ZeroVector(number_dofs);
    this->InternalForces(internal_forces, integration_method, rCurrentProcessInfo);

    rRightHandSideVector.resize(number_dofs, false);
    noalias(rRightHandSideVector) = -internal_forces;

    this->TotalStiffnessMatrix(rLeftHandSideMatrix, integration_method, rCurrentProcessInfo);

    KRATOS_CATCH("");
  }


  void MembraneCuttingPatternElement::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)

{
    TotalStiffnessMatrix(rLeftHandSideMatrix,GetGeometry().GetDefaultIntegrationMethod(),rCurrentProcessInfo);
}


void MembraneCuttingPatternElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)

{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType system_size = number_of_nodes * dimension;

    Vector internal_forces = ZeroVector(system_size);
    InternalForces(internal_forces,GetGeometry().GetDefaultIntegrationMethod(),rCurrentProcessInfo);
    rRightHandSideVector.resize(system_size);
    noalias(rRightHandSideVector) = ZeroVector(system_size);
    noalias(rRightHandSideVector) -= internal_forces;
    //CalculateAndAddBodyForce(rRightHandSideVector,rCurrentProcessInfo);
}


void MembraneCuttingPatternElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)

{
    CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);
    CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
}


void MembraneCuttingPatternElement::GetValuesVector(
    Vector& rValues,
    int Step) const
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;

    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    for (SizeType i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3>& disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        const SizeType index = i * 3;
        rValues[index] = disp[0];
        rValues[index + 1] = disp[1];
        rValues[index + 2] = disp[2];
    }
}


void MembraneCuttingPatternElement::OptimizationLeastSquare(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, double& rResponse, const ProcessInfo& rCurrentProcessInfo) {


    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType number_dofs = dimension * number_of_nodes;

    const IntegrationMethod integration_method = r_geom.GetDefaultIntegrationMethod();

    rLeftHandSideMatrix = ZeroMatrix(number_dofs);

    rRightHandSideVector = ZeroVector(number_dofs);

    rResponse = 0.0;

    this->StiffnessMatrixLeastSquare(rLeftHandSideMatrix, integration_method, rCurrentProcessInfo);

    this->InternalForcesLeastSquare(rRightHandSideVector, integration_method, rCurrentProcessInfo);

    noalias(rRightHandSideVector) = -rRightHandSideVector;

    this->ResponseFunctionLeastSquare(rResponse, integration_method, rCurrentProcessInfo);

    KRATOS_CATCH("");
  }


  void MembraneCuttingPatternElement::DerivativeJacobiDeterminant(double& rDerivDetJ, const Matrix& rShapeFunctionGradientValues,
    const SizeType DofR, const array_1d<Vector, 2>& rBaseVectors, const double DetJ)
  {
    array_1d<Vector, 2> deriv_r_base_vectors;
    this->DeriveCurrentCovariantBaseVectors(deriv_r_base_vectors, rShapeFunctionGradientValues, DofR);

    Vector cross_g1g2 = ZeroVector(3);
    MathUtils<double>::CrossProduct(cross_g1g2, rBaseVectors[0], rBaseVectors[1]);

    Vector deriv_r_cross_1 = ZeroVector(3);
    MathUtils<double>::CrossProduct(deriv_r_cross_1, deriv_r_base_vectors[0], rBaseVectors[1]);
    Vector deriv_r_cross_2 = ZeroVector(3);
    MathUtils<double>::CrossProduct(deriv_r_cross_2, rBaseVectors[0], deriv_r_base_vectors[1]);

    const Vector deriv_r_cross = deriv_r_cross_1 + deriv_r_cross_2;

    rDerivDetJ = inner_prod(deriv_r_cross, cross_g1g2) / DetJ;
  }


  void MembraneCuttingPatternElement::Derivative2JacobiDeterminant(double& rDeriv2DetJ, const Matrix& rShapeFunctionGradientValues,
    const SizeType DofR, const SizeType DofS, const array_1d<Vector, 2>& rBaseVectors, const double DetJ,
    const double DerivRDetJ, const double DerivSDetJ)
  {
    array_1d<Vector, 2> deriv_r_base_vectors;
    this->DeriveCurrentCovariantBaseVectors(deriv_r_base_vectors, rShapeFunctionGradientValues, DofR);
    array_1d<Vector, 2> deriv_s_base_vectors;
    this->DeriveCurrentCovariantBaseVectors(deriv_s_base_vectors, rShapeFunctionGradientValues, DofS);

    Vector cross_g1g2 = ZeroVector(3);
    MathUtils<double>::CrossProduct(cross_g1g2, rBaseVectors[0], rBaseVectors[1]);

    Vector deriv_r_cross_1 = ZeroVector(3);
    MathUtils<double>::CrossProduct(deriv_r_cross_1, deriv_r_base_vectors[0], rBaseVectors[1]);
    Vector deriv_r_cross_2 = ZeroVector(3);
    MathUtils<double>::CrossProduct(deriv_r_cross_2, rBaseVectors[0], deriv_r_base_vectors[1]);
    const Vector deriv_r_cross = deriv_r_cross_1 + deriv_r_cross_2;

    Vector deriv_s_cross_1 = ZeroVector(3);
    MathUtils<double>::CrossProduct(deriv_s_cross_1, deriv_s_base_vectors[0], rBaseVectors[1]);
    Vector deriv_s_cross_2 = ZeroVector(3);
    MathUtils<double>::CrossProduct(deriv_s_cross_2, rBaseVectors[0], deriv_s_base_vectors[1]);
    const Vector deriv_s_cross = deriv_s_cross_1 + deriv_s_cross_2;

    Vector cross_rs_1 = ZeroVector(3);
    MathUtils<double>::CrossProduct(cross_rs_1, deriv_r_base_vectors[0], deriv_s_base_vectors[1]);
    Vector cross_rs_2 = ZeroVector(3);
    MathUtils<double>::CrossProduct(cross_rs_2, deriv_s_base_vectors[0], deriv_r_base_vectors[1]);
    const Vector deriv_rs_cross = cross_rs_1 + cross_rs_2;

    rDeriv2DetJ = (inner_prod(deriv_rs_cross, cross_g1g2) + inner_prod(deriv_r_cross, deriv_s_cross)) / DetJ
      - (DerivRDetJ * DerivSDetJ) / DetJ;
  }


  void MembraneCuttingPatternElement::ComputeLeastSquareBaseState(
    const Matrix& rShapeFunctionGradientValues,
    array_1d<Vector, 2>& rReferenceBaseVectors,
    array_1d<Vector, 2>& rCurrentBaseVectors,
    Matrix& rReferenceMetric,
    Matrix& rCurrentContravariantMetric,
    double& rDetJReference,
    double& rDetJCurrent,
    double rC_ref[2][2][2][2],
    double rC_act[2][2][2][2],
    Matrix& rEpsilon,
    Matrix& rSigma,
    Matrix& rSigmaP,
    const double YoungModulus,
    const double PoissonRatio)
  {
    this->CovariantBaseVectors(rReferenceBaseVectors, rShapeFunctionGradientValues, ConfigurationType::Reference);
    this->CovariantBaseVectors(rCurrentBaseVectors, rShapeFunctionGradientValues, ConfigurationType::Current);

    Matrix current_metric = ZeroMatrix(2);
    this->CovariantMetric(rReferenceMetric, rReferenceBaseVectors);
    this->CovariantMetric(current_metric, rCurrentBaseVectors);
    this->ContravariantMetric(rCurrentContravariantMetric, current_metric);

    array_1d<Vector, 2> current_contravariant_base_vectors;
    this->ContraVariantBaseVectors(current_contravariant_base_vectors, rCurrentContravariantMetric, rCurrentBaseVectors);

    this->JacobiDeterminante(rDetJReference, rReferenceBaseVectors);
    this->JacobiDeterminante(rDetJCurrent, rCurrentBaseVectors);

    const double inv_det_f = rDetJCurrent / rDetJReference;

    const double lambda = YoungModulus * PoissonRatio / (1.0 - PoissonRatio * PoissonRatio);
    const double mu = YoungModulus / (2.0 * (1.0 + PoissonRatio));

    for (SizeType alpha = 0; alpha < 2; ++alpha) {
      for (SizeType beta = 0; beta < 2; ++beta) {
        for (SizeType gamma = 0; gamma < 2; ++gamma) {
          for (SizeType delta = 0; delta < 2; ++delta) {
            rC_ref[alpha][beta][gamma][delta] = (lambda * (rCurrentContravariantMetric(alpha, beta) * rCurrentContravariantMetric(gamma, delta)))
              + (mu * ((rCurrentContravariantMetric(alpha, gamma) * rCurrentContravariantMetric(beta, delta))
                + (rCurrentContravariantMetric(alpha, delta) * rCurrentContravariantMetric(beta, gamma))));

            rC_act[alpha][beta][gamma][delta] = inv_det_f * rC_ref[alpha][beta][gamma][delta];
          }
        }
      }
    }

    rEpsilon = ZeroMatrix(2);
    noalias(rEpsilon) = 0.5 * (rReferenceMetric - current_metric);

    rSigma = ZeroMatrix(2);
    for (SizeType alpha = 0; alpha < 2; ++alpha) {
      for (SizeType beta = 0; beta < 2; ++beta) {
        for (SizeType gamma = 0; gamma < 2; ++gamma) {
          for (SizeType delta = 0; delta < 2; ++delta) {
            rSigma(alpha, beta) += rC_act[alpha][beta][gamma][delta] * rEpsilon(gamma, delta);
          }
        }
      }
    }

    this->PreStress(rSigmaP, current_contravariant_base_vectors);
  }


  void MembraneCuttingPatternElement::StiffnessMatrixLeastSquare(Matrix& rStiffnessMatrix, const IntegrationMethod& ThisMethod, const ProcessInfo& rCurrentProcessInfo) {

    const auto& r_geom = GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType number_dofs = dimension * number_of_nodes;

    const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = r_geom.ShapeFunctionsLocalGradients(ThisMethod);
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geom.IntegrationPoints(ThisMethod);

    const double thickness = GetProperties()[THICKNESS];

    const double young_modulus = GetProperties()[YOUNG_MODULUS];
    const double poisson_ratio = GetProperties()[POISSON_RATIO];
    const double lambda = young_modulus * poisson_ratio / (1.0 - poisson_ratio * poisson_ratio);
    const double mu = young_modulus / (2.0 * (1.0 + poisson_ratio));

    rStiffnessMatrix = ZeroMatrix(number_dofs);

    for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number) {

      const double integration_weight_i = r_integration_points[point_number].Weight();
      const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];

      array_1d<Vector, 2> reference_base_vectors, current_base_vectors;
      Matrix reference_metric = ZeroMatrix(2);
      Matrix current_contravariant_metric = ZeroMatrix(2);
      double det_j_reference = 0.0;
      double det_j_current = 0.0;
      double C_ref[2][2][2][2];
      double C_act[2][2][2][2];
      Matrix epsilon = ZeroMatrix(2);
      Matrix sigma = ZeroMatrix(2);
      Matrix sigma_p = ZeroMatrix(2);

      this->ComputeLeastSquareBaseState(shape_functions_gradients_i, reference_base_vectors, current_base_vectors,
        reference_metric, current_contravariant_metric, det_j_reference, det_j_current, C_ref, C_act,
        epsilon, sigma, sigma_p, young_modulus, poisson_ratio);

      const Matrix delta_sigma = sigma - sigma_p;
      const double inv_det_f = det_j_current / det_j_reference;

      for (SizeType dof_r = 0; dof_r < number_dofs; ++dof_r) {

        Matrix deriv_r_current_metric = ZeroMatrix(2);
        this->DerivativeCurrentCovariantMetric(deriv_r_current_metric, shape_functions_gradients_i, dof_r, current_base_vectors);
        Matrix deriv_r_current_contra_metric = ZeroMatrix(2);
        this->DerivativeContravariantMetric(deriv_r_current_contra_metric, shape_functions_gradients_i, dof_r, current_base_vectors);

        double deriv_r_det_j_current = 0.0;
        this->DerivativeJacobiDeterminant(deriv_r_det_j_current, shape_functions_gradients_i, dof_r, current_base_vectors, det_j_current);
        const double deriv_r_inv_det_f = deriv_r_det_j_current / det_j_reference;

        const Matrix deriv_r_epsilon = -0.5 * deriv_r_current_metric;

        double deriv_r_C_ref[2][2][2][2];
        double deriv_r_C_act[2][2][2][2];
        for (SizeType alpha = 0; alpha < 2; ++alpha) {
          for (SizeType beta = 0; beta < 2; ++beta) {
            for (SizeType gamma = 0; gamma < 2; ++gamma) {
              for (SizeType delta = 0; delta < 2; ++delta) {
                deriv_r_C_ref[alpha][beta][gamma][delta] =
                  lambda * (deriv_r_current_contra_metric(alpha, beta) * current_contravariant_metric(gamma, delta)
                    + current_contravariant_metric(alpha, beta) * deriv_r_current_contra_metric(gamma, delta))
                  + mu * (deriv_r_current_contra_metric(alpha, gamma) * current_contravariant_metric(beta, delta)
                    + deriv_r_current_contra_metric(alpha, delta) * current_contravariant_metric(beta, gamma)
                    + current_contravariant_metric(alpha, gamma) * deriv_r_current_contra_metric(beta, delta)
                    + current_contravariant_metric(alpha, delta) * deriv_r_current_contra_metric(beta, gamma));

                deriv_r_C_act[alpha][beta][gamma][delta] = deriv_r_inv_det_f * C_ref[alpha][beta][gamma][delta]
                  + inv_det_f * deriv_r_C_ref[alpha][beta][gamma][delta];
              }
            }
          }
        }

        Matrix deriv_r_sigma = ZeroMatrix(2);
        for (SizeType alpha = 0; alpha < 2; ++alpha) {
          for (SizeType beta = 0; beta < 2; ++beta) {
            for (SizeType gamma = 0; gamma < 2; ++gamma) {
              for (SizeType delta = 0; delta < 2; ++delta) {
                deriv_r_sigma(alpha, beta) += deriv_r_C_act[alpha][beta][gamma][delta] * epsilon(gamma, delta)
                  + C_act[alpha][beta][gamma][delta] * deriv_r_epsilon(gamma, delta);
              }
            }
          }
        }

        for (SizeType dof_s = 0; dof_s < number_dofs; ++dof_s) {

          Matrix deriv_s_current_metric = ZeroMatrix(2);
          this->DerivativeCurrentCovariantMetric(deriv_s_current_metric, shape_functions_gradients_i, dof_s, current_base_vectors);
          Matrix deriv_s_current_contra_metric = ZeroMatrix(2);
          this->DerivativeContravariantMetric(deriv_s_current_contra_metric, shape_functions_gradients_i, dof_s, current_base_vectors);

          double deriv_s_det_j_current = 0.0;
          this->DerivativeJacobiDeterminant(deriv_s_det_j_current, shape_functions_gradients_i, dof_s, current_base_vectors, det_j_current);
          const double deriv_s_inv_det_f = deriv_s_det_j_current / det_j_reference;

          const Matrix deriv_s_epsilon = -0.5 * deriv_s_current_metric;

          double deriv_s_C_ref[2][2][2][2];
          double deriv_s_C_act[2][2][2][2];
          for (SizeType alpha = 0; alpha < 2; ++alpha) {
            for (SizeType beta = 0; beta < 2; ++beta) {
              for (SizeType gamma = 0; gamma < 2; ++gamma) {
                for (SizeType delta = 0; delta < 2; ++delta) {
                  deriv_s_C_ref[alpha][beta][gamma][delta] =
                    lambda * (deriv_s_current_contra_metric(alpha, beta) * current_contravariant_metric(gamma, delta)
                      + current_contravariant_metric(alpha, beta) * deriv_s_current_contra_metric(gamma, delta))
                    + mu * (deriv_s_current_contra_metric(alpha, gamma) * current_contravariant_metric(beta, delta)
                      + deriv_s_current_contra_metric(alpha, delta) * current_contravariant_metric(beta, gamma)
                      + current_contravariant_metric(alpha, gamma) * deriv_s_current_contra_metric(beta, delta)
                      + current_contravariant_metric(alpha, delta) * deriv_s_current_contra_metric(beta, gamma));

                  deriv_s_C_act[alpha][beta][gamma][delta] = deriv_s_inv_det_f * C_ref[alpha][beta][gamma][delta]
                    + inv_det_f * deriv_s_C_ref[alpha][beta][gamma][delta];
                }
              }
            }
          }

          Matrix deriv_s_sigma = ZeroMatrix(2);
          for (SizeType alpha = 0; alpha < 2; ++alpha) {
            for (SizeType beta = 0; beta < 2; ++beta) {
              for (SizeType gamma = 0; gamma < 2; ++gamma) {
                for (SizeType delta = 0; delta < 2; ++delta) {
                  deriv_s_sigma(alpha, beta) += deriv_s_C_act[alpha][beta][gamma][delta] * epsilon(gamma, delta)
                    + C_act[alpha][beta][gamma][delta] * deriv_s_epsilon(gamma, delta);
                }
              }
            }
          }

          Matrix deriv2_current_metric = ZeroMatrix(2);
          this->Derivative2CurrentCovariantMetric(deriv2_current_metric, shape_functions_gradients_i, dof_r, dof_s);
          Matrix deriv2_current_contra_metric = ZeroMatrix(2);
          this->Derivative2ContravariantMetric(deriv2_current_contra_metric, shape_functions_gradients_i, dof_r, dof_s, current_base_vectors);

          double deriv2_det_j_current = 0.0;
          this->Derivative2JacobiDeterminant(deriv2_det_j_current, shape_functions_gradients_i, dof_r, dof_s, current_base_vectors,
            det_j_current, deriv_r_det_j_current, deriv_s_det_j_current);
          const double deriv2_inv_det_f = deriv2_det_j_current / det_j_reference;

          const Matrix deriv2_epsilon = -0.5 * deriv2_current_metric;

          double deriv2_C_ref[2][2][2][2];
          double deriv2_C_act[2][2][2][2];
          for (SizeType alpha = 0; alpha < 2; ++alpha) {
            for (SizeType beta = 0; beta < 2; ++beta) {
              for (SizeType gamma = 0; gamma < 2; ++gamma) {
                for (SizeType delta = 0; delta < 2; ++delta) {
                  deriv2_C_ref[alpha][beta][gamma][delta] =
                    lambda * (deriv2_current_contra_metric(alpha, beta) * current_contravariant_metric(gamma, delta)
                      + deriv_r_current_contra_metric(alpha, beta) * deriv_s_current_contra_metric(gamma, delta)
                      + deriv_s_current_contra_metric(alpha, beta) * deriv_r_current_contra_metric(gamma, delta)
                      + current_contravariant_metric(alpha, beta) * deriv2_current_contra_metric(gamma, delta))
                    + mu * (deriv2_current_contra_metric(alpha, gamma) * current_contravariant_metric(beta, delta)
                      + deriv_r_current_contra_metric(alpha, gamma) * deriv_s_current_contra_metric(beta, delta)
                      + deriv_s_current_contra_metric(alpha, gamma) * deriv_r_current_contra_metric(beta, delta)
                      + current_contravariant_metric(alpha, gamma) * deriv2_current_contra_metric(beta, delta)
                      + deriv2_current_contra_metric(alpha, delta) * current_contravariant_metric(beta, gamma)
                      + deriv_r_current_contra_metric(alpha, delta) * deriv_s_current_contra_metric(beta, gamma)
                      + deriv_s_current_contra_metric(alpha, delta) * deriv_r_current_contra_metric(beta, gamma)
                      + current_contravariant_metric(alpha, delta) * deriv2_current_contra_metric(beta, gamma));

                  deriv2_C_act[alpha][beta][gamma][delta] = deriv2_inv_det_f * C_ref[alpha][beta][gamma][delta]
                    + deriv_r_inv_det_f * deriv_s_C_ref[alpha][beta][gamma][delta]
                    + deriv_s_inv_det_f * deriv_r_C_ref[alpha][beta][gamma][delta]
                    + inv_det_f * deriv2_C_ref[alpha][beta][gamma][delta];
                }
              }
            }
          }

          Matrix deriv2_sigma = ZeroMatrix(2);
          for (SizeType alpha = 0; alpha < 2; ++alpha) {
            for (SizeType beta = 0; beta < 2; ++beta) {
              for (SizeType gamma = 0; gamma < 2; ++gamma) {
                for (SizeType delta = 0; delta < 2; ++delta) {
                  deriv2_sigma(alpha, beta) += deriv2_C_act[alpha][beta][gamma][delta] * epsilon(gamma, delta)
                    + deriv_r_C_act[alpha][beta][gamma][delta] * deriv_s_epsilon(gamma, delta)
                    + deriv_s_C_act[alpha][beta][gamma][delta] * deriv_r_epsilon(gamma, delta)
                    + C_act[alpha][beta][gamma][delta] * deriv2_epsilon(gamma, delta);
                }
              }
            }
          }

          double k_elm = 0.0;
          for (SizeType m = 0; m < 2; m++) {
            for (SizeType n = 0; n < 2; n++) {
              for (SizeType l = 0; l < 2; l++) {
                for (SizeType k = 0; k < 2; k++) {
                  k_elm += ((deriv_s_sigma(m, n) * deriv_r_sigma(l, k)) + (delta_sigma(m, n) * deriv2_sigma(l, k)))
                    * reference_metric(m, l) * reference_metric(n, k);
                }
              }
            }
          }

          rStiffnessMatrix(dof_r, dof_s) += k_elm * integration_weight_i * det_j_reference * thickness;
        }
      }
    }
  }


  void MembraneCuttingPatternElement::InternalForcesLeastSquare(Vector& rInternalForces, const IntegrationMethod& ThisMethod, const ProcessInfo& rCurrentProcessInfo) {

    const auto& r_geom = GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType number_dofs = dimension * number_of_nodes;
    rInternalForces = ZeroVector(number_dofs);

    const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = r_geom.ShapeFunctionsLocalGradients(ThisMethod);
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geom.IntegrationPoints(ThisMethod);

    const double thickness = GetProperties()[THICKNESS];

    const double young_modulus = GetProperties()[YOUNG_MODULUS];
    const double poisson_ratio = GetProperties()[POISSON_RATIO];
    const double lambda = young_modulus * poisson_ratio / (1.0 - poisson_ratio * poisson_ratio);
    const double mu = young_modulus / (2.0 * (1.0 + poisson_ratio));

    for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number) {

      const double integration_weight_i = r_integration_points[point_number].Weight();
      const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];

      array_1d<Vector, 2> reference_base_vectors, current_base_vectors;
      Matrix reference_metric = ZeroMatrix(2);
      Matrix current_contravariant_metric = ZeroMatrix(2);
      double det_j_reference = 0.0;
      double det_j_current = 0.0;
      double C_ref[2][2][2][2];
      double C_act[2][2][2][2];
      Matrix epsilon = ZeroMatrix(2);
      Matrix sigma = ZeroMatrix(2);
      Matrix sigma_p = ZeroMatrix(2);

      this->ComputeLeastSquareBaseState(shape_functions_gradients_i, reference_base_vectors, current_base_vectors,
        reference_metric, current_contravariant_metric, det_j_reference, det_j_current, C_ref, C_act,
        epsilon, sigma, sigma_p, young_modulus, poisson_ratio);

      const Matrix delta_sigma = sigma - sigma_p;
      const double inv_det_f = det_j_current / det_j_reference;

      for (SizeType dof_r = 0; dof_r < number_dofs; ++dof_r)
      {
        Matrix deriv_r_current_metric = ZeroMatrix(2);
        this->DerivativeCurrentCovariantMetric(deriv_r_current_metric, shape_functions_gradients_i, dof_r, current_base_vectors);
        Matrix deriv_r_current_contra_metric = ZeroMatrix(2);
        this->DerivativeContravariantMetric(deriv_r_current_contra_metric, shape_functions_gradients_i, dof_r, current_base_vectors);

        double deriv_r_det_j_current = 0.0;
        this->DerivativeJacobiDeterminant(deriv_r_det_j_current, shape_functions_gradients_i, dof_r, current_base_vectors, det_j_current);
        const double deriv_r_inv_det_f = deriv_r_det_j_current / det_j_reference;

        const Matrix deriv_r_epsilon = -0.5 * deriv_r_current_metric;

        double deriv_r_C_ref[2][2][2][2];
        double deriv_r_C_act[2][2][2][2];
        for (SizeType alpha = 0; alpha < 2; ++alpha) {
          for (SizeType beta = 0; beta < 2; ++beta) {
            for (SizeType gamma = 0; gamma < 2; ++gamma) {
              for (SizeType delta = 0; delta < 2; ++delta) {
                deriv_r_C_ref[alpha][beta][gamma][delta] =
                  lambda * (deriv_r_current_contra_metric(alpha, beta) * current_contravariant_metric(gamma, delta)
                    + current_contravariant_metric(alpha, beta) * deriv_r_current_contra_metric(gamma, delta))
                  + mu * (deriv_r_current_contra_metric(alpha, gamma) * current_contravariant_metric(beta, delta)
                    + deriv_r_current_contra_metric(alpha, delta) * current_contravariant_metric(beta, gamma)
                    + current_contravariant_metric(alpha, gamma) * deriv_r_current_contra_metric(beta, delta)
                    + current_contravariant_metric(alpha, delta) * deriv_r_current_contra_metric(beta, gamma));

                deriv_r_C_act[alpha][beta][gamma][delta] = deriv_r_inv_det_f * C_ref[alpha][beta][gamma][delta]
                  + inv_det_f * deriv_r_C_ref[alpha][beta][gamma][delta];
              }
            }
          }
        }

        Matrix deriv_r_sigma = ZeroMatrix(2);
        for (SizeType alpha = 0; alpha < 2; ++alpha) {
          for (SizeType beta = 0; beta < 2; ++beta) {
            for (SizeType gamma = 0; gamma < 2; ++gamma) {
              for (SizeType delta = 0; delta < 2; ++delta) {
                deriv_r_sigma(alpha, beta) += deriv_r_C_act[alpha][beta][gamma][delta] * epsilon(gamma, delta)
                  + C_act[alpha][beta][gamma][delta] * deriv_r_epsilon(gamma, delta);
              }
            }
          }
        }

        double f_int = 0.0;
        for (SizeType m = 0; m < 2; m++) {
          for (SizeType n = 0; n < 2; n++) {
            for (SizeType l = 0; l < 2; l++) {
              for (SizeType k = 0; k < 2; k++) {
                f_int += delta_sigma(m, n) * deriv_r_sigma(l, k) * reference_metric(m, l) * reference_metric(n, k);
              }
            }
          }
        }

        rInternalForces[dof_r] += f_int * integration_weight_i * det_j_reference * thickness;
      }
    }
  }


  void MembraneCuttingPatternElement::ResponseFunctionLeastSquare(double& rResponseLS, const IntegrationMethod& ThisMethod, const ProcessInfo& rCurrentProcessInfo) {

    const auto& r_geom = GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();

    const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = r_geom.ShapeFunctionsLocalGradients(ThisMethod);
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geom.IntegrationPoints(ThisMethod);

    const double thickness = GetProperties()[THICKNESS];

    const double young_modulus = GetProperties()[YOUNG_MODULUS];
    const double poisson_ratio = GetProperties()[POISSON_RATIO];

    rResponseLS = 0.0;

    for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number) {

      const double integration_weight_i = r_integration_points[point_number].Weight();
      const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];

      array_1d<Vector, 2> reference_base_vectors, current_base_vectors;
      Matrix reference_metric = ZeroMatrix(2);
      Matrix current_contravariant_metric = ZeroMatrix(2);
      double det_j_reference = 0.0;
      double det_j_current = 0.0;
      double C_ref[2][2][2][2];
      double C_act[2][2][2][2];
      Matrix epsilon = ZeroMatrix(2);
      Matrix sigma = ZeroMatrix(2);
      Matrix sigma_p = ZeroMatrix(2);

      this->ComputeLeastSquareBaseState(shape_functions_gradients_i, reference_base_vectors, current_base_vectors,
        reference_metric, current_contravariant_metric, det_j_reference, det_j_current, C_ref, C_act,
        epsilon, sigma, sigma_p, young_modulus, poisson_ratio);

      const Matrix delta_sigma = sigma - sigma_p;

      double product_delta_sigma = 0.0;

      for (SizeType m = 0; m < 2; m++) {
        for (SizeType n = 0; n < 2; n++) {
          for (SizeType l = 0; l < 2; l++) {
            for (SizeType k = 0; k < 2; k++) {
              product_delta_sigma += (delta_sigma(m, n) * delta_sigma(l, k)) * reference_metric(m, l) * reference_metric(n, k);
            }
          }
        }
      }

      product_delta_sigma *= det_j_reference * integration_weight_i * thickness;

      rResponseLS += 0.5 * product_delta_sigma;

    }
  }


  void MembraneCuttingPatternElement::InternalForces(Vector& rInternalForces, const IntegrationMethod& ThisMethod, const ProcessInfo& rCurrentProcessInfo)
  {
    const auto& r_geom = GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType number_dofs = dimension * number_of_nodes;
    rInternalForces = ZeroVector(number_dofs);

    const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = r_geom.ShapeFunctionsLocalGradients(ThisMethod);
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geom.IntegrationPoints(ThisMethod);

    const double thickness = GetProperties()[THICKNESS];

    array_1d<Vector, 2> current_covariant_base_vectors;
    array_1d<Vector, 2> reference_covariant_base_vectors;
    array_1d<Vector, 2> reference_contravariant_base_vectors;

    array_1d<Vector, 2> transformed_base_vectors;

    Matrix covariant_metric_current = ZeroMatrix(3);
    Matrix covariant_metric_reference = ZeroMatrix(3);
    Matrix contravariant_metric_reference = ZeroMatrix(3);
    Matrix inplane_transformation_matrix_material = ZeroMatrix(3);
    double detJ = 0.0;
    Vector stress = ZeroVector(3);
    Vector derivative_strain = ZeroVector(3);

    for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
      // getting information for integration
      const double integration_weight_i = r_integration_points[point_number].Weight();
      const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];

      this->CovariantBaseVectors(current_covariant_base_vectors, shape_functions_gradients_i, ConfigurationType::Current);
      this->CovariantBaseVectors(reference_covariant_base_vectors, shape_functions_gradients_i, ConfigurationType::Reference);

      this->CovariantMetric(covariant_metric_current, current_covariant_base_vectors);
      this->CovariantMetric(covariant_metric_reference, reference_covariant_base_vectors);
      this->ContravariantMetric(contravariant_metric_reference, covariant_metric_reference);

      this->ContraVariantBaseVectors(reference_contravariant_base_vectors, contravariant_metric_reference, reference_covariant_base_vectors);

      this->TransformBaseVectors(transformed_base_vectors, reference_contravariant_base_vectors);

      this->InPlaneTransformationMatrix(inplane_transformation_matrix_material, transformed_base_vectors, reference_contravariant_base_vectors);


      this->JacobiDeterminante(detJ, reference_covariant_base_vectors);
      Matrix material_tangent_modulus = ZeroMatrix(dimension);
      this->MaterialResponse(stress, contravariant_metric_reference, covariant_metric_reference, covariant_metric_current,
        transformed_base_vectors, inplane_transformation_matrix_material, point_number, material_tangent_modulus,
        rCurrentProcessInfo);

      for (SizeType dof_r = 0; dof_r < number_dofs; ++dof_r)
      {
        this->DerivativeStrainGreenLagrange(derivative_strain, shape_functions_gradients_i,
          dof_r, current_covariant_base_vectors, inplane_transformation_matrix_material);
        rInternalForces[dof_r] += inner_prod(stress, derivative_strain) * detJ * integration_weight_i * thickness;
      }
    }
  }


  void MembraneCuttingPatternElement::TotalStiffnessMatrix(Matrix& rStiffnessMatrix, const IntegrationMethod& ThisMethod,
    const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_WATCH("debug")
    const auto& r_geom = GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType number_dofs = dimension * number_of_nodes;
    rStiffnessMatrix = ZeroMatrix(number_dofs);

    const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = r_geom.ShapeFunctionsLocalGradients(ThisMethod);
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geom.IntegrationPoints(ThisMethod);

    const double thickness = GetProperties()[THICKNESS];

    array_1d<Vector, 2> current_covariant_base_vectors;
    array_1d<Vector, 2> reference_covariant_base_vectors;
    array_1d<Vector, 2> reference_contravariant_base_vectors;
    array_1d<Vector, 2> transformed_base_vectors;

    Matrix covariant_metric_current = ZeroMatrix(3);
    Matrix covariant_metric_reference = ZeroMatrix(3);
    Matrix contravariant_metric_reference = ZeroMatrix(3);
    Matrix inplane_transformation_matrix_material = ZeroMatrix(3);
    double detJ = 0.0;
    double temp_stiffness_entry;
    Vector stress = ZeroVector(3);

    for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
      // getting information for integration
      const double integration_weight_i = r_integration_points[point_number].Weight();
      const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];

      this->CovariantBaseVectors(current_covariant_base_vectors, shape_functions_gradients_i, ConfigurationType::Current);
      this->CovariantBaseVectors(reference_covariant_base_vectors, shape_functions_gradients_i, ConfigurationType::Reference);

      this->CovariantMetric(covariant_metric_current, current_covariant_base_vectors);
      this->CovariantMetric(covariant_metric_reference, reference_covariant_base_vectors);
      this->ContravariantMetric(contravariant_metric_reference, covariant_metric_reference);

      this->ContraVariantBaseVectors(reference_contravariant_base_vectors, contravariant_metric_reference, reference_covariant_base_vectors);

      this->TransformBaseVectors(transformed_base_vectors, reference_contravariant_base_vectors);

      this->InPlaneTransformationMatrix(inplane_transformation_matrix_material, transformed_base_vectors, reference_contravariant_base_vectors);

      this->JacobiDeterminante(detJ, reference_covariant_base_vectors);

      Matrix material_tangent_modulus = ZeroMatrix(dimension);
      this->MaterialResponse(stress, contravariant_metric_reference, covariant_metric_reference, covariant_metric_current,
        transformed_base_vectors, inplane_transformation_matrix_material, point_number, material_tangent_modulus,
        rCurrentProcessInfo);


      for (SizeType dof_s = 0; dof_s < number_dofs; ++dof_s) {
        for (SizeType dof_r = 0; dof_r < number_dofs; ++dof_r) {

          //do not calculate symmetric entries
          if (dof_s > dof_r) {
            if (point_number == (r_integration_points.size() - 1)) rStiffnessMatrix(dof_s, dof_r) = rStiffnessMatrix(dof_r, dof_s);
          }
          else {
            temp_stiffness_entry = 0.0;
            this->MaterialStiffnessMatrixEntryIJ(temp_stiffness_entry,
              material_tangent_modulus, dof_s, dof_r, shape_functions_gradients_i,
              current_covariant_base_vectors, inplane_transformation_matrix_material);
            this->InitialStressStiffnessMatrixEntryIJ(temp_stiffness_entry,
              stress, dof_s, dof_r, shape_functions_gradients_i,
              inplane_transformation_matrix_material);
            rStiffnessMatrix(dof_s, dof_r) += temp_stiffness_entry * detJ * integration_weight_i * thickness;
          }
        }
      }
    }
  }


  void MembraneCuttingPatternElement::MaterialStiffnessMatrixEntryIJ(double& rEntryIJ,
    const Matrix& rMaterialTangentModulus, const SizeType& rPositionI,
    const SizeType& rPositionJ, const Matrix& rShapeFunctionGradientValues,
    const array_1d<Vector, 2>& rCurrentCovariantBaseVectors, const Matrix& rTransformationMatrix)
  {
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    Vector strain_derivative = ZeroVector(dimension);
    this->DerivativeStrainGreenLagrange(strain_derivative, rShapeFunctionGradientValues, rPositionI,
      rCurrentCovariantBaseVectors, rTransformationMatrix);

    Vector stress_derivative = prod(rMaterialTangentModulus, strain_derivative);

    this->DerivativeStrainGreenLagrange(strain_derivative, rShapeFunctionGradientValues, rPositionJ,
      rCurrentCovariantBaseVectors, rTransformationMatrix);

    rEntryIJ += inner_prod(stress_derivative, strain_derivative);
  }


  void MembraneCuttingPatternElement::InitialStressStiffnessMatrixEntryIJ(double& rEntryIJ,
    const Vector& rStressVector, const SizeType& rPositionI,
    const SizeType& rPositionJ, const Matrix& rShapeFunctionGradientValues,
    const Matrix& rTransformationMatrix)
  {
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    Vector strain_derivative_2 = ZeroVector(dimension);
    this->Derivative2StrainGreenLagrange(strain_derivative_2, rShapeFunctionGradientValues, rPositionI, rPositionJ,
      rTransformationMatrix);
    rEntryIJ += inner_prod(rStressVector, strain_derivative_2);
  }



  void MembraneCuttingPatternElement::StrainEulerAlmansi(Matrix& rStrain, const Matrix& rReferenceCoVariantMetric, const Matrix& rCurrentCoVariantMetric)
  {
    rStrain = ZeroMatrix(2);
    noalias(rStrain) = 0.50 * (rCurrentCoVariantMetric - rReferenceCoVariantMetric);
  }


  void MembraneCuttingPatternElement::DerivativeStrainEulerAlmansi(Matrix& rStrain, const Matrix& rShapeFunctionGradientValues, const SizeType DofR,
    const array_1d<Vector, 2> rReferenceCovariantBaseVectors)
  {

    rStrain = ZeroMatrix(2);

    Matrix reference_covariant_metric_derivative = ZeroMatrix(2);
    this->DerivativeCurrentCovariantMetric(reference_covariant_metric_derivative, rShapeFunctionGradientValues, DofR, rReferenceCovariantBaseVectors);

    rStrain = -0.50 * reference_covariant_metric_derivative;
  
  }

  void MembraneCuttingPatternElement::Derivative2StrainEulerAlmansi(Matrix& rStrain, const Matrix& rShapeFunctionGradientValues, const SizeType DofR, const SizeType DofS)
  {
    
    rStrain = ZeroMatrix(2);
    
    Matrix reference_covariant_metric_derivative2 = ZeroMatrix(2);
    this->Derivative2CurrentCovariantMetric(reference_covariant_metric_derivative2, rShapeFunctionGradientValues, DofR, DofS);

    rStrain = -0.50 * reference_covariant_metric_derivative2;

  }


  void MembraneCuttingPatternElement::StrainGreenLagrange(Vector& rStrain, const Matrix& rReferenceCoVariantMetric, const Matrix& rCurrentCoVariantMetric,
    const Matrix& rTransformationMatrix)
  {
    Matrix strain_matrix = 0.50 * (rCurrentCoVariantMetric - rReferenceCoVariantMetric);
    Vector reference_strain = MathUtils<double>::StrainTensorToVector(strain_matrix, 3);
    TransformStrains(rStrain, reference_strain, rTransformationMatrix);
  }


  void MembraneCuttingPatternElement::DerivativeStrainGreenLagrange(Vector& rStrain, const Matrix& rShapeFunctionGradientValues, const SizeType DofR,
    const array_1d<Vector, 2> rCurrentCovariantBaseVectors, const Matrix& rTransformationMatrix)
  {
    Matrix current_covariant_metric_derivative = ZeroMatrix(2);
    DerivativeCurrentCovariantMetric(current_covariant_metric_derivative, rShapeFunctionGradientValues, DofR, rCurrentCovariantBaseVectors);
    Matrix strain_matrix_derivative = 0.50 * current_covariant_metric_derivative;
    Vector reference_strain = MathUtils<double>::StrainTensorToVector(strain_matrix_derivative, 3);
    TransformStrains(rStrain, reference_strain, rTransformationMatrix);
  }


  void MembraneCuttingPatternElement::Derivative2StrainGreenLagrange(Vector& rStrain,
    const Matrix& rShapeFunctionGradientValues, const SizeType DofR, const SizeType DofS,
    const Matrix& rTransformationMatrix)
  {
    Matrix current_covariant_metric_derivative = ZeroMatrix(2);
    Derivative2CurrentCovariantMetric(current_covariant_metric_derivative, rShapeFunctionGradientValues, DofR, DofS);

    Matrix strain_matrix_derivative = 0.50 * current_covariant_metric_derivative;

    Vector reference_strain = MathUtils<double>::StrainTensorToVector(strain_matrix_derivative, 3);
    TransformStrains(rStrain, reference_strain, rTransformationMatrix);
  }


  void MembraneCuttingPatternElement::TransformStrains(Vector& rStrains,
    Vector& rReferenceStrains, const Matrix& rTransformationMatrix)
  {
    // use contravariant basevectors here
    // transform base vecs needs only G3 which is equal for co and contra if it is orthogonal
    // tranform strains needs contra-variant !
    rStrains = ZeroVector(3);
    rReferenceStrains[2] /= 2.0; // extract E12 from voigt strain vector
    noalias(rStrains) = prod(rTransformationMatrix, rReferenceStrains);
    rStrains[2] *= 2.0; // include E12 and E21 for voigt strain vector
  }


  void MembraneCuttingPatternElement::ElasticityTensorKirchhoff(double rC_reference[2][2][2][2], double rc_current[2][2][2][2], const Matrix& rShapeFunctionGradientValues, const double YoungModulus, const double PoissonRatio) {

    
     array_1d<Vector, 2> reference_covariant_base_vectors;
     array_1d<Vector, 2> reference_contravariant_base_vectors;
     array_1d<Vector, 2> current_covariant_base_vectors;
    
     Matrix covariant_metric_reference = ZeroMatrix(2);
     Matrix contravariant_metric_reference = ZeroMatrix(2);

     Matrix deformation_gradient = ZeroMatrix(3);
     double det_deformation_gradient = 0.0;

     double E = YoungModulus;
     double NU = PoissonRatio;
     const double lambda = E * NU / (1.0 - NU * NU);
     const double MU = E / (2.0 * (1.0 + NU));


     this->CovariantBaseVectors(current_covariant_base_vectors, rShapeFunctionGradientValues, ConfigurationType::Current);
     this->CovariantBaseVectors(reference_covariant_base_vectors, rShapeFunctionGradientValues, ConfigurationType::Reference);
      
     this->CovariantMetric(covariant_metric_reference, reference_covariant_base_vectors);
     this->ContravariantMetric(contravariant_metric_reference, covariant_metric_reference);

     this->ContraVariantBaseVectors(reference_contravariant_base_vectors, contravariant_metric_reference, reference_covariant_base_vectors);
      
     this->DeformationGradient(deformation_gradient, det_deformation_gradient, current_covariant_base_vectors, reference_contravariant_base_vectors);

     for (SizeType alpha = 0; alpha < 2; ++alpha) {
       for (SizeType beta = 0; beta < 2; ++beta) {
         for (SizeType gamma = 0; gamma < 2; ++gamma) {
           for (SizeType delta = 0; delta < 2; ++delta) {
             rC_reference[alpha][beta][gamma][delta] = (lambda * (contravariant_metric_reference(alpha, beta) * contravariant_metric_reference(gamma, delta))) 
               + (MU * ((contravariant_metric_reference(alpha, gamma) * contravariant_metric_reference(beta, delta)) 
                 + (contravariant_metric_reference(alpha, delta) * contravariant_metric_reference(beta, gamma))));


             rc_current[alpha][beta][gamma][delta] = (1.0 / det_deformation_gradient) * rC_reference[alpha][beta][gamma][delta];

             
           }
         }
       }
     }
         
  }


  void MembraneCuttingPatternElement::DerivativeElasticityTensorKirchhoff(double rDerivC_reference[2][2][2][2], double rDerivc_current[2][2][2][2], const Matrix& rShapeFunctionGradientValues, const SizeType DofR, const double YoungModulus, const double PoissonRatio) {


    array_1d<Vector, 2> reference_covariant_base_vectors;
    array_1d<Vector, 2> current_covariant_base_vectors;
    array_1d<Vector, 2> reference_contravariant_base_vectors;
    
    Matrix covariant_metric_reference = ZeroMatrix(2);
    Matrix contravariant_metric_reference = ZeroMatrix(2); 
    Matrix derivative_contravariant_metric_reference = ZeroMatrix(2);

    Matrix deformation_gradient = ZeroMatrix(3);
    double det_deformation_gradient = 0.0;
    double derivative_inverse_det_deformation_gradient = 0.0;

    double rC_reference[2][2][2][2] = { {{{0.0}}} };
    double rc_current[2][2][2][2] = { {{{0.0}}} };
    
    
    double E = YoungModulus;
    double NU = PoissonRatio;
    const double lambda = E * NU / (1.0 - NU * NU);
    const double MU = E / (2.0 * (1.0 + NU));

      
    this->CovariantBaseVectors(reference_covariant_base_vectors, rShapeFunctionGradientValues, ConfigurationType::Reference);
    this->CovariantBaseVectors(current_covariant_base_vectors, rShapeFunctionGradientValues, ConfigurationType::Current);
      
    this->CovariantMetric(covariant_metric_reference, reference_covariant_base_vectors);
    this->ContravariantMetric(contravariant_metric_reference, covariant_metric_reference);

    this->ContraVariantBaseVectors(reference_contravariant_base_vectors, contravariant_metric_reference, reference_covariant_base_vectors);

    this->DerivativeContravariantMetric(derivative_contravariant_metric_reference, rShapeFunctionGradientValues, DofR, reference_covariant_base_vectors);

    this->DeformationGradient(deformation_gradient, det_deformation_gradient, current_covariant_base_vectors, reference_contravariant_base_vectors);
    this->DerivativeInvDetDeformationGradient(derivative_inverse_det_deformation_gradient, rShapeFunctionGradientValues, current_covariant_base_vectors, reference_covariant_base_vectors, DofR);

    this->ElasticityTensorKirchhoff(rC_reference, rc_current, rShapeFunctionGradientValues, YoungModulus, PoissonRatio);
      
      
            
    for (SizeType alpha = 0; alpha < 2; alpha++) {
      for (SizeType beta = 0; beta < 2; beta++) {
        for (SizeType gamma = 0; gamma < 2; gamma++) {
          for (SizeType delta = 0; delta < 2; delta++)
          {
            rDerivC_reference[alpha][beta][gamma][delta] = lambda * (derivative_contravariant_metric_reference(alpha, beta) * contravariant_metric_reference(gamma, delta) 
              + contravariant_metric_reference(alpha, beta) * derivative_contravariant_metric_reference(gamma, delta))
              + MU * (derivative_contravariant_metric_reference(alpha, gamma) * contravariant_metric_reference(beta, delta) 
                + derivative_contravariant_metric_reference(alpha, delta) * contravariant_metric_reference(beta, gamma) 
                + contravariant_metric_reference(alpha, gamma) * derivative_contravariant_metric_reference(beta, delta) 
                + contravariant_metric_reference(alpha, delta) * derivative_contravariant_metric_reference(beta, gamma));

              rDerivc_current[alpha][beta][gamma][delta] = derivative_inverse_det_deformation_gradient * rC_reference[alpha][beta][gamma][delta] +
                ((1.0 / det_deformation_gradient) * rDerivC_reference[alpha][beta][gamma][delta]);
            
          }
        }
      }
    }
              
  }


  void MembraneCuttingPatternElement::Derivative2ElasticityTensorKirchhoff(double rDeriv2C_reference[2][2][2][2], double rDeriv2c_current[2][2][2][2], const Matrix& rShapeFunctionGradientValues, const SizeType DofR, const SizeType DofS, const double YoungModulus, const double PoissonRatio) {

    double E = YoungModulus;
    double NU = PoissonRatio;
    const double lambda = E * NU / (1.0 - NU * NU);
    const double MU = E / (2.0 * (1.0 + NU));

    array_1d<Vector, 2> reference_covariant_base_vectors;
    array_1d<Vector, 2> current_covariant_base_vectors;
    array_1d<Vector, 2> reference_contravariant_base_vectors;
    
    Matrix covariant_metric_reference = ZeroMatrix(2);
    Matrix contravariant_metric_reference = ZeroMatrix(2);

    Matrix derivative_contravariant_metric_reference_r = ZeroMatrix(2);
    Matrix derivative_contravariant_metric_reference_s = ZeroMatrix(2);
    Matrix derivative2_contravariant_metric_reference = ZeroMatrix(2);

    Matrix deformation_gradient = ZeroMatrix(3);
    double det_deformation_gradient = 0.0;

    double derivative_r_inverse_det_deformation_gradient = 0.0;
    double derivative_s_inverse_det_deformation_gradient = 0.0;

    double derivative2_inverse_det_deformation_gradient = 0.0;

    double rC_reference[2][2][2][2] = { {{{0.0}}} };
    double rc_current[2][2][2][2] = { {{{0.0}}} };

    double deriv_r_C_reference[2][2][2][2] = { {{{0.0}}} };
    double deriv_r_c_current[2][2][2][2] = { {{{0.0}}} };

    double deriv_s_C_reference[2][2][2][2] = { {{{0.0}}} };
    double deriv_s_c_current[2][2][2][2] = { {{{0.0}}} };

    this->CovariantBaseVectors(reference_covariant_base_vectors, rShapeFunctionGradientValues, ConfigurationType::Reference);
    this->CovariantBaseVectors(current_covariant_base_vectors, rShapeFunctionGradientValues, ConfigurationType::Current);
    
    this->CovariantMetric(covariant_metric_reference, reference_covariant_base_vectors);
    this->ContravariantMetric(contravariant_metric_reference, covariant_metric_reference);

    this->ContraVariantBaseVectors(reference_contravariant_base_vectors, contravariant_metric_reference, reference_covariant_base_vectors);

    this->DerivativeContravariantMetric(derivative_contravariant_metric_reference_r, rShapeFunctionGradientValues, DofR, reference_covariant_base_vectors);
    this->DerivativeContravariantMetric(derivative_contravariant_metric_reference_s, rShapeFunctionGradientValues, DofS, reference_covariant_base_vectors);
    this->Derivative2ContravariantMetric(derivative2_contravariant_metric_reference, rShapeFunctionGradientValues, DofR, DofS, reference_covariant_base_vectors);

    this->DeformationGradient(deformation_gradient, det_deformation_gradient, current_covariant_base_vectors, reference_contravariant_base_vectors);

    this->DerivativeInvDetDeformationGradient(derivative_r_inverse_det_deformation_gradient, rShapeFunctionGradientValues, current_covariant_base_vectors, reference_covariant_base_vectors, DofR);
    this->DerivativeInvDetDeformationGradient(derivative_s_inverse_det_deformation_gradient, rShapeFunctionGradientValues, current_covariant_base_vectors, reference_covariant_base_vectors, DofS);

    this->Derivative2InvDetDeformationGradient(derivative2_inverse_det_deformation_gradient, rShapeFunctionGradientValues, current_covariant_base_vectors, reference_covariant_base_vectors, DofR, DofS);

    this->ElasticityTensorKirchhoff(rC_reference, rc_current, rShapeFunctionGradientValues, YoungModulus, PoissonRatio);

    this->DerivativeElasticityTensorKirchhoff(deriv_r_C_reference, deriv_r_c_current, rShapeFunctionGradientValues, DofR, YoungModulus, PoissonRatio);
    this->DerivativeElasticityTensorKirchhoff(deriv_s_C_reference, deriv_s_c_current, rShapeFunctionGradientValues, DofS, YoungModulus, PoissonRatio);

    
    for (SizeType alpha = 0; alpha < 2; alpha++) {
      for (SizeType beta = 0; beta < 2; beta++) {
        for (SizeType gamma = 0; gamma < 2; gamma++) {
          for (SizeType delta = 0; delta < 2; delta++)
          {

            rDeriv2C_reference[alpha][beta][gamma][delta] =
              lambda * (derivative2_contravariant_metric_reference(alpha, beta) * contravariant_metric_reference(gamma, delta) +
                derivative_contravariant_metric_reference_r(alpha, beta) * derivative_contravariant_metric_reference_s(gamma, delta) +
                derivative_contravariant_metric_reference_s(alpha, beta) * derivative_contravariant_metric_reference_r(gamma, delta) +
                contravariant_metric_reference(alpha, beta) * derivative2_contravariant_metric_reference(gamma, delta)) +
              MU * (derivative2_contravariant_metric_reference(alpha, gamma) * contravariant_metric_reference(beta, delta) +
                derivative_contravariant_metric_reference_r(alpha, gamma) * derivative_contravariant_metric_reference_s(beta, delta) +
                derivative_contravariant_metric_reference_s(alpha, gamma) * derivative_contravariant_metric_reference_r(beta, delta) +
                contravariant_metric_reference(alpha, gamma) * derivative2_contravariant_metric_reference(beta, delta) +
                derivative2_contravariant_metric_reference(alpha, delta) * contravariant_metric_reference(beta, gamma) +
                derivative_contravariant_metric_reference_r(alpha, delta) * derivative_contravariant_metric_reference_s(beta, gamma) +
                derivative_contravariant_metric_reference_s(alpha, delta) * derivative_contravariant_metric_reference_r(beta, gamma) +
                contravariant_metric_reference(alpha, delta) * derivative2_contravariant_metric_reference(beta, gamma));


            rDeriv2c_current[alpha][beta][gamma][delta] = derivative2_inverse_det_deformation_gradient * rC_reference[alpha][beta][gamma][delta] +
              derivative_r_inverse_det_deformation_gradient * deriv_s_C_reference[alpha][beta][gamma][delta] +
              derivative_s_inverse_det_deformation_gradient * deriv_r_C_reference[alpha][beta][gamma][delta] +
              ((1.0/ det_deformation_gradient)*rDeriv2C_reference[alpha][beta][gamma][delta]);


          }
        }
      }
    }


  }


  void MembraneCuttingPatternElement::CauchyStress(Matrix& rStress, const Matrix& rShapeFunctionGradientValues, const double YoungModulus, const double PoissonRatio) {

    double rC_reference[2][2][2][2] = { {{{0.0}}} };
    double rc_current[2][2][2][2] = { {{{0.0}}} };

    array_1d<Vector, 2> reference_covariant_base_vectors;
    array_1d<Vector, 2> current_covariant_base_vectors;
    
    Matrix covariant_metric_reference = ZeroMatrix(2);
    Matrix covariant_metric_current = ZeroMatrix(2);

    Matrix e = ZeroMatrix(2);

    rStress = ZeroMatrix(2);

    this->ElasticityTensorKirchhoff(rC_reference, rc_current, rShapeFunctionGradientValues, YoungModulus, PoissonRatio);

    this->CovariantBaseVectors(current_covariant_base_vectors, rShapeFunctionGradientValues, ConfigurationType::Current);
    this->CovariantBaseVectors(reference_covariant_base_vectors, rShapeFunctionGradientValues, ConfigurationType::Reference);

    this->CovariantMetric(covariant_metric_current, current_covariant_base_vectors);
    this->CovariantMetric(covariant_metric_reference, reference_covariant_base_vectors);

    this->StrainEulerAlmansi(e, covariant_metric_reference, covariant_metric_current);


    for (SizeType alpha = 0; alpha < 2; alpha++) {
      for (SizeType beta = 0; beta < 2; beta++) {
        for (SizeType gamma = 0; gamma < 2; gamma++) {
          for (SizeType delta = 0; delta < 2; delta++) {

            rStress(alpha, beta) += rc_current[alpha][beta][gamma][delta] * e(gamma, delta);

          }
        }
      }
    }

  }


  void MembraneCuttingPatternElement::DerivativeCauchyStress(Matrix& rStress, const Matrix& rShapeFunctionGradientValues, const SizeType DofR, const double YoungModulus, const double PoissonRatio) {


    array_1d<Vector, 2> reference_covariant_base_vectors;
    array_1d<Vector, 2> current_covariant_base_vectors;

    Matrix covariant_metric_reference = ZeroMatrix(2);
    Matrix covariant_metric_current = ZeroMatrix(2);
    
    double rC_reference[2][2][2][2] = { {{{0.0}}} };
    double rc_current[2][2][2][2] = { {{{0.0}}} };

    double rDerivC_reference[2][2][2][2] = { {{{0.0}}} };
    double rDerivc_current[2][2][2][2] = { {{{0.0}}} };

    Matrix e = ZeroMatrix(2);
    Matrix derivative_e = ZeroMatrix(2);

    rStress = ZeroMatrix(2);

    this->CovariantBaseVectors(current_covariant_base_vectors, rShapeFunctionGradientValues, ConfigurationType::Current);
    this->CovariantBaseVectors(reference_covariant_base_vectors, rShapeFunctionGradientValues, ConfigurationType::Reference);

    this->CovariantMetric(covariant_metric_current, current_covariant_base_vectors);
    this->CovariantMetric(covariant_metric_reference, reference_covariant_base_vectors);

    this->ElasticityTensorKirchhoff(rC_reference, rc_current, rShapeFunctionGradientValues, YoungModulus, PoissonRatio);

    this->DerivativeElasticityTensorKirchhoff(rDerivC_reference, rDerivc_current, rShapeFunctionGradientValues, DofR, YoungModulus, PoissonRatio);

    this->StrainEulerAlmansi(e, covariant_metric_reference, covariant_metric_current);
    this->DerivativeStrainEulerAlmansi(derivative_e, rShapeFunctionGradientValues, DofR, reference_covariant_base_vectors);

    for (SizeType alpha = 0; alpha < 2; alpha++) {
      for (SizeType beta = 0; beta < 2; beta++) {
        for (SizeType gamma = 0; gamma < 2; gamma++) {
          for (SizeType delta = 0; delta < 2; delta++) {

            rStress(alpha, beta) += rDerivc_current[alpha][beta][gamma][delta] * e(gamma, delta) + rc_current[alpha][beta][gamma][delta] * derivative_e(gamma, delta);

          }
        }
      }
    }

  }


  void MembraneCuttingPatternElement::Derivative2CauchyStress(Matrix& rStress, const Matrix& rShapeFunctionGradientValues, const SizeType DofR, const SizeType DofS, const double YoungModulus, const double PoissonRatio) {

    array_1d<Vector, 2> reference_covariant_base_vectors;
    array_1d<Vector, 2> current_covariant_base_vectors;

    Matrix covariant_metric_reference = ZeroMatrix(2);
    Matrix covariant_metric_current = ZeroMatrix(2);

    double C_reference[2][2][2][2] = { {{{0.0}}} };
    double c_current[2][2][2][2] = { {{{0.0}}} };

    double deriv_r_C_reference[2][2][2][2] = { {{{0.0}}} };
    double deriv_r_c_current[2][2][2][2] = { {{{0.0}}} };

    double deriv_s_C_reference[2][2][2][2] = { {{{0.0}}} };
    double deriv_s_c_current[2][2][2][2] = { {{{0.0}}} };

    double deriv2_C_reference[2][2][2][2] = { {{{0.0}}} };
    double deriv2_c_current[2][2][2][2] = { {{{0.0}}} };

    Matrix e = ZeroMatrix(2);
    Matrix deriv_r_e = ZeroMatrix(2);
    Matrix deriv_s_e = ZeroMatrix(2);
    Matrix deriv2_e = ZeroMatrix(2);

    rStress = ZeroMatrix(2);

    this->CovariantBaseVectors(current_covariant_base_vectors, rShapeFunctionGradientValues, ConfigurationType::Current);
    this->CovariantBaseVectors(reference_covariant_base_vectors, rShapeFunctionGradientValues, ConfigurationType::Reference);

    this->CovariantMetric(covariant_metric_current, current_covariant_base_vectors);
    this->CovariantMetric(covariant_metric_reference, reference_covariant_base_vectors);

    this->ElasticityTensorKirchhoff(C_reference, c_current, rShapeFunctionGradientValues, YoungModulus, PoissonRatio);

    this->DerivativeElasticityTensorKirchhoff(deriv_r_C_reference, deriv_r_c_current, rShapeFunctionGradientValues, DofR, YoungModulus, PoissonRatio);
    this->DerivativeElasticityTensorKirchhoff(deriv_s_C_reference, deriv_s_c_current, rShapeFunctionGradientValues, DofS, YoungModulus, PoissonRatio);

    this->Derivative2ElasticityTensorKirchhoff(deriv2_C_reference, deriv2_c_current, rShapeFunctionGradientValues, DofR, DofS, YoungModulus, PoissonRatio);

    this->StrainEulerAlmansi(e, covariant_metric_reference, covariant_metric_current);
    this->DerivativeStrainEulerAlmansi(deriv_r_e, rShapeFunctionGradientValues, DofR, reference_covariant_base_vectors);
    this->DerivativeStrainEulerAlmansi(deriv_s_e, rShapeFunctionGradientValues, DofS, reference_covariant_base_vectors);
    this->Derivative2StrainEulerAlmansi(deriv2_e, rShapeFunctionGradientValues, DofR, DofS);

    for (SizeType alpha = 0; alpha < 2; alpha++) {
      for (SizeType beta = 0; beta < 2; beta++) {
        for (SizeType gamma = 0; gamma < 2; gamma++) {
          for (SizeType delta = 0; delta < 2; delta++) {

            rStress(alpha, beta) += deriv2_c_current[alpha][beta][gamma][delta] * e(gamma, delta) +
              deriv_r_c_current[alpha][beta][gamma][delta] * deriv_s_e(gamma, delta) +
              deriv_s_c_current[alpha][beta][gamma][delta] * deriv_r_e(gamma, delta) +
              c_current[alpha][beta][gamma][delta] * deriv2_e(gamma, delta);

          }
        }
      }
    }

  }


  void MembraneCuttingPatternElement::AddPreStressPk2(Vector& rStress, const array_1d<Vector, 2>& rTransformedBaseVectors) {

    Vector pre_stress = ZeroVector(3);
    if (GetProperties().Has(PRESTRESS_VECTOR)) {
      pre_stress = GetProperties()(PRESTRESS_VECTOR);

      if (Has(LOCAL_PRESTRESS_AXIS_1) && Has(LOCAL_PRESTRESS_AXIS_2)) {

        array_1d<array_1d<double, 3>, 2> local_prestress_axis;
        local_prestress_axis[0] = GetValue(LOCAL_PRESTRESS_AXIS_1) / MathUtils<double>::Norm(GetValue(LOCAL_PRESTRESS_AXIS_1));
        local_prestress_axis[1] = GetValue(LOCAL_PRESTRESS_AXIS_2) / MathUtils<double>::Norm(GetValue(LOCAL_PRESTRESS_AXIS_2));

        Matrix transformation_matrix = ZeroMatrix(3);
        InPlaneTransformationMatrix(transformation_matrix, rTransformedBaseVectors, local_prestress_axis);
        pre_stress = prod(transformation_matrix, pre_stress);

      }
      else if (Has(LOCAL_PRESTRESS_AXIS_1)) {

        Vector base_3 = ZeroVector(3);
        MathUtils<double>::UnitCrossProduct(base_3, rTransformedBaseVectors[0], rTransformedBaseVectors[1]);

        array_1d<array_1d<double, 3>, 2> local_prestress_axis;
        local_prestress_axis[0] = GetValue(LOCAL_PRESTRESS_AXIS_1) / MathUtils<double>::Norm(GetValue(LOCAL_PRESTRESS_AXIS_1));

        MathUtils<double>::UnitCrossProduct(local_prestress_axis[1], base_3, local_prestress_axis[0]);

        Matrix transformation_matrix = ZeroMatrix(3);
        InPlaneTransformationMatrix(transformation_matrix, rTransformedBaseVectors, local_prestress_axis);
        pre_stress = prod(transformation_matrix, pre_stress);
      }
    }
    noalias(rStress) += pre_stress;
  }


  void MembraneCuttingPatternElement::MaterialResponse(Vector& rStress,
    const Matrix& rReferenceContraVariantMetric, const Matrix& rReferenceCoVariantMetric, const Matrix& rCurrentCoVariantMetric,
    const array_1d<Vector, 2>& rTransformedBaseVectors, const Matrix& rTransformationMatrix, const SizeType& rIntegrationPointNumber,
    Matrix& rTangentModulus, const ProcessInfo& rCurrentProcessInfo)
  {
    Vector strain_vector = ZeroVector(3);
    noalias(rStress) = ZeroVector(3);
    StrainGreenLagrange(strain_vector, rReferenceCoVariantMetric,
      rCurrentCoVariantMetric, rTransformationMatrix);

    // do this to consider the pre-stress influence in the check of the membrane state in the claw
    Vector initial_stress = ZeroVector(3);
    if (Has(MEMBRANE_PRESTRESS)) {
      const Matrix& r_stress_input = GetValue(MEMBRANE_PRESTRESS);
      initial_stress += column(r_stress_input, rIntegrationPointNumber);
    }
    else {
      AddPreStressPk2(initial_stress, rTransformedBaseVectors);
    }
    rStress += initial_stress;

    ConstitutiveLaw::Parameters element_parameters(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    element_parameters.SetStrainVector(strain_vector);
    element_parameters.SetStressVector(rStress);
    element_parameters.SetConstitutiveMatrix(rTangentModulus);
    element_parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    element_parameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    element_parameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    mConstitutiveLawVector[rIntegrationPointNumber]->CalculateMaterialResponse(element_parameters, ConstitutiveLaw::StressMeasure_PK2);

    // do this to include the pre-stress in the actual stress state
    // rStress is reset in the claw and thus does not consider the initial stress anymore
    rStress += initial_stress;
  }


  void MembraneCuttingPatternElement::Derivative2ContravariantMetric(Matrix& rMetric, const Matrix& rShapeFunctionGradientValues, const SizeType DofR, const SizeType DofS, const array_1d<Vector, 2>& rCovariantBaseVectors) {

    Matrix covariant_metric = ZeroMatrix(2);

    Matrix derivative_covariant_metric_r = ZeroMatrix(2);
    Matrix derivative_covariant_metric_s = ZeroMatrix(2);

    Matrix derivative2_covariant_metric = ZeroMatrix(2);
    
    rMetric = ZeroMatrix(2);

    this->CovariantMetric(covariant_metric, rCovariantBaseVectors);

    double determinant_covariant_metric = (covariant_metric(1, 1) * covariant_metric(0, 0)) - (covariant_metric(1, 0) * covariant_metric(0, 1));
    double inverse_determinant_covariant_metric = 1.0 / determinant_covariant_metric;

    this->DerivativeCurrentCovariantMetric(derivative_covariant_metric_r, rShapeFunctionGradientValues, DofR, rCovariantBaseVectors);
    this->DerivativeCurrentCovariantMetric(derivative_covariant_metric_s, rShapeFunctionGradientValues, DofS, rCovariantBaseVectors);

    this->Derivative2CurrentCovariantMetric(derivative2_covariant_metric, rShapeFunctionGradientValues, DofR, DofS);

    double derivative_r_determinant_covariant_metric = derivative_covariant_metric_r(0, 0) * covariant_metric(1, 1) + covariant_metric(0, 0) * derivative_covariant_metric_r(1, 1)
      - derivative_covariant_metric_r(1, 0) * covariant_metric(0, 1) - covariant_metric(1, 0) * derivative_covariant_metric_r(0, 1);
    double derivative_s_determinant_covariant_metric = derivative_covariant_metric_s(0, 0) * covariant_metric(1, 1) + covariant_metric(0, 0) * derivative_covariant_metric_s(1, 1)
      - derivative_covariant_metric_s(1, 0) * covariant_metric(0, 1) - covariant_metric(1, 0) * derivative_covariant_metric_s(0, 1);

    double derivative2_determinant_covariant_metric = derivative2_covariant_metric(0, 0) * covariant_metric(1, 1) + derivative_covariant_metric_r(0, 0) * derivative_covariant_metric_s(1, 1)
      + derivative_covariant_metric_s(0, 0) * derivative_covariant_metric_r(1, 1) + covariant_metric(0, 0) * derivative2_covariant_metric(1, 1)
      - derivative2_covariant_metric(0, 1) * covariant_metric(1, 0) - derivative_covariant_metric_r(0, 1) * derivative_covariant_metric_s(1, 0)
      - derivative_covariant_metric_s(0, 1) * derivative_covariant_metric_r(1, 0) - covariant_metric(0, 1) * derivative2_covariant_metric(1, 0);

    double derivative_r_inverse_determinant_covariant_metric = -derivative_r_determinant_covariant_metric / (determinant_covariant_metric * determinant_covariant_metric);
    double derivative_s_inverse_determinant_covariant_metric = -derivative_s_determinant_covariant_metric / (determinant_covariant_metric * determinant_covariant_metric);

    double derivative2_inverse_determinant_covariant_metric = (-derivative2_determinant_covariant_metric * determinant_covariant_metric
      + 2 * derivative_r_determinant_covariant_metric * derivative_s_determinant_covariant_metric) / (determinant_covariant_metric * determinant_covariant_metric * determinant_covariant_metric);


    rMetric(0, 0) = derivative2_inverse_determinant_covariant_metric * covariant_metric(1, 1)
      + derivative_r_inverse_determinant_covariant_metric * derivative_covariant_metric_s(1, 1)
      + derivative_s_inverse_determinant_covariant_metric * derivative_covariant_metric_r(1, 1)
      + inverse_determinant_covariant_metric * derivative2_covariant_metric(1, 1);

    rMetric(0, 1) = -derivative2_inverse_determinant_covariant_metric * covariant_metric(0, 1)
      - derivative_r_inverse_determinant_covariant_metric * derivative_covariant_metric_s(0, 1)
      - derivative_s_inverse_determinant_covariant_metric * derivative_covariant_metric_r(0, 1)
      - inverse_determinant_covariant_metric * derivative2_covariant_metric(0, 1);

    rMetric(1, 0) = -derivative2_inverse_determinant_covariant_metric * covariant_metric(1, 0)
      - derivative_r_inverse_determinant_covariant_metric * derivative_covariant_metric_s(1, 0)
      - derivative_s_inverse_determinant_covariant_metric * derivative_covariant_metric_r(1, 0)
      - inverse_determinant_covariant_metric * derivative2_covariant_metric(1, 0);

    rMetric(1, 1) = derivative2_inverse_determinant_covariant_metric * covariant_metric(0, 0)
      + derivative_r_inverse_determinant_covariant_metric * derivative_covariant_metric_s(0, 0)
      + derivative_s_inverse_determinant_covariant_metric * derivative_covariant_metric_r(0, 0)
      + inverse_determinant_covariant_metric * derivative2_covariant_metric(0, 0);

  }


  void MembraneCuttingPatternElement::DerivativeContravariantMetric(Matrix& rMetric, const Matrix& rShapeFunctionGradientValues, const SizeType DofR, const array_1d<Vector, 2>& rCovariantBaseVectors) {

    Matrix covariant_metric = ZeroMatrix(2);
    Matrix derivative_covariant_metric = ZeroMatrix(2);

    rMetric = ZeroMatrix(2);

    this->CovariantMetric(covariant_metric, rCovariantBaseVectors);

    double determinant_covariant_metric = (covariant_metric(1, 1) * covariant_metric(0, 0)) - (covariant_metric(1, 0) * covariant_metric(0, 1));
    double inverse_determinant_covariant_metric = 1.0 / determinant_covariant_metric;

    this->DerivativeCurrentCovariantMetric(derivative_covariant_metric, rShapeFunctionGradientValues, DofR, rCovariantBaseVectors);

    double derivative_determinant_covariant_metric = derivative_covariant_metric(0, 0) * covariant_metric(1, 1) + covariant_metric(0, 0) * derivative_covariant_metric(1, 1)
      - derivative_covariant_metric(1, 0) * covariant_metric(0, 1) - covariant_metric(1, 0) * derivative_covariant_metric(0, 1);
    double derivative_inverse_determinant_covariant_metric = -derivative_determinant_covariant_metric / (determinant_covariant_metric * determinant_covariant_metric);

    rMetric(0, 0) = derivative_inverse_determinant_covariant_metric * covariant_metric(1, 1) + inverse_determinant_covariant_metric * derivative_covariant_metric(1, 1);
    rMetric(0, 1) = -derivative_inverse_determinant_covariant_metric * covariant_metric(0, 1) - inverse_determinant_covariant_metric * derivative_covariant_metric(0, 1);
    rMetric(1, 0) = -derivative_inverse_determinant_covariant_metric * covariant_metric(1, 0) - inverse_determinant_covariant_metric * derivative_covariant_metric(1, 0);
    rMetric(1, 1) = derivative_inverse_determinant_covariant_metric * covariant_metric(0, 0) + inverse_determinant_covariant_metric * derivative_covariant_metric(0, 0);

  }


  void MembraneCuttingPatternElement::ContravariantMetric(Matrix& rMetric, const Matrix& rCovariantMetric)
  {
    rMetric = ZeroMatrix(2);
    rMetric(0, 0) = rCovariantMetric(1, 1);
    rMetric(1, 1) = rCovariantMetric(0, 0);
    rMetric(0, 1) = -1.0 * rCovariantMetric(0, 1);
    rMetric(1, 0) = -1.0 * rCovariantMetric(1, 0);
    rMetric /= (rCovariantMetric(1, 1) * rCovariantMetric(0, 0)) - (rCovariantMetric(1, 0) * rCovariantMetric(0, 1));

  }


  void MembraneCuttingPatternElement::ContraVariantBaseVectors(array_1d<Vector, 2>& rBaseVectors, const Matrix& rContraVariantMetric,
    const array_1d<Vector, 2> rCovariantBaseVectors)
  {
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    rBaseVectors[0] = ZeroVector(dimension);
    rBaseVectors[1] = ZeroVector(dimension);

    rBaseVectors[0] = rContraVariantMetric(0, 0) * rCovariantBaseVectors[0] + rContraVariantMetric(0, 1) * rCovariantBaseVectors[1];
    rBaseVectors[1] = rContraVariantMetric(1, 0) * rCovariantBaseVectors[0] + rContraVariantMetric(1, 1) * rCovariantBaseVectors[1];
  }


  void MembraneCuttingPatternElement::CovariantMetric(Matrix& rMetric, const array_1d<Vector, 2>& rBaseVectorCovariant)
  {
    rMetric = ZeroMatrix(2);
    for (SizeType i = 0; i < 2; ++i) {
      for (SizeType j = 0; j < 2; ++j) {
        rMetric(i, j) = inner_prod(rBaseVectorCovariant[i], rBaseVectorCovariant[j]);
      }
    }
  }


  void MembraneCuttingPatternElement::Derivative2CurrentCovariantMetric(Matrix& rMetric,
    const Matrix& rShapeFunctionGradientValues, const SizeType DofR, const SizeType DofS)
  {
    rMetric = ZeroMatrix(2);
    array_1d<Vector, 2> derivative_covariant_base_vectors_dur;
    DeriveCurrentCovariantBaseVectors(derivative_covariant_base_vectors_dur, rShapeFunctionGradientValues, DofR);
    array_1d<Vector, 2> derivative_covariant_base_vectors_dus;
    DeriveCurrentCovariantBaseVectors(derivative_covariant_base_vectors_dus, rShapeFunctionGradientValues, DofS);

    for (SizeType i = 0; i < 2; ++i) {
      for (SizeType j = 0; j < 2; ++j) {
        rMetric(i, j) = inner_prod(derivative_covariant_base_vectors_dur[i], derivative_covariant_base_vectors_dus[j]);
        rMetric(i, j) += inner_prod(derivative_covariant_base_vectors_dus[i], derivative_covariant_base_vectors_dur[j]);
      }
    }
  }
  
  
  void MembraneCuttingPatternElement::DerivativeCurrentCovariantMetric(Matrix& rMetric,
    const Matrix& rShapeFunctionGradientValues, const SizeType DofR, const array_1d<Vector, 2> rCurrentCovariantBaseVectors)
  {
    rMetric = ZeroMatrix(2);
    array_1d<Vector, 2> derivative_covariant_base_vectors;
    DeriveCurrentCovariantBaseVectors(derivative_covariant_base_vectors, rShapeFunctionGradientValues, DofR);


    for (SizeType i = 0; i < 2; ++i) {
      for (SizeType j = 0; j < 2; ++j) {
        rMetric(i, j) = inner_prod(derivative_covariant_base_vectors[i], rCurrentCovariantBaseVectors[j]);
        rMetric(i, j) += inner_prod(derivative_covariant_base_vectors[j], rCurrentCovariantBaseVectors[i]);
      }
    }
  }


  void MembraneCuttingPatternElement::DeriveCurrentCovariantBaseVectors(array_1d<Vector, 2>& rBaseVectors,
    const Matrix& rShapeFunctionGradientValues, const SizeType DofR)
  {
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType dof_nr = DofR % dimension;
    const SizeType node_nr = (DofR - dof_nr) / dimension;
    for (SizeType i = 0; i < 2; ++i) {
      rBaseVectors[i] = ZeroVector(dimension);
      rBaseVectors[i][dof_nr] = rShapeFunctionGradientValues(node_nr, i);
    }
  }


  void MembraneCuttingPatternElement::CovariantBaseVectors(array_1d<Vector, 2>& rBaseVectors,
    const Matrix& rShapeFunctionGradientValues, const ConfigurationType& rConfiguration) const
  {
    // pass/call this ShapeFunctionsLocalGradients[pnt]
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType number_of_nodes = GetGeometry().size();
    Vector g1 = ZeroVector(dimension);
    Vector g2 = ZeroVector(dimension);

    Vector current_displacement = ZeroVector(dimension * number_of_nodes);
    if (rConfiguration == ConfigurationType::Current) GetValuesVector(current_displacement);


    for (SizeType i = 0; i < number_of_nodes; ++i) {
      g1[0] += (GetGeometry().GetPoint(i).X0() + current_displacement[i * dimension]) * rShapeFunctionGradientValues(i, 0);
      g1[1] += (GetGeometry().GetPoint(i).Y0() + current_displacement[(i * dimension) + 1]) * rShapeFunctionGradientValues(i, 0);
      g1[2] += (GetGeometry().GetPoint(i).Z0() + current_displacement[(i * dimension) + 2]) * rShapeFunctionGradientValues(i, 0);

      g2[0] += (GetGeometry().GetPoint(i).X0() + current_displacement[i * dimension]) * rShapeFunctionGradientValues(i, 1);
      g2[1] += (GetGeometry().GetPoint(i).Y0() + current_displacement[(i * dimension) + 1]) * rShapeFunctionGradientValues(i, 1);
      g2[2] += (GetGeometry().GetPoint(i).Z0() + current_displacement[(i * dimension) + 2]) * rShapeFunctionGradientValues(i, 1);
    }
    rBaseVectors[0] = g1;
    rBaseVectors[1] = g2;
  }


  void MembraneCuttingPatternElement::TransformBaseVectors(array_1d<Vector, 2>& rBaseVectors,
    const array_1d<Vector, 2>& rLocalBaseVectors) {

    // create local cartesian coordinate system aligned to global material vectors (orthotropic)
    if (Has(LOCAL_MATERIAL_AXIS_1) && Has(LOCAL_MATERIAL_AXIS_2)) {
      rBaseVectors[0] = GetValue(LOCAL_MATERIAL_AXIS_1) / MathUtils<double>::Norm(GetValue(LOCAL_MATERIAL_AXIS_1));
      rBaseVectors[1] = GetValue(LOCAL_MATERIAL_AXIS_2) / MathUtils<double>::Norm(GetValue(LOCAL_MATERIAL_AXIS_2));
    }
    else if (Has(LOCAL_MATERIAL_AXIS_1)) {
      Vector base_3 = ZeroVector(3);
      MathUtils<double>::UnitCrossProduct(base_3, rLocalBaseVectors[0], rLocalBaseVectors[1]);

      rBaseVectors[0] = GetValue(LOCAL_MATERIAL_AXIS_1) / MathUtils<double>::Norm(GetValue(LOCAL_MATERIAL_AXIS_1));
      MathUtils<double>::UnitCrossProduct(rBaseVectors[1], base_3, rBaseVectors[0]);

    }
    else {
      // create local cartesian coordinate system
      rBaseVectors[0] = ZeroVector(3);
      rBaseVectors[1] = ZeroVector(3);
      rBaseVectors[0] = rLocalBaseVectors[0] / MathUtils<double>::Norm(rLocalBaseVectors[0]);
      rBaseVectors[1] = rLocalBaseVectors[1] - (inner_prod(rLocalBaseVectors[1], rBaseVectors[0]) * rBaseVectors[0]);
      rBaseVectors[1] /= MathUtils<double>::Norm(rBaseVectors[1]);
    }
  }


  template <class T>
  void MembraneCuttingPatternElement::InPlaneTransformationMatrix(Matrix& rTransformationMatrix, const array_1d<Vector, 2>& rTransformedBaseVectors,
    const T& rLocalReferenceBaseVectors)
  {
    const double e_g_11 = inner_prod(rTransformedBaseVectors[0], rLocalReferenceBaseVectors[0]);
    const double e_g_12 = inner_prod(rTransformedBaseVectors[0], rLocalReferenceBaseVectors[1]);
    const double e_g_21 = inner_prod(rTransformedBaseVectors[1], rLocalReferenceBaseVectors[0]);
    const double e_g_22 = inner_prod(rTransformedBaseVectors[1], rLocalReferenceBaseVectors[1]);
    rTransformationMatrix = ZeroMatrix(3);
    rTransformationMatrix(0, 0) = e_g_11 * e_g_11;
    rTransformationMatrix(0, 1) = e_g_12 * e_g_12;
    rTransformationMatrix(0, 2) = 2.0 * e_g_11 * e_g_12;
    rTransformationMatrix(1, 0) = e_g_21 * e_g_21;
    rTransformationMatrix(1, 1) = e_g_22 * e_g_22;
    rTransformationMatrix(1, 2) = 2.0 * e_g_21 * e_g_22;
    rTransformationMatrix(2, 0) = e_g_11 * e_g_21;
    rTransformationMatrix(2, 1) = e_g_12 * e_g_22;
    rTransformationMatrix(2, 2) = (e_g_11 * e_g_22) + (e_g_12 * e_g_21);
  }


  void MembraneCuttingPatternElement::JacobiDeterminante(double& rDetJacobi, const array_1d<Vector, 2>& rReferenceBaseVectors) const
  {
    Vector3 g3 = ZeroVector(3);
    MathUtils<double>::CrossProduct(g3, rReferenceBaseVectors[0], rReferenceBaseVectors[1]);
    rDetJacobi = MathUtils<double>::Norm(g3);
    KRATOS_ERROR_IF(rDetJacobi < std::numeric_limits<double>::epsilon()) << "det of Jacobi smaller 0 for element with id" << Id() << std::endl;
  }


  void MembraneCuttingPatternElement::DeformationGradient(Matrix& rDeformationGradient, double& rDetDeformationGradient,
    const array_1d<Vector, 2>& rCurrentCovariantBase, const array_1d<Vector, 2>& rReferenceContraVariantBase)
  {
     //attention: this is not in the local orthonogal coordinate system

     //calculate out of plane local vectors (membrane has no thickness change so g3 and G3 are normalized)
    Vector current_cov_3 = ZeroVector(3);
    MathUtils<double>::UnitCrossProduct(current_cov_3, rCurrentCovariantBase[0], rCurrentCovariantBase[1]);

    Vector reference_contra_3 = ZeroVector(3);
    MathUtils<double>::UnitCrossProduct(reference_contra_3, rReferenceContraVariantBase[0], rReferenceContraVariantBase[1]);

    // calculate deformation gradient
    rDeformationGradient = ZeroMatrix(3);
    for (SizeType i = 0; i < 2; ++i) {
      rDeformationGradient += outer_prod(rCurrentCovariantBase[i], rReferenceContraVariantBase[i]);
    }

    // add contribution of out of plane base vectors
    rDeformationGradient += outer_prod(current_cov_3, reference_contra_3);


    // calculate det(F)
    rDetDeformationGradient = MathUtils<double>::Det(rDeformationGradient);


  }


  void MembraneCuttingPatternElement::DerivativeInvDetDeformationGradient(double& rDerivInvDetDeformationGradient, const Matrix& rShapeFunctionGradientValues, const array_1d<Vector, 2>& rCurrentCovariantBaseVectors, const array_1d<Vector, 2>& rReferenceCovariantBaseVectors, const SizeType DofR) {


    array_1d<Vector, 2> derivative_reference_covariant_base_vectors;

    Vector cross_G1G2 = ZeroVector(3);
    Vector cross_g1g2 = ZeroVector(3);

    Vector total_deriv_cross_G1G2 = ZeroVector(3);
    Vector deriv_cross_G1G2_1 = ZeroVector(3);
    Vector deriv_cross_G1G2_2 = ZeroVector(3);

    double norm_cross_g1g2 = 0.0;
    double norm_cross_G1G2 = 0.0;
    double deriv_norm_cross_G1G2 = 0.0;


    MathUtils<double>::CrossProduct(cross_G1G2, rReferenceCovariantBaseVectors[0], rReferenceCovariantBaseVectors[1]);
    norm_cross_G1G2 = MathUtils<double>::Norm(cross_G1G2);

    MathUtils<double>::CrossProduct(cross_g1g2, rCurrentCovariantBaseVectors[0], rCurrentCovariantBaseVectors[1]);
    norm_cross_g1g2 = MathUtils<double>::Norm(cross_g1g2);

    this->DeriveCurrentCovariantBaseVectors(derivative_reference_covariant_base_vectors, rShapeFunctionGradientValues, DofR);

    MathUtils<double>::CrossProduct(deriv_cross_G1G2_1, derivative_reference_covariant_base_vectors[0], rReferenceCovariantBaseVectors[1]);
    MathUtils<double>::CrossProduct(deriv_cross_G1G2_2, rReferenceCovariantBaseVectors[0], derivative_reference_covariant_base_vectors[1]);
    total_deriv_cross_G1G2 = deriv_cross_G1G2_1 + deriv_cross_G1G2_2;

    deriv_norm_cross_G1G2 = inner_prod(cross_G1G2, total_deriv_cross_G1G2) / norm_cross_G1G2;

    
    rDerivInvDetDeformationGradient = deriv_norm_cross_G1G2 / norm_cross_g1g2;


  }


  void MembraneCuttingPatternElement::Derivative2InvDetDeformationGradient(double& rDeriv2InvDetDeformationGradient, const Matrix& rShapeFunctionGradientValues, const array_1d<Vector, 2>& rCurrentCovariantBaseVectors, const array_1d<Vector, 2>& rReferenceCovariantBaseVectors, const SizeType DofR, const SizeType DofS) {

    array_1d<Vector, 2> derivative_r_reference_covariant_base_vectors;
    array_1d<Vector, 2> derivative_s_reference_covariant_base_vectors;

    Vector cross_G1G2 = ZeroVector(3);
    Vector cross_g1g2 = ZeroVector(3);

    Vector total_deriv_r_cross_G1G2 = ZeroVector(3);
    Vector deriv_cross_r_G1G2_1 = ZeroVector(3);
    Vector deriv_cross_r_G1G2_2 = ZeroVector(3);

    Vector total_deriv_s_cross_G1G2 = ZeroVector(3);
    Vector deriv_cross_s_G1G2_1 = ZeroVector(3);
    Vector deriv_cross_s_G1G2_2 = ZeroVector(3);

    Vector total_deriv2_cross_G1G2 = ZeroVector(3);
    Vector deriv2_cross_G1G2_1 = ZeroVector(3);
    Vector deriv2_cross_G1G2_2 = ZeroVector(3);

    double norm_cross_g1g2 = 0.0;
    double norm_cross_G1G2 = 0.0;

    double deriv2_norm_cross_G1G2 = 0.0;
    double deriv2_norm_cross_G1G2_1 = 0.0;
    double deriv2_norm_cross_G1G2_2 = 0.0;

    rDeriv2InvDetDeformationGradient = 0.0;

    MathUtils<double>::CrossProduct(cross_g1g2, rCurrentCovariantBaseVectors[0], rCurrentCovariantBaseVectors[1]);
    norm_cross_g1g2 = MathUtils<double>::Norm(cross_g1g2);

    MathUtils<double>::CrossProduct(cross_G1G2, rReferenceCovariantBaseVectors[0], rReferenceCovariantBaseVectors[1]);
    norm_cross_G1G2 = MathUtils<double>::Norm(cross_G1G2);

    this->DeriveCurrentCovariantBaseVectors(derivative_r_reference_covariant_base_vectors, rShapeFunctionGradientValues, DofR);
    this->DeriveCurrentCovariantBaseVectors(derivative_s_reference_covariant_base_vectors, rShapeFunctionGradientValues, DofS);

    MathUtils<double>::CrossProduct(deriv_cross_r_G1G2_1, derivative_r_reference_covariant_base_vectors[0], rReferenceCovariantBaseVectors[1]);
    MathUtils<double>::CrossProduct(deriv_cross_r_G1G2_2, rReferenceCovariantBaseVectors[0], derivative_r_reference_covariant_base_vectors[1]);
    total_deriv_r_cross_G1G2 = deriv_cross_r_G1G2_1 + deriv_cross_r_G1G2_2;

    MathUtils<double>::CrossProduct(deriv_cross_s_G1G2_1, derivative_s_reference_covariant_base_vectors[0], rReferenceCovariantBaseVectors[1]);
    MathUtils<double>::CrossProduct(deriv_cross_s_G1G2_2, rReferenceCovariantBaseVectors[0], derivative_s_reference_covariant_base_vectors[1]);
    total_deriv_s_cross_G1G2 = deriv_cross_s_G1G2_1 + deriv_cross_s_G1G2_2;

    MathUtils<double>::CrossProduct(deriv2_cross_G1G2_1, derivative_r_reference_covariant_base_vectors[0], derivative_s_reference_covariant_base_vectors[1]);
    MathUtils<double>::CrossProduct(deriv2_cross_G1G2_2, derivative_s_reference_covariant_base_vectors[0], derivative_r_reference_covariant_base_vectors[1]);
    total_deriv2_cross_G1G2 = deriv2_cross_G1G2_1 + deriv2_cross_G1G2_2;

    deriv2_norm_cross_G1G2_1 = (inner_prod(total_deriv_s_cross_G1G2, total_deriv_r_cross_G1G2) + inner_prod(cross_G1G2, total_deriv2_cross_G1G2))/ norm_cross_G1G2;
    deriv2_norm_cross_G1G2_2 = (inner_prod(cross_G1G2, total_deriv_r_cross_G1G2) * inner_prod(cross_G1G2, total_deriv_s_cross_G1G2))/ (norm_cross_G1G2 * norm_cross_G1G2 * norm_cross_G1G2);
    deriv2_norm_cross_G1G2 = deriv2_norm_cross_G1G2_1 - deriv2_norm_cross_G1G2_2;


    rDeriv2InvDetDeformationGradient = deriv2_norm_cross_G1G2 / norm_cross_g1g2;

  }


  void MembraneCuttingPatternElement::PreStress(Matrix& rPreStress, const array_1d<Vector, 2>& rTransformedBaseVectors)
  {

    rPreStress = ZeroMatrix(2);
     
    Vector rPreStress_vector = ZeroVector(3);
     if (GetProperties().Has(PRESTRESS_VECTOR)) {
       rPreStress_vector = GetProperties()(PRESTRESS_VECTOR);

       if (Has(LOCAL_PRESTRESS_AXIS_1) && Has(LOCAL_PRESTRESS_AXIS_2)) {

         array_1d<array_1d<double, 3>, 2> local_prestress_axis;
         local_prestress_axis[0] = GetValue(LOCAL_PRESTRESS_AXIS_1) / MathUtils<double>::Norm(GetValue(LOCAL_PRESTRESS_AXIS_1));
         local_prestress_axis[1] = GetValue(LOCAL_PRESTRESS_AXIS_2) / MathUtils<double>::Norm(GetValue(LOCAL_PRESTRESS_AXIS_2));

         Matrix transformation_matrix = ZeroMatrix(3);
         InPlaneTransformationMatrix(transformation_matrix, rTransformedBaseVectors, local_prestress_axis);
         rPreStress_vector = prod(transformation_matrix, rPreStress_vector);

       }
       else if (Has(LOCAL_PRESTRESS_AXIS_1)) {

         Vector base_3 = ZeroVector(3);
         MathUtils<double>::UnitCrossProduct(base_3, rTransformedBaseVectors[0], rTransformedBaseVectors[1]);

         array_1d<array_1d<double, 3>, 2> local_prestress_axis;
         local_prestress_axis[0] = GetValue(LOCAL_PRESTRESS_AXIS_1) / MathUtils<double>::Norm(GetValue(LOCAL_PRESTRESS_AXIS_1));

         MathUtils<double>::UnitCrossProduct(local_prestress_axis[1], base_3, local_prestress_axis[0]);

         Matrix transformation_matrix = ZeroMatrix(3);
         InPlaneTransformationMatrix(transformation_matrix, rTransformedBaseVectors, local_prestress_axis);
         rPreStress_vector = prod(transformation_matrix, rPreStress_vector);
       }

       rPreStress = MathUtils<double>::StressVectorToTensor(rPreStress_vector);
     }
  }


  void MembraneCuttingPatternElement::save(Serializer& rSerializer) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    rSerializer.save("mConstitutiveLawVector", mConstitutiveLawVector);
  }

  void MembraneCuttingPatternElement::load(Serializer& rSerializer)
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    rSerializer.load("mConstitutiveLawVector", mConstitutiveLawVector);
  }



}