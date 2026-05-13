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

  void MembraneCuttingPatternElement::comp_Relaxation(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) {

    // KRATOS_TRY;

    // const auto& r_geom = GetGeometry();
    // const SizeType dimension = r_geom.WorkingSpaceDimension();
    // const SizeType number_of_nodes = r_geom.size();
    // const SizeType number_dofs = dimension * number_of_nodes;

    // const IntegrationMethod integration_method = r_geom.GetDefaultIntegrationMethod();


    // Vector rRightHandSideVector = ZeroVector(number_dofs);

    // Vector InternalForces = ZeroVector(3);
    // this->InternalForces(InternalForces, integration_method, rCurrentProcessInfo);

    // rRightHandSideVector = -InternalForces;


    // Matrix rLeftHandSideMatrix = ZeroMatrix(number_dofs);
    // this->TotalStiffnessMatrix(rLeftHandSideMatrix, integration_method, rCurrentProcessInfo);


    // KRATOS_CATCH("");
  }


  

  void MembraneCuttingPatternElement::InternalForces_Least_Square(Vector& rInternalForces, const IntegrationMethod& ThisMethod, const ProcessInfo& rCurrentProcessInfo) {


    // const auto& r_geom = GetGeometry();
    // const SizeType dimension = r_geom.WorkingSpaceDimension();
    // const SizeType number_of_nodes = r_geom.size();
    // const SizeType number_dofs = dimension * number_of_nodes;
    // rInternalForces = ZeroVector(number_dofs);

    // const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = r_geom.ShapeFunctionsLocalGradients(ThisMethod);
    // const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geom.IntegrationPoints(ThisMethod);

    // const double thickness = GetProperties()[THICKNESS];

    // array_1d<Vector, 2> current_covariant_base_vectors;
    // array_1d<Vector, 2> current_contravariant_base_vectors;
    // array_1d<Vector, 2> reference_covariant_base_vectors;
    // array_1d<Vector, 2> reference_contravariant_base_vectors;

    // array_1d<Vector, 2> transformed_base_vectors;
    // array_1d<Vector, 2> transformed_base_vectors_current;

    // Matrix covariant_metric_current = ZeroMatrix(3);
    // Matrix covariant_metric_reference = ZeroMatrix(3);
    // Matrix contravariant_metric_reference = ZeroMatrix(3);
    // Matrix contravariant_metric_current = ZeroMatrix(3);
    // Matrix inplane_transformation_matrix_material = ZeroMatrix(3);
    // Matrix inplane_transformation_matrix_material_current = ZeroMatrix(3);

    // Matrix material_tangent_modulus = ZeroMatrix(dimension);
    // Matrix material_tangent_modulus_Cauchy = ZeroMatrix(dimension);
    // double detJ = 0.0;
    // Vector stress = ZeroVector(3);
    // Vector sigma_pre = ZeroVector(3);
    // Vector derivative_strain = ZeroVector(3);

    // Matrix deformation_gradient = ZeroMatrix(3);
    // double det_deformation_gradient = 0.0;
    

    // for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
    //   // getting information for integration
    //   const double integration_weight_i = r_integration_points[point_number].Weight();
    //   const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];

    //   this->CovariantBaseVectors(current_covariant_base_vectors, shape_functions_gradients_i, ConfigurationType::Current);
    //   this->CovariantBaseVectors(reference_covariant_base_vectors, shape_functions_gradients_i, ConfigurationType::Reference);

    //   this->CovariantMetric(covariant_metric_current, current_covariant_base_vectors);
    //   this->CovariantMetric(covariant_metric_reference, reference_covariant_base_vectors);
    //   this->ContravariantMetric(contravariant_metric_reference, covariant_metric_reference);
    //   this->ContravariantMetric(contravariant_metric_current, covariant_metric_current);

    //   this->ContraVariantBaseVectors(reference_contravariant_base_vectors, contravariant_metric_reference, reference_covariant_base_vectors);
    //   this->ContraVariantBaseVectors(current_contravariant_base_vectors, contravariant_metric_current, current_covariant_base_vectors);

    //   this->TransformBaseVectors(transformed_base_vectors, reference_contravariant_base_vectors);//****** transforming reference contravariant base vectors to local cartesian systems *******
    //   this->TransformBaseVectors(transformed_base_vectors_current, current_contravariant_base_vectors);

    //   this->InPlaneTransformationMatrix(inplane_transformation_matrix_material, transformed_base_vectors, reference_contravariant_base_vectors);//***** inner product betn local cartesian of ref_contra, make the transformation matrix applicable for Voigt notation *****
    //   this->InPlaneTransformationMatrix(inplane_transformation_matrix_material_current, transformed_base_vectors_current, current_contravariant_base_vectors);

    //   this->JacobiDeterminante(detJ, current_covariant_base_vectors);

    //   this->MaterialResponse(stress, contravariant_metric_reference, covariant_metric_reference, covariant_metric_current,
    //     transformed_base_vectors, inplane_transformation_matrix_material, point_number, material_tangent_modulus,
    //     rCurrentProcessInfo); // uses Green-Lagrange strain tensor to compute PK2 and C

    //   this->DeformationGradient(deformation_gradient, det_deformation_gradient, current_covariant_base_vectors, reference_contravariant_base_vectors);

    //   Matrix stress_matrix_local_cs = MathUtils<double>::StressVectorToTensor(stress);

    //   material_tangent_modulus_Cauchy = (1.0 / det_deformation_gradient) * material_tangent_modulus; //****** this operation might not be suitable for tangent modulus for vectors in Voigt notation *******

    //   //transform stresses to original bases
    //   Matrix stress_matrix = ZeroMatrix(3);
    //   for (SizeType i = 0; i < 2; ++i) {
    //     for (SizeType j = 0; j < 2; ++j) {
    //       stress_matrix += outer_prod(transformed_base_vectors[i], transformed_base_vectors[j]) * stress_matrix_local_cs(i, j);
    //     }
    //   }

    //   // calculate cauchy (this needs to be done in the original base)
    //   Matrix temp_stress_matrix = prod(deformation_gradient, stress_matrix);
    //   Matrix temp_stress_matrix_2 = prod(temp_stress_matrix, trans(deformation_gradient));
    //   Matrix cauchy_stress_matrix = temp_stress_matrix_2 / det_deformation_gradient; // Dieringer - equation 2.57, page 39

    //   // transform stresses to local orthogonal base
    //   Matrix local_stress = ZeroMatrix(2);
    //   for (SizeType i = 0; i < 2; ++i) {
    //     for (SizeType j = 0; j < 2; ++j) {
    //       local_stress(i, j) = inner_prod(transformed_base_vectors_current[i], prod(cauchy_stress_matrix, transformed_base_vectors_current[j]));
    //     }
    //   }
    //   stress = MathUtils<double>::StressTensorToVector(local_stress, 3);

    //   this->PreStress(sigma_pre, transformed_base_vectors_current); //****** is it necessary to convert pre-stress as PK2 is converted to Cauchy? ********

    //   for (SizeType dof_r = 0; dof_r < number_dofs; ++dof_r)
    //   {
    //     DerivativeStrainEulerAlmansi(derivative_strain, shape_functions_gradients_i,
    //       dof_r, reference_covariant_base_vectors, inplane_transformation_matrix_material_current);

    //     Vector stress_derivative = prod(material_tangent_modulus_Cauchy, derivative_strain); // ********* derivative of inverse det_deformation_gradient is not zero!!! **********

    //     Vector delta_sigma = stress - sigma_pre;

    //     rInternalForces[dof_r] += inner_prod(delta_sigma, stress_derivative) * detJ * integration_weight_i * thickness;
    //   }
    // }
  }


  void MembraneCuttingPatternElement::ResponseFunction_Least_Square(double& rResponseLS, const IntegrationMethod& ThisMethod, const ProcessInfo& rCurrentProcessInfo) {

    // const auto& r_geom = GetGeometry();
    // const SizeType dimension = r_geom.WorkingSpaceDimension();

    // const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = r_geom.ShapeFunctionsLocalGradients(ThisMethod);
    // const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geom.IntegrationPoints(ThisMethod);

    // const double thickness = GetProperties()[THICKNESS];

    // array_1d<Vector, 2> current_covariant_base_vectors;
    // array_1d<Vector, 2> current_contravariant_base_vectors;
    // array_1d<Vector, 2> reference_covariant_base_vectors;
    // array_1d<Vector, 2> reference_contravariant_base_vectors;

    // array_1d<Vector, 2> transformed_base_vectors;
    // array_1d<Vector, 2> transformed_base_vectors_current;

    // Matrix covariant_metric_current = ZeroMatrix(3);
    // Matrix covariant_metric_reference = ZeroMatrix(3);
    // Matrix contravariant_metric_reference = ZeroMatrix(3);
    // Matrix contravariant_metric_current = ZeroMatrix(3);
    // Matrix inplane_transformation_matrix_material = ZeroMatrix(3);
    
    // Matrix material_tangent_modulus = ZeroMatrix(dimension);
    
    // double detJ = 0.0;
    // Vector stress = ZeroVector(3);
    // Vector sigma_pre = ZeroVector(3);

    // Matrix deformation_gradient = ZeroMatrix(3);
    // double det_deformation_gradient = 0.0;

    // rResponseLS = 0.0;


    // for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
    //   // getting information for integration
    //   const double integration_weight_i = r_integration_points[point_number].Weight();
    //   const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];

    //   this->CovariantBaseVectors(current_covariant_base_vectors, shape_functions_gradients_i, ConfigurationType::Current);
    //   this->CovariantBaseVectors(reference_covariant_base_vectors, shape_functions_gradients_i, ConfigurationType::Reference);

    //   this->CovariantMetric(covariant_metric_current, current_covariant_base_vectors);
    //   this->CovariantMetric(covariant_metric_reference, reference_covariant_base_vectors);
    //   this->ContravariantMetric(contravariant_metric_reference, covariant_metric_reference);
    //   this->ContravariantMetric(contravariant_metric_current, covariant_metric_current);

    //   this->ContraVariantBaseVectors(reference_contravariant_base_vectors, contravariant_metric_reference, reference_covariant_base_vectors);
    //   this->ContraVariantBaseVectors(current_contravariant_base_vectors, contravariant_metric_current, current_covariant_base_vectors);

    //   this->TransformBaseVectors(transformed_base_vectors, reference_contravariant_base_vectors);//****** transforming reference contravariant base vectors to local cartesian systems *******
    //   this->TransformBaseVectors(transformed_base_vectors_current, current_contravariant_base_vectors);

    //   this->InPlaneTransformationMatrix(inplane_transformation_matrix_material, transformed_base_vectors, reference_contravariant_base_vectors);//***** inner product betn local cartesian of ref_contra, make the transformation matrix applicable for Voigt notation *****

    //   this->JacobiDeterminante(detJ, current_covariant_base_vectors);

    //   this->MaterialResponse(stress, contravariant_metric_reference, covariant_metric_reference, covariant_metric_current,
    //     transformed_base_vectors, inplane_transformation_matrix_material, point_number, material_tangent_modulus,
    //     rCurrentProcessInfo); // use Green-Lagrange strain tensor to compute PK2 and C

    //   this->DeformationGradient(deformation_gradient, det_deformation_gradient, current_covariant_base_vectors, reference_contravariant_base_vectors);

    //   Matrix stress_matrix_local_cs = MathUtils<double>::StressVectorToTensor(stress);

    //   //transform stresses to original bases
    //   Matrix stress_matrix = ZeroMatrix(3);
    //   for (SizeType i = 0; i < 2; ++i) {
    //     for (SizeType j = 0; j < 2; ++j) {
    //       stress_matrix += outer_prod(transformed_base_vectors[i], transformed_base_vectors[j]) * stress_matrix_local_cs(i, j);
    //     }
    //   }

    //   // calculate cauchy (this needs to be done in the original base)
    //   Matrix temp_stress_matrix = prod(deformation_gradient, stress_matrix);
    //   Matrix temp_stress_matrix_2 = prod(temp_stress_matrix, trans(deformation_gradient));
    //   Matrix cauchy_stress_matrix = temp_stress_matrix_2 / det_deformation_gradient; // Dieringer - equation 2.57, page 39

    //   // transform stresses to local orthogonal base
    //   Matrix local_stress = ZeroMatrix(2);
    //   for (SizeType i = 0; i < 2; ++i) {
    //     for (SizeType j = 0; j < 2; ++j) {
    //       local_stress(i, j) = inner_prod(transformed_base_vectors_current[i], prod(cauchy_stress_matrix, transformed_base_vectors_current[j]));
    //     }
    //   }
    //   stress = MathUtils<double>::StressTensorToVector(local_stress, 3);

    //   this->PreStress(sigma_pre, transformed_base_vectors_current);

    //   Vector delta_sigma = stress - sigma_pre;
    //   double product_delta_sigma = ((delta_sigma[0] * delta_sigma[0]) + (delta_sigma[1] * delta_sigma[1]) + (2.0 * delta_sigma[2] * delta_sigma[2])) * detJ * integration_weight_i * thickness; 

    //   rResponseLS += 0.5 * product_delta_sigma;
        
    // }
  }



  void MembraneCuttingPatternElement::StrainEulerAlmansi(Vector& rStrain, const Matrix& rReferenceCoVariantMetric, const Matrix& rCurrentCoVariantMetric,
    const Matrix& rTransformationMatrix)
  {
    // Matrix strain_matrix = 0.50 * (rCurrentCoVariantMetric - rReferenceCoVariantMetric);
    // Vector current_strain = MathUtils<double>::StrainTensorToVector(strain_matrix, 3);
    // this->TransformStrains(rStrain, current_strain, rTransformationMatrix); //transforms to local Cartesian coordinate system
  }


  void MembraneCuttingPatternElement::DerivativeStrainEulerAlmansi(Vector& rStrain, const Matrix& rShapeFunctionGradientValues, const SizeType DofR,
    const array_1d<Vector, 2> rReferenceCovariantBaseVectors, const Matrix& rTransformationMatrix)
  {
    // Matrix reference_covariant_metric_derivative = ZeroMatrix(2);
    // this->DerivativeCurrentCovariantMetric(reference_covariant_metric_derivative, rShapeFunctionGradientValues, DofR, rReferenceCovariantBaseVectors);

    // Matrix strain_matrix_derivative = -0.50 * reference_covariant_metric_derivative;
    // Vector current_strain = MathUtils<double>::StrainTensorToVector(strain_matrix_derivative, 3);
    // TransformStrains(rStrain, current_strain, rTransformationMatrix);
  }



  void MembraneCuttingPatternElement::Elasticity_Tensor_Kirchhoff(double rC_reference[2][2][2][2], double rc_current[2][2][2][2], const IntegrationMethod& ThisMethod, const double YoungModulus, const double PoissonRatio) {


    // const auto& r_geom = GetGeometry();
    // const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = r_geom.ShapeFunctionsLocalGradients(ThisMethod);
    // const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geom.IntegrationPoints(ThisMethod);
    
    // array_1d<Vector, 2> reference_covariant_base_vectors;
    // array_1d<Vector, 2> reference_contravariant_base_vectors;
    // array_1d<Vector, 2> current_covariant_base_vectors;
    
    // Matrix covariant_metric_reference = ZeroMatrix(3);
    // Matrix contravariant_metric_reference = ZeroMatrix(3);

    // Matrix deformation_gradient = ZeroMatrix(3);
    // double det_deformation_gradient = 0.0;

    // double E = YoungModulus;
    // double NU = PoissonRatio;
    // const double lambda = E * NU / (1.0 - NU * NU);
    // const double MU = E / (2.0 * (1.0 + NU));


    // for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number) {

    //   const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];

    //   this->CovariantBaseVectors(current_covariant_base_vectors, shape_functions_gradients_i, ConfigurationType::Current);
    //   this->CovariantBaseVectors(reference_covariant_base_vectors, shape_functions_gradients_i, ConfigurationType::Reference);
      
    //   this->CovariantMetric(covariant_metric_reference, reference_covariant_base_vectors);
    //   this->ContravariantMetric(contravariant_metric_reference, covariant_metric_reference);

    //   this->ContraVariantBaseVectors(reference_contravariant_base_vectors, contravariant_metric_reference, reference_covariant_base_vectors);
      
    //   this->DeformationGradient(deformation_gradient, det_deformation_gradient, current_covariant_base_vectors, reference_contravariant_base_vectors);

    //   for (SizeType alpha = 0; alpha < 2; ++alpha) {
    //     for (SizeType beta = 0; beta < 2; ++beta) {
    //       for (SizeType gamma = 0; gamma < 2; ++gamma) {
    //         for (SizeType delta = 0; delta < 2; ++delta) {
    //           rC_reference[alpha][beta][gamma][delta] =
    //             (lambda * (contravariant_metric_reference(alpha, beta) * contravariant_metric_reference(gamma, delta))) +
    //             (MU * ((contravariant_metric_reference(alpha, gamma) * contravariant_metric_reference(beta, delta)) +
    //               (contravariant_metric_reference(alpha, delta) * contravariant_metric_reference(beta, gamma))));

    //           rc_current[alpha][beta][gamma][delta] = (1 / det_deformation_gradient) * rC_reference[alpha][beta][gamma][delta];

    //         }
    //       }
    //     }
    //   }
    // }
         
  }


  void MembraneCuttingPatternElement::PreStress(Vector& rPreStress, const array_1d<Vector, 2>& rTransformedBaseVectors)//TO DO
  {

    // rPreStress = ZeroVector(3);
    // if (GetProperties().Has(PRESTRESS_VECTOR)) {
    //   rPreStress = GetProperties()(PRESTRESS_VECTOR);

    //   if (Has(LOCAL_PRESTRESS_AXIS_1) && Has(LOCAL_PRESTRESS_AXIS_2)) {

    //     array_1d<array_1d<double, 3>, 2> local_prestress_axis;
    //     local_prestress_axis[0] = GetValue(LOCAL_PRESTRESS_AXIS_1) / MathUtils<double>::Norm(GetValue(LOCAL_PRESTRESS_AXIS_1));
    //     local_prestress_axis[1] = GetValue(LOCAL_PRESTRESS_AXIS_2) / MathUtils<double>::Norm(GetValue(LOCAL_PRESTRESS_AXIS_2));

    //     Matrix transformation_matrix = ZeroMatrix(3);
    //     InPlaneTransformationMatrix(transformation_matrix, rTransformedBaseVectors, local_prestress_axis);
    //     rPreStress = prod(transformation_matrix, rPreStress);

    //   }
    //   else if (Has(LOCAL_PRESTRESS_AXIS_1)) {

    //     Vector base_3 = ZeroVector(3);
    //     MathUtils<double>::UnitCrossProduct(base_3, rTransformedBaseVectors[0], rTransformedBaseVectors[1]);

    //     array_1d<array_1d<double, 3>, 2> local_prestress_axis;
    //     local_prestress_axis[0] = GetValue(LOCAL_PRESTRESS_AXIS_1) / MathUtils<double>::Norm(GetValue(LOCAL_PRESTRESS_AXIS_1));

    //     MathUtils<double>::UnitCrossProduct(local_prestress_axis[1], base_3, local_prestress_axis[0]);

    //     Matrix transformation_matrix = ZeroMatrix(3);
    //     InPlaneTransformationMatrix(transformation_matrix, rTransformedBaseVectors, local_prestress_axis);
    //     rPreStress = prod(transformation_matrix, rPreStress);
    //   }
    // }
  }


  ////***********************************************************************************
  //// comp_Optimization_Least_Square - calculation of the optimized cutting patterning
  ////***********************************************************************************

  void MembraneCuttingPatternElement::comp_Optimization_Least_Square(
    Matrix& rLeftHandSideMatrix,
    Vector& rRightHandSideVector,
    ConstitutiveLaw::Parameters& rValues,
    const Properties& rMaterialProperties,
    double& rArea3DElem,
    double& rArea2DElem,
    double& rResponse,
    const Matrix& rFiberDirectionsRef,
    const ProcessInfo& rCurrentProcessInfo)
  {
    // KRATOS_TRY;

    


    // const auto& r_geom = GetGeometry();
    // SizeType dimension = r_geom.WorkingSpaceDimension();
    // SizeType number_of_nodes = r_geom.size();
    // SizeType number_dofs  = dimension * number_of_nodes; 

    // IntegrationMethod integration_method = r_geom.GetDefaultIntegrationMethod();

    // rArea3DElem = 0.0;
    // rArea2DElem = 0.0;
    // rResponse = 0.0;
    // rLeftHandSideMatrix = ZeroMatrix(number_dofs);
    // rRightHandSideVector.resize(number_dofs);


    // Vector InternalForces_LestSquare = ZeroVector(3);
    // this->InternalForces_Least_Square(InternalForces_LestSquare, integration_method, rCurrentProcessInfo);

    // rResponse = 0.0;
    // this->ResponseFunction_Least_Square(rResponse, integration_method, rCurrentProcessInfo);


    // //get all integration points

    
    // const GeometryType::IntegrationPointsArrayType& r_integrations_points = r_geom.IntegrationPoints(integration_method);
    
    // //const auto& N_values = r_geom.ShapeFunctionsValues(integration_method);
    // const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = r_geom.ShapeFunctionsLocalGradients(integration_method);


    // // ---------------------------------------------------------------------
    // // Material parameters for full continuum mechanics description
    // // ---------------------------------------------------------------------
    // const double thickness = GetProperties()[THICKNESS];

    // //const bool is_orthotropic = GetProperties().Has(ORTHOTROPIC_LAW); //***CHECK LATER ***
    // //const bool is_isotropic = !is_orthotropic; //***CHECK LATER ***

    // // Isotropic material parameters

    // //if (is_isotropic)//***CHECK LATER *** 
    // //{

    //   //linear_plane_stress.cpp, CalcualteElasticMatrixPlaneStresss from constitutive_law_utilities.cpp cannot be used because it calculates C for small Deformation
    // const auto& r_props = rValues.GetMaterialProperties();
    // double E = r_props[YOUNG_MODULUS];
    // double NU = r_props[POISSON_RATIO];
    // double lambda = E * NU / (1.0 - NU * NU);
    // const double MU = E / (2.0 * (1.0 + NU));
    // //}

    // // Orthotropic parameters (Muensch-Reinhardt style)

    // //if (is_orthotropic)//***CHECK LATER *** 
    // //{

    // //  //advanced_constitutive_law_utilities.cpp, why not using CalculateOrthotropicElasticMatrix?!
    // //  const Vector& r_ortho_elastic_constants = rMaterialProperties[ORTHOTROPIC_ELASTIC_CONSTANTS];
    // //  const double E1 = r_ortho_elastic_constants[0];
    // //  const double E2 = r_ortho_elastic_constants[1];
    // //  const double NU12 = r_ortho_elastic_constants[3];
    // //  //const double G = ConstitutiveLawUtilities<3>::CalculateShearModulus(r_material_properties);
    // //  double G12 = (rMaterialProperties.Has(SHEAR_MODULUS_XY)) ? rMaterialProperties[SHEAR_MODULUS_XY] : 1.0 / ((1.0 + NU21) / E1 + (1.0 + NU12) / E2);

    // //  const double NU21 = NU12 * E2 / E1;

    // //  KRATOS_ERROR_IF(NU21 > 0.5) << "The Poisson_yx is greater than 0.5." << std::endl;

    // //  double coeff_C[2][2][2][2] = { {{{0.0}}} };

    // //  coeff_C[0][0][0][0] = 1.0 / (1.0 - NU12 * NU21) * E1;
    // //  coeff_C[0][0][1][1] = 1.0 / (1.0 - NU12 * NU21) * NU12 * E1;
    // //  coeff_C[0][1][0][1] = G12;
    // //  coeff_C[0][1][1][0] = G12;
    // //  coeff_C[1][0][0][1] = G12;
    // //  coeff_C[1][0][1][0] = G12;
    // //  coeff_C[1][1][0][0] = 1.0 / (1.0 - NU12 * NU21) * NU21 * E2;
    // //  coeff_C[1][1][1][1] = 1.0 / (1.0 - NU12 * NU21) * E2;

    // //}
    
    
   

    // // ---------------------------------------------------------------------
    // // Prestress in Cartesian basis on current 3D surface
    // // ---------------------------------------------------------------------
    // double sigma11_pre = 0.0;
    // double sigma22_pre = 0.0;
    // double sigma12_pre = 0.0;


    // // ======================================================================
    // // Gauss loop
    // // ======================================================================
    // for (IndexType t = 0; t < r_integrations_points.size(); ++t)
    // {
    //   // --------------------------------------------------------------------
    //   // Get prestress information for each Gauss point
    //   // assuming constant prestress over the whole element
    //   // ---------------------------------------------------------------------

    //   Vector pre_stress = ZeroVector(3);
    //   if (GetProperties().Has(PRESTRESS_VECTOR)) {

    //     pre_stress = GetProperties()(PRESTRESS_VECTOR);

    //     sigma11_pre += pre_stress[0];
    //     sigma22_pre += pre_stress[1];
    //     sigma12_pre += pre_stress[2];
    //   }


    //   const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[t];
    //   double integration_weight_i = r_integrations_points[t].Weight();

    //   // ------------------------------------------------------------
    //   // 1. Covariant base vectors in 3D actual configuration
    //   // ------------------------------------------------------------

    //   array_1d<Vector, 2> current_covariant_base_vectors;
    //   this->CovariantBaseVectors(current_covariant_base_vectors, shape_functions_gradients_i, ConfigurationType::Current);

    //   // ------------------------------------------------------------
    //   // 2. Covariant base vectors in reference configuration
    //   // ------------------------------------------------------------
    //   array_1d<Vector, 2> reference_covariant_base_vectors;
    //   this->CovariantBaseVectors(reference_covariant_base_vectors, shape_functions_gradients_i, ConfigurationType::Reference);

    //   // ------------------------------------------------------------
    //   // 3. Covariant metric in the current confliguration 
    //   // ------------------------------------------------------------
    //   Matrix covariant_metric_current = ZeroMatrix(3);
    //   this->CovariantMetric(covariant_metric_current, current_covariant_base_vectors);

    //   // ------------------------------------------------------------
    //   // 4. Covariant metric in the reference configuration covariant
    //   // ------------------------------------------------------------
    //   Matrix covariant_metric_reference = ZeroMatrix(3);
    //   this->CovariantMetric(covariant_metric_reference, reference_covariant_base_vectors);

    //   // ------------------------------------------------------------
    //   // 5. Contravariant metric in the 3D actual configuration
    //   // ------------------------------------------------------------
    //   Matrix contravariant_metric_current = ZeroMatrix(3);
    //   this->ContravariantMetric(contravariant_metric_current, covariant_metric_current);

    //   // ------------------------------------------------------------
    //   // 6. Contravariant metric in reference configuration 
    //   // ------------------------------------------------------------
    //   Matrix contravariant_metric_reference = ZeroMatrix(3);
    //   this->ContravariantMetric(contravariant_metric_reference, covariant_metric_reference);

    //   // ---------------------------------------------------------------------
    //   // 7. Determinant of covariant metric in reference configuration
    //   // ---------------------------------------------------------------------
    //   double det_covariant_metric_reference = (covariant_metric_reference(0, 0) * covariant_metric_reference(1, 1)) - (covariant_metric_reference(1, 0) * covariant_metric_reference(0, 1));
    //   double det_contravariant_metric_reference = 1.0 / det_covariant_metric_reference;

    //   // ------------------------------------------------------------
    //   // 8. Contravariant base vectors in 3D actual configuration
    //   // ------------------------------------------------------------
    //   array_1d<Vector, 2> current_contravariant_base_vectors;
    //   this->ContraVariantBaseVectors(current_contravariant_base_vectors, contravariant_metric_current, current_covariant_base_vectors);

    //   // ------------------------------------------------------------
    //   // 9. Contravariant base vectors in reference configuration
    //   // ------------------------------------------------------------
    //   array_1d<Vector, 2> reference_contravariant_base_vectors;
    //   this->ContraVariantBaseVectors(reference_contravariant_base_vectors, contravariant_metric_reference, reference_covariant_base_vectors);

    //   // ------------------------------------------------------------
    //   // 10. Transformation matrix for orthotropic material
    //   // ***** rFiberDirectionsReference has to be defined in CuttingPatternStrategy ******
    //   // ------------------------------------------------------------
    //   //Matrix transformation_matrix_orthotropic = ZeroMatrix(2);
    //   //this->comp_Transformation_Matrix(transformation_matrix_orthotropic, reference_contravariant_base_vectors, rFiberDirectionsReference);

    //   // ---------------------------------------------------------------------------------
    //   // 11. calculation of local cartesian coordinate system in the current configuration
    //   // ---------------------------------------------------------------------------------
    //   array_1d<Vector, 2> current_Euclidean_base_vectors;
    //   this->TransformBaseVectors(current_Euclidean_base_vectors, current_covariant_base_vectors);


    //   // ---------------------------------------------------------------------------------
    //   // 12. Transformation between current contravariant and cartesian
    //   // ---------------------------------------------------------------------------------
    //   Matrix transformation_matrix_current_contra = ZeroMatrix(2);
    //   this->comp_Transformation_Matrix(transformation_matrix_current_contra, current_contravariant_base_vectors, current_Euclidean_base_vectors);
    //   const double g1e1 = transformation_matrix_current_contra(0, 0);
    //   const double g1e2 = transformation_matrix_current_contra(0, 1);
    //   const double g2e1 = transformation_matrix_current_contra(1, 0);
    //   const double g2e2 = transformation_matrix_current_contra(1, 1);

    //   // ------------------------------------------------------------
    //   // 13. Elemental area in current configuration
    //   // ------------------------------------------------------------
    //   array_1d<double, 3> cross_g1g2_cov;
    //   cross_g1g2_cov.clear();
    //   MathUtils<double>::CrossProduct(cross_g1g2_cov, current_covariant_base_vectors[0], current_covariant_base_vectors[1]);
    //   double VectorNorm_g1g2_cov = MathUtils<double>::Norm(cross_g1g2_cov);

    //   double da_current = VectorNorm_g1g2_cov;
    //   rArea3DElem += da_current * integration_weight_i;

    //   // ------------------------------------------------------------
    //   // 14. Elemental area in reference configuration
    //   // ------------------------------------------------------------
    //   array_1d<double, 3> cross_G1G2_cov;
    //   MathUtils<double>::CrossProduct(cross_G1G2_cov, reference_covariant_base_vectors[0], reference_covariant_base_vectors[1]);
    //   double VectorNorm_G1G2_cov = MathUtils<double>::Norm(cross_G1G2_cov);

    //   double dA_reference = VectorNorm_G1G2_cov;
    //   rArea2DElem += dA_reference * integration_weight_i;

    //   // ------------------------------------------------------------
    //   // 15. Determinant of deformation gradient and it's inverse
    //   // ------------------------------------------------------------
    //   double DetF = da_current / dA_reference;
    //   double InverseDetF = 1.0 / DetF;

    //   // ------------------------------------------------------------
    //   // 16. Euler-Almansi-Strain-Tensor
    //   // ------------------------------------------------------------
    //   Matrix e = ZeroMatrix(2);
    //   e(0, 0) = 0.5 * (covariant_metric_current(0, 0) - covariant_metric_reference(0, 0));
    //   e(0, 1) = 0.5 * (covariant_metric_current(0, 1) - covariant_metric_reference(0, 1));
    //   e(1, 0) = 0.5 * (covariant_metric_current(1, 0) - covariant_metric_reference(1, 0));
    //   e(1, 1) = 0.5 * (covariant_metric_current(1, 1) - covariant_metric_reference(1, 1));

    //   // ------------------------------------------------------------
    //   // 17. Elastic Cauchy stress in current configuration
    //   // ------------------------------------------------------------
    //   Matrix sigma = ZeroMatrix(2);

    //   double C_reference[2][2][2][2] = { {{{0.0}}} };
    //   double c_current[2][2][2][2] = { {{{0.0}}} };

    //   for (SizeType alpha = 0; alpha < 2; ++alpha) {
    //     for (SizeType beta = 0; beta < 2; ++beta) {
    //       for (SizeType gamma = 0; gamma < 2; ++gamma) {
    //         for (SizeType delta = 0; delta < 2; ++delta) {

    //           //if (is_isotropic) //*********** need to check whether 'is_isotropic' would really work ************
    //           //{
    //             C_reference[alpha][beta][gamma][delta] =
    //               (lambda * (contravariant_metric_reference(alpha, beta) * contravariant_metric_reference(gamma, delta))) +
    //               (MU * ((contravariant_metric_reference(alpha, gamma) * contravariant_metric_reference(beta, delta)) +
    //                 (contravariant_metric_reference(alpha, delta) * contravariant_metric_reference(beta, gamma))));
    //           //}
    //           /*else if (is_orthotropic) {
    //             for (SizeType eps = 0; eps < 2; ++eps) {
    //               for (SizeType zet = 0; zet < 2; ++zet) {
    //                 for (SizeType eta = 0; eta < 2; ++eta) {
    //                   for (SizeType the = 0; the < 2; ++the) {
    //                     C_reference[alpha][beta][gamma][delta] +=
    //                       coeff_C[eps][zet][eta][the] *
    //                       transformation_matrix_orthotropic(alpha, eps) *
    //                       transformation_matrix_orthotropic(beta, zet) *
    //                       transformation_matrix_orthotropic(gamma, eta) *
    //                       transformation_matrix_orthotropic(delta, the);
    //                   }
    //                 }
    //               }
    //             }
    //           }*/

    //           c_current[alpha][beta][gamma][delta] = InverseDetF * C_reference[alpha][beta][gamma][delta];

    //           sigma(alpha, beta) += c_current[alpha][beta][gamma][delta] * e(gamma, delta);
    //         }
    //       }
    //     }
    //   }

    //   // ------------------------------------------------------------------------------------------------------------------------
    //   // 18. Prestress transformed from Cartesian to covariant basis (contravariant coefficient) in the current configuration
    //   // ------------------------------------------------------------------------------------------------------------------------
    //   Matrix sigma_pre_current_contra = ZeroMatrix(2);
    //   sigma_pre_current_contra(0, 0) = sigma11_pre * (g1e1 * g1e1) + sigma12_pre * (g1e1 * g1e2)
    //     + sigma12_pre * (g1e2 * g1e1) + sigma22_pre * (g1e2 * g1e2);

    //   sigma_pre_current_contra(0, 1) = sigma11_pre * (g1e1 * g2e1) + sigma12_pre * (g1e1 * g2e2)
    //     + sigma12_pre * (g1e2 * g2e1) + sigma22_pre * (g1e2 * g2e2);

    //   sigma_pre_current_contra(1, 0) = sigma_pre_current_contra(0, 1);

    //   sigma_pre_current_contra(1, 1) = sigma11_pre * (g2e1 * g2e1) + sigma12_pre * (g2e1 * g2e2)
    //     + sigma12_pre * (g2e2 * g2e1) + sigma22_pre * (g2e2 * g2e2);


    //   // =================================================================
    //   // 19. Loop over design DOFs X_r
    //   // =================================================================

    //   for (SizeType r = 0; r < number_dofs; ++r)
    //   {

    //     // -------------------------------------------------------------------------------------
    //     // 1. first derivative of covariant base vectors in reference configuration w.r.t DOFs
    //     // -------------------------------------------------------------------------------------
    //     array_1d<Vector, 2> ref_cov_base_vectors_rDeriv;
    //     this->DeriveCurrentCovariantBaseVectors(ref_cov_base_vectors_rDeriv, shape_functions_gradients_i, r);

    //     // -------------------------------------------------------------------------------------
    //     // 2. first derivative of covariant metric in reference configuration w.r.t DOFs
    //     // -------------------------------------------------------------------------------------
    //     Matrix ref_cov_metric_rDeriv = ZeroMatrix(2);
    //     this->DerivativeCurrentCovariantMetric(ref_cov_metric_rDeriv, shape_functions_gradients_i, r, reference_covariant_base_vectors);

    //     // -------------------------------------------------------------------------------------
    //     // 3. first derivative of Euler Almansi strain tensor w.r.t DOFs
    //     // -------------------------------------------------------------------------------------
    //     Matrix e_rDeriv = ZeroMatrix(2);
    //     e_rDeriv(0, 0) = -0.5 * ref_cov_metric_rDeriv(0, 0);
    //     e_rDeriv(0, 1) = -0.5 * ref_cov_metric_rDeriv(0, 1);
    //     e_rDeriv(1, 0) = -0.5 * ref_cov_metric_rDeriv(1, 0);
    //     e_rDeriv(1, 1) = -0.5 * ref_cov_metric_rDeriv(1, 1);

    //     // -------------------------------------------------------------------------------------
    //     // 4. first derivative of determinant of reference covariant metric w.r.t DOFs
    //     // -------------------------------------------------------------------------------------
    //     double rDeriv_det_cov_metric_ref = ref_cov_metric_rDeriv(0, 0) * covariant_metric_reference(1, 1) + covariant_metric_reference(0, 0) * ref_cov_metric_rDeriv(1, 1)
    //       - ref_cov_metric_rDeriv(1, 0) * covariant_metric_reference(0, 1) - covariant_metric_reference(1, 0) * ref_cov_metric_rDeriv(0, 1);

    //     // -------------------------------------------------------------------------------------
    //     // 5. first derivative of determinant of reference contravariant metric w.r.t DOFs
    //     // -------------------------------------------------------------------------------------
    //     double rDeriv_det_contra_metric_ref = -rDeriv_det_cov_metric_ref / (det_covariant_metric_reference * det_covariant_metric_reference);

    //     // -------------------------------------------------------------------------------------
    //     // 6. first derivative of contravariant metric in reference configuration w.r.t DOFs
    //     // -------------------------------------------------------------------------------------
    //     Matrix ref_contra_metric_rDeriv = ZeroMatrix(2);

    //     ref_contra_metric_rDeriv(0, 0) = rDeriv_det_contra_metric_ref * covariant_metric_reference(1, 1) + det_contravariant_metric_reference * ref_cov_metric_rDeriv(1, 1);
    //     ref_contra_metric_rDeriv(0, 1) = -rDeriv_det_contra_metric_ref * covariant_metric_reference(0, 1) - det_contravariant_metric_reference * ref_cov_metric_rDeriv(0, 1);
    //     ref_contra_metric_rDeriv(1, 0) = -rDeriv_det_contra_metric_ref * covariant_metric_reference(1, 0) - det_contravariant_metric_reference * ref_cov_metric_rDeriv(1, 0);
    //     ref_contra_metric_rDeriv(1, 1) = rDeriv_det_contra_metric_ref * covariant_metric_reference(0, 0) + det_contravariant_metric_reference * ref_cov_metric_rDeriv(0, 0);

    //     // ----------------------------------------------------------------------------------------
    //     // 7. first derivative of contravariant base vectors in reference configuration w.r.t DOFs
    //     // ----------------------------------------------------------------------------------------
    //     array_1d<Vector, 2> ref_contra_base_vectors_rDeriv;

    //     for (SizeType alpha = 0; alpha < 2; ++alpha) {
    //       ref_contra_base_vectors_rDeriv[alpha] = ZeroVector(3);

    //       for (SizeType beta = 0; beta < 2; ++beta) {
    //         ref_contra_base_vectors_rDeriv[alpha] +=
    //           ref_contra_metric_rDeriv(alpha, beta) * reference_covariant_base_vectors[beta]
    //           + contravariant_metric_reference(alpha, beta) * ref_cov_base_vectors_rDeriv[beta];
    //       }
    //     }

    //     // ---------------------------------------------------------------------------------------
    //     // 8. first derivative of elemental area in reference configuration w.r.t DOFs
    //     // ---------------------------------------------------------------------------------------
    //     double dA_reference_rDeriv = 0.0;
    //     this->DerivativeElementalArea(dA_reference_rDeriv, reference_covariant_base_vectors, ref_cov_base_vectors_rDeriv);

    //     // ---------------------------------------------------------------------------------------
    //     // 9. first derivative of inverse of deformation gradient w.r.t DOFs
    //     // ---------------------------------------------------------------------------------------
    //     double InverseDetF_rDeriv = dA_reference_rDeriv / da_current;


    //     // ---------------------------------------------------------------------------------------
    //     // 10. first derivative of transformation matrix for orthotropic material w.r.t DOFs
    //     // ---------------------------------------------------------------------------------------
    //     //Matrix transformation_matrix_orthotropic_rDeriv = ZeroMatrix(2);
    //     //this->comp_Transformation_Matrix(transformation_matrix_orthotropic_rDeriv, ref_contra_base_vectors_rDeriv, rFiberDirectionsReference);

    //     // ---------------------------------------------------------------------------------------
    //     // 11. first derivative of Cauchy stress tensor w.r.t DOFs
    //     // ---------------------------------------------------------------------------------------
    //     Matrix sigma_rDeriv = ZeroMatrix(2);


    //     double C_reference_rDeriv[2][2][2][2] = { {{{0.0}}} };
    //     double c_current_rDeriv[2][2][2][2] = { {{{0.0}}} };

    //     for (SizeType alpha = 0; alpha < 2; ++alpha) {
    //       for (SizeType beta = 0; beta < 2; ++beta) {
    //         for (SizeType gamma = 0; gamma < 2; ++gamma) {
    //           for (SizeType delta = 0; delta < 2; ++delta) {

    //             //if (is_isotropic) //************************ check whether 'is_isotropic works' *************************
    //             //{
    //               C_reference_rDeriv[alpha][beta][gamma][delta] =
    //                 lambda * (ref_contra_metric_rDeriv(alpha, beta) * contravariant_metric_reference(gamma, delta) +
    //                   contravariant_metric_reference(alpha, beta) * ref_contra_metric_rDeriv(gamma, delta))
    //                 + MU * (ref_contra_metric_rDeriv(alpha, gamma) * contravariant_metric_reference(beta, delta) +
    //                   ref_contra_metric_rDeriv(alpha, delta) * contravariant_metric_reference(beta, gamma) +
    //                   contravariant_metric_reference(alpha, gamma) * ref_contra_metric_rDeriv(beta, delta) +
    //                   contravariant_metric_reference(alpha, delta) * ref_contra_metric_rDeriv(beta, gamma));
    //             /*}*/
    //             //else if (is_orthotropic) //************************ check whether 'is_orthotropic works' *************************
    //             //{
    //             //  for (SizeType eps = 0; eps < 2; ++eps) {
    //             //    for (SizeType zet = 0; zet < 2; ++zet) {
    //             //      for (SizeType eta = 0; eta < 2; ++eta) {
    //             //        for (SizeType the = 0; the < 2; ++the) {
    //             //          C_reference_rDeriv[alpha][beta][gamma][delta] +=
    //             //            coeff_C[eps][zet][eta][the] *
    //             //            (
    //             //              transformation_matrix_orthotropic_rDeriv(alpha, eps) * transformation_matrix_orthotropic(beta, zet) * transformation_matrix_orthotropic(gamma, eta) * transformation_matrix_orthotropic(delta, the) +
    //             //              transformation_matrix_orthotropic(alpha, eps) * transformation_matrix_orthotropic_rDeriv(beta, zet) * transformation_matrix_orthotropic(gamma, eta) * transformation_matrix_orthotropic(delta, the) +
    //             //              transformation_matrix_orthotropic(alpha, eps) * transformation_matrix_orthotropic(beta, zet) * transformation_matrix_orthotropic_rDeriv(gamma, eta) * transformation_matrix_orthotropic(delta, the) +
    //             //              transformation_matrix_orthotropic(alpha, eps) * transformation_matrix_orthotropic(beta, zet) * transformation_matrix_orthotropic(gamma, eta) * transformation_matrix_orthotropic_rDeriv(delta, the)
    //             //              );
    //             //        }
    //             //      }
    //             //    }
    //             //  }
    //             //}

    //             c_current_rDeriv[alpha][beta][gamma][delta] =
    //               InverseDetF_rDeriv * C_reference[alpha][beta][gamma][delta] +
    //               InverseDetF * C_reference_rDeriv[alpha][beta][gamma][delta];

    //             sigma_rDeriv(alpha, beta) +=
    //               c_current_rDeriv[alpha][beta][gamma][delta] * e(gamma, delta) +
    //               c_current[alpha][beta][gamma][delta] * e_rDeriv(gamma, delta);
    //           }
    //         }
    //       }
    //     }
    //     // =================================================================
    //     // 19. Loop over design DOFs X_r
    //     // =================================================================
    //     for (SizeType s = 0; s < number_dofs; ++s)
    //     {

    //       // -------------------------------------------------------------------------------------
    //       // 1. first derivative of covariant base vectors in reference configuration w.r.t DOFs
    //       // -------------------------------------------------------------------------------------
    //       array_1d<Vector, 2> ref_cov_base_vectors_sDeriv;
    //       this->DeriveCurrentCovariantBaseVectors(ref_cov_base_vectors_sDeriv, shape_functions_gradients_i, s);


    //       // -------------------------------------------------------------------------------------
    //       // 2. first derivative of covariant metric in reference configuration w.r.t DOFs
    //       // -------------------------------------------------------------------------------------
    //       Matrix ref_cov_metric_sDeriv = ZeroMatrix(2);
    //       this->DerivativeCurrentCovariantMetric(ref_cov_metric_sDeriv, shape_functions_gradients_i, s, reference_covariant_base_vectors);

    //       //second derivative of covariant metric in reference configuration w.r.t DOFs       
    //       Matrix ref_cov_metric_rsDeriv = ZeroMatrix(2);
    //       this->Derivative2CurrentCovariantMetric(ref_cov_metric_rsDeriv, shape_functions_gradients_i, r, s);


    //       // -------------------------------------------------------------------------------------
    //       // 3. first derivative of Euler Almansi strain tensor w.r.t DOFs
    //       // -------------------------------------------------------------------------------------
    //       Matrix e_sDeriv = ZeroMatrix(2);
    //       e_sDeriv(0, 0) = -0.5 * ref_cov_metric_sDeriv(0, 0);
    //       e_sDeriv(0, 1) = -0.5 * ref_cov_metric_sDeriv(0, 1);
    //       e_sDeriv(1, 0) = -0.5 * ref_cov_metric_sDeriv(1, 0);
    //       e_sDeriv(1, 1) = -0.5 * ref_cov_metric_sDeriv(1, 1);

    //       // second derivative of Euler Almansi strain tensor w.r.t DOFs
    //       Matrix e_rsDeriv = ZeroMatrix(2);
    //       e_rsDeriv(0, 0) = -0.5 * ref_cov_metric_rsDeriv(0, 0);
    //       e_rsDeriv(0, 1) = -0.5 * ref_cov_metric_rsDeriv(0, 1);
    //       e_rsDeriv(1, 0) = -0.5 * ref_cov_metric_rsDeriv(1, 0);
    //       e_rsDeriv(1, 1) = -0.5 * ref_cov_metric_rsDeriv(1, 1);


    //       // -------------------------------------------------------------------------------------
    //       // 4. first derivative of determinant of reference covariant metric w.r.t DOFs
    //       // -------------------------------------------------------------------------------------
    //       double sDeriv_det_cov_metric_ref = ref_cov_metric_sDeriv(0, 0) * covariant_metric_reference(1, 1) + covariant_metric_reference(0, 0) * ref_cov_metric_sDeriv(1, 1)
    //         - ref_cov_metric_sDeriv(1, 0) * covariant_metric_reference(0, 1) - covariant_metric_reference(1, 0) * ref_cov_metric_sDeriv(0, 1);

    //       // second derivative of determinant of reference covariant metric w.r.t DOFs
    //       double rsDeriv_det_cov_metric_ref = ref_cov_metric_rsDeriv(0, 0) * covariant_metric_reference(1, 1) + ref_cov_metric_rDeriv(0, 0) * ref_cov_metric_sDeriv(1, 1)
    //         + ref_cov_metric_sDeriv(0, 0) * ref_cov_metric_rDeriv(1, 1) + covariant_metric_reference(0, 0) * ref_cov_metric_rsDeriv(1, 1)
    //         - ref_cov_metric_rsDeriv(0, 1) * covariant_metric_reference(1, 0) - ref_cov_metric_rDeriv(0, 1) * ref_cov_metric_sDeriv(1, 0)
    //         - ref_cov_metric_sDeriv(0, 1) * ref_cov_metric_rDeriv(1, 0) - covariant_metric_reference(0, 1) * ref_cov_metric_rsDeriv(1, 0);


    //       // -------------------------------------------------------------------------------------
    //       // 5. first derivative of determinant of reference contravariant metric w.r.t DOFs
    //       // -------------------------------------------------------------------------------------
    //       double sDeriv_det_contra_metric_ref = -sDeriv_det_cov_metric_ref / (det_covariant_metric_reference * det_covariant_metric_reference);

    //       // second derivative of determinant of reference contravariant metric w.r.t DOFs
    //       double rsDeriv_det_contra_metric_ref = (-rsDeriv_det_cov_metric_ref * det_covariant_metric_reference
    //         + 2 * rDeriv_det_cov_metric_ref * sDeriv_det_cov_metric_ref) / (det_covariant_metric_reference * det_covariant_metric_reference * det_covariant_metric_reference);


    //       // -------------------------------------------------------------------------------------
    //       // 6. first derivative of contravariant metric in reference configuration w.r.t DOFs
    //       // -------------------------------------------------------------------------------------
    //       Matrix ref_contra_metric_sDeriv = ZeroMatrix(2);

    //       ref_contra_metric_sDeriv(0, 0) = sDeriv_det_contra_metric_ref * covariant_metric_reference(1, 1) + det_contravariant_metric_reference * ref_cov_metric_sDeriv(1, 1);
    //       ref_contra_metric_sDeriv(0, 1) = -sDeriv_det_contra_metric_ref * covariant_metric_reference(0, 1) - det_contravariant_metric_reference * ref_cov_metric_sDeriv(0, 1);
    //       ref_contra_metric_sDeriv(1, 0) = -sDeriv_det_contra_metric_ref * covariant_metric_reference(1, 0) - det_contravariant_metric_reference * ref_cov_metric_sDeriv(1, 0);
    //       ref_contra_metric_sDeriv(1, 1) = sDeriv_det_contra_metric_ref * covariant_metric_reference(0, 0) + det_contravariant_metric_reference * ref_cov_metric_sDeriv(0, 0);

    //       // second derivative of contravariant metric in reference configuration w.r.t DOFs
    //       Matrix ref_contra_metric_rsDeriv = ZeroMatrix(2);

    //       ref_contra_metric_rsDeriv(0, 0) = rsDeriv_det_contra_metric_ref * covariant_metric_reference(1, 1)
    //         + rDeriv_det_contra_metric_ref * ref_cov_metric_sDeriv(1, 1)
    //         + sDeriv_det_contra_metric_ref * ref_cov_metric_rDeriv(1, 1)
    //         + det_contravariant_metric_reference * ref_cov_metric_rsDeriv(1, 1);

    //       ref_contra_metric_rsDeriv(0, 1) = -rsDeriv_det_contra_metric_ref * covariant_metric_reference(0, 1)
    //         - rDeriv_det_contra_metric_ref * ref_cov_metric_sDeriv(0, 1)
    //         - sDeriv_det_contra_metric_ref * ref_cov_metric_rDeriv(0, 1)
    //         - det_contravariant_metric_reference * ref_cov_metric_rsDeriv(0, 1);

    //       ref_contra_metric_rsDeriv(1, 0) = -rsDeriv_det_contra_metric_ref * covariant_metric_reference(1, 0)
    //         - rDeriv_det_contra_metric_ref * ref_cov_metric_sDeriv(1, 0)
    //         - sDeriv_det_contra_metric_ref * ref_cov_metric_rDeriv(1, 0)
    //         - det_contravariant_metric_reference * ref_cov_metric_rsDeriv(1, 0);

    //       ref_contra_metric_rsDeriv(1, 1) = rsDeriv_det_contra_metric_ref * covariant_metric_reference(0, 0)
    //         + rDeriv_det_contra_metric_ref * ref_cov_metric_sDeriv(0, 0)
    //         + sDeriv_det_contra_metric_ref * ref_cov_metric_sDeriv(0, 0)
    //         + det_contravariant_metric_reference * ref_cov_metric_rsDeriv(0, 0);


    //       // ----------------------------------------------------------------------------------------
    //       // 7. first derivative of contravariant base vectors in reference configuration w.r.t DOFs
    //       // ----------------------------------------------------------------------------------------
    //       array_1d<Vector, 2> ref_contra_base_vectors_sDeriv;

    //       for (SizeType alpha = 0; alpha < 2; ++alpha) {
    //         ref_contra_base_vectors_sDeriv[alpha] = ZeroVector(3);

    //         for (SizeType beta = 0; beta < 2; ++beta) {
    //           ref_contra_base_vectors_sDeriv[alpha] +=
    //             ref_contra_metric_sDeriv(alpha, beta) * reference_covariant_base_vectors[beta]
    //             + contravariant_metric_reference(alpha, beta) * ref_cov_base_vectors_sDeriv[beta];
    //         }
    //       }

    //       // second derivative of contravariant base vectors in reference configuration w.r.t DOFs
    //       array_1d<Vector, 2> ref_contra_base_vectors_rsDeriv;


    //       // ---------------------------------------------------------------------------------------
    //       // 8. first derivative of elemental area in reference configuration w.r.t DOFs
    //       // ---------------------------------------------------------------------------------------
    //       double dA_reference_sDeriv = 0.0;
    //       this->DerivativeElementalArea(dA_reference_sDeriv, reference_covariant_base_vectors, ref_cov_base_vectors_sDeriv);

    //       // second derivative of elemental area in reference configuration w.r.t DOFs
    //       double dA_reference_rsDeriv = 0.0;
    //       this->Derivative2ElementalArea(dA_reference_rsDeriv, reference_covariant_base_vectors, ref_cov_base_vectors_rDeriv, ref_cov_base_vectors_sDeriv);


    //       // ---------------------------------------------------------------------------------------
    //       // 9. first derivative of inverse of determinant of deformation gradient w.r.t DOFs
    //       // ---------------------------------------------------------------------------------------
    //       double InverseDetF_sDeriv = dA_reference_sDeriv / da_current;

    //       // second derivative of inverse of determinant of deformation gradient w.r.t DOFs
    //       double InverseDetF_rsDeriv = dA_reference_rsDeriv / da_current;


    //       // ---------------------------------------------------------------------------------------
    //       // 10. first derivative of transformation matrix for orthotropic material w.r.t DOFs
    //       // ---------------------------------------------------------------------------------------
    //       //Matrix transformation_matrix_orthotropic_sDeriv = ZeroMatrix(2);
    //       //this->comp_Transformation_Matrix(transformation_matrix_orthotropic_sDeriv, ref_contra_base_vectors_sDeriv, rFiberDirectionsReference);

    //       // second derivative of transformation matrix for orthotropic material w.r.t DOFs
    //       // ************************ Check for accuracy **************************
    //       //Matrix transformation_matrix_orthotropic_rsDeriv = ZeroMatrix(2);
    //       //this->comp_Transformation_Matrix(transformation_matrix_orthotropic_rsDeriv, ref_contra_base_vectors_rsDeriv, rFiberDirectionsReference);


    //       // ---------------------------------------------------------------------------------------
    //       // 11. first derivative of Cauchy stress tensor w.r.t DOFs
    //       // ---------------------------------------------------------------------------------------
    //       Matrix sigma_sDeriv = ZeroMatrix(2);
    //       Matrix sigma_rsDeriv = ZeroMatrix(2);

    //       double C_reference_sDeriv[2][2][2][2] = { {{{0.0}}} };
    //       double c_current_sDeriv[2][2][2][2] = { {{{0.0}}} };

    //       double C_reference_rsDeriv[2][2][2][2] = { {{{0.0}}} };
    //       double c_current_rsDeriv[2][2][2][2] = { {{{0.0}}} };

    //       for (SizeType alpha = 0; alpha < 2; ++alpha) {
    //         for (SizeType beta = 0; beta < 2; ++beta) {
    //           for (SizeType gamma = 0; gamma < 2; ++gamma) {
    //             for (SizeType delta = 0; delta < 2; ++delta) {

    //               //if (is_isotropic) //************************ check whether 'is_isotropic works' *************************
    //               //{
    //                 C_reference_sDeriv[alpha][beta][gamma][delta] =
    //                   lambda * (ref_contra_metric_sDeriv(alpha, beta) * contravariant_metric_reference(gamma, delta) +
    //                     contravariant_metric_reference(alpha, beta) * ref_contra_metric_sDeriv(gamma, delta))
    //                   + MU * (ref_contra_metric_sDeriv(alpha, gamma) * contravariant_metric_reference(beta, delta) +
    //                     ref_contra_metric_sDeriv(alpha, delta) * contravariant_metric_reference(beta, gamma) +
    //                     contravariant_metric_reference(alpha, gamma) * ref_contra_metric_sDeriv(beta, delta) +
    //                     contravariant_metric_reference(alpha, delta) * ref_contra_metric_sDeriv(beta, gamma));

    //                 C_reference_rsDeriv[alpha][beta][gamma][delta] =
    //                   lambda * (ref_contra_metric_rsDeriv(alpha, beta) * contravariant_metric_reference(gamma, delta) +
    //                     ref_contra_metric_rDeriv(alpha, beta) * ref_contra_metric_sDeriv(gamma, delta) +
    //                     ref_contra_metric_sDeriv(alpha, beta) * ref_contra_metric_rDeriv(gamma, delta) +
    //                     contravariant_metric_reference(alpha, beta) * ref_contra_metric_rsDeriv(gamma, delta)) +
    //                   MU * (ref_contra_metric_rsDeriv(alpha, gamma) * contravariant_metric_reference(beta, delta) +
    //                     ref_contra_metric_rDeriv(alpha, gamma) * ref_contra_metric_sDeriv(beta, delta) +
    //                     ref_contra_metric_sDeriv(alpha, gamma) * ref_contra_metric_rDeriv(beta, delta) +
    //                     contravariant_metric_reference(alpha, gamma) * ref_contra_metric_rsDeriv(beta, delta) +
    //                     ref_contra_metric_rsDeriv(alpha, delta) * contravariant_metric_reference(beta, gamma) +
    //                     ref_contra_metric_rDeriv(alpha, delta) * ref_contra_metric_sDeriv(beta, gamma) +
    //                     ref_contra_metric_sDeriv(alpha, delta) * ref_contra_metric_rDeriv(beta, gamma) +
    //                     contravariant_metric_reference(alpha, delta) * ref_contra_metric_rsDeriv(beta, gamma));
    //               /*}*/
    //               //else if (is_orthotropic) //************************ check whether 'is_orthotropic works' *************************
    //               //{
    //               //  for (SizeType eps = 0; eps < 2; ++eps) {
    //               //    for (SizeType zet = 0; zet < 2; ++zet) {
    //               //      for (SizeType eta = 0; eta < 2; ++eta) {
    //               //        for (SizeType the = 0; the < 2; ++the) {
    //               //          C_reference_sDeriv[alpha][beta][gamma][delta] +=
    //               //            coeff_C[eps][zet][eta][the] *
    //               //            (
    //               //              transformation_matrix_orthotropic_sDeriv(alpha, eps) * transformation_matrix_orthotropic(beta, zet) * transformation_matrix_orthotropic(gamma, eta) * transformation_matrix_orthotropic(delta, the) +
    //               //              transformation_matrix_orthotropic(alpha, eps) * transformation_matrix_orthotropic_sDeriv(beta, zet) * transformation_matrix_orthotropic(gamma, eta) * transformation_matrix_orthotropic(delta, the) +
    //               //              transformation_matrix_orthotropic(alpha, eps) * transformation_matrix_orthotropic(beta, zet) * transformation_matrix_orthotropic_sDeriv(gamma, eta) * transformation_matrix_orthotropic(delta, the) +
    //               //              transformation_matrix_orthotropic(alpha, eps) * transformation_matrix_orthotropic(beta, zet) * transformation_matrix_orthotropic(gamma, eta) * transformation_matrix_orthotropic_sDeriv(delta, the)
    //               //              );

    //               //          C_reference_rsDeriv[alpha][beta][gamma][delta] +=
    //               //            coeff_C[eps][zet][eta][the] *
    //               //            (
          
    //               //              transformation_matrix_orthotropic_rsDeriv(alpha, eps) * transformation_matrix_orthotropic(beta, zet) * transformation_matrix_orthotropic(gamma, eta) * transformation_matrix_orthotropic(delta, the) + 
    //               //              transformation_matrix_orthotropic(alpha, eps) * transformation_matrix_orthotropic_rsDeriv(beta, zet) * transformation_matrix_orthotropic(gamma, eta) * transformation_matrix_orthotropic(delta, the) + 
    //               //              transformation_matrix_orthotropic(alpha, eps) * transformation_matrix_orthotropic(beta, zet) * transformation_matrix_orthotropic_rsDeriv(gamma, eta) * transformation_matrix_orthotropic(delta, the) + 
    //               //              transformation_matrix_orthotropic(alpha, eps) * transformation_matrix_orthotropic(beta, zet) * transformation_matrix_orthotropic(gamma, eta) * transformation_matrix_orthotropic_rsDeriv(delta, the) + 
    //               //              transformation_matrix_orthotropic_rDeriv(alpha, eps) * transformation_matrix_orthotropic_sDeriv(beta, zet) * transformation_matrix_orthotropic(gamma, eta) * transformation_matrix_orthotropic(delta, the) + 
    //               //              transformation_matrix_orthotropic_sDeriv(alpha, eps) * transformation_matrix_orthotropic_rDeriv(beta, zet) * transformation_matrix_orthotropic(gamma, eta) * transformation_matrix_orthotropic(delta, the) + 
    //               //              transformation_matrix_orthotropic_rDeriv(alpha, eps) * transformation_matrix_orthotropic(beta, zet) * transformation_matrix_orthotropic_sDeriv(gamma, eta) * transformation_matrix_orthotropic(delta, the) + 
    //               //              transformation_matrix_orthotropic_sDeriv(alpha, eps) * transformation_matrix_orthotropic(beta, zet) * transformation_matrix_orthotropic_rDeriv(gamma, eta) * transformation_matrix_orthotropic(delta, the) + 
    //               //              transformation_matrix_orthotropic_rDeriv(alpha, eps) * transformation_matrix_orthotropic(beta, zet) * transformation_matrix_orthotropic(gamma, eta) * transformation_matrix_orthotropic_sDeriv(delta, the) + 
    //               //              transformation_matrix_orthotropic_sDeriv(alpha, eps) * transformation_matrix_orthotropic(beta, zet) * transformation_matrix_orthotropic(gamma, eta) * transformation_matrix_orthotropic_rDeriv(delta, the) + 
    //               //              transformation_matrix_orthotropic(alpha, eps) * transformation_matrix_orthotropic_rDeriv(beta, zet) * transformation_matrix_orthotropic_sDeriv(gamma, eta) * transformation_matrix_orthotropic(delta, the) + 
    //               //              transformation_matrix_orthotropic(alpha, eps) * transformation_matrix_orthotropic_sDeriv(beta, zet) * transformation_matrix_orthotropic_rDeriv(gamma, eta) * transformation_matrix_orthotropic(delta, the) + 
    //               //              transformation_matrix_orthotropic(alpha, eps) * transformation_matrix_orthotropic_rDeriv(beta, zet) * transformation_matrix_orthotropic(gamma, eta) * transformation_matrix_orthotropic_sDeriv(delta, the) + 
    //               //              transformation_matrix_orthotropic(alpha, eps) * transformation_matrix_orthotropic_sDeriv(beta, zet) * transformation_matrix_orthotropic(gamma, eta) * transformation_matrix_orthotropic_rDeriv(delta, the) + 
    //               //              transformation_matrix_orthotropic(alpha, eps) * transformation_matrix_orthotropic(beta, zet) * transformation_matrix_orthotropic_rDeriv(gamma, eta) * transformation_matrix_orthotropic_sDeriv(delta, the) + 
    //               //              transformation_matrix_orthotropic(alpha, eps) * transformation_matrix_orthotropic(beta, zet) * transformation_matrix_orthotropic_sDeriv(gamma, eta) * transformation_matrix_orthotropic_rDeriv(delta, the)
    //               //              );
    //               //         

    //               //        }
    //               //      }
    //               //    }
    //               //  }
    //               //}

    //               c_current_sDeriv[alpha][beta][gamma][delta] = InverseDetF_sDeriv * C_reference[alpha][beta][gamma][delta] + InverseDetF * C_reference_sDeriv[alpha][beta][gamma][delta];

    //               c_current_rsDeriv[alpha][beta][gamma][delta] = InverseDetF_rsDeriv * C_reference[alpha][beta][gamma][delta] +
    //                 InverseDetF_rDeriv * C_reference_sDeriv[alpha][beta][gamma][delta] +
    //                 InverseDetF_sDeriv * C_reference_rDeriv[alpha][beta][gamma][delta] +
    //                 InverseDetF * C_reference_rsDeriv[alpha][beta][gamma][delta];


    //               sigma_sDeriv(alpha, beta) +=
    //                 c_current_sDeriv[alpha][beta][gamma][delta] * e(gamma, delta) +
    //                 c_current[alpha][beta][gamma][delta] * e_sDeriv(gamma, delta);

    //               sigma_rsDeriv(alpha, beta) += c_current_rsDeriv[alpha][beta][gamma][delta] * e(gamma, delta) +
    //                 c_current_rDeriv[alpha][beta][gamma][delta] * e_sDeriv(gamma, delta) +
    //                 c_current_sDeriv[alpha][beta][gamma][delta] * e_rDeriv(gamma, delta) +
    //                 c_current[alpha][beta][gamma][delta] * e_rsDeriv(gamma, delta);
    //             }
    //           }
    //         }
    //       }


    //       // =================================================================
    //       // elemental stiffness matrix (Dieringer - equation 5.15)
    //       // =================================================================
    //       double k_elm = 0.0;
    //       for (SizeType m = 0; m < 2; m++){
    //         for (SizeType n = 0; n < 2; n++) {
    //           for (SizeType l = 0; l < 2; l++) {
    //             for (SizeType k = 0; k < 2; k++) {
    //               k_elm += ((sigma_sDeriv(m, n) * sigma_rDeriv(l, k)) + ((sigma(m, n) - sigma_pre_current_contra(m, n)) * sigma_rsDeriv(l, k))) * covariant_metric_current(m, l) * covariant_metric_current(n, k);
    //             }
    //           }
    //         }
    //       }

    //       rLeftHandSideMatrix(r, s) += (k_elm * integration_weight_i * da_current);

    //     }

    //   }


    //   // =================================================================
    //   // calculation of response function (Dieringer - equation 5.12)
    //   // =================================================================

    //   /*Matrix delta_sigma = ZeroMatrix(2);
    //   delta_sigma = sigma - sigma_pre_current_contra;

    //   double delta_sigma_2 = 0.0;

    //   for (SizeType m = 0; m < 2; m++) {
    //     for (SizeType n = 0; n < 2; n++) {
    //       for (SizeType l = 0; l < 2; l++) {
    //         for (SizeType k = 0; k < 2; k++) {
    //           delta_sigma_2 += delta_sigma(m, n) * delta_sigma(l, k) * covariant_metric_current(m, l) * covariant_metric_current(n, k);
    //         }
    //       }
    //     }
    //   }

    //   delta_sigma_2 = delta_sigma_2 * integration_weight_i * da_current;

    //   rResponse += 0.5 * delta_sigma_2;*/
    // }

    // KRATOS_CATCH("");
  }



  
  void MembraneCuttingPatternElement::comp_Transformation_Matrix(Matrix& rTransMat, const array_1d<Vector, 2>& rBaseVector1, const array_1d<Vector, 2>& rBaseVector2)
  {
    // rTransMat = ZeroMatrix(2);
    // rTransMat(0, 0) = inner_prod(rBaseVector1[0], rBaseVector2[0]);
    // rTransMat(0, 1) = inner_prod(rBaseVector1[0], rBaseVector2[1]);
    // rTransMat(1, 0) = inner_prod(rBaseVector1[1], rBaseVector2[0]);
    // rTransMat(1, 1) = inner_prod(rBaseVector1[1], rBaseVector2[1]);
  }



  void MembraneCuttingPatternElement::DerivativeElementalArea(double& rDerivElementalArea, const array_1d<Vector, 2>& rCovariantBaseVectors, const array_1d<Vector, 2>& rCovariantBaseVectorsDerivative)
  {
    // Vector cross_G1G2 = ZeroVector(3);
    // MathUtils<double>::CrossProduct(cross_G1G2, rCovariantBaseVectors[0], rCovariantBaseVectors[1]);

    // const double dA = MathUtils<double>::Norm(cross_G1G2);

    // Vector cross_dG1_G2 = ZeroVector(3);
    // MathUtils<double>::CrossProduct(cross_dG1_G2, rCovariantBaseVectorsDerivative[0], rCovariantBaseVectors[1]);

    // Vector cross_G1_dG2 = ZeroVector(3);
    // MathUtils<double>::CrossProduct(cross_G1_dG2, rCovariantBaseVectors[0], rCovariantBaseVectorsDerivative[1]);

    // Vector cross_derivative = cross_dG1_G2 + cross_G1_dG2;

    // rDerivElementalArea = inner_prod(cross_derivative, cross_G1G2) / dA;
  }



  void MembraneCuttingPatternElement::Derivative2ElementalArea(double& rDeriv2ElementalArea, const array_1d<Vector, 2>& rCovariantBaseVectors, const array_1d<Vector, 2>& rCovariantBaseVectorsDerivativeR, const array_1d<Vector, 2>& rCovariantBaseVectorsDerivativeS)
  {
    // Vector cross_G1G2 = ZeroVector(3);
    // MathUtils<double>::CrossProduct(cross_G1G2, rCovariantBaseVectors[0], rCovariantBaseVectors[1]);

    // const double dA = MathUtils<double>::Norm(cross_G1G2);

    // Vector cross_rG1_G2 = ZeroVector(3);
    // MathUtils<double>::CrossProduct(cross_rG1_G2, rCovariantBaseVectorsDerivativeR[0], rCovariantBaseVectors[1]);

    // Vector cross_G1_rG2 = ZeroVector(3);
    // MathUtils<double>::CrossProduct(cross_G1_rG2, rCovariantBaseVectors[0], rCovariantBaseVectorsDerivativeR[1]);

    // const Vector cross_rDeriv = cross_rG1_G2 + cross_G1_rG2;

    // Vector cross_sG1_G2 = ZeroVector(3);
    // MathUtils<double>::CrossProduct(cross_sG1_G2, rCovariantBaseVectorsDerivativeS[0], rCovariantBaseVectors[1]);

    // Vector cross_G1_sG2 = ZeroVector(3);
    // MathUtils<double>::CrossProduct(cross_G1_sG2, rCovariantBaseVectors[0], rCovariantBaseVectorsDerivativeS[1]);

    // const Vector cross_sDeriv = cross_sG1_G2 + cross_G1_sG2;

    // Vector cross_rG1_sG2 = ZeroVector(3);
    // MathUtils<double>::CrossProduct(cross_rG1_sG2, rCovariantBaseVectorsDerivativeR[0], rCovariantBaseVectorsDerivativeS[1]);

    // Vector cross_sG1_rG2 = ZeroVector(3);
    // MathUtils<double>::CrossProduct(cross_sG1_rG2, rCovariantBaseVectorsDerivativeS[0], rCovariantBaseVectorsDerivativeR[1]);

    // const Vector cross_rsDeriv = cross_rG1_sG2 + cross_sG1_rG2;

    // const double cross_r_dot_cross = inner_prod(cross_rDeriv, cross_G1G2);
    // const double cross_s_dot_cross = inner_prod(cross_sDeriv, cross_G1G2);

    // rDeriv2ElementalArea =
    //   (inner_prod(cross_rsDeriv, cross_G1G2) + inner_prod(cross_rDeriv, cross_sDeriv)) / dA
    //   - (cross_r_dot_cross * cross_s_dot_cross) / (dA * dA * dA);
  }

}