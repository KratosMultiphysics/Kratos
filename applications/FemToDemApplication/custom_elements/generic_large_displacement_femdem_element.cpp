//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

// System includes

// External includes

// Project includes
#include "generic_large_displacement_femdem_element.hpp"
#include "utilities/geometry_utilities.h"
#include "fem_to_dem_application_variables.h"

namespace Kratos
{


//***********************DEFAULT CONSTRUCTOR******************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
GenericLargeDisplacementFemDemElement<TDim, TyieldSurf>::GenericLargeDisplacementFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry)
	: GenericSmallStrainFemDemElement(NewId, pGeometry)
{
	//DO NOT ADD DOFS HERE!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
GenericLargeDisplacementFemDemElement<TDim, TyieldSurf>::GenericLargeDisplacementFemDemElement(
    IndexType NewId, 
    GeometryType::Pointer pGeometry, 
    PropertiesType::Pointer pProperties
    )
	: GenericSmallStrainFemDemElement(NewId, pGeometry, pProperties)
{
	// BY DEFAULT, THE GEOMETRY WILL DEFINE THE INTEGRATION METHOD
	mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();

	// Each component == Each edge
	if (mNonConvergedThresholds.size() != NumberOfEdges)
    	mNonConvergedThresholds.resize(NumberOfEdges);
	noalias(mNonConvergedThresholds) = ZeroVector(NumberOfEdges);   // Equivalent stress

	if (mThresholds.size() != NumberOfEdges)
    	mThresholds.resize(NumberOfEdges);
	noalias(mThresholds) = ZeroVector(NumberOfEdges); // Stress mThreshold on edge

	if (mNonConvergedDamages.size() != NumberOfEdges)
    	mNonConvergedDamages.resize(NumberOfEdges);
	noalias(mDamages) = ZeroVector(NumberOfEdges); // Converged mDamage on each edge

	if (mNonConvergedDamages.size() != NumberOfEdges)
    	mNonConvergedDamages.resize(NumberOfEdges);
	noalias(mNonConvergedDamages) = ZeroVector(NumberOfEdges); // mDamages on edges of "i" iteration
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

template<unsigned int TDim, unsigned int TyieldSurf>
GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::GenericLargeDisplacementFemDemElement(GenericLargeDisplacementFemDemElement const& rOther)
	: GenericSmallStrainFemDemElement(rOther)
{
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

template<unsigned int TDim, unsigned int TyieldSurf>
GenericLargeDisplacementFemDemElement<TDim, TyieldSurf> &GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::operator=(GenericLargeDisplacementFemDemElement const& rOther)
{
	GenericSmallStrainFemDemElement::operator=(rOther);
	return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::Create(
    IndexType NewId, 
    NodesArrayType const& rThisNodes, 
    PropertiesType::Pointer pProperties
    ) const
{
	return Element::Pointer(new GenericLargeDisplacementFemDemElement(NewId, GetGeometry().Create(rThisNodes), pProperties));
}

//************************************CLONE*******************************************
//************************************************************************************

template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
{
	GenericLargeDisplacementFemDemElement NewElement(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
	return Element::Pointer(new GenericLargeDisplacementFemDemElement(NewElement));
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

template<unsigned int TDim, unsigned int TyieldSurf>
GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::~GenericLargeDisplacementFemDemElement()
{
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::InitializeNonLinearIteration(
    ProcessInfo& rCurrentProcessInfo
    )
{

}


/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateLocalSystem(
	MatrixType& rLeftHandSideMatrix,
	VectorType& rRightHandSideVector,
	ProcessInfo& rCurrentProcessInfo
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, 
    ProcessInfo& rCurrentProcessInfo
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    const SizeType number_of_nodes = this->GetGeometry().size();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();

    // Kinematic variables
    Matrix B, F, DN_DX, InvJ0, J, J0;
    double detJ0;

    const SizeType mat_size = number_of_nodes * dimension;
    B.resize(VoigtSize, dimension * number_of_nodes);

    // Resizing as needed the RHS
    if (rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size, false);
    rRightHandSideVector = ZeroVector(mat_size); //resetting RHS

    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    Matrix delta_position(number_of_nodes, dimension);
    noalias(delta_position) = ZeroMatrix(number_of_nodes, dimension);
    delta_position = this->CalculateDeltaPosition(delta_position);

    // Loop over Gauss Points
    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {

		J = this->GetGeometry().Jacobian(J, point_number, mThisIntegrationMethod);
        detJ0 = this->CalculateDerivativesOnReferenceConfiguration(J0, InvJ0, DN_DX, point_number, mThisIntegrationMethod);

        const double integration_weigth = integration_points[point_number].Weight() * detJ0;
        const Matrix &Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);
		Vector N = row(Ncontainer, point_number);

        Vector VolumeForce = ZeroVector(dimension);
		VolumeForce = this->CalculateVolumeForce(VolumeForce, N);
		// Taking into account Volume Force into de RHS
		for (unsigned int i = 0; i < number_of_nodes; i++) {
			int index = dimension * i;
			for (unsigned int j = 0; j < dimension; j++) {
				rRightHandSideVector[index + j] += integration_weigth * N[i] * VolumeForce[j];
			}
		}

        GeometryUtils::DeformationGradient(J, InvJ0, F);
        this->CalculateB(B, F, DN_DX);

        const double damage_element = this->CalculateElementalDamage(mNonConvergedDamages);
		const Vector& r_stress_vector = this->GetValue(STRESS_VECTOR);
        const Vector& r_integrated_stress_vector = (1.0 - damage_element) * r_stress_vector;
        this->CalculateAndAddInternalForcesVector(rRightHandSideVector, B, r_integrated_stress_vector, integration_weigth);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void GenericLargeDisplacementFemDemElement<2,0>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB2D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<3,0>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB3D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<2,1>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB2D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<3,1>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB3D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<2,2>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB2D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<3,2>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB3D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<2,3>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB2D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<3,3>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB3D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<2,4>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB2D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<3,4>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB3D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<2,5>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB2D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<3,5>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB3D(rB, rF, rDN_DX);
}


template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateB2D(
    Matrix& rB, 
    const Matrix& rF, 
    const Matrix& rDN_DX
    )
{
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
	const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    for (SizeType i = 0; i < number_of_nodes; i++) {
        unsigned int index = 2 * i;
        rB( 0, index + 0 ) = rF( 0, 0 ) * rDN_DX( i, 0 );
        rB( 0, index + 1 ) = rF( 1, 0 ) * rDN_DX( i, 0 );
        rB( 1, index + 0 ) = rF( 0, 1 ) * rDN_DX( i, 1 );
        rB( 1, index + 1 ) = rF( 1, 1 ) * rDN_DX( i, 1 );
        rB( 2, index + 0 ) = rF( 0, 0 ) * rDN_DX( i, 1 ) + rF( 0, 1 ) * rDN_DX( i, 0 );
        rB( 2, index + 1 ) = rF( 1, 0 ) * rDN_DX( i, 1 ) + rF( 1, 1 ) * rDN_DX( i, 0 );

    }
}
template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateB3D(
    Matrix& rB, 
    const Matrix& rF, 
    const Matrix& rDN_DX
    )
{
	const unsigned int number_of_nodes = GetGeometry().PointsNumber();
	const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    for (SizeType i = 0; i < number_of_nodes; ++i) {
        const auto index = dimension * i;
        rB(0, index + 0) = rF(0, 0) * rDN_DX(i, 0);
        rB(0, index + 1) = rF(1, 0) * rDN_DX(i, 0);
        rB(0, index + 2) = rF(2, 0) * rDN_DX(i, 0);
        rB(1, index + 0) = rF(0, 1) * rDN_DX(i, 1);
        rB(1, index + 1) = rF(1, 1) * rDN_DX(i, 1);
        rB(1, index + 2) = rF(2, 1) * rDN_DX(i, 1);
        rB(2, index + 0) = rF(0, 2) * rDN_DX(i, 2);
        rB(2, index + 1) = rF(1, 2) * rDN_DX(i, 2);
        rB(2, index + 2) = rF(2, 2) * rDN_DX(i, 2);
        rB(3, index + 0) = rF(0, 0) * rDN_DX(i, 1) + rF(0, 1) * rDN_DX(i, 0);
        rB(3, index + 1) = rF(1, 0) * rDN_DX(i, 1) + rF(1, 1) * rDN_DX(i, 0);
        rB(3, index + 2) = rF(2, 0) * rDN_DX(i, 1) + rF(2, 1) * rDN_DX(i, 0);
        rB(4, index + 0) = rF(0, 1) * rDN_DX(i, 2) + rF(0, 2) * rDN_DX(i, 1);
        rB(4, index + 1) = rF(1, 1) * rDN_DX(i, 2) + rF(1, 2) * rDN_DX(i, 1);
        rB(4, index + 2) = rF(2, 1) * rDN_DX(i, 2) + rF(2, 2) * rDN_DX(i, 1);
        rB(5, index + 0) = rF(0, 2) * rDN_DX(i, 0) + rF(0, 0) * rDN_DX(i, 2);
        rB(5, index + 1) = rF(1, 2) * rDN_DX(i, 0) + rF(1, 0) * rDN_DX(i, 2);
        rB(5, index + 2) = rF(2, 2) * rDN_DX(i, 0) + rF(2, 0) * rDN_DX(i, 2);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateGreenLagrangeStrainVector(
    Vector& rStrainVector,
    const Matrix& rF
    )
{
    Matrix strain_tensor;
    strain_tensor.resize(TDim, TDim);
    Matrix identity = identity_matrix<double>(TDim);
    noalias(strain_tensor) = 0.5 * (prod(trans(rF), rF) - identity);
    rStrainVector = MathUtils<double>::StrainTensorToVector(strain_tensor, rStrainVector.size());
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateStressVectorPredictor(
    Vector& rStressVector, 
    const Matrix& rConstitutiveMAtrix, 
    const Vector& rStrainVector)
{
    noalias(rStressVector) = prod(rConstitutiveMAtrix, rStrainVector);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateAndAddInternalForcesVector(
    Vector& rRightHandSideVector, 
    const Matrix& rB, 
    const Vector& rStressVector, 
    const double IntegrationWeight
    )
{
    noalias(rRightHandSideVector) -= IntegrationWeight * prod(trans(rB), rStressVector);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateGeometricK(
    MatrixType& rLeftHandSideMatrix,
    const Matrix& rDN_DX,
    const Vector& rStressVector,
    const double IntegrationWeight
    )
{
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    Matrix stress_tensor = MathUtils<double>::StressVectorToTensor(rStressVector);
    Matrix reduced_Kg = prod(rDN_DX, IntegrationWeight * Matrix(prod(stress_tensor, trans(rDN_DX))));
    MathUtils<double>::ExpandAndAddReducedMatrix(rLeftHandSideMatrix, reduced_Kg, dimension);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateAndAddMaterialK(
    MatrixType& rLeftHandSideMatrix,
    const Matrix& B,
    const Matrix& D,
    const double IntegrationWeight,
	const double Damage
    )
{
    // Secant Constitutive Tensor
    noalias(rLeftHandSideMatrix) += (1.0 - Damage) * IntegrationWeight * prod(trans(B), Matrix(prod(D, B)));
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
double GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateDerivativesOnReferenceConfiguration(
    Matrix& rJ0,
    Matrix& rInvJ0,
    Matrix& rDN_DX,
    const IndexType PointNumber,
    IntegrationMethod ThisIntegrationMethod
    )
{
    GeometryType &r_geom = GetGeometry();
    GeometryUtils::JacobianOnInitialConfiguration(r_geom, r_geom.IntegrationPoints(ThisIntegrationMethod)[PointNumber], rJ0);
    double detJ0;
    MathUtils<double>::InvertMatrix(rJ0, rInvJ0, detJ0);
    const Matrix &rDN_De = GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod)[PointNumber];
    GeometryUtils::ShapeFunctionsGradients(rDN_De, rInvJ0, rDN_DX);
    return detJ0;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void  GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateTangentTensor(
    Matrix& rTangentTensor,
    const Vector& rStrainVectorGP,
    const Vector& rStressVectorGP,
    const Matrix& rDeformationGradientGP,
    const Matrix& rElasticMatrix,
    ConstitutiveLaw::Parameters& rValues
	)
{
	const double number_components = rStrainVectorGP.size();
	rTangentTensor.resize(number_components, number_components);
	Vector perturbed_stress, perturbed_strain;
	perturbed_strain.resize(number_components);
	perturbed_stress.resize(number_components);
    Matrix perturbed_deformation_gradient;
    const double size_1 = rDeformationGradientGP.size1();
    const double size_2 = rDeformationGradientGP.size2();
	
	for (unsigned int i_component = 0; i_component < size_1; i_component++) {
	    for (unsigned int j_component = i_component; j_component < size_2; j_component++) {
            double perturbation;
            const int component_voigt_index = this->CalculateVoigtIndex(number_components, i_component, j_component);
            this->CalculatePerturbation(rStrainVectorGP, perturbation, component_voigt_index);
            this->PerturbateDeformationGradient(perturbed_deformation_gradient, rDeformationGradientGP, perturbation, i_component, j_component);
            this->CalculateGreenLagrangeStrainVector(perturbed_strain, perturbed_deformation_gradient);
            this->IntegratePerturbedStrain(perturbed_stress, perturbed_strain, rElasticMatrix, rValues);
            const Vector& r_delta_stress = perturbed_stress - rStressVectorGP;
            this->AssignComponentsToTangentTensor(rTangentTensor, r_delta_stress, perturbation, component_voigt_index);
        }
	}
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::PerturbateDeformationGradient(
    Matrix& rPerturbedDeformationGradient,
    const Matrix& rDeformationGradientGP,
    const double Perturbation,
    const int ComponentI,
    const int ComponentJ
    )
{
    rPerturbedDeformationGradient = rDeformationGradientGP;
    if (ComponentI == ComponentJ) {
        rPerturbedDeformationGradient(ComponentI, ComponentJ) += Perturbation;
    } else {
        rPerturbedDeformationGradient(ComponentI, ComponentJ) += 0.5 * Perturbation;
        rPerturbedDeformationGradient(ComponentJ, ComponentI) += 0.5 * Perturbation;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
int GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateVoigtIndex(
    const SizeType VoigtSize,
    const int ComponentI,
    const int ComponentJ
    )
{
    if (VoigtSize == 6) {
        switch(ComponentI) {
            case 0:
                switch(ComponentJ) {
                    case 0:
                        return 0;
                    case 1:
                        return 3;
                    case 2:
                        return 5;
                    default:
                        return 0;
                }
            case 1:
                switch(ComponentJ) {
                    case 0:
                        return 3;
                    case 1:
                        return 1;
                    case 2:
                        return 4;
                    default:
                        return 0;
                }
            case 2:
                switch(ComponentJ) {
                    case 0:
                        return 5;
                    case 1:
                        return 4;
                    case 2:
                        return 2;
                    default:
                        return 0;
                }
            default:
                return 0;
        }
    } else {
        switch(ComponentI) {
            case 0:
                switch(ComponentJ) {
                    case 0:
                        return 0;
                    case 1:
                        return 2;
                    default:
                        return 0;
                }
            case 1:
                switch(ComponentJ) {
                    case 0:
                        return 2;
                    case 1:
                        return 1;
                    default:
                        return 0;
                }
            default:
                return 0;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class GenericLargeDisplacementFemDemElement<2,0>;
template class GenericLargeDisplacementFemDemElement<2,1>;
template class GenericLargeDisplacementFemDemElement<2,2>;
template class GenericLargeDisplacementFemDemElement<2,3>;
template class GenericLargeDisplacementFemDemElement<2,4>;
template class GenericLargeDisplacementFemDemElement<2,5>;
template class GenericLargeDisplacementFemDemElement<3,0>;
template class GenericLargeDisplacementFemDemElement<3,1>;
template class GenericLargeDisplacementFemDemElement<3,2>;
template class GenericLargeDisplacementFemDemElement<3,3>;
template class GenericLargeDisplacementFemDemElement<3,4>;
template class GenericLargeDisplacementFemDemElement<3,5>;

} // namespace Kratos