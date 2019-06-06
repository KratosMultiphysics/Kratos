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

#if !defined(KRATOS_GENERIC_LARGE_DISPLACEMENT_FEMDEM_ELEMENT_H_INCLUDED)
#define KRATOS_GENERIC_LARGE_DISPLACEMENT_FEMDEM_ELEMENT_H_INCLUDED

// System includes


// External include

// Project includes

#include "custom_elements/generic_small_strain_femdem_element.hpp"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{
///@name Kratos Globals
///@{
///@}

///@name Type Definitions
///@{
///@}

///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class GenericLargeDisplacementFemDemElement
 * @ingroup FemToDemApplication
 * @brief Small Displacement element for the 2D and 3D cases
 * @author Alejandro Cornejo
 */
template<unsigned int TDim, unsigned int TyieldSurf>
class GenericLargeDisplacementFemDemElement 
    : public GenericSmallStrainFemDemElement<TDim, TyieldSurf>
{
public:

    /// Counted pointer of GenericLargeDisplacementFemDemElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GenericLargeDisplacementFemDemElement);


    ///@}
    ///@name Life Cycle
    ///@{

	/// Default constructors
	GenericLargeDisplacementFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry);
	GenericLargeDisplacementFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

	///Copy constructor
	GenericLargeDisplacementFemDemElement(GenericLargeDisplacementFemDemElement const &rOther);

	/// Destructor.
	virtual ~GenericLargeDisplacementFemDemElement();

	/// Assignment operator.
	GenericLargeDisplacementFemDemElement &operator=(GenericLargeDisplacementFemDemElement const &rOther);
	Element::Pointer Create(IndexType NewId, NodesArrayType const &ThisNodes, PropertiesType::Pointer pProperties) const;
	Element::Pointer Clone(IndexType NewId, NodesArrayType const &ThisNodes) const;

	GenericLargeDisplacementFemDemElement()
	{
	}

    /**
     * this is called for non-linear analysis at the beginning of the iteration process
     */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);

    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix the elemental left hand side matrix
     * @param rRightHandSideVector the elemental right hand side
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side matrix only
     * @param rLeftHandSideMatrix the elemental left hand side matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector the elemental right hand side vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    /**
     * this computes the deformation matrix B
     */
    void CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX);
    void CalculateB2D(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX);
    void CalculateB3D(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX);

    /**
     * this computes the Green-Lagrange Strain vector from F
     */
    void CalculateGreenLagrangeStrainVector(Vector& rStrainVector, const Matrix& rF);

    /**
     * this computes stress predictor S = C:E
     */
    void CalculateStressVectorPredictor(Vector& rStressVector, const Matrix& rConstitutiveMAtrix, const Vector& rStrainVector);

    /**
     * this adds the internal forces
     */
    void CalculateAndAddInternalForcesVector(Vector& rRightHandSideVector, const Matrix& rB, const Vector& rStressVector, const double IntegrationWeight);

    /**
     * this adds geometric contribution to the LHS
     */
    void CalculateGeometricK(MatrixType& rLeftHandSideMatrix, const Matrix& rDN_DX, const Vector& rStressVector, const double IntegrationWeight);

    /**
     * this adds material contribution to the LHS when secant
     */
    void CalculateAndAddMaterialK(MatrixType& rLeftHandSideMatrix,const Matrix& B, const Matrix& D, const double IntegrationWeight, const double Damage);

    /**
     * this computes the derivatives of the kinematics in the ref conf
     */
    double CalculateDerivativesOnReferenceConfiguration(Matrix& rJ0, Matrix& rInvJ0, Matrix& rDN_DX, const IndexType PointNumber, IntegrationMethod ThisIntegrationMethod);

    /**
     * this computes constituive tangent tensor via perturbations
     */
    void CalculateTangentTensor(Matrix& rTangentTensor, const Vector& rStrainVectorGP, const Vector& rStressVectorGP, const Matrix& rDeformationGradientGP, const Matrix& rElasticMatrix, ConstitutiveLaw::Parameters& rValues);

    /**
     * this perturbates F
     */
    void PerturbateDeformationGradient(Matrix& rPerturbedDeformationGradient, const Matrix& rDeformationGradientGP, const double Perturbation, const int ComponentI, const int ComponentJ);

    /**
     * this gets the voigt index for a set of components
     */
    int CalculateVoigtIndex(const SizeType VoigtSize, const int ComponentI, const int ComponentJ);

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

}; // Class GenericLargeDisplacementFemDemElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif