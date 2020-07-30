// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Klaus B. Sautter
//

#if !defined(MEMBRANE_ELEMENT_3D_H_INCLUDED )
#define  MEMBRANE_ELEMENT_3D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"

namespace Kratos
{

  class MembraneElement
    : public Element
  {
  public:

    // Counted pointer of MembraneElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(MembraneElement);

    // Constructor using an array of nodes
    MembraneElement(IndexType NewId, GeometryType::Pointer pGeometry);

    // Constructor using an array of nodes with properties
    MembraneElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    // Destructor
    ~MembraneElement() = default;



    enum class VoigtType {
      Strain,
      Stress
    };

    enum class ConfigurationType {
      Current,
      Reference
    };

    // Name Operations

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param pGeom The pointer to the geometry of the element
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param ThisNodes The array containing nodes
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    void EquationIdVector(
      EquationIdVectorType& rResult,
      ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(
      DofsVectorType& ElementalDofList,
      ProcessInfo& rCurrentProcessInfo) override;

    void Initialize() override;

    void CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(
      VectorType& rRightHandSideVector,
      ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(
      MatrixType& rLeftHandSideMatrix,
      VectorType& rRightHandSideVector,
      ProcessInfo& rCurrentProcessInfo) override;


    void GetValuesVector(
      Vector& rValues,
      int Step = 0) override;

    void GetFirstDerivativesVector(
      Vector& rValues,
      int Step = 0) override;

    void GetSecondDerivativesVector(
      Vector& rValues,
      int Step = 0) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
      const Variable<array_1d<double, 3>>& rVariable,
      std::vector<array_1d<double, 3>>& rOutput,
      const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(
      const Variable<array_1d<double, 3>>& rVariable,
      std::vector<array_1d<double, 3>>& rOutput,
      const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMassMatrix(MatrixType& rMassMatrix,ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLumpedMassVector(VectorType& rMassVector);

    void AddExplicitContribution(
      const VectorType& rRHSVector, const Variable<VectorType>& rRHSVariable,
      Variable<array_1d<double, 3>>& rDestinationVariable,
      const ProcessInfo& rCurrentProcessInfo) override;

    void AddExplicitContribution(
      const VectorType& rRHSVector,
      const Variable<VectorType>& rRHSVariable,
      Variable<double >& rDestinationVariable,
      const ProcessInfo& rCurrentProcessInfo) override;

    void Calculate(const Variable<Matrix>& rVariable,
      Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo) override;


    void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) override;

  private:
     /**
     * @brief Calculates the covariant base vectors
     * @param rBaseVectors The base vectors to be calculated
     * @param rShapeFunctionGradientValues Current shapefunction gradients
     * @param Configuration Reference/Current
     */
    void CovariantBaseVectors(array_1d<Vector,2>& rBaseVectors,
     const Matrix& rShapeFunctionGradientValues, const ConfigurationType& rConfiguration);

      /**
     * @brief Calculates the covariant metric
     * @param rMetric The metric to be calculated
     * @param rBaseVectorCovariant Covariant base vectors
     */
    void CovariantMetric(Matrix& rMetric,const array_1d<Vector,2>& rBaseVectorCovariant);

     /**
     * @brief Calculates the contra variant base vectors
     * @param rBaseVectors The base vectors to be calculated
     * @param rContraVariantMetric Contra variant metric
     * @param rCovariantBaseVectors Covariant base vectors
     */
    void ContraVariantBaseVectors(array_1d<Vector,2>& rBaseVectors,const Matrix& rContraVariantMetric,
      const array_1d<Vector,2> rCovariantBaseVectors);

      /**
     * @brief Calculates the contra variant metric
     * @param rMetric The metric to be calculated
     * @param rCovariantMetric Covariant metric
     */
    void ContravariantMetric(Matrix& rMetric,const Matrix& rCovariantMetric);


      /**
     * @brief Calculates 1st derivative of the current covariant base vectors
     * @param rBaseVectors The derived base bectors
     * @param rShapeFunctionGradientValues Current shapefunction gradients
     * @param DofR current degree of freedom 1
     */
    void DeriveCurrentCovariantBaseVectors(array_1d<Vector,2>& rBaseVectors,
     const Matrix& rShapeFunctionGradientValues, const SizeType DofR);


      /**
     * @brief Calculates 2nd derivative of the current covariant metric
     * @param rMetric The derived metric
     * @param rShapeFunctionGradientValues Current shapefunction gradients
     * @param DofR current degree of freedom 1
     * @param DofS current degree of freedom 2
     */
    void Derivative2CurrentCovariantMetric(Matrix& rMetric,
      const Matrix& rShapeFunctionGradientValues, const SizeType DofR, const SizeType DofS);


      /**
     * @brief Calculates the determinant of the Jacobian
     * @param rDetJacobi The determinant of the Jacobian
     * @param rReferenceBaseVectors Reference base vectors
     */
    void JacobiDeterminante(double& rDetJacobi, const array_1d<Vector,2>& rReferenceBaseVectors);


      /**
     * @brief Calculates 2nd derivative of the green lagrange strain
     * @param rStrain The derived strain
     * @param rShapeFunctionGradientValues Current shapefunction gradients
     * @param DofR current degree of freedom 1
     * @param DofS current degree of freedom 2
     * @param rTransformationMatrix local coordinate system transformation
     */
    void Derivative2StrainGreenLagrange(Vector& rStrain,
      const Matrix& rShapeFunctionGradientValues, const SizeType DofR, const SizeType DofS,
      const Matrix& rTransformationMatrix);



      /**
     * @brief Calculates 1st derivative of the green lagrange strain
     * @param rStrain The derived strain
     * @param rShapeFunctionGradientValues Current shapefunction gradients
     * @param DofR current degree of freedom 1
     * @param rCurrentCovariantBaseVectors current covariant base vectors
     * @param rTransformationMatrix local coordinate system transformation
     */
    void DerivativeStrainGreenLagrange(Vector& rStrain, const Matrix& rShapeFunctionGradientValues, const SizeType DofR,
      const array_1d<Vector,2> rCurrentCovariantBaseVectors, const Matrix& rTransformationMatrix);


      /**
     * @brief Calculates green lagrange strain
     * @param rStrain The strain
     * @param rReferenceCoVariantMetric reference covariant metric
     * @param rCurrentCoVariantMetric current covariant metric
     * @param rTransformationMatrix local coordinate system transformation
     */
    void StrainGreenLagrange(Vector& rStrain, const Matrix& rReferenceCoVariantMetric,const Matrix& rCurrentCoVariantMetric,
       const Matrix& rTransformationMatrix);

      /**
     * @brief Calculates the piola-kirchhoff-2 stress
     * @param rStress The stress
     * @param rReferenceContraVariantMetric reference contra variant metric
     * @param rReferenceCoVariantMetric reference covariant metric
     * @param rCurrentCoVariantMetric current covariant metric
     * @param rTransformedBaseVectors local coordinate system
     * @param rTransformationMatrix local coordinate system transformation
     * @param rIntegrationPointNumber current integration point number
     */
    void MaterialResponse(Vector& rStress,
      const Matrix& rReferenceContraVariantMetric,const Matrix& rReferenceCoVariantMetric,const Matrix& rCurrentCoVariantMetric,
      const array_1d<Vector,2>& rTransformedBaseVectors, const Matrix& rTransformationMatrix, const SizeType& rIntegrationPointNumber,
      Matrix& rTangentModulus);


      /**
     * @brief Adds pre-stress to a given stress vector
     * @param rStress The stress
     * @param rTransformedBaseVectors local coordinate system
     */
    void AddPreStressPk2(Vector& rStress, const array_1d<Vector,2>& rTransformedBaseVectors);

      /**
     * @brief Calculates 1st derivative of the current covariant metric
     * @param rMetric The derived metric
     * @param rShapeFunctionGradientValues Current shapefunction gradients
     * @param DofR current degree of freedom 1
     * @param rCurrentCovariantBaseVectors current covariant base vectors
     */
    void DerivativeCurrentCovariantMetric(Matrix& rMetric,
      const Matrix& rShapeFunctionGradientValues, const SizeType DofR, const array_1d<Vector,2> rCurrentCovariantBaseVectors);


      /**
     * @brief Calculates the internal forces
     * @param rInternalForces The internal forces
     * @param ThisMethod numerical integration method
     */
    void InternalForces(Vector& rInternalForces,const IntegrationMethod& ThisMethod);


      /**
     * @brief Calculates the stiffness matrix
     * @param rStiffnessMatrix The stiffness matrix
     * @param ThisMethod numerical integration method
     */
    void TotalStiffnessMatrix(Matrix& rStiffnessMatrix,const IntegrationMethod& ThisMethod);


      /**
     * @brief Calculates initial stress part of the total stiffness matrix
     * @param rEntryIJ the matrix entry to be calculated
     * @param rStressVector Current stress
     * @param rPositionI current degree of freedom 1
     * @param rPositionJ current degree of freedom 2
     * @param rShapeFunctionGradientValues Current shapefunction gradients
     * @param rTransformationMatrix local coordinate system transformation
     */
    void InitialStressStiffnessMatrixEntryIJ(double& rEntryIJ,
      const Vector& rStressVector,
      const SizeType& rPositionI, const SizeType& rPositionJ, const Matrix& rShapeFunctionGradientValues,
      const Matrix& rTransformationMatrix);



      /**
     * @brief Calculates material part of the total stiffness matrix
     * @param rEntryIJ the matrix entry to be calculated
     * @param rMaterialTangentModulus material tangent modulus
     * @param rPositionI current degree of freedom 1
     * @param rPositionJ current degree of freedom 2
     * @param rShapeFunctionGradientValues Current shapefunction gradients
     * @param rCurrentCovariantBaseVectors current covariant base vectors
     * @param rTransformationMatrix local coordinate system transformation
     */
    void MaterialStiffnessMatrixEntryIJ(double& rEntryIJ,
      const Matrix& rMaterialTangentModulus,
      const SizeType& rPositionI, const SizeType& rPositionJ, const Matrix& rShapeFunctionGradientValues,
      const array_1d<Vector,2>& rCurrentCovariantBaseVectors,const Matrix& rTransformationMatrix);



      /**
     * @brief Transform strains to a given local coordinate system
     * @param rStrains the strains in the given local coordinate system
     * @param rReferenceStrains the strains in the reference coordinate system
     * @param rTransformationMatrix local coordinate system transformation
     */
    void TransformStrains(Vector& rStrains, Vector& rReferenceStrains, const Matrix& rTransformationMatrix);


      /**
     * @brief Creates a given local coordinate system
     * @param rBaseVectors the new local coordinate system
     * @param rLocalBaseVectors the base coordinate system
     */
    void TransformBaseVectors(array_1d<Vector,2>& rBaseVectors,
     const array_1d<Vector,2>& rLocalBaseVectors);


      /**
     * @brief Creates the transformation matrix for a given local coordinate system
     * @param rTransformationMatrix the transformation matrix
     * @param rTransformedBaseVectors local coordinate system
     * @param rLocalReferenceBaseVectors the base coordinate system
     */
    template <class T>
    void InPlaneTransformationMatrix(Matrix& rTransformationMatrix, const array_1d<Vector,2>& rTransformedBaseVectors,
      const T& rLocalReferenceBaseVectors);


      /**
     * @brief Calculates the principal vectors
     * @param rPrincipalVector the principal vectors
     * @param rNonPrincipalVector reference state
     */
    void PrincipalVector(Vector& rPrincipalVector, const Vector& rNonPrincipalVector);


    void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Vector >& rVariable,
        std::vector< Vector >& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

    void DeformationGradient(Matrix& rDeformationGradient, double& rDetDeformationGradient,
       const array_1d<Vector,2>& rCurrentCovariantBase, const array_1d<Vector,2>& rReferenceContraVariantBase);


    void CalculateAndAddBodyForce(VectorType& rRightHandSideVector);

    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector; /// The vector containing the constitutive laws
    double mReferenceArea = 0.0;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization
    MembraneElement() = default;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}

  };	// class MembraneElement.

}	// namespace Kratos.

#endif // KRATOS_MEMBRANE_ELEMENT_3D_H_INCLUDED  defined
