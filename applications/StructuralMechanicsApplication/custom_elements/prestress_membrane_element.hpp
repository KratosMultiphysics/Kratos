// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Long Chen
//

#if !defined(KRATOS_MEMBRANE_ELEMENT_3D_H_INCLUDED )
#define  KRATOS_MEMBRANE_ELEMENT_3D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "custom_utilities/structural_mechanics_math_utilities.hpp"

namespace Kratos
{

/** \brief PrestressMembraneElement
* This is a 3D nonlinear isoparametric membrane element, which deals with large displacements.
* Its functionalities for Formfinding (using Updated Reference Strategy) and Cutting Pattern are
* implemented based on the dissertations of Roland Wuechner, Johannes Linhard, and Falko Dieringer
* at TUM
*/
  class PrestressMembraneElement
    : public Element
  {
  public:

    // Counted pointer of MembraneElement
    KRATOS_CLASS_POINTER_DEFINITION(PrestressMembraneElement);

    // Constructor using an array of nodes
    PrestressMembraneElement(IndexType NewId, GeometryType::Pointer pGeometry);

    // Constructor using an array of nodes with properties
    PrestressMembraneElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    // Destructor
    ~PrestressMembraneElement() override;


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

    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(
      VectorType& rRightHandSideVector,
      ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(
      MatrixType& rLeftHandSideMatrix,
      VectorType& rRightHandSideVector,
      ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
      const Variable<Matrix>& rVariable,
      std::vector<Matrix>& Output,
      const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMassMatrix(
      MatrixType& rMassMatrix,
      ProcessInfo& rCurrentProcessInfo) override;

    void CalculateDampingMatrix(
      MatrixType& rDampingMatrix,
      ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep(
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

    void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
      std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

  protected:


  private:
    ///@name Static Member Variables
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;
    Vector mDetJ0;
    double mTotalDomainInitialSize;

    std::vector< Matrix >              mGVector;
    std::vector< array_1d<double, 3> > mGab0;

    bool mAnisotropicPrestress;                               // determines if isotropic or anisotropic prestress is applied
    std::vector< array_1d<double, 3> > mG1Initial;            // Base vector 1 in initial reference configuration
    std::vector< array_1d<double, 3> > mG2Initial;            // Base vector 2 in initial reference configuration
    std::vector< array_1d<double, 3> > mG3Initial;            // Base vector 3 in initial reference configuration


    void CalculateAll(
      MatrixType& rLeftHandSideMatrix,
      VectorType& rRightHandSideVector,
      const ProcessInfo& rCurrentProcessInfo,
      const bool CalculateStiffnessMatrixFlag,
      const bool CalculateResidualVectorFlag);

    void CalculateAndAddKm(
      Matrix& rK,
      Matrix& rB,
      Matrix& rD,
      const double& rWeight);


    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    void CalculateAndAddNonlinearKm(
        Matrix& rK,
        Matrix& rB11,
        Matrix& rB22,
        Matrix& rB12,
        Vector& rSD,
        const double& rWeight);

    void CalculateQ(
      BoundedMatrix<double, 3, 3>& rQ,
      const unsigned int PointNumber);

    void CalculateB(
        Matrix& B,
        const BoundedMatrix<double, 3, 3>& Q,
        const Matrix& DN_De,
        const array_1d<double, 3>& g1,
        const array_1d<double, 3>& g2);

    void CalculateStrain(
        Vector& rStrainVector,
        array_1d<double, 3>& rgab,
        array_1d<double, 3>& rGab);

    void CalculateAndAdd_BodyForce(
      const Vector& rN,
      const ProcessInfo& rCurrentProcessInfo,
      array_1d<double, 3>& BodyForce,
      VectorType& rRightHandSideVector,
      const double& rWeight);

    void CalculateAndAdd_PressureForce(
      VectorType& rResidualVector,
      const Vector& N,
      const array_1d<double, 3>& rv3,
      const double& rPressure,
      const double& rWeight,
      const ProcessInfo& rCurrentProcessInfo);

    void CalculateMetricDeformed(const unsigned int PointNumber,
        Matrix DN_De,
        array_1d<double, 3>& rgab,
        array_1d<double, 3>& rg1,
        array_1d<double, 3>& rg2);


    void CalculateSecondVariationStrain(Matrix DN_De,
        Matrix& Strain_locCartesian11,
        Matrix& Strain_locCartesian22,
        Matrix& Strain_locCartesian12,
        BoundedMatrix<double, 3, 3>& Q);

    void InitializeFormfinding(const unsigned int rIntegrationPointSize);

    void ProjectPrestress(const unsigned int PointNumber);

    void UpdatePrestress(const unsigned int PointNumber);

    void ComputePrestress(const unsigned int rIntegrationPointSize);

    void ComputeBaseVectors(const GeometryType::IntegrationPointsArrayType& rIntegrationPoints);

    void InitializeMaterial(const unsigned int NumberIntegrationPoints);

    void ComputeContravariantBaseVectors(
                        array_1d<double, 3>& rG1Contra,
                        array_1d<double, 3>& rG2Contra,
                        const unsigned int PointNumber);

    void ComputeRelevantCoSys(const unsigned int PointNumber,
             array_1d<double, 3>& rg1,array_1d<double, 3>& rg2,array_1d<double, 3>& rg3, array_1d<double, 3>& rgab,
             array_1d<double, 3>& rG3,
             array_1d<double, 3>& rE1Tot, array_1d<double, 3>& rE2Tot,array_1d<double, 3>& rE3Tot,
             array_1d<double, 3>& rE1,array_1d<double, 3>& rE2,array_1d<double, 3>& rE3,
             array_1d<double, 3>& rBaseRefContraTot1,array_1d<double, 3>& rBaseRefContraTot2);

    void ComputeEigenvaluesDeformationGradient(const unsigned int PointNumber,
                    BoundedMatrix<double,3,3>& rOrigin, BoundedMatrix<double,3,3>& rTarget, BoundedMatrix<double,3,3>& rTensor,
                    const array_1d<double, 3>& rBaseRefContraTot1, const array_1d<double, 3>& rBaseRefContraTot2,
                    const array_1d<double, 3>& rE1Tot, const array_1d<double, 3>& rE2Tot, const array_1d<double, 3>& rE3Tot,
                    const array_1d<double, 3>& rgab,
                    double& rLambda1, double& rLambda2);

    void ComputeEigenvectorsDeformationGradient(const unsigned int PointNumber,
                                BoundedMatrix<double,3,3>& rTensor, BoundedMatrix<double,3,3>& rOrigin,
                                const BoundedMatrix<double,3,3>& rDeformationGradientTotal,
                                const array_1d<double, 3>& rE1Tot, const array_1d<double, 3>& rE2Tot,
                                const double Lambda1, const double Lambda2,
                                BoundedMatrix<double,3,3>& rNAct);

    void ModifyPrestress(const unsigned int PointNumber,
                    BoundedMatrix<double,3,3>& rOrigin, BoundedMatrix<double,3,3>& rTarget,BoundedMatrix<double,3,3>& rTensor,
                    const array_1d<double, 3>& rE1, const array_1d<double, 3>& rE2, const array_1d<double, 3>& rE3, const array_1d<double, 3>& rG3,
                    const array_1d<double, 3>& rg1, const array_1d<double, 3>& rg2, const array_1d<double, 3>& rg3, const BoundedMatrix<double,3,3>& rNAct,
                    const double Lambda1, const double Lambda2);

    const Matrix CalculateDeformationGradient(const unsigned int PointNumber);

    int  Check(const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization
    PrestressMembraneElement() {}

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}

  };	// class MembraneElement.

}	// namespace Kratos.

#endif // KRATOS_MEMBRANE_ELEMENT_3D_H_INCLUDED  defined
