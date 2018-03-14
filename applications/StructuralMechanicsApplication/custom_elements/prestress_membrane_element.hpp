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

    Element::Pointer Create(
      IndexType NewId,
      NodesArrayType const& ThisNodes,
      PropertiesType::Pointer pProperties) const override;

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
      Vector& values,
      int Step = 0) override;

    void GetFirstDerivativesVector(
      Vector& values,
      int Step = 0) override;

    void GetSecondDerivativesVector(
      Vector& values,
      int Step = 0) override;

    void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
      std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

  protected:


  private:
    ///@name Static Member Variables
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;
    Vector mDetJ0;
    double mTotalDomainInitialSize;

    std::vector< array_1d<double, 3> > mStrainsVector;	      //container of Strain // TODO is this needed?
    std::vector< array_1d<double, 6> > mStressesVector;	      //container of Stress // TODO is this needed?
    std::vector< array_1d<double, 6> > mCauchyStressesVector;	//container of Stress // TODO is this needed?

    std::vector< Matrix >              mGVector;
    std::vector< array_1d<double, 3> > mGab0;
    std::vector< array_1d<double, 3> > mG1;                   // Base vector 1 in updated reference configuration
    std::vector< array_1d<double, 3> > mG2;                   // Base vector 2 in updated reference configuration

    // Using this variable is a potential bug if the element is not used in formfinding!
    // In the future this should be a Processinfo Variable (e.g. FROMFINDING_STEP), which
    // is set by the Fromfinding Strategy
    // The element can then check if this Var is set and use it accordingly
    unsigned int mStep;                                       // Simulation step for formfinding

    bool mAnisotropicPrestress;                               // determines if isotropic or anisotropic prestress is applied
    std::vector< array_1d<double, 3> > mG1Initial;            // Base vector 1 in initial reference configuration
    std::vector< array_1d<double, 3> > mG2Initial;            // Base vector 2 in initial reference configuration
    std::vector< array_1d<double, 3> > mG3Initial;            // Base vector 2 in initial reference configuration


    void CalculateAll(
      MatrixType& rLeftHandSideMatrix,
      VectorType& rRightHandSideVector,
      const ProcessInfo& rCurrentProcessInfo,
      bool CalculateStiffnessMatrixFlag,
      bool CalculateResidualVectorFlag);

    void CalculateAndAddKm(
      Matrix& K,
      Matrix& msB,
      Matrix& msD,
      double weight);

    void InitializeNonLinearIteration();

    void CalculateAndAddNonlinearKm(
        Matrix& K,
        Matrix& B11,
        Matrix& B22,
        Matrix& B12,
        Vector& SD,
        double weight);

    //void CalculateAndAddKg(
    //  Matrix& K,
    //  boost::numeric::ublas::bounded_matrix<double, 3, 3>& msQ,
    //  const Matrix& DN_De,
    //  Vector& msStressVector,
    //  double weight);

    //void CalculateAndSubKp(
    //  Matrix& K,
    //  array_1d<double, 3>& ge,
    //  array_1d<double, 3>& gn,
    //  const Matrix& DN_De,
    //  const Vector& N,
    //  double pressure,
    //  double weight);

    void ClearNodalForces();

    //void AddExplicitContribution(
    //  const VectorType& rRHSVector,
    //  const Variable<VectorType>& rRHSVariable,
    //  Variable<array_1d<double, 3> >& rDestinationVariable,
    //  const ProcessInfo& rCurrentProcessInfo);


    //void MakeCrossMatrix(
    //  boost::numeric::ublas::bounded_matrix<double, 3, 3>& M,
    //  array_1d<double, 3>& U);


    void CalculateQ(
      boost::numeric::ublas::bounded_matrix<double, 3, 3>& msQ,
      Matrix& msG);

    void CalculateB(
        Matrix& B,
        boost::numeric::ublas::bounded_matrix<double, 3, 3>& Q,
        const Matrix& DN_De,
        const array_1d<double, 3>& g1,
        const array_1d<double, 3>& g2);

    //void CalculateJ(
    //  boost::numeric::ublas::bounded_matrix<double, 2, 2>& j,
    //  array_1d<double, 3>& ge,
    //  array_1d<double, 3>& gn,
    //  array_1d<double, 3>& v3);

    void CalculateStrain(
        Vector& StrainVector,
        array_1d<double, 3>& gab,
        array_1d<double, 3>& gab_ref);

    void CalculateAndAdd_BodyForce(
      const Vector& N,
      const ProcessInfo& rCurrentProcessInfo,
      array_1d<double, 3>& BodyForce,
      VectorType& rRightHandSideVector,
      double weight);

    void CalculateAndAdd_PressureForce(
      VectorType& residualvector,
      const Vector& N,
      const array_1d<double, 3>& v3,
      double pressure,
      double weight,
      const ProcessInfo& rCurrentProcessInfo);

    void CalculateMetricDeformed(const unsigned int& PointNumber,
        Matrix DN_De,
        array_1d<double, 3>& gab,
        array_1d<double, 3>& g1,
        array_1d<double, 3>& g2);


    void CalculateSecondVariationStrain(Matrix DN_De,
        Matrix& Strain_locCartesian11,
        Matrix& Strain_locCartesian22,
        Matrix& Strain_locCartesian12,
        boost::numeric::ublas::bounded_matrix<double, 3, 3>& Q,
        array_1d<double, 3>& g1,
        array_1d<double, 3>& g2);

    void CalculateMembraneElasticityTensor(
        Matrix& D
        );

    void TransformPrestress(const unsigned int PointNumber);

    void UpdatePrestress(const unsigned int PointNumber);

    void PrestressComputation(const unsigned int PointNumber);

    void ComputeBaseVectors(const GeometryType::IntegrationPointsArrayType& rIntegrationPoints);

    void InitializeMaterial(const unsigned int NumberIntegrationPoints);

    void ComputeContravariantBaseVectors(
                        array_1d<double, 3>& rG1Contra, 
                        array_1d<double, 3>& rG2Contra,
                        const unsigned int& rPointNumber);

    void ComputeRelevantCoSys(const unsigned int PointNumber,
             array_1d<double, 3>& rg1,array_1d<double, 3>& rg2,array_1d<double, 3>& rg3, array_1d<double, 3>& rgab,
             array_1d<double, 3>& rG3,
             array_1d<double, 3>& rE1Tot, array_1d<double, 3>& rE2Tot,array_1d<double, 3>& rE3Tot,
             array_1d<double, 3>& rE1,array_1d<double, 3>& rE2,array_1d<double, 3>& rE3,
             array_1d<double, 3>& rBaseRefContraTot1,array_1d<double, 3>& rBaseRefContraTot2);

    void ComputeEigenvaluesDeformationGradient(const unsigned int PointNumber,
                    bounded_matrix<double,3,3>& rOrigin, bounded_matrix<double,3,3>& rTarget, bounded_matrix<double,3,3>& rTensor,
                    const array_1d<double, 3>& rBaseRefContraTot1, const array_1d<double, 3>& rBaseRefContraTot2,
                    const array_1d<double, 3>& rE1Tot, const array_1d<double, 3>& rE2Tot, const array_1d<double, 3>& rE3Tot,
                    const array_1d<double, 3>& rgab,
                    double& rLambda1, double& rLambda2);

    void ComputeEigenvectorsDeformationGradient(const unsigned int PointNumber,
                                bounded_matrix<double,3,3>& rTensor, bounded_matrix<double,3,3>& rOrigin,
                                const bounded_matrix<double,3,3>& rDeformationGradientTotal,
                                const array_1d<double, 3>& rE1Tot, const array_1d<double, 3>& rE2Tot,
                                const double Lambda1, const double Lambda2,
                                bounded_matrix<double,3,3>& rNAct);

    void ModifyPrestress(const unsigned int PointNumber,
                    bounded_matrix<double,3,3>& rOrigin, bounded_matrix<double,3,3>& rTarget,bounded_matrix<double,3,3>& rTensor,
                    const array_1d<double, 3>& rE1, const array_1d<double, 3>& rE2, const array_1d<double, 3>& rE3, const array_1d<double, 3>& rG3,
                    const array_1d<double, 3>& rg1, const array_1d<double, 3>& rg2, const array_1d<double, 3>& rg3, const bounded_matrix<double,3,3>& rNAct,
                    const double Lambda1, const double Lambda2);

    const Matrix CalculateDeformationGradient(const unsigned int& rPointNumber);

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
