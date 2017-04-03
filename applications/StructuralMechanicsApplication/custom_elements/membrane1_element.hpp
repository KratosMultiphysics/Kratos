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

#if !defined(KRATOS_MEMBRANE1_ELEMENT_H_INCLUDED )
#define  KRATOS_MEMBRANE1_ELEMENT_H_INCLUDED



// System includes

// External includes

// Project includes
#include "includes/element.h"


namespace Kratos
{

/** \brief MembraneElement3D
* This is a 3D nonlinear isoparametric membrane element, which deals with large displacements
* Its functionalities for Formfinding (Updated Reference Strategy) and Cutting Pattern are
* implemented based on the dissertations of Roland Wuechner, Johannes Linhard, and Falko Dieringer
* at TUM
*/
  class Membrane1Element
    : public Element
  {
  public:

    // Counted pointer of MembraneElement
    KRATOS_CLASS_POINTER_DEFINITION(Membrane1Element);

    // Constructor using an array of nodes
    Membrane1Element(IndexType NewId, GeometryType::Pointer pGeometry);

    // Constructor using an array of nodes with properties
    Membrane1Element(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    // Destructor
    virtual ~Membrane1Element();


    // Name Operations

    Element::Pointer Create(
      IndexType NewId,
      NodesArrayType const& ThisNodes,
      PropertiesType::Pointer pProperties) const;

    void EquationIdVector(
      EquationIdVectorType& rResult,
      ProcessInfo& rCurrentProcessInfo);

    void GetDofList(
      DofsVectorType& ElementalDofList,
      ProcessInfo& rCurrentProcessInfo);

    void Initialize();

    void CalculateRightHandSide(
      VectorType& rRightHandSideVector,
      ProcessInfo& rCurrentProcessInfo);

    void CalculateLocalSystem(
      MatrixType& rLeftHandSideMatrix,
      VectorType& rRightHandSideVector,
      ProcessInfo& rCurrentProcessInfo);

    void CalculateOnIntegrationPoints(
      const Variable<Matrix>& rVariable,
      std::vector<Matrix>& Output,
      const ProcessInfo& rCurrentProcessInfo);

    void CalculateMassMatrix(
      MatrixType& rMassMatrix,
      ProcessInfo& rCurrentProcessInfo);

    void CalculateDampingMatrix(
      MatrixType& rDampingMatrix,
      ProcessInfo& rCurrentProcessInfo);

    void FinalizeSolutionStep(
      ProcessInfo& rCurrentProcessInfo);

    void GetValuesVector(
      Vector& values,
      int Step = 0);

    void GetFirstDerivativesVector(
      Vector& values,
      int Step = 0);

    void GetSecondDerivativesVector(
      Vector& values,
      int Step = 0);

    void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
      std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);


  protected:


  private:
    ///@name Static Member Variables


    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;
    Geometry< Point<3, double> >::Pointer  mpReferenceGeometry;

    Vector mDetJ0;

    double mTotalDomainInitialSize;
    double mdensity;
    double mThickness0; //thickness in the reference configuration
    double mThickness;									// thickness in actual configuration

    std::vector< array_1d<double, 3> > mStrainsVector;	//container of Strain
    std::vector< array_1d<double, 6> > mStressesVector;	//container of Stress
    std::vector< array_1d<double, 6> > mCauchyStressesVector;	//container of Stress


    std::vector< array_1d<double, 3> >  mV1;
    std::vector< array_1d<double, 3> >  mV2;
    std::vector< Matrix >              mG_Vector;

    std::vector< array_1d<double, 3> > mGab0;
    std::vector< array_1d<double, 3> > mG1;
    std::vector< array_1d<double, 3> > mG2;

    array_1d<double, 3> mPreStress;         // Pre-Stress which cannot be read through mpda file. this is a temporary solution



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

    void CrossProduct(
      array_1d<double, 3>& cross,
      array_1d<double, 3>& a,
      array_1d<double, 3>& b);

    //void SubtractMatrix(
    //  MatrixType& Destination,
    //  boost::numeric::ublas::bounded_matrix<double, 3, 3>& InputMatrix,
    //  int InitialRow,
    //  int InitialCol);

    //void ExpandReducedMatrix(
    //  Matrix& Destination,
    //  Matrix& ReducedMatrix);

    void CalculateQ(
      boost::numeric::ublas::bounded_matrix<double, 3, 3>& msQ,
      Matrix& msG);

    void Membrane1Element::CalculateB(
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

    //// this function transforms the local stress (with 3 components)
    //// to the global one (with 6 components)
    //void Calculate_GlobalStressVector(
    //  array_1d<double, 6>& GlobalVector,
    //  Vector& LocalStressVector,
    //  array_1d<double, 3>& v1,
    //  array_1d<double, 3>& v2);

    ////auxiliary function needed in the calculation of output stresses
    //inline array_1d<double, 6> VoigtTensorComponents(
    //  array_1d<double, 3>& a,
    //  array_1d<double, 3>& b);

    void CalculateMetricDeformed(unsigned int& PointNumber,
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


    int  Check(const ProcessInfo& rCurrentProcessInfo);

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization
    Membrane1Element() {}

    void save(Serializer& rSerializer) const
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
        rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
      rSerializer.save("ReferenceGeometry", mpReferenceGeometry);
      rSerializer.save("DetJ0", mDetJ0);
      rSerializer.save("TotalDomainInitialSize", mTotalDomainInitialSize);
      rSerializer.save("density", mdensity);
      rSerializer.save("Thickness0", mThickness0);
      rSerializer.save("Thickness", mThickness);
      rSerializer.save("StrainsVector", mStrainsVector);
      rSerializer.save("StressesVector", mStressesVector);
      rSerializer.save("CauchyStressesVector", mCauchyStressesVector);
      rSerializer.save("Thickness0", mThickness0);
      rSerializer.save("V1", mV1);
      rSerializer.save("V2", mV2);
      rSerializer.save("G_Vector", mG_Vector);
    }

    void load(Serializer& rSerializer)
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
        rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
      rSerializer.load("ReferenceGeometry", mpReferenceGeometry);
      rSerializer.load("DetJ0", mDetJ0);
      rSerializer.load("TotalDomainInitialSize", mTotalDomainInitialSize);
      rSerializer.load("density", mdensity);
      rSerializer.load("Thickness0", mThickness0);
      rSerializer.load("Thickness", mThickness);
      rSerializer.load("StrainsVector", mStrainsVector);
      rSerializer.load("StressesVector", mStressesVector);
      rSerializer.load("CauchyStressesVector", mCauchyStressesVector);
      rSerializer.load("Thickness0", mThickness0);
      rSerializer.load("V1", mV1);
      rSerializer.load("V2", mV2);
      rSerializer.load("G_Vector", mG_Vector);
    }

    ///@}

  };	// class MembraneElement.

}	// namespace Kratos.

#endif // KRATOS_MEMBRANE1_ELEMENT_H_INCLUDED  defined 
