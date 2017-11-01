// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_MEMBRANE_ELEMENT_H_INCLUDED )
#define  KRATOS_MEMBRANE_ELEMENT_H_INCLUDED



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
    KRATOS_CLASS_POINTER_DEFINITION( MembraneElement );

    // Constructor using an array of nodes
    MembraneElement(IndexType NewId, GeometryType::Pointer pGeometry);

    // Constructor using an array of nodes with properties
    MembraneElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    // Destructor
    ~MembraneElement() override;


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

    void Initialize () override;

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
    Geometry<Point>::Pointer  mpReferenceGeometry;

    Vector mDetJ0;

    double mTotalDomainInitialSize;
    double mdensity;
    double mThickness0; //thickness in the reference configuration


    Vector mThickness;									//container of thickness
    std::vector< array_1d<double,3> > mStrainsVector;	//container of Strain
    std::vector< array_1d<double,6> > mStressesVector;	//container of Stress
    std::vector< array_1d<double,6> > mCauchyStressesVector;	//container of Stress


    std::vector< array_1d<double,3> >  mV1;
    std::vector< array_1d<double,3> >  mV2;
    std::vector< Matrix >              mG_Vector;

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

    void CalculateAndAddKg(
        Matrix& K,
        boost::numeric::ublas::bounded_matrix<double,3,3>& msQ,
        const Matrix& DN_De,
        Vector& msStressVector,
        double weight);

    void CalculateAndSubKp(
        Matrix& K,
        array_1d<double,3>& ge,
        array_1d<double,3>& gn,
        const Matrix& DN_De,
        const Vector& N,
        double pressure,
        double weight);

    void ClearNodalForces();

    void AddExplicitContribution(
        const VectorType& rRHSVector,
        const Variable<VectorType>& rRHSVariable,
        Variable<array_1d<double,3> >& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo) override;


    void MakeCrossMatrix(
        boost::numeric::ublas::bounded_matrix<double,3,3>& M,
        array_1d<double,3>& U);

    void CrossProduct(
        array_1d<double,3>& cross,
        array_1d<double,3>& a,
        array_1d<double,3>& b);

    void SubtractMatrix(
        MatrixType& Destination,
        boost::numeric::ublas::bounded_matrix<double,3,3>& InputMatrix,
        int InitialRow,
        int InitialCol);

    void ExpandReducedMatrix(
        Matrix& Destination,
        Matrix& ReducedMatrix);

    void CalculateQ(
        boost::numeric::ublas::bounded_matrix<double,3,3>& msQ,
        Matrix& msG);

    void CalculateB(
        Matrix& msB,
        boost::numeric::ublas::bounded_matrix<double,3,3>& msQ,
        const Matrix& DN_De,
        array_1d<double,3>& ge,
        array_1d<double,3>& gn);

    void CalculateJ(
        boost::numeric::ublas::bounded_matrix<double,2,2>& j,
        array_1d<double,3>& ge,
        array_1d<double,3>& gn,
        array_1d<double,3>& v3);

    void CalculateStrain(
        Vector& StrainVector,
        boost::numeric::ublas::bounded_matrix<double,2,2>& C);

    void CalculateAndAdd_BodyForce(
        const Vector& N,
        const ProcessInfo& rCurrentProcessInfo,
        array_1d<double,3>& BodyForce,
        VectorType& rRightHandSideVector,
        double weight);

    void CalculateAndAdd_PressureForce(
        VectorType& residualvector,
        const Vector& N,
        const array_1d<double,3>& v3,
        double pressure,
        double weight,
        const ProcessInfo& rCurrentProcessInfo);

    // this function transforms the local stress (with 3 components)
    // to the global one (with 6 components)
    void Calculate_GlobalStressVector(
        array_1d<double,6>& GlobalVector,
        Vector& LocalStressVector,
        array_1d<double,3>& v1,
        array_1d<double,3>& v2);

    //auxiliary function needed in the calculation of output stresses
    inline array_1d<double,6> VoigtTensorComponents(
        array_1d<double,3>& a,
        array_1d<double,3>& b);

    int  Check( const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization
    MembraneElement() {}

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
        rSerializer.save("ConstitutiveLawVector",mConstitutiveLawVector);
        rSerializer.save("ReferenceGeometry",mpReferenceGeometry);
        rSerializer.save("DetJ0",mDetJ0);
        rSerializer.save("TotalDomainInitialSize",mTotalDomainInitialSize);
        rSerializer.save("density",mdensity);
        rSerializer.save("Thickness0",mThickness0);
        rSerializer.save("Thickness",mThickness);
        rSerializer.save("StrainsVector",mStrainsVector);
        rSerializer.save("StressesVector",mStressesVector);
        rSerializer.save("CauchyStressesVector",mCauchyStressesVector);
        rSerializer.save("Thickness0",mThickness0);
        rSerializer.save("V1",mV1);
        rSerializer.save("V2",mV2);
        rSerializer.save("G_Vector",mG_Vector);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
        rSerializer.load("ConstitutiveLawVector",mConstitutiveLawVector);
        rSerializer.load("ReferenceGeometry",mpReferenceGeometry);
        rSerializer.load("DetJ0",mDetJ0);
        rSerializer.load("TotalDomainInitialSize",mTotalDomainInitialSize);
        rSerializer.load("density",mdensity);
        rSerializer.load("Thickness0",mThickness0);
        rSerializer.load("Thickness",mThickness);
        rSerializer.load("StrainsVector",mStrainsVector);
        rSerializer.load("StressesVector",mStressesVector);
        rSerializer.load("CauchyStressesVector",mCauchyStressesVector);
        rSerializer.load("Thickness0",mThickness0);
        rSerializer.load("V1",mV1);
        rSerializer.load("V2",mV2);
        rSerializer.load("G_Vector",mG_Vector);
    }

    ///@}

};	// class MembraneElement.

}	// namespace Kratos.

#endif // KRATOS_MEMBRANE_ELEMENT_H_INCLUDED  defined 
