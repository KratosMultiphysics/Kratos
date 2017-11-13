#if !defined(KRATOS_MESHLESS_SHELL_ELEMENT_H_INCLUDED )
#define  KRATOS_MESHLESS_SHELL_ELEMENT_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

#include "custom_elements/meshless_base_element.h"

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

/// Short class definition.
/** Detail class definition.
*/
class MeshlessShellElement
    : public MeshlessBaseElement
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of MeshlessShellElement
    KRATOS_CLASS_POINTER_DEFINITION(MeshlessShellElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
	MeshlessShellElement(IndexType NewId, GeometryType::Pointer pGeometry);
	MeshlessShellElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~MeshlessShellElement();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations

	Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;
	

	void EquationIdVector(
		EquationIdVectorType& rResult,
		ProcessInfo& rCurrentProcessInfo) override;

	void GetDofList(
		DofsVectorType& ElementalDofList,
		ProcessInfo& rCurrentProcessInfo) override;

	void Initialize();

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

	void FinalizeSolutionStep(
		ProcessInfo& rCurrentProcessInfo) override;

	//void GetValuesVector(
	//	Vector& values,
	//	int Step = 0);

	//void GetFirstDerivativesVector(
	//	Vector& values,
	//	int Step = 0);

	//void GetSecondDerivativesVector(
	//	Vector& values,
	//	int Step = 0);

	void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
		std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo) override;



protected:


private:
	///@name Static Member Variables


	ConstitutiveLaw::Pointer mConstitutiveLawVector;
	Geometry< Point >::Pointer  mpReferenceGeometry;

	double mDetJ0;

	double mTotalDomainInitialSize;
	double mdensity;
	double mThickness0; //thickness in the reference configuration


	double mThickness;									//container of thickness
	array_1d<double, 3> mStrainsVector;	//container of Strain
	array_1d<double, 6> mStressesVector;	//container of Stress
	array_1d<double, 6> mCauchyStressesVector;	//container of Stress


	array_1d<double, 3> mV1;
	array_1d<double, 3> mV2;
	Matrix              mG_Vector;

	array_1d<double, 3> mGab0;
	array_1d<double, 3> mCurvature0;


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


	//void ClearNodalForces();


	void CalculateQ(
		boost::numeric::ublas::bounded_matrix<double, 3, 3>& msQ,
		Matrix& msG);

	void CalculateMetricDeformed(Matrix DN_De,
		array_1d<double, 3>& gab,
		array_1d<double, 3>& curvature_coefficient,
		array_1d<double, 3>& g1,
		array_1d<double, 3>& g2);

	void CalculateBMembrane(
		Matrix& B,
		boost::numeric::ublas::bounded_matrix<double, 3, 3>& Q,
		const Matrix& DN_De,
		const array_1d<double, 3>& g1,
		const array_1d<double, 3>& g2);

	void CalculateBCurvature(
		Matrix& B,
		boost::numeric::ublas::bounded_matrix<double, 3, 3>& Q,
		const Matrix& DN_De,
		const array_1d<double, 3>& g1,
		const array_1d<double, 3>& g2);

	void CalculateSecondVariationStrainCurvature(Matrix DN_De,
		Matrix& Strain_curvature11,
		Matrix& Strain_curvature22,
		Matrix& Strain_curvature12,
		Matrix& Curvature_curvature11,
		Matrix& Curvature_curvature22,
		Matrix& Curvature_curvature12,
		boost::numeric::ublas::bounded_matrix<double, 3, 3>& Q,
		array_1d<double, 3>& g1,
		array_1d<double, 3>& g2);

	void CalculateStrain(
		Vector& StrainVector,
		array_1d<double, 3>& gab,
		array_1d<double, 3>& gab_ref);

	void CalculateCurvature(
		Vector& CurvatureVector,
		array_1d<double, 3>& bv,
		array_1d<double, 3>& bv_ref);

	int  Check(const ProcessInfo& rCurrentProcessInfo);

	///@}
	///@name Serialization
	///@{

	friend class Serializer;

	// A private default constructor necessary for serialization
	MeshlessShellElement() {}

	void save(Serializer& rSerializer) const override
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
		rSerializer.save("Gab0", mGab0);
		rSerializer.save("Curvature0", mCurvature0);
	}

	void load(Serializer& rSerializer) override
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
		rSerializer.load("Gab0", mGab0);
		rSerializer.load("Curvature0", mCurvature0);
	}

	///@}

};	 // Class MeshlessShellElement

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    MeshlessShellElement& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const MeshlessShellElement& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_MESHLESS_MESHLESS_SHELL_ELEMENT_H_INCLUDED  defined 