#if !defined(KRATOS_MESHLESS_BASE_SURFACE_ELEMENT_H_INCLUDED )
#define  KRATOS_MESHLESS_BASE_SURFACE_ELEMENT_H_INCLUDED


// System includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "custom_elements/meshless_base_element.h"

namespace Kratos
{
///@name Kratos Classes
///@{
/// Short class definition.
/** Meshless shell element deals as base class for thin walled structures.
KRATOS_API(IGA_STRUCTURAL_MECHANICS_APPLICATION)
*/
class  MeshlessBaseSurfaceElement
    : public MeshlessBaseElement
{
protected:
	/**
	* Internal variables used for metric transformation
	*/
	struct MetricVariables
	{
		Vector gab; // covariant metric
		Vector gab_con; // contravariant metric
		Vector curvature; //
		Matrix J; //Jacobian
		double  detJ;
		Vector g1; //base vector 1
		Vector g2; //base vector 2
		Vector g3; //base vector 3
		double dA; //differential area
		Matrix H; //Hessian
		Matrix Q; //Transformation matrix Q from contravariant to cartesian basis
		Matrix T; //Transformation matrix T from contravariant to local cartesian basis

		/**
		* The default constructor
		* @param Dimension: The size of working space dimension
		*/
		MetricVariables(const unsigned int& Dimension)
		{
			gab = ZeroVector(Dimension);
			gab_con = ZeroVector(Dimension);

			curvature = ZeroVector(Dimension);

			J = ZeroMatrix(Dimension, Dimension);
			detJ = 1.0;

			g1 = ZeroVector(Dimension);
			g2 = ZeroVector(Dimension);
			g3 = ZeroVector(Dimension);

			dA = 1.0;

			Matrix H = ZeroMatrix(3, 3);
			Matrix Q = ZeroMatrix(3, 3);
			Matrix T = ZeroMatrix(3, 3);
		}
	};

	/**
	* Internal variables used in the constitutive equations
	*/
	struct ConstitutiveVariables
	{
		Vector StrainVector;
		Vector StrainCurvatureVector;
		Vector StressVector;
		Vector StressCurvatureVector;
		Matrix DMembrane;
		Matrix DCurvature;

		/**
		* The default constructor
		* @param StrainSize: The size of the strain vector in Voigt notation
		*/
		ConstitutiveVariables(const unsigned int& StrainSize)
		{
			StrainVector = ZeroVector(StrainSize);
			StrainCurvatureVector = ZeroVector(StrainSize);
			StressVector = ZeroVector(StrainSize);
			StressCurvatureVector = ZeroVector(StrainSize);
			DMembrane = ZeroMatrix(StrainSize, StrainSize);
			DCurvature = ZeroMatrix(StrainSize, StrainSize);
		}
	};

public:
    ///@name Type Definitions
    ///@{
    /// Counted pointer of MeshlessBaseSurfaceElement
    KRATOS_CLASS_POINTER_DEFINITION(MeshlessBaseSurfaceElement);
    ///@}
    ///@name Life Cycle
    ///@{
    /// Default constructor.
	 // Constructor using an array of nodes
	MeshlessBaseSurfaceElement(IndexType NewId, GeometryType::Pointer pGeometry)
		: MeshlessBaseElement(NewId, pGeometry)
	{};
	 // Constructor using an array of nodes with properties
	MeshlessBaseSurfaceElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: MeshlessBaseElement(NewId, pGeometry, pProperties)
	{};

	MeshlessBaseSurfaceElement() : MeshlessBaseElement()
	{};

    /// Destructor.
	virtual ~MeshlessBaseSurfaceElement() override
	{};

	Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
	{
		return boost::make_shared< MeshlessBaseSurfaceElement >(NewId, GetGeometry().Create(ThisNodes), pProperties);
	};

    ///@}
    ///@name Operations
	///@{

	/**
	* Called to initialize the element.
	* Must be called before any calculation is done
	*/
	void Initialize() override;


	///@}
protected:

	ConstitutiveLaw::Pointer mConstitutiveLaw;
	MetricVariables mInitialMetric = MetricVariables(3);

	///@name Operations
	///@{
	/**
	* It initializes the material
	*/
	virtual void InitializeMaterial();

	/**
	* Called at the end of eahc solution step
	* @param rCurrentProcessInfo: the current process info instance
	*/
	void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

	/**
	* This functions calculates both the RHS and the LHS
	* @param rLeftHandSideMatrix: The LHS
	* @param rRightHandSideVector: The RHS
	* @param rCurrentProcessInfo: The current process info instance
	* @param CalculateStiffnessMatrixFlag: The flag to set if compute the LHS
	* @param CalculateResidualVectorFlag: The flag to set if compute the RHS
	*/
	virtual void CalculateAll(
		MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo,
		const bool CalculateStiffnessMatrixFlag,
		const bool CalculateResidualVectorFlag
	);

	/**
	* This function provides a more general interface to the element.
	* It is designed so that rLHSvariables and rRHSvariables are passed to the element thus telling what is the desired output
	* @param rLeftHandSideMatrices: container with the output left hand side matrices
	* @param rLHSVariables: paramter describing the expected LHSs
	* @param rRightHandSideVectors: container for the desired RHS output
	* @param rRHSVariables: parameter describing the expected RHSs
	*/
	void CalculateLocalSystem(
		MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo ) override;

	/**
	* This is called during the assembling process in order to calculate the elemental right hand side vector only
	* @param rRightHandSideVector: the elemental right hand side vector
	* @param rCurrentProcessInfo: the current process info instance
	*/
	void CalculateRightHandSide(
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo ) override;

	/**
	* Calculation of the Material Stiffness Matrix. Km = B^T * D *B
	*/
	void CalculateAndAddKm(
		MatrixType& rLeftHandSideMatrix,
		const Matrix& B,
		const Matrix& D,
		const double IntegrationWeight );

	void CalculateAndAddNonlinearKm(
		Matrix& rLeftHandSideMatrix,
		const Matrix& B11,
		const Matrix& B22,
		const Matrix& B12,
		const Vector& SD,
		double IntegrationWeight);

	/**
	* Is called to compute the respective metric depending on the deformation.
	* @param metric: the current metric
	*/
	virtual void CalculateMetric(
		MetricVariables& metric);

	/**
	* This functions updates the constitutive variables
	* @param rActual
	Metric: The actual metric
	* @param rThisConstitutiveVariables: The constitutive variables to be calculated
	* @param rValues: The CL parameters
	* @param ThisStressMeasure: The stress measure considered
	*/
	virtual void CalculateConstitutiveVariables(
		MetricVariables& rActualMetric,
		ConstitutiveVariables& rThisConstitutiveVariables,
		ConstitutiveLaw::Parameters& rValues,
		const ConstitutiveLaw::StressMeasure ThisStressMeasure
	);

	void CalculateStrain(
		Vector& StrainVector,
		Vector& gab,
		Vector& gab0);

	void CalculateCurvature(
		Vector& CurvatureVector,
		Vector& bv,
		Vector& bv_ref);

	void CalculateBMembrane(
		Matrix& rB,
		const MetricVariables& metric);

	void CalculateBCurvature(
		Matrix& rB,
		const MetricVariables& metric);

	void CalculateSecondVariationStrainCurvature(
		Matrix& Strain_in_Q_coordinates11,
		Matrix& Strain_in_Q_coordinates22,
		Matrix& Strain_in_Q_coordinates12,
		Matrix& Curvature_in_Q_coordinates11,
		Matrix& Curvature_in_Q_coordinates22,
		Matrix& Curvature_in_Q_coordinates12,
		const MetricVariables& rMetric);
	///@}

private:
	///@name Operations
	///@{

	///@}

	///@name Static Member Variables

	///@}
	///@name Serialization
	///@{
	friend class Serializer;

	void save(Serializer& rSerializer) const override
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
		rSerializer.save("ConstitutiveLaw", mConstitutiveLaw);
	}
	void load(Serializer& rSerializer) override
	{
		KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
		rSerializer.load("ConstitutiveLaw", mConstitutiveLaw);
	}
	///@}
};	 // Class MeshlessBaseSurfaceElement
///@}
}  // namespace Kratos.

#endif // KRATOS_MESHLESS_BASE_SURFACE_ELEMENT_H_INCLUDED  defined 