#if !defined(KRATOS_CURVE_BASE_DISCRETE_ELEMENT_H_INCLUDED )
#define  KRATOS_CURVE_BASE_DISCRETE_ELEMENT_H_INCLUDED


// System includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

// External includes
//#include "boost/smart_ptr.hpp"

// Project includes
#include "custom_elements/meshless_base_element.h"

namespace Kratos
{
///@name Kratos Classes
///@{
/// Short class definition.
/** CurveBaseDiscreteElement deals as base class for curve structured element formulations.
*/
class  CurveBaseDiscreteElement
    : public MeshlessBaseElement
{
protected:

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "CurveBaseDiscreteElement #" << Id();
        return buffer.str();
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CurveBaseDiscreteElement #" << Id();
    }

public:
    ///@name Type Definitions
    ///@{
    /// Counted pointer of CurveBaseDiscreteElement
    KRATOS_CLASS_POINTER_DEFINITION(CurveBaseDiscreteElement);
    ///@}
    ///@name Life Cycle
    ///@{
    /// Default constructor.
	 // Constructor using an array of nodes
	CurveBaseDiscreteElement(IndexType NewId, GeometryType::Pointer pGeometry)
		: MeshlessBaseElement(NewId, pGeometry)
	{};
	 // Constructor using an array of nodes with properties
	CurveBaseDiscreteElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: MeshlessBaseElement(NewId, pGeometry, pProperties)
	{};

	CurveBaseDiscreteElement() : MeshlessBaseElement()
	{};

    /// Destructor.
	virtual ~CurveBaseDiscreteElement() override
	{};

	Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
	{
		KRATOS_ERROR << "Trying to create a \"CurveBaseDiscreteElement\"" << std::endl;
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
	///@name Static Member Variables
	///@{
	ConstitutiveLaw::Pointer mConstitutiveLaw;
	Vector mBaseVector0;
	///@}
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
	virtual void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

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
	///@}

	void GetBaseVector(
		Vector& rBaseVector, 
		const Matrix& rDN_De);

	/**
	* MappingGeometricToParameter calculates the J-tilde for the mapping from
	Geometry to Parameter Space. This paramater is needed for all condition
	integrations on edges.
	*
	* @param DN_De derivatives of shape functions.
	* @param JGeometricToParameter Mapping parameter for Geometric Space to
	Parameter Space
	*
	* @see Jacobian
	*/
	void GetBoundaryEdgeBaseVector(const Matrix& DN_De,
		const array_1d<double, 2>& Tangents,
		Vector& rBaseVector);

	void Get1stVariationsAxialStrain(
		Vector& rEpsilon1stVariationDoF,
		const Vector& rBaseVector,
		const int& rNumberOfDoFs, 
		const Matrix& rDN_De);

	void Get2ndVariationsAxialStrain(
		Matrix& rEpsilon2ndVariationDoF,
		const int& rNumberOfDoFs, 
		const Matrix& rDN_De);

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
};	 // Class CurveBaseDiscreteElement
///@}
}  // namespace Kratos.

#endif // KRATOS_CURVE_BASE_DISCRETE_ELEMENT_H_INCLUDED  defined