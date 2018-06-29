#if !defined(KRATOS_TRUSS_DISCRETE_ELEMENT_H_INCLUDED )
#define  KRATOS_TRUSS_DISCRETE_ELEMENT_H_INCLUDED


// System includes
#include "includes/define.h"
#include "includes/element.h"
//#include "includes/ublas_interface.h"
#include "includes/variables.h"

// External includes
//#include "boost/smart_ptr.hpp"

// Project includes
#include "custom_elements/meshless_base_element.h"
#include "custom_elements/curve_base_discrete_element.h"

namespace Kratos
{
///@name Kratos Classes
///@{
/// Short class definition.
/** Truss element.
*/
class TrussDiscreteElement
    : public CurveBaseDiscreteElement
{
public:
    ///@name Type Definitions
    ///@{
    /// Counted pointer of TrussDiscreteElement
    KRATOS_CLASS_POINTER_DEFINITION(TrussDiscreteElement);
    ///@}
    ///@name Life Cycle
    ///@{
    /// Default constructor.
	// Constructor using an array of nodes
	TrussDiscreteElement(IndexType NewId, GeometryType::Pointer pGeometry)
		: CurveBaseDiscreteElement(NewId, pGeometry)
	{};
	// Constructor using an array of nodes with properties
	TrussDiscreteElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
		: CurveBaseDiscreteElement(NewId, pGeometry, pProperties)
	{};

	// default constructor necessary for serialization
	TrussDiscreteElement() : CurveBaseDiscreteElement() {};

	/// Destructor.
	virtual ~TrussDiscreteElement() override
	{};

	Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
	{
		return Kratos::make_shared< TrussDiscreteElement >(NewId, GetGeometry().Create(ThisNodes), pProperties);
	};

    ///@}
    ///@name Operations
	///@{

	/**
	* This functions calculates both the RHS and the LHS
	* @param rLeftHandSideMatrix: The LHS
	* @param rRightHandSideVector: The RHS
	* @param rCurrentProcessInfo: The current process info instance
	* @param CalculateStiffnessMatrixFlag: The flag to set if compute the LHS
	* @param CalculateResidualVectorFlag: The flag to set if compute the RHS
	*/
	void CalculateAll(
		MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo,
		const bool CalculateStiffnessMatrixFlag,
		const bool CalculateResidualVectorFlag
	) override;

	/**
	* Calculate a double Variable on the Element Constitutive Law
	* @param rVariable: The variable we want to get
	* @param rOutput: The values obtained int the integration points
	* @param rCurrentProcessInfo: the current process info instance
	*/
	void CalculateOnIntegrationPoints(
		const Variable<double>& rVariable,
		std::vector<double>& rOutput,
		const ProcessInfo& rCurrentProcessInfo
	) override;

	/**
	* Calculate a Vector Variable on the Element Constitutive Law
	* @param rVariable: The variable we want to get
	* @param rOutput: The values obtained int the integration points
	* @param rCurrentProcessInfo: the current process info instance
	*/
	void CalculateOnIntegrationPoints(
		const Variable<Vector>& rVariable,
		std::vector<Vector>& rValues,
		const ProcessInfo& rCurrentProcessInfo
	) override;

	/**
	* Sets on rResult the ID's of the element degrees of freedom
	* @param rResult: The vector containing the equation id
	* @param rCurrentProcessInfo: The current process info instance
	*/
	void EquationIdVector(
		EquationIdVectorType& rResult,
		ProcessInfo& rCurrentProcessInfo
	) override;

	/**
	* Sets on rElementalDofList the degrees of freedom of the considered element geometry
	* @param rElementalDofList: The vector containing the dof of the element
	* @param rCurrentProcessInfo: The current process info instance
	*/
	void GetDofList(
		DofsVectorType& rElementalDofList,
		ProcessInfo& rCurrentProcessInfo
	) override;

	///@}

protected:
	/**
	* It initializes the material
	*/
	virtual void InitializeMaterial();
	/**
	* Called at the end of eahc solution step
	* @param rCurrentProcessInfo: the current process info instance
	*/
	void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

private:
	///@name Static Member Variables
	///@{
	///@}
	///@name Operations
	///@{


	///@}
	///@name Serialization
	///@{

	friend class Serializer;

	void save(Serializer& rSerializer) const override
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
	}

	void load(Serializer& rSerializer) override
	{
		KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
	}

	///@}

};	 // Class TrussDiscreteElement
///@}

}  // namespace Kratos.

#endif // KRATOS_MESHLESS_TRUSS_DISCRETE_ELEMENT_H_INCLUDED  defined 