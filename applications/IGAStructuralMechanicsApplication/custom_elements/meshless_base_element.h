#if !defined(KRATOS_MESHLESS_BASE_ELEMENT_H_INCLUDED )
#define  KRATOS_MESHLESS_BASE_ELEMENT_H_INCLUDED


// System includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
// External includes
#include "boost/smart_ptr.hpp"
// Project includes


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
class MeshlessBaseElement
    : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of MeshlessBaseElement
    KRATOS_CLASS_POINTER_DEFINITION(MeshlessBaseElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
	MeshlessBaseElement(IndexType NewId, GeometryType::Pointer pGeometry);
	MeshlessBaseElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~MeshlessBaseElement();
    ///@}
    ///@name Operations
    ///@{
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;
    ///@}

protected:
	// A protected default constructor necessary for serialization
	// As this class is base class the constructor has to be protected.
	MeshlessBaseElement() : Element()
	{
	};
    ///@name Protected Operations
    ///@{
	void GetGeometryData(
		double& integration_weight,
		Vector& N,
		Matrix& DN_De);

	void Jacobian(
		const Matrix& DN_De,
		Matrix& Jacobian);

	void Hessian(Matrix& Hessian,
		const Matrix& DDN_DDe);

	void ComputeGlobalDerivatives(
		const Matrix& DN_De,
		Matrix& DN_DX,
		Matrix& Jacobian);

	void CrossProduct(
		array_1d<double, 3>& cross,
		const array_1d<double, 3>& a,
		const array_1d<double, 3>& b);

  void CrossProduct2(
    Vector& cross,
    const Vector& a,
    const Vector& b);

	array_1d<double, 3> CrossProduct(
		const array_1d<double, 3>& a,
		const array_1d<double, 3>& b);
    ///@}

private:
    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }
    ///@}
}; // Class MeshlessBaseElement
}  // namespace Kratos.

#endif // KRATOS_MESHLESS_BASE_ELEMENT_H_INCLUDED  defined 


