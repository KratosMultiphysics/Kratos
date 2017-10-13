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
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;
    ///@}
    ///@name Access
    ///@{
    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{
    ///@}
    ///@name Friends
    ///@{
    ///@}

protected:
	// A protected default constructor necessary for serialization
	// As this class is base class the constructor has to be protected.
	MeshlessBaseElement() : Element()
	{
	}
    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
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
		const Matrix DN_De);

	void ComputeGlobalDerivatives(
		const Matrix& DN_De,
		Matrix& DN_DX,
		Matrix& Jacobian);

	void CrossProduct(
		array_1d<double, 3>& cross,
		const array_1d<double, 3>& a,
		const array_1d<double, 3>& b);

	array_1d<double, 3> CrossProduct(
		const array_1d<double, 3>& a,
		const array_1d<double, 3>& b);
    ///@}
    ///@name Protected  Access
    ///@{
    ///@}
    ///@name Protected Inquiry
    ///@{
    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

private:
    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{
    ///@}
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
    ///@name Private Operators
    ///@{
    ///@}
    ///@name Private Operations
    ///@{
    ///@}
    ///@name Private  Access
    ///@{
    ///@}
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //MeshlessBaseElement& operator=(const MeshlessBaseElement& rOther);

    /// Copy constructor.
    //MeshlessBaseElement(const MeshlessBaseElement& rOther);

    ///@}

}; // Class MeshlessBaseElement

///@}
///@name Type Definitions
///@{MeshlessBaseElement
///@}
///@name Input and output
///@{
/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    MeshlessBaseElement& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const MeshlessBaseElement& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_MESHLESS_BASE_ELEMENT_H_INCLUDED  defined 


