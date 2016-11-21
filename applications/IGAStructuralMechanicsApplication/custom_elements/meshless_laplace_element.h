#if !defined(KRATOS_MESHLESS_LAPLACE_ELEMENT_H_INCLUDED )
#define  KRATOS_MESHLESS_LAPLACE_ELEMENT_H_INCLUDED



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
class MeshlessLaplaceElement
    : public MeshlessBaseElement
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of MembraneElement
    KRATOS_CLASS_POINTER_DEFINITION(MeshlessLaplaceElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
	MeshlessLaplaceElement(IndexType NewId, GeometryType::Pointer pGeometry);
	MeshlessLaplaceElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~MeshlessLaplaceElement();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

	Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;
	
	void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

	void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) override;

	void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

	void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

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

	    // A private default constructor necessar	r serializa	
	MeshlessLaplaceElement() : MeshlessBaseElement()
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

	/*void GetGeometryData(
		double& integration_weight,
		Vector& N,
		Matrix& DN_De
		);

	void ComputeGlobalDerivatives(
		const Matrix& DN_De,
		Matrix& DN_DX,
		Matrix& Jacobian);*/


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



    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    virtual void load(Serializer& rSerializer)
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
    //MembraneElement& operator=(const MembraneElement& rOther);

    /// Copy constructor.
    //MembraneElement(const MembraneElement& rOther);


    ///@}

}; // Class MeshlessLaplaceElement

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    MeshlessLaplaceElement& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const MeshlessLaplaceElement& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_MESHLESS_LAPLACE_ELEMENT_H_INCLUDED  defined 