//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#if !defined(KRATOS_EULER_PRIMITIVE_VAR_ELEM_H_INCLUDED)
#define  KRATOS_EULER_PRIMITIVE_VAR_ELEM_H_INCLUDED 

// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h" 
#include "includes/serializer.h"

namespace Kratos
{

  template< unsigned int TNumNodes >
  class EulerPrimVarElement : public Element
  {
  public:
     
    /// Counted pointer of EulerPrimVarElement
    KRATOS_CLASS_POINTER_DEFINITION( EulerPrimVarElement );

//----------------------------------------------------------------------

    /// Default constructor.
    EulerPrimVarElement()
	: Element()
	{}
	
    EulerPrimVarElement(IndexType NewId, GeometryType::Pointer pGeometry)
	: Element(NewId, pGeometry)
	{}

	EulerPrimVarElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
	: Element(NewId, pGeometry, pProperties)
	{}

    /// Destructor.
    virtual ~ EulerPrimVarElement() {};

//----------------------------------------------------------------------

	Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
	{
		KRATOS_TRY
		return Element::Pointer(new EulerPrimVarElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
		KRATOS_CATCH("")
	}


    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    void GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo);

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);

//----------------------------------------------------------------------

  protected:

    void CalculateGeometry(boost::numeric::ublas::bounded_matrix<double, TNumNodes, 2>& rDN_DX, double& rArea);
    
    double ComputeElemSize(boost::numeric::ublas::bounded_matrix<double, TNumNodes, 2>& rDN_DX);
    
    void GetNodalValues(array_1d<double, TNumNodes*3>& rdepth, array_1d<double, TNumNodes*3>& runkn, array_1d<double, TNumNodes*3>& rprev);
    
    void GetElementValues(boost::numeric::ublas::bounded_matrix<double, TNumNodes, 2>& rDN_DX, array_1d<double, TNumNodes*3>& r_nodal_var, double& rheight, boost::numeric::ublas::bounded_matrix<double,1, 2>& rvel);
    
//----------------------------------------------------------------------

  private:

    friend class Serializer;


  }; // Class EulerPrimVarElement

}  // namespace Kratos.

#endif // KRATOS_EULER_PRIMITIVE_VAR_ELEM_H_INCLUDED  defined
