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

#if !defined(KRATOS_EULER_CONSERVED_VAR_ELEM_H_INCLUDED)
#define  KRATOS_EULER_CONSERVED_VAR_ELEM_H_INCLUDED 

// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h" 

namespace Kratos
{

  template< unsigned int TNumNodes >
  class EulerConsVarElement : public Element
  {
  public:
     
    /// Counted pointer of EulerConsVarElement
    KRATOS_CLASS_POINTER_DEFINITION( EulerConsVarElement );

//----------------------------------------------------------------------

    /// Default constructor.
    EulerConsVarElement()
    : Element()
    {}
    
    EulerConsVarElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {}
    
    EulerConsVarElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~ EulerConsVarElement() {};

//----------------------------------------------------------------------

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
    {
        KRATOS_TRY
        return Element::Pointer(new EulerConsVarElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
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
    
    void GetNodalValues(array_1d<double, TNumNodes*3>& rdepth, array_1d<double, TNumNodes*3>& runkn, array_1d<double, TNumNodes*3>& rproj);
    
    void GetElementValues(boost::numeric::ublas::bounded_matrix<double, TNumNodes, 2>& rDN_DX, array_1d<double, TNumNodes*3>& rNodalVar, array_1d<double,2>& rMom, double& rHeight, array_1d<double,2>& rVel, boost::numeric::ublas::bounded_matrix<double,2,2>& rDivU);

    void CalculateLumpedMassMatrix(boost::numeric::ublas::bounded_matrix<double, TNumNodes*3, TNumNodes*3>& rMassMatrix);

  private:

    friend class Serializer;


  }; // Class EulerConsVarElement

}  // namespace Kratos.

#endif // KRATOS_EULER_CONSERVED_VAR_ELEM_H_INCLUDED  defined
