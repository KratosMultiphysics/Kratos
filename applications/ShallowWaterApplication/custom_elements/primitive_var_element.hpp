//   
//   Project Name:        Kratos       
//   Last modified by:    Miguel Masó Sotomayor
//   Date:                June 28th 2017
//   Revision:            1.3
//
//

#if !defined(KRATOS_PRIMITIVE_VAR_ELEM_H_INCLUDED)
#define  KRATOS_PRIMITIVE_VAR_ELEM_H_INCLUDED 

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
  class PrimitiveVarElement : public Element
  {
  public:
     
    /// Counted pointer of PrimitiveVarElement
    KRATOS_CLASS_POINTER_DEFINITION( PrimitiveVarElement );

//----------------------------------------------------------------------

    /// Default constructor.
    PrimitiveVarElement() :
        Element()
    {}
    
    PrimitiveVarElement(IndexType NewId, GeometryType::Pointer pGeometry) :
        Element(NewId, pGeometry)
    {}

    PrimitiveVarElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
        Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~ PrimitiveVarElement() {};

//----------------------------------------------------------------------

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        KRATOS_TRY
        return Element::Pointer(new PrimitiveVarElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
        KRATOS_CATCH("")
    }


    int Check(const ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    void GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo);

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);

//----------------------------------------------------------------------

  protected:

    void CalculateGeometry(boost::numeric::ublas::bounded_matrix<double, TNumNodes, 2>& rDN_DX, double& rArea);
    
    double ComputeElemSize(const boost::numeric::ublas::bounded_matrix<double, TNumNodes, 2>& rDN_DX);
    
    void GetNodalValues(array_1d<double, TNumNodes*3>& rDepth, array_1d<double,TNumNodes*3>& rRain, array_1d<double,TNumNodes*3>& rUnkn, array_1d<double, TNumNodes*3>& rProj);
    
    void GetElementValues(const boost::numeric::ublas::bounded_matrix<double,TNumNodes, 2>& rDN_DX, const array_1d<double,TNumNodes*3>& rNodalVar, array_1d<double,2>& rVelocity, double& rHeight, array_1d<double,2>& rHeightGrad);
    
    void ComputeStabilizationParameters(const double& rHeight, const array_1d<double,2>& rHeightGrad, const double& rElemSize, double& rTauU, double& rTauH, double& rKdc);
    
    double mGravity;
    double mHeightUnitConvert;

//----------------------------------------------------------------------

  private:

    friend class Serializer;


  }; // Class PrimitiveVarElement

}  // namespace Kratos.

#endif // KRATOS_PRIMITIVE_VAR_ELEM_H_INCLUDED  defined
