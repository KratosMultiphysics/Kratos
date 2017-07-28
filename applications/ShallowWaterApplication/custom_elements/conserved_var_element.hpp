//   
//   Project Name:        Kratos       
//   Last modified by:    Miguel Masó Sotomayor
//   Date:                July 3rd 2017
//   Revision:            1.2
//
//

#if !defined(KRATOS_CONSERVED_VAR_ELEM_H_INCLUDED)
#define  KRATOS_CONSERVED_VAR_ELEM_H_INCLUDED 

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
  class ConservedVarElement : public Element
  {
    public:
     
    /// Counted pointer of ConservedVarElement
    KRATOS_CLASS_POINTER_DEFINITION( ConservedVarElement );


    /// Default constructor.
    ConservedVarElement(IndexType NewId, GeometryType::Pointer pGeometry);
    ConservedVarElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~ ConservedVarElement();


    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    void GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo);

    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);



    protected:

    void CalculateConsistentMassMatrix(boost::numeric::ublas::bounded_matrix<double,9,9>& rMassMatrix);

    void CalculateLumpedMassMatrix(boost::numeric::ublas::bounded_matrix<double,9,9>& rMassMatrix);

    private:

    friend class Serializer;

    // A private default constructor necessary for serialization
    ConservedVarElement() : Element()
    {
    }
       
       
  }; // Class ConservedVarElement
}  // namespace Kratos.

#endif // KRATOS_CONSERVED_VAR_ELEM_H_INCLUDED  defined
