//   
//   Project Name:        Kratos       
//   Last modified by:    $Author:  Miguel Masó Sotomayor $
//   Date:                $Date:             june 14 2017 $
//   Revision:            $Revision:                  1.1 $
//
//

#if !defined(KRATOS_NON_CONSERVATIVE_DC_ELEM_H_INCLUDED)
#define  KRATOS_NON_CONSERVATIVE_DC_ELEM_H_INCLUDED 

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

  class NonConservativeDC : public Element
  {
    public:
     
    /// Counted pointer of NonConservativeDC
    KRATOS_CLASS_POINTER_DEFINITION(NonConservativeDC);


    /// Default constructor.
    NonConservativeDC(IndexType NewId, GeometryType::Pointer pGeometry);
    NonConservativeDC(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~ NonConservativeDC();


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
    NonConservativeDC() : Element()
    {
    }
       
       
  }; // Class NonConservativeDC
}  // namespace Kratos.

#endif // KRATOS_NON_CONSERVATIVE_DC_ELEM_H_INCLUDED  defined
