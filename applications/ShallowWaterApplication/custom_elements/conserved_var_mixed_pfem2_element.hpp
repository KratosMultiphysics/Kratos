//   
//   Project Name:        Kratos       
//   Last modified by:    Miguel Masó Sotomayor
//   Date:                June 28th 2017
//   Revision:            1.3
//
//

#if !defined(KRATOS_CONSERVED_VAR_MIXED_PFEM2_ELEM_H_INCLUDED)
#define  KRATOS_CONSERVED_VAR_MIXED_PFEM2_ELEM_H_INCLUDED 

// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h" 
#include "includes/serializer.h"
#include "custom_elements/primitive_var_element.hpp"

namespace Kratos
{

  template< unsigned int TNumNodes >
  class ConservedVarMixedPfem2Element : public PrimitiveVarElement<TNumNodes>
  {
  public:
     
    /// Counted pointer of ConservedVarMixedPfem2Element
    KRATOS_CLASS_POINTER_DEFINITION( ConservedVarMixedPfem2Element );

    typedef PrimitiveVarElement<TNumNodes>                             BaseType;

    typedef typename BaseType::IndexType                              IndexType;

    typedef typename BaseType::ElementVariables                ElementVariables;

    typedef typename BaseType::EquationIdVectorType        EquationIdVectorType;

    typedef typename BaseType::DofsVectorType                    DofsVectorType;

    typedef typename BaseType::VectorType                            VectorType;

    typedef typename BaseType::MatrixType                            MatrixType;

    typedef typename BaseType::GeometryType                        GeometryType;

    typedef typename BaseType::GeometryType::Pointer        GeometryPointerType;

    typedef typename BaseType::NodesArrayType                    NodesArrayType;

    typedef typename BaseType::PropertiesType                    PropertiesType;

    typedef typename BaseType::PropertiesType::Pointer    PropertiesPointerType;

//----------------------------------------------------------------------

    /// Default constructor.
    ConservedVarMixedPfem2Element() :
        BaseType()
    {}
    
    ConservedVarMixedPfem2Element(IndexType NewId, GeometryPointerType pGeometry) :
        BaseType(NewId, pGeometry)
    {}

    ConservedVarMixedPfem2Element(IndexType NewId, GeometryPointerType pGeometry, PropertiesPointerType pProperties) :
        BaseType(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~ ConservedVarMixedPfem2Element() override {};

//----------------------------------------------------------------------

    Element::Pointer Create(IndexType NewId, NodesArrayType const& rThisNodes, PropertiesPointerType pProperties) const override
    {
        KRATOS_TRY
        return boost::make_shared< ConservedVarMixedPfem2Element < TNumNodes > >(NewId, this->GetGeometry().Create(rThisNodes), pProperties);
        KRATOS_CATCH("")
    }

    Element::Pointer Create(IndexType NewId, GeometryPointerType pGeom, PropertiesPointerType pProperties) const override
    {
        KRATOS_TRY
        return boost::make_shared< ConservedVarMixedPfem2Element < TNumNodes > >(NewId, pGeom, pProperties);
        KRATOS_CATCH("")
    }

    //~ Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesPointerType pProperties) const override
    //~ {
        //~ KRATOS_TRY
        //~ return Element::Pointer(new ConservedVarMixedPfem2Element<TNumNodes>(NewId, GetGeometry().Create(ThisNodes), pProperties));
        //~ KRATOS_CATCH("")
    //~ }

//----------------------------------------------------------------------

    int Check( const ProcessInfo& rCurrentProcessInfo ) override;

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

//----------------------------------------------------------------------

  protected:

    void GetNodalValues(ElementVariables& rVariables);

    void GetElementValues(const boost::numeric::ublas::bounded_matrix<double,TNumNodes, 2>& rDN_DX, ElementVariables& rVariables);
    
    void ComputeAuxMatrices(
            const boost::numeric::ublas::bounded_matrix<double,TNumNodes, TNumNodes>& rNcontainer,
            const boost::numeric::ublas::bounded_matrix<double,TNumNodes,2>& rDN_DX,
            const ElementVariables& rVariables,
            boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rMassMatrixScalar,
            boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rMassMatrixVector,
            boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rScalarGrad,
            boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rVectorDiv,
            boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rNonLinearTerm,
            boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rScalarDiff,
            boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rVectorDiff );


//----------------------------------------------------------------------

  private:

    friend class Serializer;


  }; // Class ConservedVarMixedPfem2Element

}  // namespace Kratos.

#endif // KRATOS_CONSERVED_VAR_MIXED_PFEM2_ELEM_H_INCLUDED  defined
