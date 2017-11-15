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

    struct ElementVariables
    {
        double dt_inv;
        double lumping_factor;
        double dyn_tau;
        double gravity;
        double manning2;
        double height_units;

        //~ double height;
        //~ array_1d<double,2> velocity;
        //~ array_1d<double,2> height_grad;
        double scalar;
        array_1d<double,2> vector;
        array_1d<double,2> scalar_grad;
        boost::numeric::ublas::bounded_matrix<double,2,2> vector_grad;
        double vector_div;
        
        array_1d<double, TNumNodes*3> depth;
        array_1d<double, TNumNodes*3> rain;
        array_1d<double, TNumNodes*3> unknown;
        array_1d<double, TNumNodes*3> proj_unk;
    };

    void InitializeElement(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo);

    void CalculateGeometry(boost::numeric::ublas::bounded_matrix<double, TNumNodes, 2>& rDN_DX, double& rArea);
    
    double ComputeElemSize(const boost::numeric::ublas::bounded_matrix<double, TNumNodes, 2>& rDN_DX);
    
    void GetNodalValues(ElementVariables& rVariables);
    
    void GetElementValues(const boost::numeric::ublas::bounded_matrix<double,TNumNodes, 2>& rDN_DX, ElementVariables& rVariables);
    
    void ComputeStabilizationParameters(const ElementVariables& rVariables,
                                        const double& rElemSize,
                                        double& rTauU,
                                        double& rTauH,
                                        double& rKdc);
    
    void ComputeAuxMatrices(
            const boost::numeric::ublas::bounded_matrix<double,TNumNodes, TNumNodes>& rNcontainer,
            const boost::numeric::ublas::bounded_matrix<double,TNumNodes,2>& rDN_DX,
            const ElementVariables& rVariables,
            boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rMassMatrixScalar,
            boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rMassMatrixVector,
            boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rScalarGrad,
            boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rVectorDiv,
            boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rScalarDiff,
            boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rVectorDiff );

//----------------------------------------------------------------------

  private:

    friend class Serializer;


  }; // Class PrimitiveVarElement

}  // namespace Kratos.

#endif // KRATOS_PRIMITIVE_VAR_ELEM_H_INCLUDED  defined
