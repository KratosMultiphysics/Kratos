//
//   Project Name:        Kratos
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_FRACTIONAL_STEP_PFEM_2_3D_ELEM_H_INCLUDED)
#define  KRATOS_FRACTIONAL_STEP_PFEM_2_3D_ELEM_H_INCLUDED

// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"


namespace Kratos
{

  class FractionalStepPFEM23D
	  : public Element
   {
   public:

     /// Counted pointer of PFEM22D
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(FractionalStepPFEM23D);
    ///base type: an IndexedObject that automatically has a unique number
    ///typedef IndexedObject BaseType;
    ///Element from which it is derived
    ///typedef VMS<TDim, TNumNodes> ElementBaseType;
    ///definition of node type (default is: Node<3>)

    //typedef Node < 3 > NodeType;
    /**
     * Properties are used to store any parameters
     * related to the constitutive law
     */

    typedef Properties PropertiesType;
    ///definition of the geometry type with given NodeType
    typedef Geometry<NodeType> GeometryType;
    ///definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    //typedef typename ElementBaseType::MatrixType MatrixType;
    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef std::vector<std::size_t> EquationIdVectorType;
    typedef std::vector< Dof<double>::Pointer > DofsVectorType;
    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    /// Default constructor.
     FractionalStepPFEM23D(IndexType NewId, GeometryType::Pointer pGeometry);
     FractionalStepPFEM23D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

     /// Destructor.
     virtual ~ FractionalStepPFEM23D() override;


     Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;

     void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

     void AddExplicitContribution(const ProcessInfo& CurrentProcessInfo) override;

     void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

     void GetDofList(DofsVectorType& ElementalDofList,const ProcessInfo& CurrentProcessInfo) const override;

     void InitializeSolutionStep(const ProcessInfo& CurrentProcessInfo) override;

    virtual void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo) override;

    virtual void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
            std::vector<array_1d<double, 3 > >& rValues,
            const ProcessInfo& rCurrentProcessInfo) override;

   protected:

       	void CalculateLocalPressureSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo);

       	void PressureEquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const;

        void GetPressureDofList(DofsVectorType& ElementalDofList,const ProcessInfo& CurrentProcessInfo) const;

        void CalculateViscousRHS(const ProcessInfo& CurrentProcessInfo);

       	void CalculatePressureProjection(const ProcessInfo& CurrentProcessInfo);

        void AddViscousTerm(BoundedMatrix<double, (3-1)*6, (3-1)*6 >& rDampMatrix,
                         BoundedMatrix<double, (3+1), 3 >& rShapeDeriv,
                         const double Weight);



template<class T>
bool InvertMatrix(const T& input, T& inverse)  ;


   private:
	friend class Serializer;

       FractionalStepPFEM23D() : Element()
       {
       }



   }; // Class PFEM22D
}  // namespace Kratos.

#endif // KRATOS_MONOLITIC_PFEM2_2D_ELEM_H_INCLUDED  defined
