//
//   Project Name:        Kratos
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_MONOLITHIC_3FLUID_PFEM2_3D_ELEM_H_INCLUDED)
#define  KRATOS_MONOLITHIC_3FLUID_PFEM2_3D_ELEM_H_INCLUDED

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

class Monolithic3FluidPFEM23D : public Element
{
public:

    /// Counted pointer of PFEM22D
    KRATOS_CLASS_POINTER_DEFINITION(Monolithic3FluidPFEM23D);
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
    typedef VectorMap<IndexType, DataValueContainer> SolutionStepsElementalDataContainerType;

    /// Default constructor.
     Monolithic3FluidPFEM23D(IndexType NewId, GeometryType::Pointer pGeometry);
     Monolithic3FluidPFEM23D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

     /// Destructor.
     virtual ~ Monolithic3FluidPFEM23D() override;


     Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;

     void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

     //void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

     void AddExplicitContribution(ProcessInfo& CurrentProcessInfo) override;

     void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

     void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo) override;

     void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo) override;


protected:

    void CalculateViscousRHS(ProcessInfo& CurrentProcessInfo);

   	void CalculatePressureProjection(ProcessInfo& CurrentProcessInfo);

    void AddViscousTerm(MatrixType& rDampMatrix,
        const boost::numeric::ublas::bounded_matrix<double, 4, 3>& rShapeDeriv,
        const double Weight);

    void AddViscousTerm(boost::numeric::ublas::bounded_matrix<double, 21, 21 > & output,
        boost::numeric::ublas::bounded_matrix<double, 4, 3 >& rShapeDeriv,
        array_1d<double,4>&  distances,
        std::vector< Matrix >& gauss_gradients,
        array_1d<double,6>&  viscosities,
        array_1d<double,6>&  signs,
        array_1d<double,6>&  volumes ,
        const unsigned int ndivisions);

    template<class T>
    bool InvertMatrix(const T& input, T& inverse)  ;


private:
    friend class Serializer;

    Monolithic3FluidPFEM23D() : Element()
    {
    }


}; // Class PFEM22D
}  // namespace Kratos.

#endif // KRATOS_MONOLITHIC_PFEM2_2D_ELEM_H_INCLUDED  defined
