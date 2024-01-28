//
//   Project Name:        Kratos
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

#ifndef KRATOS_MONOLITHIC_MODIFIED_PFEM2_3D_ELEM_H_INCLUDED
#define KRATOS_MONOLITHIC_MODIFIED_PFEM2_3D_ELEM_H_INCLUDED

// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"


namespace Kratos
{

class MonolithicModifiedPFEM23D : public Element
{
public:

     /// Counted pointer of MonolithicModifiedPFEM23D
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(MonolithicModifiedPFEM23D);
    ///base type: an IndexedObject that automatically has a unique number
    ///typedef IndexedObject BaseType;
    ///Element from which it is derived
    ///typedef VMS<TDim, TNumNodes> ElementBaseType;
    ///definition of node type (default is: Node)

    //typedef Node NodeType;
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
    MonolithicModifiedPFEM23D(IndexType NewId = 0) :
        Element(NewId)
    {}
    MonolithicModifiedPFEM23D(IndexType NewId, const NodesArrayType& ThisNodes) :
        Element(NewId, ThisNodes)
    {}

    /// Default constructor.
    MonolithicModifiedPFEM23D(IndexType NewId, GeometryType::Pointer pGeometry);

    MonolithicModifiedPFEM23D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~ MonolithicModifiedPFEM23D() override;


    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void AddExplicitContribution(const ProcessInfo& CurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& ElementalDofList,const ProcessInfo& CurrentProcessInfo) const override;

    void InitializeSolutionStep(const ProcessInfo& CurrentProcessInfo) override;


protected:

    void CalculatePressureProjection(const ProcessInfo& CurrentProcessInfo);

    virtual void AddViscousTerm(MatrixType& rDampMatrix,
        const BoundedMatrix<double, 4, 3>& rShapeDeriv,
        double& Viscosity,const double Area);

    void AddViscousTerm(BoundedMatrix<double, 21, 21 > & output,
        BoundedMatrix<double, (4), 3 >& rShapeDeriv,
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

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

    MonolithicModifiedPFEM23D & operator=(MonolithicModifiedPFEM23D const& rOther);

    /// Copy constructor.
    MonolithicModifiedPFEM23D(MonolithicModifiedPFEM23D const& rOther);

}; // Class MonolithicModifiedPFEM23D

}  // namespace Kratos.

#endif // KRATOS_MONOLITHIC_MODIFIED_PFEM2_3D_ELEM_H_INCLUDED  defined
