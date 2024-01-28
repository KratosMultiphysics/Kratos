//
//   Project Name:        Kratos
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

#ifndef KRATOS_MONOLITHIC_MODIFIED_PFEM2_2D_ELEM_H_INCLUDED
#define KRATOS_MONOLITHIC_MODIFIED_PFEM2_2D_ELEM_H_INCLUDED

// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "custom_utilities/pfem_particle_fluidonly.h"



namespace Kratos
{

class MonolithicModifiedPFEM22D : public Element
{
public:

     /// Counted pointer of PFEM22D
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(MonolithicModifiedPFEM22D);
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
	typedef PointerVector< PFEM_Particle_Fluid, PFEM_Particle_Fluid*, std::vector<PFEM_Particle_Fluid*> > ParticlePointerVector;

    /// Default constructor.
    MonolithicModifiedPFEM22D(IndexType NewId = 0) :
        Element(NewId)
    {}

    MonolithicModifiedPFEM22D(IndexType NewId, const NodesArrayType& ThisNodes) :
        Element(NewId, ThisNodes)
    {}

    MonolithicModifiedPFEM22D(IndexType NewId, GeometryType::Pointer pGeometry);

    MonolithicModifiedPFEM22D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~ MonolithicModifiedPFEM22D() override;


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
        const BoundedMatrix<double, 3, 2>& rShapeDeriv,
        double& Viscosity,const double Area);

    void AddViscousTerm(BoundedMatrix<double, 13, 13 > & output,
        BoundedMatrix<double, (2+1), 2 >& rShapeDeriv,
        array_1d<double,3>&  distances,
        std::vector< Matrix >& gauss_gradients,
        array_1d<double,3>&  viscosities,
        array_1d<double,3>&  signs,
        array_1d<double,3>&  volumes ,
        const unsigned int ndivisions);

    void AddViscousTerm(BoundedMatrix<double, 12, 12 > & output,
        BoundedMatrix<double, (2+1), 2 >& rShapeDeriv,
        array_1d<double,3>&  distances,
        std::vector< Matrix >& gauss_gradients,
        array_1d<double,3>&  viscosities,
        array_1d<double,3>&  signs,
        array_1d<double,3>&  volumes ,
        const unsigned int ndivisions);

    void  AddViscousTerm(MatrixType& rDampMatrix,
        std::vector< Matrix > & gauss_gradients_discontinuous,
        array_1d<double,3>&  volumes,
        array_1d<double,3>&  viscosities,
        const int ndivisions);

    template<class T>
    bool InvertMatrix(const T& input, T& inverse);


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

    MonolithicModifiedPFEM22D & operator=(MonolithicModifiedPFEM22D const& rOther);

    /// Copy constructor.
    MonolithicModifiedPFEM22D(MonolithicModifiedPFEM22D const& rOther);

}; // Class MODIFIEDPFEM22D

}  // namespace Kratos.

#endif // KRATOS_MONOLITHIC_MODIFIED_PFEM2_2D_ELEM_H_INCLUDED  defined
