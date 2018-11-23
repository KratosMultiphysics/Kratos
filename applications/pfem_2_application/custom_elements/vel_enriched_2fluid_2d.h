//
//   Project Name:        Kratos
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_VELOCITY_ENRICHED_PFEM2_2D_ELEM_H_INCLUDED)
#define  KRATOS_VELOCITY_ENRICHED_PFEM2_2D_ELEM_H_INCLUDED

// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "custom_utilities/pfem_particle_fluidonly.h"



namespace Kratos
{

class VelocityEnrichedPFEM22D : public Element
{
public:

    /// Counted pointer of PFEM22D
    KRATOS_CLASS_POINTER_DEFINITION(VelocityEnrichedPFEM22D);
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
	typedef PointerVector< PFEM_Particle_Fluid, PFEM_Particle_Fluid*, std::vector<PFEM_Particle_Fluid*> > ParticlePointerVector;

    /// Default constructor.
    VelocityEnrichedPFEM22D(IndexType NewId = 0) :
        Element(NewId)
    {}
    VelocityEnrichedPFEM22D(IndexType NewId, const NodesArrayType& ThisNodes) :
        Element(NewId, ThisNodes)
    {}

     VelocityEnrichedPFEM22D(IndexType NewId, GeometryType::Pointer pGeometry);
     VelocityEnrichedPFEM22D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

     /// Destructor.
     virtual ~ VelocityEnrichedPFEM22D() override;


     Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;

     void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

     //void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

     void AddExplicitContribution(ProcessInfo& CurrentProcessInfo) override;

     void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

     void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo) override;

     void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo) override;



protected:

    void CalculatePressureProjection(ProcessInfo& CurrentProcessInfo);


    virtual void AddViscousTerm(MatrixType& rDampMatrix,
        const BoundedMatrix<double, 3, 2>& rShapeDeriv,
        double Viscosity,const double Area);

    void AddViscousTerm(BoundedMatrix<double, 17, 17 > & output,
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

    VelocityEnrichedPFEM22D & operator=(VelocityEnrichedPFEM22D const& rOther);

    /// Copy constructor.
    VelocityEnrichedPFEM22D(VelocityEnrichedPFEM22D const& rOther);

    void CalculateInterfaceNormal(
        BoundedMatrix<double, 3, 2 >& rPoints,
        array_1d<double,3>&  rDistances,
        array_1d<double,2>&  normal,
        double & interface_area,
        array_1d<double,3>&  Ninterface,
        BoundedMatrix<double, 2, 2 >& rInterfacePoints);

    inline void CalculatePosition(
        const bounded_matrix<double, 3, 3 > & coordinates,
        const double xc, const double yc, const double zc, array_1d<double, 3 > & N);

    inline double CalculateVol(
        const double x0, const double y0,const double x1, const double y1, const double x2, const double y2 );


}; // Class PFEM22D
}  // namespace Kratos.

#endif // KRATOS_MONOLITIC_PFEM2_2D_ELEM_H_INCLUDED  defined
