//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//



#if !defined(KRATOS_MONOLITHIC_MODIFIED_PFEM2_2D_H_INCLUDED )
#define  KRATOS_MONOLITHIC_MODIFIED_PFEM2_2D_H_INCLUDED

// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "custom_utilities/pfem_particle_fluidonly.h"
#include "includes/serializer.h"
#include "utilities/geometry_utilities.h"
#include "includes/cfd_variables.h"

#include "utilities/divide_triangle_2d_3.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"



namespace Kratos
{

  class MonolithicModifiedPFEM22D
      : public Element
   {
   public:

     /// Counted pointer of PFEM22D
    KRATOS_CLASS_POINTER_DEFINITION(MonolithicModifiedPFEM22D);
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
    typedef PointerVector< PFEM_Particle_Fluid, PFEM_Particle_Fluid*, std::vector<PFEM_Particle_Fluid*> > ParticlePointerVector;

    const int TDim = 2;
    const SizeType NumNodes = TDim+1;
    const SizeType LocalSize = (TDim+1)*(TDim+1);

    struct element_data
    {
        bounded_matrix<double,3, 2> v, vn, f, f_nitsche, f_n, f_nitsche_n; 
        array_1d<double,3> p, rho;
        
        bounded_matrix<double, 3, 2 > DN_DX;
        array_1d<double, 3 > N;
        
        Matrix C;
        Vector stress;
        
        double bdf0;
        double bdf1;
        double bdf2;
        double h;
        double dyn_tau_coeff;
    };
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

     void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

     //void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);


     void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

     void GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo) const override;




protected:

    void ComputeElementAsWATER(bounded_matrix<double,9, 9>& lhs_local,
                              array_1d<double,9>& rhs_local,
                              Matrix& rLeftHandSideMatrix,
                              Vector& rRightHandSideVector,
                              const double& Volume,
                              element_data& data,
                              const ProcessInfo& rCurrentProcessInfo);

    void ComputeElementAsAIR(bounded_matrix<double,9, 9>& lhs_local,
                              array_1d<double,9>& rhs_local,
                              Matrix& rLeftHandSideMatrix,
                              Vector& rRightHandSideVector,
                              const double& Volume,
                              element_data& data,
                              const ProcessInfo& rCurrentProcessInfo);     

    void ComputeElementAsMIXED(bounded_matrix<double,18, 18>& lhs_local, array_1d<double,18>& rhs_local,
                               Matrix& rLeftHandSideMatrix,
                               Vector& rRightHandSideVector,
                               const double& Volume,
                               element_data& data,
                               const ProcessInfo& rCurrentProcessInfo,
                               array_1d<double,3>& distances
                              );                        

    void ComputeConstitutiveResponse(element_data& data, const double rho, const double mu,  const ProcessInfo& rCurrentProcessInfo);
    
    
    void ComputeGaussPointLHSandRHSContribution(bounded_matrix<double,9,9>& lhs, array_1d<double,9>& rhs, const element_data& data, const ProcessInfo& rCurrentProcessInfo, const double& weight);
    
    void ComputeGaussPointLHSandRHSContribution_NitscheDomainTerms(bounded_matrix<double,18,18>& lhs, array_1d<double,18>& rhs, const element_data& data, const ProcessInfo& rCurrentProcessInfo, const double& weight, boost::numeric::ublas::bounded_matrix<double,(18),(18)>& Tdublicate);
    
    void ComputeGaussPointLHSandRHSContribution_NitscheBoundaryTerms(bounded_matrix<double,18,18>& lhs, array_1d<double,18>& rhs, const element_data& data, const ProcessInfo& rCurrentProcessInfo, const double& weight, array_1d<double,3>& normal, boost::numeric::ublas::bounded_matrix<double,(18),(18)>& Tdublicate, boost::numeric::ublas::bounded_matrix<double,(18),(18)>& Tdublicate_neighbor, boost::numeric::ublas::bounded_matrix<double,(6),(6)>& Tself, boost::numeric::ublas::bounded_matrix<double,(6),(6)>& Tneighbor,const Matrix& C_self,const Matrix& C_neighbor);
    


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

    MonolithicModifiedPFEM22D & operator=(MonolithicModifiedPFEM22D const& rOther);

    /// Copy constructor.
    MonolithicModifiedPFEM22D(MonolithicModifiedPFEM22D const& rOther);



   }; // Class PFEM22D
}  // namespace Kratos.

#endif // KRATOS_MONOLITHIC_MODIFIED_PFEM2_2D_H_INCLUDED  defined 
