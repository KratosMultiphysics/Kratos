//
//   Project Name:        Kratos
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

#ifndef KRATOS_MONOLITHIC_PFEM2_2D_MODIFIED_ELEM_H_INCLUDED
#define KRATOS_MONOLITHIC_PFEM2_2D_MODIFIED_ELEM_H_INCLUDED

// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "custom_utilities/pfem_particle_fluidonly.h"
#include "includes/serializer.h"
#include "utilities/geometry_utilities.h"
#include "includes/cfd_variables.h"

#include "utilities/divide_triangle_2d_3.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"



namespace Kratos
{

class MonolithicPFEM22DModified : public Element
{
public:

     /// Counted pointer of PFEM22D
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(MonolithicPFEM22DModified);
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
    MonolithicPFEM22DModified(IndexType NewId = 0) :
        Element(NewId)
    {}

    MonolithicPFEM22DModified(IndexType NewId, const NodesArrayType& ThisNodes) :
        Element(NewId, ThisNodes)
    {}

    MonolithicPFEM22DModified(IndexType NewId, GeometryType::Pointer pGeometry);

    MonolithicPFEM22DModified(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~ MonolithicPFEM22DModified() override;


    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    // void AddExplicitContribution(const ProcessInfo& CurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& ElementalDofList,const ProcessInfo& CurrentProcessInfo) const override;

    void InitializeSolutionStep(const ProcessInfo& CurrentProcessInfo) override;

    


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

    // void ComputeElementAsMIXED(bounded_matrix<double,18, 18>& lhs_local, array_1d<double,18>& rhs_local,
    //                            Matrix& rLeftHandSideMatrix,
    //                            Vector& rRightHandSideVector,
    //                            const double& Volume,
    //                            element_data& data,
    //                            const ProcessInfo& rCurrentProcessInfo,
    //                            array_1d<double,3>& distances
    //                           );                        

    void ComputeConstitutiveResponse(element_data& data, const double rho, const double mu,  const ProcessInfo& rCurrentProcessInfo);
    
    
    void ComputeGaussPointLHSandRHSContribution(bounded_matrix<double,9,9>& lhs, array_1d<double,9>& rhs, const element_data& data, const ProcessInfo& rCurrentProcessInfo, const double& weight);
    
    void ComputeGaussPointLHSandRHSContribution_NitscheDomainTerms(bounded_matrix<double,18,18>& lhs, array_1d<double,18>& rhs, const element_data& data, const ProcessInfo& rCurrentProcessInfo, const double& weight, boost::numeric::ublas::bounded_matrix<double,(18),(18)>& Tdublicate);
    
    void ComputeGaussPointLHSandRHSContribution_NitscheBoundaryTerms(bounded_matrix<double,18,18>& lhs, array_1d<double,18>& rhs, const element_data& data, const ProcessInfo& rCurrentProcessInfo, const double& weight, array_1d<double,3>& normal, boost::numeric::ublas::bounded_matrix<double,(18),(18)>& Tdublicate, boost::numeric::ublas::bounded_matrix<double,(18),(18)>& Tdublicate_neighbor, boost::numeric::ublas::bounded_matrix<double,(6),(6)>& Tself, boost::numeric::ublas::bounded_matrix<double,(6),(6)>& Tneighbor,const Matrix& C_self,const Matrix& C_neighbor);

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

    array_1d<double,(2)> utilde_n = ZeroVector(2); //u tilde at the previous time step

    MonolithicPFEM22DModified & operator=(MonolithicPFEM22DModified const& rOther);

    /// Copy constructor.
    MonolithicPFEM22DModified(MonolithicPFEM22DModified const& rOther);

}; // Class PFEM22D

}  // namespace Kratos.

#endif // KRATOS_MONOLITHIC_PFEM2_2D_MODIFIED_ELEM_H_INCLUDED  defined
