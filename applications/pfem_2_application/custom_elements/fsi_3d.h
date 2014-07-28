//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_MONOLITIC_NANNO_3D_ELEM_H_INCLUDED)
#define  KRATOS_MONOLITIC_NANNO_3D_ELEM_H_INCLUDED 

// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h" 
#include "custom_utilities/pfem_particle.h"



namespace Kratos
{

  class FsiPFEM23D
	  : public Element
   {
   public:
     
     /// Counted pointer of PFEM22D
    KRATOS_CLASS_POINTER_DEFINITION(FsiPFEM23D);
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
     FsiPFEM23D(IndexType NewId, GeometryType::Pointer pGeometry);
     FsiPFEM23D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

     /// Destructor.
     virtual ~ FsiPFEM23D();


     Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

     void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
     
     //void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
     
     void AddExplicitContribution(ProcessInfo& CurrentProcessInfo);
          
     void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

     void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);

     void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

     virtual void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo);
            
     void Calculate(const Variable<Vector> &rVariable,
                                     Vector &rOutput,
                                     const ProcessInfo &rCurrentProcessInfo);        
      
   protected:                                   
                                
        void AddElasticityTerm(MatrixType& rDampMatrix,
                                       const boost::numeric::ublas::bounded_matrix<double, 4, 3>& rShapeDeriv,
                                       const double Weight);   
        
        void AddViscousTerm(MatrixType& rDampMatrix,
                                       const boost::numeric::ublas::bounded_matrix<double, 4, 3>& rShapeDeriv,
                                       const double Weight);   
                                       
        void CalculatePressureProjection(ProcessInfo& CurrentProcessInfo);
                                                                      
        void UpdateStressesToNewConfiguration(array_1d<double,6>& OutputVector,
									   double& pressure,
                                       const boost::numeric::ublas::bounded_matrix<double, 4, 3>& rShapeDeriv,
                                       const array_1d<double,12>&  previous_step_vel,
                                       const double ShearModulus,
                                       const double delta_t);           
       
       inline bool CalculatePosition(const bounded_matrix<double, 4, 3 > & coordinates,
                                         const double xc, const double yc, const double zc,
                                         array_1d<double, 4 > & N
                                        );
                                        
       inline bool CalculatePosition(Geometry<Node < 3 > >&geom,
										const double xc, const double yc, const double zc,
										array_1d<double, 4 > & N
										);                                                       
                                        
       inline double CalculateVol(const double x0, const double y0, const double z0,
                                      const double x1, const double y1, const double z1,
                                      const double x2, const double y2, const double z2,
                                      const double x3, const double y3, const double z3
                                     );   
       
       void UpdateParticlePressure(
						 PFEM_Particle & pparticle,
						 Geometry< Node<3> >& geom);
                                     
       void UpdateParticleStresses(
						 PFEM_Particle & pparticle,
						 Geometry< Node<3> >& geom,
						 const boost::numeric::ublas::bounded_matrix<double, 12, 1 >& velocities,
						 const array_1d<double,4> & distances,
						 const boost::numeric::ublas::bounded_matrix<double, 12, 6 >& B_matrix,
						 const double delta_t
						 );                              
                                     
       void TestParticleWithDruckerPrager(
						 PFEM_Particle & pparticle
						 );
                                     
       void CalculateReducedTimeStepForPositiveJacobian(const boost::numeric::ublas::bounded_matrix<double, 4, 3>& rShapeDeriv,
														const boost::numeric::ublas::bounded_matrix<double, 12, 1 >& velocities,
														const double& delta_t, 
														double& reduced_delta_t);
														
	
                              
                                          
template<class T>                                     
bool InvertMatrix(const T& input, T& inverse)  ;                                      

   
   private:
	friend class Serializer;

       FsiPFEM23D() : Element()
       {
       }


       
   }; // Class PFEM22D
}  // namespace Kratos.

#endif // KRATOS_MONOLITIC_PFEM2_3D_ELEM_H_INCLUDED  defined
