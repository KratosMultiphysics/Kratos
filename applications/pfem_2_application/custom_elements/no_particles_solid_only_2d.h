//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_NO_PARTICLES_MONOLITHIC_SOLID_ONLY_2D_ELEM_H_INCLUDED)
#define  KRATOS_NO_PARTICLES_MONOLITHIC_SOLID_ONLY_2D_ELEM_H_INCLUDED 

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

  class NoParticlesSolidOnlyPFEM22D
	  : public Element
   {
   public:
     
     /// Counted pointer of PFEM22D
    KRATOS_CLASS_POINTER_DEFINITION(NoParticlesSolidOnlyPFEM22D);
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
     NoParticlesSolidOnlyPFEM22D(IndexType NewId, GeometryType::Pointer pGeometry);
     NoParticlesSolidOnlyPFEM22D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

     /// Destructor.
     virtual ~ NoParticlesSolidOnlyPFEM22D();


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
            
    virtual void GetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
            std::vector<array_1d<double, 3 > >& rValues,
            const ProcessInfo& rCurrentProcessInfo);
            
     void Calculate(const Variable<Vector> &rVariable,
                                     Vector &rOutput,
                                     const ProcessInfo &rCurrentProcessInfo);
      
   protected:
   
        void CalculateLocalFractionalVelocitySystem(MatrixType& rLeftHandSideMatrix,
                                                    VectorType& rRightHandSideVector,
                                                    ProcessInfo& rCurrentProcessInfo);

        void CalculateLocalPressureSystem(MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          ProcessInfo& rCurrentProcessInfo);
        
        //void CalculateLocalFinalVelocitySystem(MatrixType& rLeftHandSideMatrix,
          //                                      VectorType& rRightHandSideVector,
            //                                    ProcessInfo& rCurrentProcessInfo);    
                                                
        void CalculateViscousRHS(ProcessInfo& CurrentProcessInfo);
       
       	void CalculatePressureProjection(ProcessInfo& CurrentProcessInfo);
                                                                                                                     
                                
        void AddElasticityTerm(MatrixType& rDampMatrix,
                                       const boost::numeric::ublas::bounded_matrix<double, 3, 2>& rShapeDeriv,
                                       const double Weight);   
        
        void AddViscousTerm(MatrixType& rDampMatrix,
                                       const boost::numeric::ublas::bounded_matrix<double, 3, 2>& rShapeDeriv,
                                       const double Weight);   
                                                                      
        void UpdateStressesToNewConfiguration(array_1d<double,3>& OutputVector,
									   double& pressure,
                                       const boost::numeric::ublas::bounded_matrix<double, 3, 2>& rShapeDeriv,
                                       const array_1d<double,6>&  previous_step_vel,
                                       const array_1d<double,6>&  previous_step_accel,
                                       const double ShearModulus,
                                       const double delta_t);           
                                       
		void AddViscousTerm(MatrixType& rDampMatrix,
                                       const boost::numeric::ublas::bounded_matrix<double, 3, 2>& rShapeDeriv,
                                       const double viscosity_air,
                                       const double viscosity_water,
                                       array_1d<double,3>  volumes,
                                       array_1d<double,3>  signs,
                                       Matrix Ngauss,
                                       const double Area);
                                       
       void AddViscousTerm(boost::numeric::ublas::bounded_matrix<double, (2-1)*6, (2-1)*6 >& rDampMatrix,
                         boost::numeric::ublas::bounded_matrix<double, (2+1), 2 >& rShapeDeriv,
                         const double Weight);
                         
                         
       void AddViscousTerm(boost::numeric::ublas::bounded_matrix<double, 8, 8 > & LocalAxisExtendedDampMatrix,
						  boost::numeric::ublas::bounded_matrix<double, (2+1), 2 >& rShapeDeriv,
                          std::vector< Matrix >& gauss_gradients, 
						  array_1d<double,3>&  viscosities,
						  array_1d<double,3>&  signs,
						  array_1d<double,3>&  volumes ,
						  const unsigned int ndivisions);
						  
						  
		void AddViscousTerm(boost::numeric::ublas::bounded_matrix<double, 13, 13 > & output,
						  boost::numeric::ublas::bounded_matrix<double, (2+1), 2 >& rShapeDeriv,
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
                                       
       void CalculateInterfaceNormal(boost::numeric::ublas::bounded_matrix<double, 3, 2 >& rPoints, array_1d<double,3>&  distances, array_1d<double,2>&  normal, double & interface_area, array_1d<double,3>&  Ninterface);

       bool invert44(boost::numeric::ublas::bounded_matrix<double, 4, 4 > & m , boost::numeric::ublas::bounded_matrix<double, 4, 4 > & inverse );
       bool invert33(boost::numeric::ublas::bounded_matrix<double, 3, 3 > & m , boost::numeric::ublas::bounded_matrix<double, 3, 3 > & result );
       
       inline bool CalculatePosition(const bounded_matrix<double, 3, 3 > & coordinates,
                                         const double xc, const double yc, const double zc,
                                         array_1d<double, 3 > & N
                                        );
                                        
       inline bool CalculatePosition(Geometry<Node < 3 > >&geom,
                const double xc, const double yc, const double zc,
                array_1d<double, 3 > & N
                );                           
                                        
        inline double CalculateVol(const double x0, const double y0,
                                      const double x1, const double y1,
                                      const double x2, const double y2
                                     );     
                                     
                                     
                                     
         void UpdateParticlePressure(
						 PFEM_Particle & pparticle,
						 Geometry< Node<3> >& geom
						 );

	   void UpdateElementStresses(
						 Vector & element_stresses,
						 const double mu,
						 Geometry< Node<3> >& geom,
						 const boost::numeric::ublas::bounded_matrix<double, 6, 1 >& velocities,
						 const boost::numeric::ublas::bounded_matrix<double, 6, 3 >& B_matrix,
						 const double delta_t
						 );

	   void TestElementWithDruckerPrager(Geometry< Node<3> >& geom);
						 
	   void CalculateReducedTimeStepForPositiveJacobian(const boost::numeric::ublas::bounded_matrix<double, 3, 2>& rShapeDeriv,const boost::numeric::ublas::bounded_matrix<double, 6, 1 >& velocities, const double& delta_t, double& reduced_delta_t);
   
   
template<class T>                                     
bool InvertMatrix(const T& input, T& inverse)  ;                                      

   
   private:
	friend class Serializer;

       NoParticlesSolidOnlyPFEM22D() : Element()
       {
       }


       
   }; // Class PFEM22D
}  // namespace Kratos.

#endif // KRATOS_NO_PARTICLES_MONOLITIC_PFEM2_2D_ELEM_H_INCLUDED  defined
