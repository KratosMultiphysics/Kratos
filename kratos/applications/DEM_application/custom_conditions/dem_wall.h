//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

#if !defined(KRATOS_DEM_WALL_H_INCLUDED )
#define  KRATOS_DEM_WALL_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/variables.h"
//#include "../custom_elements/spheric_particle.h"

namespace Kratos
{
class SphericParticle;
class DEMWall : public Condition
{
public:

    // Counted pointer of DEMWall
    KRATOS_CLASS_POINTER_DEFINITION( DEMWall );
	
	
	typedef WeakPointerVector<Element> ParticleWeakVectorType; 
	typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
	typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;
	
	typedef WeakPointerVector<Condition> ConditionWeakVectorType; 
	typedef WeakPointerVector<Condition >::iterator ConditionWeakIteratorType;


    // Constructor void
    DEMWall();

    // Constructor using an array of nodes
    DEMWall( IndexType NewId, GeometryType::Pointer pGeometry );

    // Constructor using an array of nodes with properties
    DEMWall( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    virtual ~DEMWall();


    // Name Operations

    virtual Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties ) const;


    virtual void Initialize();
    virtual void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& r_process_info );		
    virtual void CalculateElasticForces(VectorType& rRightHandSideVector, ProcessInfo& r_process_info );
    virtual void Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& r_process_info);
    virtual void InitializeSolutionStep(ProcessInfo& r_process_info);  
    virtual void FinalizeSolutionStep(ProcessInfo& r_process_info);          
    virtual void CalculateNormal(array_1d<double, 3>& rnormal);   
    virtual void AddExplicitContribution(const VectorType& rRHS,
                                 const Variable<VectorType>& rRHSVariable,
                                 Variable<array_1d<double,3> >& rDestinationVariable,
                                 const ProcessInfo& r_process_info);
    
    virtual void GetDeltaDisplacement( array_1d<double, 3> & delta_displacement, int inode);
    /*
    double mTgOfFrictionAngle;
    double mYoungModulus;
    double mPoissonRatio;
    */
    
    double GetYoung();
    double GetPoisson();
    double GetTgOfFrictionAngle();
    
    std::vector<SphericParticle*> mNeighbourSphericParticles;
    
    
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param r_process_info
     */



protected:
  
    

private:
    ///@name Static Member Variables

    /// privat variables


    // privat name Operations

	

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
    }


}; // class DEMWall.

} // namespace Kratos.

#endif // KRATOS_DEM_WALL_H_INCLUDED  defined 
 
