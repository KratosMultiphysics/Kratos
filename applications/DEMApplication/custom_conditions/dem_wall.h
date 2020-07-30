//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

#if !defined(KRATOS_DEM_WALL_H_INCLUDED )
#define  KRATOS_DEM_WALL_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/condition.h"
//#include "includes/variables.h"
//#include "../custom_elements/spheric_particle.h"

namespace Kratos
{
class SphericParticle;
class KRATOS_API(DEM_APPLICATION) DEMWall : public Condition
{
public:

    // Counted pointer of DEMWall
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( DEMWall );


	typedef GlobalPointersVector<Element> ParticleWeakVectorType;
	typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
	typedef GlobalPointersVector<Element >::iterator ParticleWeakIteratorType;

	typedef GlobalPointersVector<Condition> ConditionWeakVectorType;
	typedef GlobalPointersVector<Condition >::iterator ConditionWeakIteratorType;


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
        PropertiesType::Pointer pProperties ) const override;


    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;
    virtual void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& r_process_info ) override;
    virtual void CalculateElasticForces(VectorType& rRightHandSideVector, ProcessInfo& r_process_info );
    virtual void InitializeSolutionStep(ProcessInfo& r_process_info) override;
    virtual void FinalizeSolutionStep(ProcessInfo& r_process_info) override;
    virtual void CalculateNormal(array_1d<double, 3>& rnormal);
    virtual void AddExplicitContribution(const VectorType& rRHS,
                                 const Variable<VectorType>& rRHSVariable,
                                 Variable<array_1d<double,3> >& rDestinationVariable,
                                 const ProcessInfo& r_process_info) override;

    virtual void GetDeltaDisplacement( array_1d<double, 3> & delta_displacement, int inode);
    virtual void ComputeConditionRelativeData(int rigid_neighbour_index,
                                              SphericParticle* const particle,
                                              double LocalCoordSystem[3][3],
                                              double& DistPToB,
                                              array_1d<double, 4>& Weight,
                                              array_1d<double, 3>& wall_delta_disp_at_contact_point,
                                              array_1d<double, 3>& wall_velocity_at_contact_point,
                                              int& ContactType){
        KRATOS_ERROR << "Base class DemWall method ComputeConditionRelativeData was called!" << std::endl;
    }
    virtual bool IsPhantom(){return false;}
    virtual int CheckSide(SphericParticle* p_particle){return 1.0;}

    /*
    double mTgOfFrictionAngle;
    double mYoungModulus;
    double mPoissonRatio;
    */

    double GetYoung();
    double GetPoisson();
    double GetTgOfFrictionAngle();

    std::vector<SphericParticle*> mNeighbourSphericParticles;
    std::vector<array_1d <double, 3> > mRightHandSideVector;

    virtual void GetRightHadSideVector(std::vector<array_1d <double, 3> >& rRightHandSideVector);
    virtual void SetRightHadSideVector(const std::vector<array_1d <double, 3> >& rRightHandSideVector);
    virtual void AddToRightHadSideVector(const std::vector<array_1d <double, 3> >& rRightHandSideVector);
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param r_process_info
     */
    std::vector<SphericParticle*>& GetVectorOfGluedParticles() {
        return mVectorOfGluedParticles;
    }


protected:



private:
    ///@name Static Member Variables
    std::vector<SphericParticle*> mVectorOfGluedParticles;
    /// private variables


    // private name Operations


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
        //rSerializer.save("mRightHandSideVector", mRightHandSideVector);
    }

    virtual void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
        //rSerializer.load("mRightHandSideVector", mRightHandSideVector);
    }


}; // class DEMWall.

} // namespace Kratos.

#endif // KRATOS_DEM_WALL_H_INCLUDED  defined

