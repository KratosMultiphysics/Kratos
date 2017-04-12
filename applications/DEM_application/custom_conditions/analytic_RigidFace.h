// Author: Guillermo Casas (gcasas@cimne.upc.edu)

#if !defined(KRATOS_ANALYTIC_AnalyticRigidFace3D_H_INCLUDED )
#define  KRATOS_ANALYTIC_AnalyticRigidFace3D_H_INCLUDED

// External includes

// Project includes
#include "RigidFace.h"

namespace Kratos
{

class AnalyticRigidFace3D : public RigidFace3D
{
public:

    // Counted pointer of AnalyticRigidFace3D
    KRATOS_CLASS_POINTER_DEFINITION( AnalyticRigidFace3D );
	
    typedef RigidFace3D BaseType;
    typedef WeakPointerVector<Element> ParticleWeakVectorType; 
    typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
    typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;

    typedef WeakPointerVector<Condition> ConditionWeakVectorType; 
    typedef WeakPointerVector<Condition >::iterator ConditionWeakIteratorType;


    // Constructor void
    AnalyticRigidFace3D();

    // Constructor using an array of nodes
    AnalyticRigidFace3D( IndexType NewId, GeometryType::Pointer pGeometry );

    // Constructor using an array of nodes with properties
    AnalyticRigidFace3D( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    virtual ~AnalyticRigidFace3D();

    Condition::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const override;

    void InitializeSolutionStep(ProcessInfo& r_process_info);

    void ComputeConditionRelativeData(int rigid_neighbour_index,
                                      SphericParticle* const particle,
                                      double LocalCoordSystem[3][3],
                                      double& DistPToB,
                                      array_1d<double, 4>& Weight,
                                      array_1d<double, 3>& wall_delta_disp_at_contact_point,
                                      array_1d<double, 3>& wall_velocity_at_contact_point,
                                      int& ContactType) override;

    int CheckSide(SphericParticle* p_particle) override;
    bool IsPhantom() override {return true;}

    int GetNumberOfCrossings();
    std::vector<int> GetSignedCollidingIds();
    std::vector<double> GetCollidingNormalRelativeVelocity();
    std::vector<double> GetCollidingTangentialRelativeVelocity();

private:

    unsigned int mNumberOfCrossingSpheres;
    std::vector<int> mContactingNeighbourSignedIds;
    std::vector<int> mOldContactingNeighbourSignedIds;
    std::vector<int> mCrossers;
    std::vector<double> mCollidingNormalVelocities;
    std::vector<double> mCollidingTangentialVelocities;

    void TestForNewCrosserAndPushBack(SphericParticle* p_particle);

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, DEMWall );
    }

    virtual void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, DEMWall );
    }

}; // class AnalyticRigidFace3D.

} // namespace Kratos.

#endif // KRATOS_ANALYTIC_AnalyticRigidFace3D_H_INCLUDED  defined
