// Authors: Miquel Santasusana msantasusana@cimne.upc.edu,
//           Guillermo Casas (gcasas@cimne.upc.edu)

#if !defined(KRATOS_RIGIDFACE3D_H_INCLUDED)
#define  KRATOS_RIGIDFACE3D_H_INCLUDED

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/variables.h"
#include "dem_wall.h"

namespace Kratos

{

class KRATOS_API(DEM_APPLICATION) RigidFace3D : public DEMWall
{
public:

    // Counted pointer of RigidFace3D
    KRATOS_CLASS_POINTER_DEFINITION( RigidFace3D );
	
	
    typedef WeakPointerVector<Element> ParticleWeakVectorType; 
    typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
    typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;

    typedef WeakPointerVector<Condition> ConditionWeakVectorType; 
    typedef WeakPointerVector<Condition >::iterator ConditionWeakIteratorType;


    // Constructor void
    RigidFace3D();

    // Constructor using an array of nodes
    RigidFace3D( IndexType NewId, GeometryType::Pointer pGeometry );

    // Constructor using an array of nodes with properties
    RigidFace3D( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    virtual ~RigidFace3D();

    class FaceDataBuffer
    {
    public:
        FaceDataBuffer(RigidFace3D* p_this_condition): mpThisCondition(p_this_condition)
        {}

        virtual ~FaceDataBuffer(){}

    void SetCurrentNeighbour(SphericParticle* p_neighbour)
    {
        mpNeighbourParticle = p_neighbour;
    }

    double mIndentation;
    double mLocalRelVel[3];
    RigidFace3D* mpThisCondition;
    SphericParticle* mpNeighbourParticle;
    };

    Condition::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const override;

    void Initialize() override;
    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& r_process_info ) override;
    void CalculateElasticForces(VectorType& rElasticForces, ProcessInfo& r_process_info) override;
    void CalculateNormal(array_1d<double, 3>& rnormal) override;
    void Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& r_process_info) override;
    void FinalizeSolutionStep(ProcessInfo& r_process_info) override;
    int CheckSide(SphericParticle* p_particle) override;
    virtual bool CheckProjectionFallsInside(SphericParticle *p_particle);
    void ComputeConditionRelativeData(int rigid_neighbour_index,
                                      SphericParticle* const particle,
                                      double LocalCoordSystem[3][3],
                                      double& DistPToB,
                                      array_1d<double, 4>& Weight,
                                      array_1d<double, 3>& wall_delta_disp_at_contact_point,
                                      array_1d<double, 3>& wall_velocity_at_contact_point,
                                      int& ContactType) override;
    
    array_1d<double, 3> GetVelocity();

protected:
  
private:

    friend class Serializer;

    template<typename T>
    inline int Sign(T x)
    {
        return (T(0) < x) - (x < T(0));
    }

    virtual void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, DEMWall );
    }

    virtual void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, DEMWall );
    }

}; // class RigidFace3D.

} // namespace Kratos.

#endif // KRATOS_RIGIDFACE3D_H_INCLUDED  defined 
