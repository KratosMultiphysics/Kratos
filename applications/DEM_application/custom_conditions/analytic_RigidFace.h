// Author: Guillermo Casas (gcasas@cimne.upc.edu)

#if !defined(KRATOS_ANALYTIC_AnalyticRigidFace3D_H_INCLUDED )
#define  KRATOS_ANALYTIC_AnalyticRigidFace3D_H_INCLUDED

// External includes

// Project includes
#include "RigidFace.h"

namespace Kratos
{

class KRATOS_API(DEM_APPLICATION) AnalyticRigidFace3D : public RigidFace3D
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

    void InitializeSolutionStep(ProcessInfo& r_process_info) override;

    int CheckSide(SphericParticle* p_particle) override;
    bool IsPhantom() override {return true;}

    int GetNumberThroughput();
    std::vector<int> GetSignedCollidingIds();
    int AreThereNewCrossings();
    std::vector<double> GetCollidingNormalRelativeVelocity();
    std::vector<double> GetCollidingTangentialRelativeVelocity();
    std::vector<double> GetMasses();
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Analytic" << RigidFace3D::Info();
        return buffer.str();
    }


private:

    unsigned int mNumberOfCrossingSpheres;
    int mNumberThroughput;
    std::vector<int> mContactingNeighbourSignedIds;
    std::vector<int> mOldContactingNeighbourSignedIds;
    std::vector<int> mCrossers;
    std::vector<double> mCollidingNormalVelocities;
    std::vector<double> mCollidingTangentialVelocities;
    std::vector<double> mMasses;

    template<class Tvalue>
    bool IsInside(const Tvalue& value, const std::vector<Tvalue>& my_vector)
    {
        auto beginning = std::begin(my_vector);
        auto end       = std::end(my_vector);
        const bool is_inside = (end != std::find(beginning, end, value));
        return is_inside;
    }

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
