//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

#if !defined(KRATOS_RIGIDFACE3D_H_INCLUDED )
#define  KRATOS_RIGIDFACE3D_H_INCLUDED

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/variables.h"
#include "dem_wall.h"

namespace Kratos
{

class RigidFace3D : public DEMWall
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

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties ) const;

    void Initialize();
    void CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& r_process_info );
		
    void CalculateElasticForces(VectorType& rElasticForces, ProcessInfo& r_process_info);
    void CalculateNormal(array_1d<double, 3>& rnormal);
    void Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& r_process_info);
    void FinalizeSolutionStep(ProcessInfo& r_process_info);
    

protected:
  
private:

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, DEMWall );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, DEMWall );
    }

}; // class RigidFace3D.

} // namespace Kratos.

#endif // KRATOS_RIGIDFACE3D_H_INCLUDED  defined 
