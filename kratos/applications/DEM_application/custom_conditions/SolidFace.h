//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

#if !defined(KRATOS_SOLIDFACE3D_H_INCLUDED )
#define  KRATOS_SOLIDFACE3D_H_INCLUDED

// External includes
//#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/condition.h"
//#include "includes/variables.h"
#include "dem_wall.h"

namespace Kratos
{

class SolidFace3D : public DEMWall
{
public:

    // Counted pointer of SolidFace3D
    KRATOS_CLASS_POINTER_DEFINITION( SolidFace3D );
	
	
    typedef WeakPointerVector<Element> ParticleWeakVectorType; 
    typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
    typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;

    typedef WeakPointerVector<Condition> ConditionWeakVectorType; 
    typedef WeakPointerVector<Condition >::iterator ConditionWeakIteratorType;


    // Constructor void
    SolidFace3D();

    // Constructor using an array of nodes
    SolidFace3D( IndexType NewId, GeometryType::Pointer pGeometry );

    // Constructor using an array of nodes with properties
    SolidFace3D( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    virtual ~SolidFace3D();

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties ) const;

    void Initialize();
    virtual void CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& r_process_info );
		
    void CalculateNormal(array_1d<double, 3>& rnormal);
    void FinalizeSolutionStep(ProcessInfo& r_process_info);
    
    void GetDeltaDisplacement( array_1d<double, 3> & delta_displacement, int inode);

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

}; // class SolidFace3D.

} // namespace Kratos.

#endif // KRATOS_SOLIDFACE3D_H_INCLUDED  defined 
