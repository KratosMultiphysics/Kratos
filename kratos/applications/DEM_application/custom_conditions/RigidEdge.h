//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

#if !defined(KRATOS_RIGIDEDGE_H_INCLUDED )
#define  KRATOS_RIGIDEDGE_H_INCLUDED

// System includes

// External includes
//#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/variables.h"
#include "dem_wall.h"

namespace Kratos
{

class RigidEdge3D : public DEMWall
{
public:
    // Counted pointer of RigidEdge3D
    KRATOS_CLASS_POINTER_DEFINITION(RigidEdge3D);
	
	typedef WeakPointerVector<Element> ParticleWeakVectorType; 
	typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
	typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;
	
	typedef WeakPointerVector<Condition> ConditionWeakVectorType; 
	typedef WeakPointerVector<Condition >::iterator ConditionWeakIteratorType;
	

    /**
     * Default constructor.
     */
    RigidEdge3D( IndexType NewId, GeometryType::Pointer pGeometry);

    RigidEdge3D( IndexType NewId, GeometryType::Pointer pGeometry,
                           PropertiesType::Pointer pProperties
                         );


    RigidEdge3D( IndexType NewId, GeometryType::Pointer pGeometry,
                           PropertiesType::Pointer pProperties,
                           Condition::Pointer Master,
                           Condition::Pointer Slave,
                           Point<3>& MasterContactLocalPoint,
                           Point<3>& SlaveContactLocalPoint,
                           int SlaveIntegrationPointIndex
                         );
    /**
     * Destructor.
     */
    virtual ~RigidEdge3D();

    /**
     * Operations.
     */



    Condition::Pointer Create( IndexType NewId,
                               NodesArrayType const& ThisNodes,
                               PropertiesType::Pointer pProperties) const;


    void Initialize();
    void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                 ProcessInfo& r_process_info);
    
    void CalculateElasticForces(VectorType& rElasticForces, ProcessInfo& r_process_info);
    void CalculateNormal(array_1d<double, 3>& rnormal);
    void Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& r_process_info);
    void FinalizeSolutionStep(ProcessInfo& r_process_info);


    /**
     * Turn back information as a string.
     * (DEACTIVATED)
     */
    //std::string Info();

    /**
     * Print information about this object.
     * (DEACTIVATED)
     */
    //virtual void PrintInfo(std::ostream& rOStream) const;

    /**
     * Print object's data.
     * (DEACTIVATED)
     */
    //virtual void PrintData(std::ostream& rOStream) const;

protected:


private:
    
	
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization
    RigidEdge3D() {};

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMWall );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMWall );
    }

}; // Class ContactLink3DExplicit
}  // namespace Kratos.

#endif // KRATOS_RIGIDEDGE_H_INCLUDED  defined 
