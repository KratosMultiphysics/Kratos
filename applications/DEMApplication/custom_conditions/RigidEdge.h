//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

#if !defined(KRATOS_RIGIDEDGE_H_INCLUDED )
#define  KRATOS_RIGIDEDGE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/variables.h"
#include "dem_wall.h"
#include "custom_strategies/schemes/glued_to_wall_scheme.h"

namespace Kratos
{

class KRATOS_API(DEM_APPLICATION) RigidEdge2D : public DEMWall
{
public:
    // Counted pointer of RigidEdge2D
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(RigidEdge2D);

	typedef GlobalPointersVector<Element> ParticleWeakVectorType;
	typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
	typedef GlobalPointersVector<Element >::iterator ParticleWeakIteratorType;

	typedef GlobalPointersVector<Condition> ConditionWeakVectorType;
	typedef GlobalPointersVector<Condition >::iterator ConditionWeakIteratorType;


    /**
     * Default constructor.
     */
    RigidEdge2D( IndexType NewId, GeometryType::Pointer pGeometry);

    RigidEdge2D( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);


    RigidEdge2D( IndexType NewId, GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties,
                        Condition::Pointer Master,
                        Condition::Pointer Slave,
                        Point& MasterContactLocalPoint,
                        Point& SlaveContactLocalPoint,
                        int SlaveIntegrationPointIndex
                        );
    /**
     * Destructor.
     */
    virtual ~RigidEdge2D();


    Condition::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;
    void CalculateNormal(array_1d<double, 3>& rnormal) override;
    void Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& r_process_info) override;
    void FinalizeSolutionStep(const ProcessInfo& r_process_info) override;
    void ComputeConditionRelativeData(int rigid_neighbour_index,
                                    SphericParticle* const particle,
                                    double LocalCoordSystem[3][3],
                                    double& DistPToB,
                                    array_1d<double, 4>& Weight,
                                    array_1d<double, 3>& edge_delta_disp_at_contact_point,
                                    array_1d<double, 3>& edge_velocity_at_contact_point,
                                    int& ContactType) override;

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
    RigidEdge2D() {};

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMWall );
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMWall );
    }

}; // Class ContactLink3DExplicit
}  // namespace Kratos.

#endif // KRATOS_RIGIDEDGE_H_INCLUDED  defined
