//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//***********************************************************/


#if !defined(KRATOS_MAPPING_CONDITION_H_INCLUDED )
#define  KRATOS_MAPPING_CONDITION_H_INCLUDED



// System includes


// External includes


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
class KRATOS_API(DEM_APPLICATION) MAPcond : public Condition
{
public:

    // Counted pointer of MAPcond
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( MAPcond );


	typedef GlobalPointersVector<Element> ParticleWeakVectorType;
	typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
	typedef GlobalPointersVector<Element >::iterator ParticleWeakIteratorType;

	typedef GlobalPointersVector<Condition> ConditionWeakVectorType;
	typedef GlobalPointersVector<Condition >::iterator ConditionWeakIteratorType;


    // Constructor void
    MAPcond();

    // Constructor using an array of nodes
    MAPcond( IndexType NewId, GeometryType::Pointer pGeometry );

    // Constructor using an array of nodes with properties
    MAPcond( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    virtual ~MAPcond();


    // Name Operations

    virtual Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const override;

    virtual void Initialize(const ProcessInfo& rCurrentProcessInfo) override;
    virtual void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& r_process_info ) override;
    virtual void AddExplicitContribution(const VectorType& rRHS, const Variable<VectorType>& rRHSVariable, Variable<array_1d<double,3> >& rDestinationVariable, const ProcessInfo& r_process_info) override;

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

    virtual void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
    }

    virtual void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
    }


}; // class MAPcond.

} // namespace Kratos.

#endif // KRATOS_MAPPING_CONDITION_H_INCLUDED  defined

