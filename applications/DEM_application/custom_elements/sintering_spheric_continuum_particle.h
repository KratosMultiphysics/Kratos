//   
//   Project Name:        ThermalD       
//   Last Modified by:    $Author: Ferran Arrufat $
//   Date:                $Date: 2015-02-01  $
//   Revision:            $Revision: 1.0.0.0 $
//
//


#if !defined(KRATOS_SINTERING_SPHERIC_CONTINUUM_PARTICLE_H_INCLUDED )
#define  KRATOS_SINTERING_SPHERIC_CONTINUUM_PARTICLE_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "spheric_particle.h"
#include "spheric_continuum_particle.h"
#include "thermal_spheric_particle.h"
#include "Particle_Contact_Element.h"
#include "containers/vector_component_adaptor.h"
#include "discrete_element.h"


namespace Kratos
{
	class KRATOS_API(DEM_APPLICATION) SinteringSphericContinuumParticle : public ThermalSphericParticle<SphericContinuumParticle>
	{
	public:

		/// Pointer definition of SinteringSphericContinuumParticle
		KRATOS_CLASS_POINTER_DEFINITION(SinteringSphericContinuumParticle);

		typedef WeakPointerVector<Element> ParticleWeakVectorType;
		typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
		typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;

		typedef Node <3> NodeType;
		typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
		typedef std::size_t IndexType;
		typedef Geometry<Node < 3 > > GeometryType;
		typedef Properties PropertiesType; 

		/// Default constructor. 
		SinteringSphericContinuumParticle() {};
		SinteringSphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry):ThermalSphericParticle<SphericContinuumParticle>(NewId, pGeometry){}
		SinteringSphericContinuumParticle(IndexType NewId, NodesArrayType const& ThisNodes):ThermalSphericParticle<SphericContinuumParticle>(NewId, ThisNodes){}
		SinteringSphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):ThermalSphericParticle<SphericContinuumParticle>(NewId, pGeometry, pProperties){}


		Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
			return Element::Pointer(new SinteringSphericContinuumParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
		}
		

	/// Destructor.
	virtual ~SinteringSphericContinuumParticle() {};

	void Initialize(const ProcessInfo& r_process_info);

	void UpdatingNeighboursVector(ProcessInfo& r_process_info);

	void SetInitialSinteringSphereContacts(ProcessInfo& r_process_info);

	void InitializeForceComputation(ProcessInfo& r_process_info);

	void AddUpForcesAndProject(double OldCoordSystem[3][3],
		double LocalCoordSystem[3][3],
		double LocalContactForce[3],
		double LocalElasticContactForce[3],
		double GlobalContactForce[3],
		double GlobalElasticContactForce[3],
		double ViscoDampingLocalContactForce[3],
		const double cohesive_force,
		double sinter_driv_force,
		array_1d<double, 3>& r_elastic_force,
		array_1d<double, 3>& r_contact_force,
		const unsigned int i_neighbour_count);

	void ComputeBallToBallContactForce(array_1d<double, 3>& rElasticForce,
		array_1d<double, 3 > & rContactForce,
		array_1d<double, 3>& rInitialRotaMoment,
		ProcessInfo& r_process_info,
		const double dt,
		const bool multi_stage_RHS);


	double sintering_displ;
	double sinter_driv_force;
	std::vector<double> mOldNeighbourSinteringDisplacement; // initialization of a container of sintering displ - old
	std::vector<double> mActualNeighbourSinteringDisplacement; // initialization of a container of sintering displ - actual
	};


}  // namespace Kratos.

#endif // KRATOS_THERMAL_SPHERIC_PARTICLE_H_INCLUDED  defined 


