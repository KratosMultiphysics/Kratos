//   Last Modified by:    $Author: Salva Latorre
//   Date:                $Date: 20151209
//   Revision:            $Revision: 1.2 $
//

#if !defined(KRATOS_DEM_APPLICATION_H_INCLUDED)
#define KRATOS_DEM_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "custom_elements/cylinder_particle.h"
#include "custom_elements/cylinder_continuum_particle.h"
#include "custom_elements/spheric_particle.h"
#include "custom_elements/spheric_continuum_particle.h"
#include "custom_elements/nanoparticle.h"
#include "custom_elements/analytic_spheric_particle.h"
#include "custom_elements/ice_continuum_particle.h"
#include "custom_elements/Particle_Contact_Element.h"
#include "custom_elements/cluster3D.h"
#include "custom_elements/thermal_spheric_particle.h"
#include "custom_elements/sintering_spheric_continuum_particle.h"
#include "custom_elements/bonding_spheric_continuum_particle.h"
#include "custom_elements/custom_clusters/cubecluster3D.h"
#include "custom_elements/custom_clusters/pillcluster3D.h"
#include "custom_elements/custom_clusters/ellipsoidcluster3D.h"
#include "custom_elements/custom_clusters/cuboidcluster3D.h"
#include "custom_elements/custom_clusters/capsulecluster3D.h"
#include "custom_elements/custom_clusters/singlespherecluster3D.h"
#include "custom_elements/custom_clusters/beadcluster3D.h"
#include "custom_conditions/mapping_condition.h"
#include "custom_conditions/SolidFace.h"
#include "custom_conditions/RigidFace.h"
#include "custom_conditions/analytic_RigidFace.h"
#include "custom_conditions/RigidEdge.h"

namespace Kratos
{

class KratosDEMApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(KratosDEMApplication);

    /// Default constructor.
    KratosDEMApplication();

    /// Destructor.
    virtual ~KratosDEMApplication() {}

    virtual void Register() override;

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        return "KratosDEMApplication";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
    }

protected:

private:

    const CylinderParticle mCylinderParticle2D;
    const CylinderContinuumParticle mCylinderContinuumParticle2D;
    const SphericParticle mSphericParticle3D;
    const NanoParticle mNanoParticle3D;
    const AnalyticSphericParticle mAnalyticSphericParticle3D;
    const SphericContinuumParticle mSphericContinuumParticle3D; 
    const IceContinuumParticle mIceContinuumParticle3D; 
    const ThermalSphericParticle<SphericContinuumParticle> mThermalSphericContinuumParticle3D;  
    const ThermalSphericParticle<SphericParticle> mThermalSphericParticle3D;  
    const SinteringSphericContinuumParticle mSinteringSphericContinuumParticle3D;
    const BondingSphericContinuumParticle mBondingSphericContinuumParticle3D;
    const ParticleContactElement mParticleContactElement;
    const SolidFace3D  mSolidFace3D3N;
    const SolidFace3D  mSolidFace3D4N;
    const RigidFace3D  mRigidFace3D3N;
    const AnalyticRigidFace3D  mAnalyticRigidFace3D3N;
    const RigidFace3D  mRigidFace3D4N;
    const RigidEdge3D  mRigidEdge3D2N;
    const Cluster3D  mCluster3D;
    const CubeCluster3D  mCubeCluster3D;
    const PillCluster3D  mPillCluster3D;
    const EllipsoidCluster3D  mEllipsoidCluster3D;
    const CuboidCluster3D  mCuboidCluster3D;
    const CapsuleCluster3D  mCapsuleCluster3D;
    const SingleSphereCluster3D  mSingleSphereCluster3D;
    const BeadCluster3D  mBeadCluster3D;
    const MAPcond    mMapCon3D3N;

    // static const ApplicationCondition  msApplicationCondition;
    KratosDEMApplication& operator=(KratosDEMApplication const& rOther);

    /// Copy constructor.
    KratosDEMApplication(KratosDEMApplication const& rOther);

}; // Class KratosDEMApplication

}  // namespace Kratos.

#endif // KRATOS_DEM_APPLICATION_H_INCLUDED  defined 


