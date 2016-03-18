//   Last Modified by:    $Author: Salva Latorre
//   Date:                $Date: 20151209
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_DEM_APPLICATION_H_INCLUDED)
#define KRATOS_DEM_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/kratos_flags.h"
#include "containers/flags.h"
#include "includes/variables.h"
#include "custom_elements/cylinder_particle.h"
#include "custom_elements/cylinder_continuum_particle.h"
#include "custom_elements/spheric_particle.h"
#include "custom_elements/spheric_continuum_particle.h"
#include "custom_elements/nanoparticle.h"
#include "custom_elements/Particle_Contact_Element.h"
#include "custom_elements/cluster3D.h"
#include "custom_elements/thermal_spheric_particle.h"
#include "custom_elements/custom_clusters/linecluster3D.h"
#include "custom_elements/custom_clusters/cubecluster3D.h"
#include "custom_elements/custom_clusters/pillcluster3D.h"
#include "custom_elements/custom_clusters/ellipsoidcluster3D.h"
#include "custom_elements/custom_clusters/ringcluster3D.h"
#include "custom_elements/custom_clusters/cuboidcluster3D.h"
#include "custom_elements/custom_clusters/cornkernelcluster3D.h"
#include "custom_elements/custom_clusters/corn3cluster3D.h"
#include "custom_elements/custom_clusters/soybeancluster3D.h"
#include "custom_elements/custom_clusters/rock1cluster3D.h"
#include "custom_elements/custom_clusters/rock2cluster3D.h"
#include "custom_elements/custom_clusters/wheat5cluster3D.h"
#include "custom_elements/custom_clusters/ballast1cluster3D.h"
#include "custom_elements/custom_clusters/ballast2cluster3D.h"
#include "custom_elements/custom_clusters/ballast3cluster3D.h"
#include "custom_elements/custom_clusters/ballast4cluster3D.h"
#include "custom_elements/custom_clusters/capsulecluster3D.h"
#include "custom_elements/custom_clusters/ballast5cluster3D.h"
#include "custom_elements/custom_clusters/ballast6cluster3D.h"
#include "custom_conditions/mapping_condition.h"
#include "custom_conditions/SolidFace.h"
#include "custom_conditions/RigidFace.h"
#include "custom_conditions/RigidEdge.h"
#include "custom_constitutive/DEM_Dempack_CL.h"
#include "custom_constitutive/DEM_Dempack_2D_CL.h"
#include "custom_constitutive/DEM_Dempack_torque_CL.h"

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

    virtual void Register();

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "KratosDEMApplication";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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
    const SphericContinuumParticle mSphericContinuumParticle3D; 
    const ThermalSphericParticle mThermalSphericContinuumParticle3D;  
    const ParticleContactElement mParticleContactElement;
    const VariablesList mVariablesList;
    const SolidFace3D  mSolidFace3D3N;
    const SolidFace3D  mSolidFace3D4N;
    const RigidFace3D  mRigidFace3D3N;
    const RigidFace3D  mRigidFace3D4N;
    const RigidEdge3D  mRigidEdge3D2N;
    const Cluster3D  mCluster3D;
    const LineCluster3D  mLineCluster3D;
    const CubeCluster3D  mCubeCluster3D;
    const PillCluster3D  mPillCluster3D;
    const EllipsoidCluster3D  mEllipsoidCluster3D;
    const RingCluster3D  mRingCluster3D;
    const CuboidCluster3D  mCuboidCluster3D;
    const CornKernelCluster3D  mCornKernelCluster3D;
    const Corn3Cluster3D  mCorn3Cluster3D;
    const SoyBeanCluster3D  mSoyBeanCluster3D;
    const Rock1Cluster3D  mRock1Cluster3D;
    const Rock2Cluster3D  mRock2Cluster3D;
    const Wheat5Cluster3D  mWheat5Cluster3D;
    const Ballast1Cluster3D  mBallast1Cluster3D;
    const Ballast2Cluster3D  mBallast2Cluster3D;
    const Ballast3Cluster3D  mBallast3Cluster3D;
    const Ballast4Cluster3D  mBallast4Cluster3D;
    const CapsuleCluster3D  mCapsuleCluster3D;
    const Ballast5Cluster3D  mBallast5Cluster3D;
    const Ballast6Cluster3D  mBallast6Cluster3D;
    const MAPcond    mMapCon3D3N;

    //       static const ApplicationCondition  msApplicationCondition;
    KratosDEMApplication& operator=(KratosDEMApplication const& rOther);

    /// Copy constructor.
    KratosDEMApplication(KratosDEMApplication const& rOther);

}; // Class KratosDEMApplication

}  // namespace Kratos.

#endif // KRATOS_DEM_APPLICATION_H_INCLUDED  defined 


