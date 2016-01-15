//
//   Project Name:        Kratos
//   Last Modified by:    $Author:  G.Casas$
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_SWIMMING_DEM_APPLICATION_H_INCLUDED )
#define  KRATOS_SWIMMING_DEM_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "includes/dem_variables.h"  //TODO: must be removed eventually
#include "includes/cfd_variables.h"  //TODO: must be removed eventually
#include "includes/legacy_structural_app_vars.h"  //TODO: must be removed eventually
#include "custom_elements/monolithic_dem_coupled.h"
#include "custom_elements/monolithic_dem_coupled_weak.h"
#include "custom_elements/spheric_swimming_particle.h"
#include "../DEM_application/custom_elements/spheric_particle.h"
#include "../DEM_application/custom_elements/nanoparticle.h"

namespace Kratos
{
  //KRATOS_DEFINE_VARIABLE(int,TRACK_SUBSCALES)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(AVERAGED_FLUID_VELOCITY) 
  

class KratosSwimmingDEMApplication : public KratosApplication
{
public:

    /// Pointer definition of KratosSwimmingDEMApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosSwimmingDEMApplication);

    /// Default constructor.
    KratosSwimmingDEMApplication();

    /// Destructor.
    virtual ~KratosSwimmingDEMApplication() {}


    virtual void Register();

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "KratosSwimmingDEMApplication";
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
        //KRATOS_WATCH("in KratosSwimmingDEMApplication");
        //KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
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
    ///@name Static Member Variables
    ///@{
    /// 2D instance of the VMS element
    const MonolithicDEMCoupled<2> mMonolithicDEMCoupled2D;
    /// 3D instance of the VMS element
    const MonolithicDEMCoupled<3> mMonolithicDEMCoupled3D;

    /// 2D instance of the VMS element
    const MonolithicDEMCoupledWeak<2> mMonolithicDEMCoupledWeak2D;
    /// 3D instance of the VMS element
    const MonolithicDEMCoupledWeak<3> mMonolithicDEMCoupledWeak3D;

    /// swimming derivation of spheric basic DEM element (SphericParticle)
    const SphericSwimmingParticle<SphericParticle> mSphericSwimmingParticle3D;
    const SphericSwimmingParticle<NanoParticle> mSwimmingNanoParticle3D;

    //const DEM_FEM_Particle mDEM_FEM_Particle2D;
    const VariablesList mVariablesList;

    //       static const ApplicationCondition  msApplicationCondition;

    ///@}
    ///@name Member Variables
    ///@{
// 		const Elem2D   mElem2D;
// 		const Elem3D   mElem3D;

    /// Assignment operator.
    KratosSwimmingDEMApplication& operator=(KratosSwimmingDEMApplication const& rOther);

    /// Copy constructor.
    KratosSwimmingDEMApplication(KratosSwimmingDEMApplication const& rOther);

}; // Class KratosSwimmingDEMApplication

}  // namespace Kratos.

#endif // KRATOS_SWIMMING_DEM_APPLICATION_H_INCLUDED  defined 


