//
//   Project Name:        Kratos
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KratosPFEM2Application_H_INCLUDED )
#define  KratosPFEM2Application_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "custom_utilities/pfem_particle.h"
#include "custom_utilities/pfem_particle_fluidonly.h"

#include "custom_elements/fractional_step_pfem_2_2d.h" //including the file for the element
#include "custom_elements/fractional_step_pfem_2_3d.h" //including the file for the element
#include "custom_elements/monolithic_2fluid_2d.h" //including the file for the element
#include "custom_elements/monolithic_2fluid_3d.h" //including the file for the element
#include "custom_elements/nonewtonian_2fluid_2d.h" //including the file for the element
#include "custom_elements/nonewtonian_2fluid_3d.h" //including the file for the element
#include "custom_elements/monolithic_2fluid_2d_partintegration.h" //including the file for the element
#include "custom_elements/monolithic_2fluid_3d_partintegration.h" //including the file for the element
//#include "custom_elements/vel_enriched_2fluid_2d.h"
//#include "custom_elements/vel_enriched_2fluid_2d_nopressure.h"
#include "custom_elements/qfluid_2d.h"
#include "custom_elements/qfluid_3d.h"
#include "custom_conditions/fixed_velocity_2d.h" //the condition
#include "custom_conditions/fixed_velocity_3d.h" //the condition
#include "custom_conditions/fixed_pressure_2d.h" //the condition
#include "custom_conditions/fixed_pressure_3d.h" //the condition
#include "custom_conditions/autoslip_inlet_3d.h" //the condition


namespace Kratos
{

    ///@name Kratos Globals
    ///@{

    // Variables definition KratosPFEM2Application
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, PRESS_GRADIENT_JUMP)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, PRESS_DISCONTINUITY)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application, double, INV_LAPLACIAN_ENRICH)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ENRICH_RHS)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, G_VALUE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, GRADIENT_DISCONTINUITY)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, PREVIOUS_ITERATION_PRESSURE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, FIRST_ITERATION_PRESSURE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, VELOCITY_OVER_ELEM_SIZE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, MEAN_SIZE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, MEAN_VELOCITY_DIFFERENCE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, SPECIFIC_HEAT_CAPACITY_WATER)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, SPECIFIC_HEAT_CAPACITY_AIR)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, DELTA_TEMPERATURE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, AVAILABLE_AIR_VOLUME)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, AVAILABLE_UNBURNED_AIR_VOLUME)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, OXYGEN_FRACTION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, CORRECTED_DISTANCE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, SOLID_PRESSURE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, SOLID_YP)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, WATER_DISTANCE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, VOLUMETRIC_STRAIN)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ELASTIC_PRESSURE)

    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,bool, USEFUL_ELEMENT_FOR_COMBUSTION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,Vector, ENRICH_LHS_ROW_3D)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,Vector, WATER_GAUSS_POINT)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, WATER_VOLUME)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,Vector, ELEMENT_MEAN_STRESS)
    typedef PointerVector< PFEM_Particle, PFEM_Particle*, std::vector<PFEM_Particle*> > ParticlePointerVector;
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,ParticlePointerVector , PARTICLE_POINTERS)
    typedef PointerVector< PFEM_Particle_Fluid, PFEM_Particle_Fluid*, std::vector<PFEM_Particle_Fluid*> > FluidParticlePointerVector;
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,FluidParticlePointerVector , FLUID_PARTICLE_POINTERS)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,int, NUMBER_OF_PARTICLES)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,int, NUMBER_OF_PARTICLES_AUX)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,int, NUMBER_OF_WATER_PARTICLES)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,int, NUMBER_OF_FLUID_PARTICLES)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,int, PARTICLE_POINTERS_OFFSET)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,int, WATER_PARTICLE_POINTERS_OFFSET)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,int, USE_PRESS_PROJ)
    //KRATOS_DEFINE_VARIABLE(double, IS_AIR)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,ENRICH_LHS_ROW)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,ENRICH_PRESS_PROJ_NEGATIVE)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,ENRICH_PRESS_PROJ_POSITIVE)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,SURFACE_NORMAL)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,SURFACE_COORDINATES)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,PRESS_PROJ_NO_RO)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,DELTA_VELOCITY)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,WATER_VELOCITY)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,WATER_MESH_VELOCITY)

    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,PROJECTED_VELOCITY)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, VOLUME_CORRECTION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, INLET_VELOCITY)


    //	KRATOS_DEFINE_VARIABLE(double, NODAL_AREA)


    ///@}
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name  Functions
    ///@{

    ///@}
    ///@name Kratos Classes
    ///@{

    /// Short class definition.
    /** Detail class definition.
	*/
    class KratosPFEM2Application : public KratosApplication
    {
    public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosPFEM2Application
    KRATOS_CLASS_POINTER_DEFINITION(KratosPFEM2Application);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosPFEM2Application();

    /// Destructor.
    virtual ~KratosPFEM2Application() override {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual void Register() override;

    int a;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

      /// Turn back information as a string.
    virtual std::string Info() const override
    {
        return "KratosPFEM2Application";
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
        KRATOS_WATCH("in my application");
        KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

    protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

    private:
    ///@name Static Member Variables
    ///@{



    ///@}
    ///@name Member Variables
    ///@{

    const FractionalStepPFEM22D   mFractionalStepPFEM22D;
    const FractionalStepPFEM23D   mFractionalStepPFEM23D;
    const MonolithicPFEM22D   mMonolithicPFEM22D;
    const MonolithicPFEM23D   mMonolithicPFEM23D;
    const NoNewtonianMonolithicPFEM22D   mNoNewtonianMonolithicPFEM22D;
    const NoNewtonianMonolithicPFEM23D   mNoNewtonianMonolithicPFEM23D;

    const MonolithicAutoSlipPFEM22D   mMonolithicAutoSlipPFEM22D;
    const MonolithicAutoSlipPFEM23D   mMonolithicAutoSlipPFEM23D;

    //const VelocityEnrichedPFEM22D   mVelocityEnrichedPFEM22D;
    //const VelocityEnrichedPFEM22DNoPressure   mVelocityEnrichedPFEM22DNoPressure;

    const QFluid2D mQFluid2D;
    const QFluid3D mQFluid3D;

    const FixedVelocity2D   mFixedVelocity2D;
    const FixedVelocity3D   mFixedVelocity3D;
    const FixedPressure2D   mFixedPressure2D;
    const FixedPressure3D   mFixedPressure3D;
    const MonolithicAutoSlipInlet3D   mMonolithicAutoSlipInlet3D;


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KratosPFEM2Application& operator=(KratosPFEM2Application const& rOther);

    /// Copy constructor.
    KratosPFEM2Application(KratosPFEM2Application const& rOther);

    ///@}

    }; // Class KratosPFEM2Application

  ///@}


  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  ///@}


}  // namespace Kratos.

#endif // KratosPFEM2Application_H_INCLUDED  defined
