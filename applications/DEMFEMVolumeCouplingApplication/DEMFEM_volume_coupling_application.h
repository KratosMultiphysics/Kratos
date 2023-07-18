/*
 * Author: Miguel Angel Celigueta
 *
 *  maceli@cimne.upc.edu
 */

// #if !defined (KRATOS_DEM_STRUCTURES_COUPLING_APPLICATION_H_INCLUDED)
// #define KRATOS_DEM_STRUCTURES_COUPLING_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"

/* CONDITIONS */
#include "custom_elements/volume_coupling_element.h"


namespace Kratos {

  KRATOS_DEFINE_APPLICATION_VARIABLE(STRUCTURAL_MECHANICS_APPLICATION, double, NODAL_COUPLING_WEIGHT)
  KRATOS_DEFINE_APPLICATION_VARIABLE(DEM_APPLICATION, double, PARTICLE_COUPLING_WEIGHT)
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(DEM_APPLICATION, DISPLACEMENT_MULTIPLIED_MASS)
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(DEM_APPLICATION, DEMFEM_VOLUME_COUPLING_FORCE)

  KRATOS_CREATE_VARIABLE(double, NODAL_COUPLING_WEIGHT )
  KRATOS_CREATE_VARIABLE(double, PARTICLE_COUPLING_WEIGHT ) // if i use same variable NODAL_COUPLING_WEIGHT for particles and nodes, then do i need create and register it twice?
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_MULTIPLIED_MASS)
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DEMFEM_VOLUME_COUPLING_FORCE)

  KRATOS_REGISTER_VARIABLE(double, NODAL_COUPLING_WEIGHT )
  KRATOS_REGISTER_VARIABLE(double, PARTICLE_COUPLING_WEIGHT )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_MULTIPLIED_MASS)
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DEMFEM_VOLUME_COUPLING_FORCE)



  

    class KRATOS_API(DEMFEM_VOLUME_COUPLING_APPLICATION) FEMDEMVolumeCouplingApplication : public KratosApplication {  // doubt about KRATOS_API(DEMFEM_VOLUME_COUPLING_APPLICATION)

        public:

            // ///@name Type Definitions
            // ///@{
            // KRATOS_CLASS_POINTER_DEFINITION(KratosDemStructuresCouplingApplication);

            // /// Default constructor.
            // KratosDemStructuresCouplingApplication();

            // /// Destructor.
            // virtual ~KratosDemStructuresCouplingApplication() {}

            // virtual void Register() override;

            // /// Turn back information as a string.
            // virtual std::string Info() const override {
            //     return "KratosDemStructuresCouplingApplication";
            // }

            // /// Print information about this object.
            // virtual void PrintInfo(std::ostream& rOStream) const override {
            //     rOStream << Info();
            //     PrintData(rOStream);
            // }

            // ///// Print object's data.
            // virtual void PrintData(std::ostream& rOStream) const override {
            //     rOStream << "Variables:" << std::endl;
            //     KratosComponents<VariableData>().PrintData(rOStream);
            //     rOStream << std::endl;
            //     rOStream << "Elements:" << std::endl;
            //     KratosComponents<Element>().PrintData(rOStream);
            //     rOStream << std::endl;
            //     rOStream << "Conditions:" << std::endl;
            //     KratosComponents<Condition>().PrintData(rOStream);
            // }

        protected:

        private:

            // 
            
            const VolumeCouplingElement mVolumeCouplingElement2D3N;
            const VolumeCouplingElement mVolumeCouplingElement2D4N;
            const VolumeCouplingElement mVolumeCouplingElement2D6N;
            const VolumeCouplingElement mVolumeCouplingElement2D8N;
            const VolumeCouplingElement mVolumeCouplingElement2D9N;
            const VolumeCouplingElement mVolumeCouplingElement2D10N;
            const VolumeCouplingElement mVolumeCouplingElement2D15N;
            const VolumeCouplingElement mVolumeCouplingElement3D4N;
            const VolumeCouplingElement mVolumeCouplingElement3D5N;
            const VolumeCouplingElement mVolumeCouplingElement3D6N;
            const VolumeCouplingElement mVolumeCouplingElement3D8N;
            const VolumeCouplingElement mVolumeCouplingElement3D10N;
            const VolumeCouplingElement mVolumeCouplingElement3D13N;
            const VolumeCouplingElement mVolumeCouplingElement3D15N;
            const VolumeCouplingElement mVolumeCouplingElement3D20N;
            const VolumeCouplingElement mVolumeCouplingElement3D27N;

    }; // Class KratosDemStructuresCouplingApplication

}  // namespace Kratos

// #endif //KRATOS_DEM_STRUCTURES_COUPLING_APPLICATION_H_INCLUDED defined


