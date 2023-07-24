/*
 * Author: Miguel Angel Celigueta
 *
 *  maceli@cimne.upc.edu
 */

#if !defined (KRATOS_DEMFEM_VOLUME_COUPLING_APPLICATION_H_INCLUDED)
#define KRATOS_DEMFEM_VOLUME_COUPLING_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "DEM_application_variables.h"
#include "structural_mechanics_application_variables.h"

/* CONDITIONS */
#include "custom_elements/volume_coupling_element.h"
#include "custom_elements/volume_coupling_particle.h"

namespace Kratos {

  KRATOS_DEFINE_APPLICATION_VARIABLE(DEMFEM_VOLUME_COUPLING_APPLICATION, double, NODAL_COUPLING_WEIGHT)
  KRATOS_DEFINE_APPLICATION_VARIABLE(DEMFEM_VOLUME_COUPLING_APPLICATION, double, PARTICLE_COUPLING_WEIGHT)
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(DEMFEM_VOLUME_COUPLING_APPLICATION, DISPLACEMENT_MULTIPLIED_MASS)
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(DEMFEM_VOLUME_COUPLING_APPLICATION, DEMFEM_VOLUME_COUPLING_FORCE)





  

    class KRATOS_API(DEMFEM_VOLUME_COUPLING_APPLICATION) DEMFEMVolumeCouplingApplication : public KratosApplication {  // doubt about KRATOS_API(DEMFEM_VOLUME_COUPLING_APPLICATION)

        public:

            // ///@name Type Definitions
            // ///@{
            KRATOS_CLASS_POINTER_DEFINITION(DEMFEMVolumeCouplingApplication);

            /// Default constructor.
            DEMFEMVolumeCouplingApplication();

            /// Destructor.
            virtual ~DEMFEMVolumeCouplingApplication() {}

            virtual void Register() override;

            /// Turn back information as a string.
            virtual std::string Info() const override {
                return "DEMFEMVolumeCouplingApplication";
            }

            /// Print information about this object.
            virtual void PrintInfo(std::ostream& rOStream) const override {
                rOStream << Info();
                PrintData(rOStream);
            }

            ///// Print object's data.
            virtual void PrintData(std::ostream& rOStream) const override {
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

            const VolumeCouplingParticle mVolumeCouplingParticle3D;

    }; // Class KratosDemStructuresCouplingApplication

}  // namespace Kratos

 #endif //KRATOS_DEMFEM_VOLUME_COUPLING_APPLICATION_H_INCLUDED defined


