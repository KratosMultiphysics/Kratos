/*
 * Author: Miguel Angel Celigueta
 *
 *  maceli@cimne.upc.edu
 */

#if !defined (KRATOS_DEM_STRUCTURES_COUPLING_APPLICATION_H_INCLUDED)
#define KRATOS_DEM_STRUCTURES_COUPLING_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"

/* CONDITIONS */
#include "custom_conditions/surface_load_from_DEM_condition_3d.h"

namespace Kratos {

    class KRATOS_API(DEM_STRUCTURES_COUPLING_APPLICATION) KratosDemStructuresCouplingApplication : public KratosApplication {

        public:

            ///@name Type Definitions
            ///@{
            KRATOS_CLASS_POINTER_DEFINITION(KratosDemStructuresCouplingApplication);

            /// Default constructor.
            KratosDemStructuresCouplingApplication();

            /// Destructor.
            virtual ~KratosDemStructuresCouplingApplication() {}

            virtual void Register() override;

            /// Turn back information as a string.
            virtual std::string Info() const override {
                return "KratosDemStructuresCouplingApplication";
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

            // Surface load from DEM
            const SurfaceLoadFromDEMCondition3D mSurfaceLoadFromDEMCondition3D3N;

            KratosDemStructuresCouplingApplication& operator=(KratosDemStructuresCouplingApplication const& rOther);

            /// Copy constructor
            KratosDemStructuresCouplingApplication(KratosDemStructuresCouplingApplication const& rOther);

    }; // Class KratosDemStructuresCouplingApplication

}  // namespace Kratos

#endif //KRATOS_DEM_STRUCTURES_COUPLING_APPLICATION_H_INCLUDED defined


