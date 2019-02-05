// Last Modified by: Salva Latorre
//

#if !defined (KRATOS_NVIDIAFLEX_APPLICATION_H_INCLUDED)
#define KRATOS_NVIDIAFLEX_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"

namespace Kratos {

    class KRATOS_API(NVIDIAFLEX_APPLICATION) KratosNvidiaFlexApplication : public KratosApplication {

        public:

            ///@name Type Definitions
            ///@{
            KRATOS_CLASS_POINTER_DEFINITION(KratosNvidiaFlexApplication);

            /// Default constructor.
            KratosNvidiaFlexApplication();

            /// Destructor.
            virtual ~KratosNvidiaFlexApplication() {}

            virtual void Register() override;

            /// Turn back information as a string.
            virtual std::string Info() const override {
                return "KratosNvidiaFlexApplication";
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

            KratosNvidiaFlexApplication& operator=(KratosNvidiaFlexApplication const& rOther);

            /// Copy constructor
            KratosNvidiaFlexApplication(KratosNvidiaFlexApplication const& rOther);

    }; // Class KratosNvidiaFlexApplication

}  // namespace Kratos

#endif //KRATOS_NVIDIAFLEX_APPLICATION_H_INCLUDED defined


