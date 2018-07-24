#ifndef FLEX_WRAPPER_H
#define FLEX_WRAPPER_H

#include "includes/define.h"
#include "custom_external_libraries/NvFlex.h"
#include "custom_external_libraries/NvFlexExt.h"
#include "custom_external_libraries/NvFlexDevice.h"

namespace Kratos {

    class KRATOS_API(NVIDIAFLEX_APPLICATION) FlexWrapper {

        public:

        KRATOS_CLASS_POINTER_DEFINITION(FlexWrapper);
        
        FlexWrapper();
        
        /// Destructor
        virtual ~FlexWrapper();

        void FlexWrapperFunction();

        /// Turn back information as a string
        virtual std::string Info() const;

        /// Print information about this object
        virtual void PrintInfo(std::ostream& rOStream) const;

        /// Print object's data
        virtual void PrintData(std::ostream& rOStream) const;

        protected:

        private:

        /// Assignment operator
        FlexWrapper & operator=(FlexWrapper const& rOther);

        }; // Class FlexWrapper
} // namespace Kratos

#endif // FLEX_WRAPPER_H
