#pragma once
// Author: Miguel Angel Celigueta maceli@cimne.upc.edu

#include <pybind11/pybind11.h>
#include "includes/define.h"
#include "includes/model_part.h"


namespace Kratos {

    class AuxiliaryUtilities {

    public:

        KRATOS_CLASS_POINTER_DEFINITION(AuxiliaryUtilities);

        AuxiliaryUtilities() {};
        virtual ~AuxiliaryUtilities() {};

    };
}
