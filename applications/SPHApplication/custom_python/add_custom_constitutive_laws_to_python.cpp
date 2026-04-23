

// Project includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"


// Laws includes
#include "custom_constitutive/volumetric_linear_elastic_2D_law.h"

namespace Kratos::Python {
void AddCustomConstitutiveLawsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_< VolumetricLinearElastic2DLaw, typename VolumetricLinearElastic2DLaw::Pointer, ConstitutiveLaw >
    (m, "VolumetricLinearElastic2DLaw").def(py::init<>() )
    ;
}

}