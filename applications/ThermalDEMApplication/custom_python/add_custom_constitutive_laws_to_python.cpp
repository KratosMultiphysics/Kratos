//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_constitutive/heat_exchange_mechanism.h"
#include "custom_constitutive/heat_generation_mechanism.h"
#include "custom_constitutive/conduction/direct_conduction_model.h"
#include "custom_constitutive/conduction/direct_conduction_bob_complete.h"
#include "custom_constitutive/conduction/direct_conduction_bob_modified.h"
#include "custom_constitutive/conduction/direct_conduction_bob_simple.h"
#include "custom_constitutive/conduction/direct_conduction_collision.h"
#include "custom_constitutive/conduction/direct_conduction_pipe.h"
#include "custom_constitutive/conduction/indirect_conduction_model.h"
#include "custom_constitutive/conduction/indirect_conduction_surround_layer.h"
#include "custom_constitutive/conduction/indirect_conduction_vargas.h"
#include "custom_constitutive/conduction/indirect_conduction_voronoi_a.h"
#include "custom_constitutive/conduction/indirect_conduction_voronoi_b.h"
#include "custom_constitutive/convection/convection_model.h"
#include "custom_constitutive/convection/nusselt_gunn.h"
#include "custom_constitutive/convection/nusselt_hanz_marshall.h"
#include "custom_constitutive/convection/nusselt_li_mason.h"
#include "custom_constitutive/convection/nusselt_whitaker.h"
#include "custom_constitutive/radiation/radiation_model.h"
#include "custom_constitutive/radiation/radiation_continuum_krause.h"
#include "custom_constitutive/radiation/radiation_continuum_zhou.h"
#include "custom_constitutive/generation/generation_model.h"
#include "custom_constitutive/generation/generation_dissipation.h"
#include "custom_constitutive/real_contact/real_contact_model.h"
#include "custom_constitutive/real_contact/real_contact_lu.h"
#include "custom_constitutive/real_contact/real_contact_zhou.h"
#include "custom_constitutive/real_contact/real_contact_morris_area.h"
#include "custom_constitutive/real_contact/real_contact_morris_area_time.h"
#include "custom_constitutive/real_contact/real_contact_rangel_area.h"
#include "custom_constitutive/real_contact/real_contact_rangel_area_time.h"
#include "custom_constitutive/sintering_continuum.h"
#include "custom_constitutive/DEM_KDEM_CL.h"

namespace Kratos
{
  namespace Python
  {
    namespace py = pybind11;

    void AddCustomConstitutiveLawsToPython(pybind11::module& m) {

      py::class_<HeatExchangeMechanism, HeatExchangeMechanism::Pointer>(m, "HeatExchangeMechanism")
        .def(py::init<>())
        .def("SetHeatExchangeMechanismInProperties", &HeatExchangeMechanism::SetHeatExchangeMechanismInProperties);

      py::class_<Variable<HeatExchangeMechanism::Pointer>, Variable<HeatExchangeMechanism::Pointer>::Pointer>(m, "HeatExchangeMechanismPointerVariable")
        .def("__str__", &Variable<HeatExchangeMechanism::Pointer>::Info);

      py::class_<Variable<HeatExchangeMechanism*>, Variable<HeatExchangeMechanism*>::Pointer>(m, "HeatExchangeMechanismRawPointerVariable")
        .def("__str__", &Variable<HeatExchangeMechanism*>::Info);


      py::class_<HeatGenerationMechanism, HeatGenerationMechanism::Pointer>(m, "HeatGenerationMechanism")
        .def(py::init<>())
        .def("SetHeatGenerationMechanismInProperties", &HeatGenerationMechanism::SetHeatGenerationMechanismInProperties);

      py::class_<Variable<HeatGenerationMechanism::Pointer>, Variable<HeatGenerationMechanism::Pointer>::Pointer>(m, "HeatGenerationMechanismPointerVariable")
        .def("__str__", &Variable<HeatGenerationMechanism::Pointer>::Info);

      py::class_<Variable<HeatGenerationMechanism*>, Variable<HeatGenerationMechanism*>::Pointer>(m, "HeatGenerationMechanismRawPointerVariable")
        .def("__str__", &Variable<HeatGenerationMechanism*>::Info);


      py::class_<RealContactModel, RealContactModel::Pointer>(m, "RealContactModel")
        .def(py::init<>())
        .def("SetRealContactModelInProperties", &RealContactModel::SetRealContactModelInProperties);

      py::class_<Variable<RealContactModel::Pointer>, Variable<RealContactModel::Pointer>::Pointer>(m, "RealContactModelPointerVariable")
        .def("__str__", &Variable<RealContactModel::Pointer>::Info);

      py::class_<Variable<RealContactModel*>, Variable<RealContactModel*>::Pointer>(m, "RealContactModelRawPointerVariable")
        .def("__str__", &Variable<RealContactModel*>::Info);


      // Conduction -----------------------------------------------------------------------------------------------------------------
      py::class_<DirectConductionModel, DirectConductionModel::Pointer, HeatExchangeMechanism>(m, "DirectConductionModel")
        .def(py::init<>());

      py::class_<DirectConductionBOBComplete, DirectConductionBOBComplete::Pointer, DirectConductionModel>(m, "DirectConductionBOBComplete")
        .def(py::init<>());

      py::class_<DirectConductionBOBModified, DirectConductionBOBModified::Pointer, DirectConductionBOBComplete>(m, "DirectConductionBOBModified")
        .def(py::init<>());

      py::class_<DirectConductionBOBSimple, DirectConductionBOBSimple::Pointer, DirectConductionModel>(m, "DirectConductionBOBSimple")
        .def(py::init<>());

      py::class_<DirectConductionCollision, DirectConductionCollision::Pointer, DirectConductionModel>(m, "DirectConductionCollision")
        .def(py::init<>());

      py::class_<DirectConductionPipe, DirectConductionPipe::Pointer, DirectConductionModel>(m, "DirectConductionPipe")
        .def(py::init<>());

      py::class_<IndirectConductionModel, IndirectConductionModel::Pointer, HeatExchangeMechanism>(m, "IndirectConductionModel")
        .def(py::init<>());

      py::class_<IndirectConductionSurroundLayer, IndirectConductionSurroundLayer::Pointer, IndirectConductionModel>(m, "IndirectConductionSurroundLayer")
        .def(py::init<>());

      py::class_<IndirectConductionVargas, IndirectConductionVargas::Pointer, IndirectConductionModel>(m, "IndirectConductionVargas")
        .def(py::init<>());

      py::class_<IndirectConductionVoronoiA, IndirectConductionVoronoiA::Pointer, IndirectConductionModel>(m, "IndirectConductionVoronoiA")
        .def(py::init<>());

      py::class_<IndirectConductionVoronoiB, IndirectConductionVoronoiB::Pointer, IndirectConductionModel>(m, "IndirectConductionVoronoiB")
        .def(py::init<>());


      // Convection -----------------------------------------------------------------------------------------------------------------
      py::class_<ConvectionModel, ConvectionModel::Pointer, HeatExchangeMechanism>(m, "ConvectionModel")
        .def(py::init<>());

      py::class_<NusseltGunn, NusseltGunn::Pointer, ConvectionModel>(m, "NusseltGunn")
        .def(py::init<>());

      py::class_<NusseltHanzMarshall, NusseltHanzMarshall::Pointer, ConvectionModel>(m, "NusseltHanzMarshall")
        .def(py::init<>());

      py::class_<NusseltLiMason, NusseltLiMason::Pointer, ConvectionModel>(m, "NusseltLiMason")
        .def(py::init<>());

      py::class_<NusseltWhitaker, NusseltWhitaker::Pointer, ConvectionModel>(m, "NusseltWhitaker")
        .def(py::init<>());


      // Radiation ------------------------------------------------------------------------------------------------------------------
      py::class_<RadiationModel, RadiationModel::Pointer, HeatExchangeMechanism>(m, "RadiationModel")
        .def(py::init<>());

      py::class_<RadiationContinuumKrause, RadiationContinuumKrause::Pointer, RadiationModel>(m, "RadiationContinuumKrause")
        .def(py::init<>());

      py::class_<RadiationContinuumZhou, RadiationContinuumZhou::Pointer, RadiationModel>(m, "RadiationContinuumZhou")
        .def(py::init<>());


      // Generation -------------------------------------------------------------------------------------------------------------------
      py::class_<GenerationModel, GenerationModel::Pointer, HeatGenerationMechanism>(m, "GenerationModel")
        .def(py::init<>());

      py::class_<GenerationDissipation, GenerationDissipation::Pointer, GenerationModel>(m, "GenerationDissipation")
        .def(py::init<>());


      // Real contact ---------------------------------------------------------------------------------------------------------------
      py::class_<RealContactLu, RealContactLu::Pointer, RealContactModel>(m, "RealContactLu")
        .def(py::init<>());

      py::class_<RealContactZhou, RealContactZhou::Pointer, RealContactModel>(m, "RealContactZhou")
        .def(py::init<>());

      py::class_<RealContactMorrisArea, RealContactMorrisArea::Pointer, RealContactModel>(m, "RealContactMorrisArea")
        .def(py::init<>());

      py::class_<RealContactMorrisAreaTime, RealContactMorrisAreaTime::Pointer, RealContactModel>(m, "RealContactMorrisAreaTime")
        .def(py::init<>());

      py::class_<RealContactRangelArea, RealContactRangelArea::Pointer, RealContactModel>(m, "RealContactRangelArea")
        .def(py::init<>());

      py::class_<RealContactRangelAreaTime, RealContactRangelAreaTime::Pointer, RealContactModel>(m, "RealContactRangelAreaTime")
        .def(py::init<>());

      // Sintering ---------------------------------------------------------------------------------------------------------------
      py::class_<SinteringContinuum, SinteringContinuum::Pointer, DEM_KDEM>(m, "SinteringContinuum")
        .def(py::init<>())
        ;
    }

  } // namespace Python
} // namespace Kratos