// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#include "swimming_DEM_application.h"
#include "auton_hunt_prudhomme_virtual_mass_force_law.h"

namespace Kratos {

    VirtualMassForceLaw::Pointer AutonHuntPrudhommeVirtualMassForceLaw::Clone() const {
        AutonHuntPrudhommeVirtualMassForceLaw::Pointer p_clone(new AutonHuntPrudhommeVirtualMassForceLaw(*this));
        return p_clone;
    }

    AutonHuntPrudhommeVirtualMassForceLaw::AutonHuntPrudhommeVirtualMassForceLaw(Parameters r_parameters)
    {
        Parameters default_parameters( R"(
            {
                "name":"AutonHuntPrudhommeVirtualMassForceLaw",
                "do_apply_faxen_corrections": false
            }  )" );

        r_parameters.ValidateAndAssignDefaults(default_parameters);
        mDoApplyFaxenCorrections = r_parameters["do_apply_faxen_corrections"].GetBool();
    }

    void AutonHuntPrudhommeVirtualMassForceLaw::Initialize(const ProcessInfo& r_process_info){}

    std::string AutonHuntPrudhommeVirtualMassForceLaw::GetTypeOfLaw() {
        std::string type_of_law = "Auton Hunt and Prud'Homme virtual mass force law";
        return type_of_law;
    }

    double AutonHuntPrudhommeVirtualMassForceLaw::GetVirtualMassCoefficient(Geometry<Node >& r_geometry,
                                                                         const array_1d<double, 3>& minus_slip_acc)
    {
        return 0.5;
    }

    void AutonHuntPrudhommeVirtualMassForceLaw::ComputeForce(Geometry<Node >& r_geometry,
                                                const double fluid_density,
                                                const double displaced_volume,
                                                array_1d<double, 3>& virtual_mass_force,
                                                const ProcessInfo& r_current_process_info)
    {
        const array_1d<double, 3>& fluid_acc = r_geometry[0].FastGetSolutionStepValue(FLUID_ACCEL_PROJECTED);
        const double radius = r_geometry[0].FastGetSolutionStepValue(RADIUS);

        array_1d<double, 3> minus_slip_acc = fluid_acc; // the particle acceleration is assumed to be treated implicitly through the added_mass
        const double fluid_mass = displaced_volume * fluid_density;
        const double virtual_mass_coeff = this->GetVirtualMassCoefficient(r_geometry, minus_slip_acc);
        mLastVirtualMassAddedMass = virtual_mass_coeff * fluid_mass;

        if (mDoApplyFaxenCorrections) {
            const array_1d<double, 3>& fluid_vel_laplacian_rate = r_geometry[0].FastGetSolutionStepValue(FLUID_VEL_LAPL_RATE_PROJECTED);
            noalias(minus_slip_acc) -= 0.1 * radius * radius * fluid_vel_laplacian_rate; // add Faxen term
        }

        noalias(virtual_mass_force) = fluid_mass * virtual_mass_coeff * minus_slip_acc;
        std::cout << "fluid_accel = " << minus_slip_acc << ", " << "f_v = " << virtual_mass_force << " (faxen terms = " << mDoApplyFaxenCorrections << ")" << std::endl;
    }

} // namespace Kratos
