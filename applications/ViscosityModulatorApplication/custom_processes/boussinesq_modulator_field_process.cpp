#include "boussinesq_modulator_field_process.h"

#include "viscosity_modulator_application.h"
#include "includes/variables.h"

namespace Kratos
{
    /* Public functions *******************************************************/
    BoussinesqModulatorFieldProcess::BoussinesqModulatorFieldProcess(
       ModelPart& rModelPart) :mrModelPart(rModelPart)
    {}

    BoussinesqModulatorFieldProcess::BoussinesqModulatorFieldProcess(
        ModelPart& rModelPart,
        Parameters& rParameters):
    mrModelPart(rModelPart),
    mrGravity(array_1d<double,3>(3,0.0))
    {
        // Gravity
        KRATOS_ERROR_IF_NOT( rParameters.Has("gravity") ) <<
        "In Boussinesq Force Process: \'gravity\' not found in parameters." << std::endl;

        Parameters GravityParam = rParameters.GetValue("gravity");

        KRATOS_ERROR_IF_NOT( GravityParam.IsArray() ) <<
        "In Boussinesq Force Process: Given \'gravity\' parameter is not an array." << std::endl;

        KRATOS_ERROR_IF_NOT( GravityParam.size() == 3) <<
        "In Boussinesq Force Process: Given \'gravity\' parameter is not a size 3 array." << std::endl;

        for (int i = 0; i < 3; i++)
        {
            mrGravity[i] = GravityParam.GetArrayItem(i).GetDouble();
        }

        // Base fluid density
        KRATOS_ERROR_IF_NOT( rParameters.Has("base_fluid_density") ) <<
        "In Boussinesq Force Process: \'base_fluid_density\' not found in parameters." << std::endl;

        Parameters BaseFluidDensityParam = rParameters.GetValue("base_fluid_density");

        KRATOS_ERROR_IF_NOT( BaseFluidDensityParam.IsDouble() ) <<
        "In Boussinesq Force Process: Given \'base_fluid_density\' parameter is not a double." << std::endl;

        mRho0 = BaseFluidDensityParam.GetDouble();

        // Max density
        KRATOS_ERROR_IF_NOT( rParameters.Has("max_density") ) <<
        "In Boussinesq Force Process: \'max_density\' not found in parameters." << std::endl;

        Parameters MaxDensityParam = rParameters.GetValue("max_density");

        KRATOS_ERROR_IF_NOT( MaxDensityParam.IsDouble() ) <<
        "In Boussinesq Force Process: Given \'max_density\' parameter is not a double." << std::endl;

        mRhoMax = MaxDensityParam.GetDouble();

        if(!rParameters.Has("modify_pressure"))
        {
            mModifyPressure = false;
        } else {
            mModifyPressure = rParameters["modify_pressure"].GetBool();
        }
        if(!rParameters.Has("modify_density"))
        {
            mModifyDensity = false;
        } else {
            mModifyDensity = rParameters["modify_density"].GetBool();
        }
        if(!rParameters.Has("r0"))
        {
            mrR0 = ZeroVector(3);
        } else {
            for (int i = 0; i < 3; i++)
            {
                mrR0[i] = rParameters["r0"].GetArrayItem(i).GetDouble();
            }
        }
    }

    BoussinesqModulatorFieldProcess::~BoussinesqModulatorFieldProcess()
    {
    }

    void BoussinesqModulatorFieldProcess::Execute()
    {
        this->AssignBoussinesqForce();
    }

    void BoussinesqModulatorFieldProcess::ExecuteInitialize()
    {
        this->ValidateModelPart();
    }

    void BoussinesqModulatorFieldProcess::ExecuteInitializeSolutionStep()
    {
        // Remove hydrostatic contribution
        if(mModifyPressure) 
        {
            int num_nodes = mrModelPart.NumberOfNodes();
            #pragma omp parallel for firstprivate(num_nodes)
            for (int i = 0; i < num_nodes; ++i)
            {
                ModelPart::NodeIterator iNode = mrModelPart.NodesBegin() + i;

                double gravity_field_potential = 0.0;
                Vector node_coords = iNode->Coordinates();
                for(unsigned d = 0; d < node_coords.size(); d++)
                {
                    gravity_field_potential += mrGravity[d] * (node_coords[d] - mrR0[d]);
                }

                iNode->FastGetSolutionStepValue(PRESSURE) -= gravity_field_potential;
            }
        }
        // this->AssignBoussinesqForce();
    }

    void BoussinesqModulatorFieldProcess::ExecuteFinalizeSolutionStep()
    {
        // Add hydrostatic contribution
        if(mModifyPressure) 
        {
            int num_nodes = mrModelPart.NumberOfNodes();
            #pragma omp parallel for firstprivate(num_nodes)
            for (int i = 0; i < num_nodes; ++i)
            {
                ModelPart::NodeIterator iNode = mrModelPart.NodesBegin() + i;

                double gravity_field_potential = 0.0;
                Vector node_coords = iNode->Coordinates();
                for(unsigned d = 0; d < node_coords.size(); d++)
                {
                    gravity_field_potential += mrGravity[d] * (node_coords[d] - mrR0[d]);
                }

                iNode->FastGetSolutionStepValue(PRESSURE) += gravity_field_potential;
                iNode->FastGetSolutionStepValue(BODY_FORCE) += mrGravity;
            }
        }
        // this->AssignBoussinesqForce();
    }

    std::string BoussinesqModulatorFieldProcess::Info() const
    {
        return "BoussinesqModulatorFieldProcess";
    }

    void BoussinesqModulatorFieldProcess::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "BoussinesqModulatorFieldProcess";
    }

    void BoussinesqModulatorFieldProcess::PrintData(std::ostream& rOStream) const
    {
    }

    /* Protected functions ****************************************************/

    void BoussinesqModulatorFieldProcess::ValidateModelPart()
    {
        ViscosityModulatorSettings::Pointer p_settings = mrModelPart.GetProcessInfo()[VISCOSITY_MODULATOR_SETTINGS];
        const Variable<double>& CONCENTRATION_VAR = p_settings->GetUnknownVariable();

        // Nodal variables
        KRATOS_ERROR_IF_NOT( mrModelPart.GetNodalSolutionStepVariablesList().Has(CONCENTRATION_VAR) ) <<
        CONCENTRATION_VAR.Name() << " variable is not added to the ModelPart nodal data." << std::endl;

        KRATOS_ERROR_IF_NOT( mrModelPart.GetNodalSolutionStepVariablesList().Has(BODY_FORCE) ) <<
        "\'BODY_FORCE\' variable is not added to the ModelPart nodal data." << std::endl;
    }

    void BoussinesqModulatorFieldProcess::AssignBoussinesqForce()
    {
        // Read the unkown variable
        ViscosityModulatorSettings::Pointer p_settings = mrModelPart.GetProcessInfo()[VISCOSITY_MODULATOR_SETTINGS];
        const Variable<double>& CONCENTRATION_VAR = p_settings->GetUnknownVariable();
        const Variable<double>& DENSITY_VAR = p_settings->GetDensityVariable();

        int num_nodes = mrModelPart.NumberOfNodes();
        #pragma omp parallel for firstprivate(num_nodes)
        for (int i = 0; i < num_nodes; ++i)
        {
            ModelPart::NodeIterator iNode = mrModelPart.NodesBegin() + i;

            // Change the body force (and pressure)
            double phi = iNode->FastGetSolutionStepValue(CONCENTRATION_VAR);
            double delta_rho = mRhoMax - mRho0;

            // double gravity_field_potential = 0.0;
            // Vector node_coords = iNode->Coordinates();
            // for(unsigned d = 0; d < node_coords.size(); d++)
            // {
            //     gravity_field_potential += (mrGravity[d] * (node_coords[d] - mrR0[d]));
            // }

            array_1d<double,3> delta_rho_gravity_term = delta_rho / mRho0 * phi * mrGravity;
            if(mModifyPressure)
            {
                iNode->FastGetSolutionStepValue(BODY_FORCE) = delta_rho_gravity_term;
                // iNode->FastGetSolutionStepValue(PRESSURE) -= (1.0 * mRho0 / delta_rho) * gravity_field_potential;

            } else {
                iNode->FastGetSolutionStepValue(BODY_FORCE) += mrGravity + delta_rho_gravity_term;
            }

            // Change density
            if(mModifyDensity)
                iNode->FastGetSolutionStepValue(DENSITY_VAR) = mRho0 + delta_rho * phi;
        }
    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const BoussinesqModulatorFieldProcess& rThis)
    {
        rThis.PrintData(rOStream);
        return rOStream;
    }

}
