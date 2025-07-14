#include "boussinesq_concentration_field_process.h"

#include "convection_diffusion_application.h"

namespace Kratos
{
    /* Public functions *******************************************************/
    BoussinesqConcentrationFieldProcess::BoussinesqConcentrationFieldProcess(
       ModelPart& rModelPart) :mrModelPart(rModelPart)
    {}

    BoussinesqConcentrationFieldProcess::BoussinesqConcentrationFieldProcess(
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

        // Solid particles density
        KRATOS_ERROR_IF_NOT( rParameters.Has("particles_density") ) <<
        "In Boussinesq Force Process: \'particles_density\' not found in parameters." << std::endl;

        Parameters SolidDensityParam = rParameters.GetValue("particles_density");

        KRATOS_ERROR_IF_NOT( SolidDensityParam.IsDouble() ) <<
        "In Boussinesq Force Process: Given \'base_fluid_density\' parameter is not a double." << std::endl;

        mRhoP = SolidDensityParam.GetDouble();

        if(!rParameters.Has("modify_pressure"))
        {
            mModifyPressure = false;
        } else {
            mModifyPressure = rParameters["modify_pressure"].GetBool();
        }
    }

    BoussinesqConcentrationFieldProcess::~BoussinesqConcentrationFieldProcess()
    {
    }

    void BoussinesqConcentrationFieldProcess::Execute()
    {
        this->AssignBoussinesqForce();
    }

    void BoussinesqConcentrationFieldProcess::ExecuteInitialize()
    {
        this->ValidateModelPart();
    }

    void BoussinesqConcentrationFieldProcess::ExecuteInitializeSolutionStep()
    {
        this->AssignBoussinesqForce();
    }

    std::string BoussinesqConcentrationFieldProcess::Info() const
    {
        return "BoussinesqConcentrationFieldProcess";
    }

    void BoussinesqConcentrationFieldProcess::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "BoussinesqConcentrationFieldProcess";
    }

    void BoussinesqConcentrationFieldProcess::PrintData(std::ostream& rOStream) const
    {
    }

    /* Protected functions ****************************************************/

    void BoussinesqConcentrationFieldProcess::ValidateModelPart()
    {
        // Nodal variables
        KRATOS_ERROR_IF_NOT( mrModelPart.GetNodalSolutionStepVariablesList().Has(TEMPERATURE) ) <<
        "\'TEMPERATURE\' variable is not added to the ModelPart nodal data." << std::endl;

        KRATOS_ERROR_IF_NOT( mrModelPart.GetNodalSolutionStepVariablesList().Has(BODY_FORCE) ) <<
        "\'BODY_FORCE\' variable is not added to the ModelPart nodal data." << std::endl;
    }

    void BoussinesqConcentrationFieldProcess::AssignBoussinesqForce()
    {
        // Read the unkown variable
        ConvectionDiffusionSettings::Pointer p_settings = mrModelPart.GetProcessInfo()[CONVECTION_DIFFUSION_SETTINGS];
        const Variable<double>& CONCENTRATION_VAR = p_settings->GetUnknownVariable();
        const Variable<double>& DENSITY_VAR = p_settings->GetDensityVariable();

        int num_nodes = mrModelPart.NumberOfNodes();
        #pragma omp parallel for firstprivate(num_nodes)
        for (int i = 0; i < num_nodes; ++i)
        {
            ModelPart::NodeIterator iNode = mrModelPart.NodesBegin() + i;

            // Change the body force (and pressure)
            double phi = iNode->FastGetSolutionStepValue(CONCENTRATION_VAR);
            double delta_rho = mRhoP - mRho0;

            if(mModifyPressure)
            {
                array_1d<double, 3> r_vec = iNode->Coordinates();
                double gravity_field_potential = r_vec[0] * mrGravity[0] + r_vec[1] * mrGravity[1] + r_vec[2] * mrGravity[2];

                iNode->FastGetSolutionStepValue(BODY_FORCE) = (delta_rho / mRho0) * phi * mrGravity;
                iNode->FastGetSolutionStepValue(PRESSURE) -= mRho0 * gravity_field_potential;
            } else {
                iNode->FastGetSolutionStepValue(BODY_FORCE) += (1. + (delta_rho / mRho0) * phi) * mrGravity;
            }

            // Change density
            iNode->FastGetSolutionStepValue(DENSITY_VAR) = mRho0 + delta_rho * phi;
        }
    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const BoussinesqConcentrationFieldProcess& rThis)
    {
        rThis.PrintData(rOStream);
        return rOStream;
    }

}
