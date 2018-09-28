#include "boussinesq_force_process.h"

namespace Kratos
{
    /* Public functions *******************************************************/
    BoussinesqForceProcess::BoussinesqForceProcess(
        ModelPart& rModelPart,
        Parameters& rParameters):
    mrModelPart(rModelPart),
    mrGravity(array_1d<double,3>(3,0.0))
    {
        // Read settings from parameters
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

        if ( rParameters.Has("thermal_expansion_coefficient") )
        {
            Parameters ThermalExpansionParam = rParameters.GetValue("thermal_expansion_coefficient");
            KRATOS_ERROR_IF_NOT( ThermalExpansionParam.IsDouble() ) <<
            "In Boussinesq Force Process: Given \'thermal_expansion_coefficient\' parameter is not a double." << std::endl;

            mThermalExpansionCoefficient = ThermalExpansionParam.GetDouble();

            KRATOS_ERROR_IF( mThermalExpansionCoefficient <= 0.0 ) <<
            "In Boussinesq Force Process: Incorrect value for \'thermal_expansion_coefficient\' parameter:" << std::endl <<
            "Expected a positive double, got " << mThermalExpansionCoefficient << std::endl;

            mUseAmbientTemperature = false;
        }
        else
        {
            mThermalExpansionCoefficient = 0.0;
            mUseAmbientTemperature = true;
        }
    }

    BoussinesqForceProcess::~BoussinesqForceProcess()
    {

    }

    void BoussinesqForceProcess::Execute()
    {
        this->AssignBoussinesqForce();
    }

    void BoussinesqForceProcess::ExecuteInitialize()
    {
        this->ValidateModelPart();
    }

    void BoussinesqForceProcess::ExecuteInitializeSolutionStep()
    {
        this->AssignBoussinesqForce();
    }

    std::string BoussinesqForceProcess::Info() const
    {
        return "BoussinesqForceProcess";
    }

    void BoussinesqForceProcess::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "BoussinesqForceProcess";
    }

    void BoussinesqForceProcess::PrintData(std::ostream& rOStream) const
    {
    }

    /* Protected functions ****************************************************/

    void BoussinesqForceProcess::ValidateModelPart()
    {
        // Variable registration
        KRATOS_ERROR_IF( AMBIENT_TEMPERATURE.Key() == 0 ) <<
        "\'AMBIENT_TEMPERATURE\' variable is not registered in Kratos." << std::endl;

        KRATOS_ERROR_IF( TEMPERATURE.Key() == 0 ) <<
        "\'TEMPERATURE\' variable is not registered in Kratos." << std::endl;

        KRATOS_ERROR_IF( BODY_FORCE.Key() == 0 ) <<
        "\'BODY_FORCE\' variable is not registered in Kratos." << std::endl;

        // Nodal variables
        KRATOS_ERROR_IF_NOT( mrModelPart.GetNodalSolutionStepVariablesList().Has(TEMPERATURE) ) <<
        "\'TEMPERATURE\' variable is not added to the ModelPart nodal data." << std::endl;

        KRATOS_ERROR_IF_NOT( mrModelPart.GetNodalSolutionStepVariablesList().Has(BODY_FORCE) ) <<
        "\'BODY_FORCE\' variable is not added to the ModelPart nodal data." << std::endl;

        // Variables in ProcessInfo
        KRATOS_ERROR_IF_NOT( mrModelPart.GetProcessInfo().Has(AMBIENT_TEMPERATURE) ) <<
        "In Boussinesq Force Process: \'AMBIENT_TEMPERATURE\' not given in ProcessInfo." << std::endl;
    }

    void BoussinesqForceProcess::AssignBoussinesqForce()
    {
        const double ambient_temperature = mrModelPart.GetProcessInfo().GetValue(AMBIENT_TEMPERATURE);
        KRATOS_ERROR_IF( ambient_temperature <= 0.0 ) <<
        "In Boussinesq Force Process: \'AMBIENT_TEMPERATURE\' obtained from ProcessInfo is incorrect." << std::endl <<
        "Expected a positive double, got " << ambient_temperature << std::endl;

        // Note: the default value of 1/AMBIENT_TEMPERATURE is the usual assumption for perfect gases.
        const double alpha = (mUseAmbientTemperature) ? 1.0 / ambient_temperature : mThermalExpansionCoefficient;

        int num_nodes = mrModelPart.NumberOfNodes();
        #pragma omp parallel for firstprivate(num_nodes,ambient_temperature)
        for (int i = 0; i < num_nodes; ++i)
        {
            ModelPart::NodeIterator iNode = mrModelPart.NodesBegin() + i;
            double temperature = iNode->FastGetSolutionStepValue(TEMPERATURE);

            iNode->FastGetSolutionStepValue(BODY_FORCE) = (1. - alpha*(temperature-ambient_temperature))*mrGravity;
        }

    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const BoussinesqForceProcess& rThis)
    {
        rThis.PrintData(rOStream);
        return rOStream;
    }

}
