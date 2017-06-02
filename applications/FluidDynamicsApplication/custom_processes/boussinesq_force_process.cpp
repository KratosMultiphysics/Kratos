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
        KRATOS_ERROR_IF_NOT( rParameters.Has("Gravity") ) <<
        "Boussinesq Force Process Error: \'Gravity\' not found in parameters." << std::endl;

        Parameters rGravityParam = rParameters.GetValue("Gravity");

        KRATOS_ERROR_IF_NOT( rGravityParam.IsArray() ) <<
        "Boussinesq Force Process Error: Given \'Gravity\' parameter is not an array." << std::endl;

        KRATOS_ERROR_IF_NOT( rGravityParam.size() == 3) <<
        "Boussinesq Force Process Error: Given \'Gravity\' parameter is not a size 3 array." << std::endl;

        for (int i = 0; i < 3; i++)
        {
            mrGravity[i] = rGravityParam.GetArrayItem(i).GetDouble();
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
        "Boussinesq Force Process Error: \'AMBIENT_TEMPERATURE\' not given in ProcessInfo." << std::endl;
    }

    void BoussinesqForceProcess::AssignBoussinesqForce()
    {
        const double AmbientTemperature = mrModelPart.GetProcessInfo().GetValue(AMBIENT_TEMPERATURE);
        int NumNodes = mrModelPart.NumberOfNodes();
        #pragma omp parallel for firstprivate(NumNodes,AmbientTemperature)
        for (int i = 0; i < NumNodes; ++i)
        {
            ModelPart::NodeIterator iNode = mrModelPart.NodesBegin() + i;
            double Temperature = iNode->FastGetSolutionStepValue(TEMPERATURE);

            iNode->FastGetSolutionStepValue(BODY_FORCE) = (1. - (Temperature-AmbientTemperature)/Temperature)*mrGravity;
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
