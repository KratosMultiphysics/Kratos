#include "boussinesq_force_process.h"

namespace Kratos
{
    /* Public functions *******************************************************/
    BoussinesqForceProcess::BoussinesqForceProcess(
        ModelPart& rModelPart,
        Parameters& rParameters):
    mrModelPart(rModelPart),
    mrParameters(rParameters),
    mrGravity(array_1d<double,3>(3,0.0))
    {

    }

    BoussinesqForceProcess::~BoussinesqForceProcess()
    {

    }

    void BoussinesqForceProcess::Execute()
    {
        this->AssignBoussinesqForce();
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
