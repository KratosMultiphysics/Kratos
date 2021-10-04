#include "tait_equation_process.h"

namespace Kratos
{
    /* Public functions *******************************************************/
    TaitEquationProcess::TaitEquationProcess(
        ModelPart& rModelPart,
        Parameters& rParameters):
    mrModelPart(rModelPart)
    {
        // Read settings from parameters
        KRATOS_ERROR_IF_NOT( rParameters.Has("rho_0") ) <<
        "In Tait Equation Process: \'rho_0\' not found in parameters." << std::endl;

        Parameters rhoParam = rParameters.GetValue("rho_0");
        mrho_0 = rhoParam.GetDouble();

        KRATOS_ERROR_IF( mrho_0 <= 0.0 ) <<
            "In Tait Equation Process: Incorrect value for \'rho_0\' parameter:" << std::endl <<
            "Expected a positive double, got " << mrho_0 << std::endl;

        KRATOS_ERROR_IF_NOT( rParameters.Has("p_0") ) <<
        "In Tait Equation Process: \'p_0\' not found in parameters." << std::endl;

        Parameters PressureParam = rParameters.GetValue("p_0");
        mp_0 = PressureParam.GetDouble();
        
        KRATOS_ERROR_IF( mp_0 <= 0.0 ) <<
            "In Tait Equation Process: Incorrect value for \'p_0\' parameter:" << std::endl <<
            "Expected a positive double, got " << mp_0 << std::endl;

        KRATOS_ERROR_IF_NOT( rParameters.Has("k_0") ) <<
        "In Tait Equation Process: \'k_0\' not found in parameters." << std::endl;

        Parameters kParam = rParameters.GetValue("k_0");
        mk_0 = kParam.GetDouble();

        KRATOS_ERROR_IF_NOT( rParameters.Has("theta") ) <<
        "In Tait Equation Process: \'theta\' not found in parameters." << std::endl;

        Parameters ThetaParam = rParameters.GetValue("theta");
        mtheta = ThetaParam.GetDouble();
    }

    TaitEquationProcess::~TaitEquationProcess()
    {

    }

    void TaitEquationProcess::Execute()
    {
        this->AssignTaitEquation();
    }

    void TaitEquationProcess::ExecuteInitialize()
    {
        this->ValidateModelPart();
    }

    void TaitEquationProcess::ExecuteInitializeSolutionStep()
    {
        this->AssignTaitEquation();
    }

    std::string TaitEquationProcess::Info() const
    {
        return "TaitEquationProcess";
    }

    void TaitEquationProcess::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "TaitEquationProcess";
    }

    void TaitEquationProcess::PrintData(std::ostream& rOStream) const
    {
    }

    /* Protected functions ****************************************************/

    void TaitEquationProcess::ValidateModelPart()
    {
        // Nodal variables
        KRATOS_ERROR_IF_NOT( mrModelPart.GetNodalSolutionStepVariablesList().Has(DENSITY) ) <<
        "\'DENSITY\' variable is not added to the ModelPart nodal data." << std::endl;

        KRATOS_ERROR_IF_NOT( mrModelPart.GetNodalSolutionStepVariablesList().Has(SOUND_VELOCITY) ) <<
        "\'SOUND_VELOCITY\' variable is not added to the ModelPart nodal data." << std::endl;

        KRATOS_ERROR_IF_NOT( mrModelPart.GetNodalSolutionStepVariablesList().Has(PRESSURE) ) <<
        "\'PRESSURE\' variable is not added to the ModelPart nodal data." << std::endl;
    }

    void TaitEquationProcess::AssignTaitEquation()
    {
        int num_nodes = mrModelPart.NumberOfNodes();
        // #pragma omp parallel for firstprivate(num_nodes,ambient_temperature)
        for (int i = 0; i < num_nodes; ++i)
        {
            ModelPart::NodeIterator iNode = mrModelPart.NodesBegin() + i;
            double pressure = iNode->FastGetSolutionStepValue(PRESSURE);
            double mod_rho = mrho_0*std::pow((pressure - mp_0)/mk_0 + 1,1/mtheta);
            double modified_c = std::pow(mk_0*mtheta*std::pow(mod_rho/mrho_0,mtheta - 1)/mrho_0,0.5); 
            iNode->FastGetSolutionStepValue(DENSITY,0) = mod_rho;
            iNode->SetValue(SOUND_VELOCITY,modified_c);
        }

    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const TaitEquationProcess& rThis)
    {
        rThis.PrintData(rOStream);
        return rOStream;
    }

}
