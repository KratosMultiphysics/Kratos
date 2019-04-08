//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//


#if !defined(KRATOS_CONTROL_MODULE_FEM_DEM_PROCESS )
#define  KRATOS_CONTROL_MODULE_FEM_DEM_PROCESS

#include "custom_processes/control_module_process.h"

#include "dem_structures_coupling_application_variables.h"

namespace Kratos
{

class ControlModuleFemDemProcess : public ControlModuleProcess
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ControlModuleFemDemProcess);

    /// Defining a table with double argument and result type as table type.
    typedef Table<double,double> TableType;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ControlModuleFemDemProcess(
        ModelPart& rModelPart,
        Parameters rParameters
        ) : ControlModuleProcess(rModelPart,rParameters) {}

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ControlModuleFemDemProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ControlModuleFemDemProcess algorithms.
    void Execute() override
    {
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY;

        const int NNodes = static_cast<int>(mrModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();
        typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > ComponentType;
        ComponentType VarComponent = KratosComponents< ComponentType >::Get(mVariableName);

        #pragma omp parallel for
        for(int i = 0; i<NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator it = it_begin + i;

            it->Fix(VarComponent);
            it->FastGetSolutionStepValue(VarComponent) = 0.0;
        }

        // TODO: change this variable name
        mrModelPart[CHARGING_VELOCITY] = mVelocity;

        KRATOS_CATCH("");
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;

        const int NNodes = static_cast<int>(mrModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();
        typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > ComponentType;
        ComponentType VarComponent = KratosComponents< ComponentType >::Get(mVariableName);
        const double DeltaTime = mrModelPart.GetProcessInfo()[DELTA_TIME];
        // TODO: change this variable name
        mVelocity = mrModelPart[CHARGING_VELOCITY];

        #pragma omp parallel for
        for(int i = 0; i<NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator it = it_begin + i;

            it->FastGetSolutionStepValue(VarComponent) += mVelocity * DeltaTime;
        }

        KRATOS_CATCH("");
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    void ExecuteFinalizeSolutionStep() override
    {
        // This must do nothing
    }


    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ControlModuleFemDemProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ControlModuleFemDemProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    ControlModuleFemDemProcess& operator=(ControlModuleFemDemProcess const& rOther);

    /// Copy constructor.
    //ControlModuleFemDemProcess(ControlModuleFemDemProcess const& rOther);

}; // Class ControlModuleFemDemProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ControlModuleFemDemProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ControlModuleFemDemProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_CONTROL_MODULE_FEM_DEM_PROCESS defined */
