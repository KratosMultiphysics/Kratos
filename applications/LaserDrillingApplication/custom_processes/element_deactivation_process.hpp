//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  Ignasi de Pouplana
//


#if !defined(KRATOS_ELEMENT_DEACTIVATION_PROCESS )
#define  KRATOS_ELEMENT_DEACTIVATION_PROCESS

#include "includes/table.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "laserdrilling_application_variables.h"

namespace Kratos
{

class ElementDeactivationProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ElementDeactivationProcess);

    /// Defining a table with double argument and result type as table type.
    typedef Table<double,double> TableType;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ElementDeactivationProcess(ModelPart& model_part,
                                Parameters rParameters
                                ) : Process(Flags()) , mr_model_part(model_part)
    {
        KRATOS_TRY

        Parameters default_parameters( R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "thermal_energy_per_volume_threshold": 1.0e20,
                "thermal_counter_threshold": 0,
                "alpha_threshold": 0.8,
                "decomposition_law" : "Prout-Tompkins",
                "decomposition_law_reference_temperature": 400.0,
                "decomposition_law_constant_1": 1e-7,
                "decomposition_law_constant_2": 0.007
            }  )" );

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mthermal_energy_per_volume_threshold = rParameters["thermal_energy_per_volume_threshold"].GetDouble();
        mthermal_counter_threshold = rParameters["thermal_counter_threshold"].GetInt();
        m_alpha_threshold = rParameters["alpha_threshold"].GetDouble();
        m_decomposition_law = rParameters["decomposition_law"].GetString();
        m_decomposition_law_reference_temperature = rParameters["decomposition_law_reference_temperature"].GetDouble();
        m_decomposition_law_constant_1 = rParameters["decomposition_law_constant_1"].GetDouble();
        m_decomposition_law_constant_2 = rParameters["decomposition_law_constant_2"].GetDouble();
        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ElementDeactivationProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ElementDeactivationProcess algorithms.
    void Execute() override
    {
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;

        const ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();

        int NElems = static_cast<int>(mr_model_part.Elements().size());
        ModelPart::ElementsContainerType::iterator el_begin = mr_model_part.ElementsBegin();
        #pragma omp parallel for
        for (int i = 0; i < NElems; i++)
        {
            ModelPart::ElementsContainerType::iterator itElem = el_begin + i;
            Element::GeometryType& rGeom = itElem->GetGeometry();
            // const unsigned int NumNodes = rGeom.PointsNumber();
            GeometryData::IntegrationMethod MyIntegrationMethod = itElem->GetIntegrationMethod();
            unsigned int NumGPoints = rGeom.IntegrationPoints(MyIntegrationMethod).size();
            std::vector<double> ThermalEnergyPerVolumeVector(NumGPoints); // All components in this vector contain the same elemental thermal energy
            std::vector<double> Temperature(NumGPoints);
            std::vector<double> ThermalDecomposition(NumGPoints);
            std::vector<double> DecomposedElementalVolume(NumGPoints);
            std::vector<double> ElementalVolume(NumGPoints);
            //std::vector<double> EnergyPerVolumeVector(NumGPoints);

            itElem->CalculateOnIntegrationPoints(THERMAL_ENERGY_PER_VOLUME, ThermalEnergyPerVolumeVector, CurrentProcessInfo);
            itElem->CalculateOnIntegrationPoints(TEMPERATURE, Temperature, CurrentProcessInfo);
            //itElem->CalculateOnIntegrationPoints(THERMAL_DECOMPOSITION, ThermalDecomposition, CurrentProcessInfo);
            itElem->CalculateOnIntegrationPoints(DECOMPOSED_ELEMENTAL_VOLUME, DecomposedElementalVolume, CurrentProcessInfo);
            itElem->CalculateOnIntegrationPoints(ELEMENTAL_VOLUME, ElementalVolume, CurrentProcessInfo);
            //itElem->CalculateOnIntegrationPoints(ENERGY_PER_VOLUME, EnergyPerVolumeVector, CurrentProcessInfo);
        }

        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ElementDeactivationProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ElementDeactivationProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    ModelPart& mr_model_part;
    double mthermal_energy_per_volume_threshold;
    int mthermal_counter_threshold;
    double m_alpha_threshold;
    std::string m_decomposition_law;
    double m_decomposition_law_reference_temperature;
    double m_decomposition_law_constant_1;
    double m_decomposition_law_constant_2;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    ElementDeactivationProcess& operator=(ElementDeactivationProcess const& rOther);

    /// Copy constructor.
    //ElementDeactivationProcess(ElementDeactivationProcess const& rOther);

}; // Class ElementDeactivationProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ElementDeactivationProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ElementDeactivationProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_ELEMENT_DEACTIVATION_PROCESS defined */
