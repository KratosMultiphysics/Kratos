//
//   Project Name:        KratosPoromechanicsApplication $
//   Last modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:               June 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_PERIODIC_INTERFACE_PROCESS )
#define  KRATOS_PERIODIC_INTERFACE_PROCESS

#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "custom_utilities/solid_mechanics_math_utilities.hpp" //from SolidMechanics
#include "utilities/math_utils.h"

#include "poromechanics_application_variables.h"

namespace Kratos
{

class PeriodicInterfaceProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(PeriodicInterfaceProcess);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    PeriodicInterfaceProcess(ModelPart& model_part,
                                Parameters rParameters
                                ) : Process(Flags()) , mr_model_part(model_part)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "dimension": 2,
                "stress_limit": 100.0e6
            }  )" );

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mDimension = rParameters["dimension"].GetInt();
        mStressLimit = rParameters["stress_limit"].GetDouble();

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~PeriodicInterfaceProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the PeriodicInterfaceProcess algorithms.
    void Execute() override
    {
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY;

        int NCons = static_cast<int>(mr_model_part.Conditions().size());
        ModelPart::ConditionsContainerType::iterator con_begin = mr_model_part.ConditionsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NCons; i++)
        {
            ModelPart::ConditionsContainerType::iterator itCond = con_begin + i;
            Condition::GeometryType& rGeom = itCond->GetGeometry();

            itCond->Set(PERIODIC,true);

            rGeom[0].FastGetSolutionStepValue(PERIODIC_PAIR_INDEX) = rGeom[1].Id();
            rGeom[1].FastGetSolutionStepValue(PERIODIC_PAIR_INDEX) = rGeom[0].Id();
        }

        int NElems = static_cast<int>(mr_model_part.Elements().size());
        ModelPart::ElementsContainerType::iterator el_begin = mr_model_part.ElementsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NElems; i++)
        {
            ModelPart::ElementsContainerType::iterator itElem = el_begin + i;
            itElem->Set(ACTIVE,false);
        }

        KRATOS_CATCH("");
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    void ExecuteFinalizeSolutionStep() override
    {
        KRATOS_TRY;

        int NCons = static_cast<int>(mr_model_part.Conditions().size());
        ModelPart::ConditionsContainerType::iterator con_begin = mr_model_part.ConditionsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NCons; i++)
        {
            ModelPart::ConditionsContainerType::iterator itCond = con_begin + i;
            Condition::GeometryType& rGeom = itCond->GetGeometry();

            Matrix NodalStressMatrix(mDimension,mDimension);
            noalias(NodalStressMatrix) = 0.5 * ( rGeom[0].FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR)
                                                + rGeom[1].FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR) );
            Vector PrincipalStresses(mDimension);
            if(mDimension == 2)
            {
                PrincipalStresses[0] = 0.5*(NodalStressMatrix(0,0)+NodalStressMatrix(1,1)) +
                                    sqrt(0.25*(NodalStressMatrix(0,0)-NodalStressMatrix(1,1))*(NodalStressMatrix(0,0)-NodalStressMatrix(1,1)) +
                                            NodalStressMatrix(0,1)*NodalStressMatrix(0,1));
                PrincipalStresses[1] = 0.5*(NodalStressMatrix(0,0)+NodalStressMatrix(1,1)) -
                                    sqrt(0.25*(NodalStressMatrix(0,0)-NodalStressMatrix(1,1))*(NodalStressMatrix(0,0)-NodalStressMatrix(1,1)) +
                                            NodalStressMatrix(0,1)*NodalStressMatrix(0,1));
            }
            else
            {
                noalias(PrincipalStresses) = SolidMechanicsMathUtilities<double>::EigenValuesDirectMethod(NodalStressMatrix);
            }

            // Check whether the principal stress S1 at the node is higher than the prescribed limit to activate the joints
            if (PrincipalStresses[0] >= mStressLimit)
            {
                itCond->Set(PERIODIC,false);
                rGeom[0].FastGetSolutionStepValue(PERIODIC_PAIR_INDEX) = 0;
                rGeom[1].FastGetSolutionStepValue(PERIODIC_PAIR_INDEX) = 0;

                WeakPointerVector<Element>& rE = rGeom[0].GetValue(NEIGHBOUR_ELEMENTS);
                for(unsigned int ie = 0; ie < rE.size(); ie++)
                {
                    #pragma omp critical
                    {
                        rE[ie].Set(ACTIVE,true);
                    }
                }
            }
        }

        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "PeriodicInterfaceProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PeriodicInterfaceProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    ModelPart& mr_model_part;
    int mDimension;
    double mStressLimit;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    PeriodicInterfaceProcess& operator=(PeriodicInterfaceProcess const& rOther);

    /// Copy constructor.
    //PeriodicInterfaceProcess(PeriodicInterfaceProcess const& rOther);

}; // Class PeriodicInterfaceProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  PeriodicInterfaceProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const PeriodicInterfaceProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_PERIODIC_INTERFACE_PROCESS defined */
