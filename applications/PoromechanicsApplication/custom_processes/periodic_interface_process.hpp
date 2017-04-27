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
#include "custom_utilities/comparison_utilities.hpp"//from SolidMechanics
#include "utilities/math_utils.h"

#include "poromechanics_application_variables.h"

namespace Kratos
{

class PeriodicInterfaceProcess : public Process
{

// protected:

    // struct PeriodicBar
    // {
    //     Condition::Pointer pCondition;
    //     std::vector<Element::Pointer> NeighbourElements;
    // };

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


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
                "mesh_id": 0,
                "dimension": 2
            }  )" );
        
        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);
        
        mmesh_id = rParameters["mesh_id"].GetInt();

        if(rParameters["dimension"].GetInt()==2)
            mVoigtSize = 3;
        else
            mVoigtSize = 6;

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------
    
    /// Destructor
    virtual ~PeriodicInterfaceProcess() {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the PeriodicInterfaceProcess algorithms.
    void Execute()
    {
    }
    
    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize()
    {
        KRATOS_TRY;
        
        int NCons = static_cast<int>(mr_model_part.Conditions().size());
        ModelPart::ConditionsContainerType::iterator con_begin = mr_model_part.ConditionsBegin();
        
        #pragma omp parallel for
        for(int i = 0; i < NCons; i++)
        {
            ModelPart::ConditionsContainerType::iterator itCond = con_begin + i;
            itCond->Set(PERIODIC,true);
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

    /// this function will be executed at every time step AFTER writing the output
    void ExecuteAfterOutputStep()
    {
        KRATOS_TRY;
        
        int NCons = static_cast<int>(mr_model_part.Conditions().size());
        ModelPart::ConditionsContainerType::iterator con_begin = mr_model_part.ConditionsBegin();
        
        #pragma omp parallel for
        for(int i = 0; i < NCons; i++)
        {
            ModelPart::ConditionsContainerType::iterator itCond = con_begin + i;
            const Condition::GeometryType& Geom = itCond->GetGeometry();
            
            const Matrix& NodalStressMatrix = Geom[0].FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
            Vector NodalStressVector(mVoigtSize);
            noalias(NodalStressVector) = MathUtils<double>::StressTensorToVector( NodalStressMatrix, mVoigtSize );
            ComparisonUtilities EquivalentStress;
            double VonMisesStress = EquivalentStress.CalculateVonMises(NodalStressVector);
            // TODO: seguir amb el limit de tensio definit per cada grup
            if (VonMisesStress > threshold)
            {
                WeakPointerVector<Element>& rE Geom[0].GetValue(NEIGHBOUR_ELEMENTS);
                for(unsigned int ie = 0; ie < rE.size(); ie++)
                {

                    #pragma omp critical
                    {
                        
                    }
                }

            }
            
            
            
        }
        
        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    std::string Info() const
    {
        return "PeriodicInterfaceProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "PeriodicInterfaceProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    ModelPart& mr_model_part;
    std::size_t mmesh_id;
    int mVoigtSize;


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
