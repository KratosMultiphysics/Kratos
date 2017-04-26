// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_ALM_VARIABLES_CALCULATION_PROCESS)
#define KRATOS_ALM_VARIABLES_CALCULATION_PROCESS

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "contact_structural_mechanics_application_variables.h"
#include "includes/model_part.h"
#include "geometries/point.h"
#include "custom_utilities/contact_utilities.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
    
///@}
///@name  Enum's
///@{
    
///@}
///@name  Functions
///@{
    
/// Short class definition.
// This process calculates the variables related with the ALM 
/** Detail class definition.
*/
class ALMVariablesCalculationProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ALMVariablesCalculationProcess
    KRATOS_CLASS_POINTER_DEFINITION(ALMVariablesCalculationProcess);
    
    // General type definitions
    typedef Node<3>                                          NodeType;
    typedef Point<3>                                        PointType;
    typedef Geometry<NodeType>                           GeometryType;
    typedef Geometry<PointType>                     GeometryPointType;
    typedef ModelPart::NodesContainerType              NodesArrayType;
    typedef ModelPart::ConditionsContainerType    ConditionsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ALMVariablesCalculationProcess(
        ModelPart& rThisModelPart,
        Variable<double>& rNodalLengthVariable = NODAL_H,
        const double FactorStiffness = 10.0
        ):mrThisModelPart(rThisModelPart), 
          mrNodalLengthVariable(rNodalLengthVariable),
          mFactorStiffness(FactorStiffness)
    {
        KRATOS_TRY;
        
        if (rThisModelPart.GetNodalSolutionStepVariablesList().Has( rNodalLengthVariable ) == false ) // TODO: Ask Riccardo if it is necessary to use GetNodalSolutionStepVariablesList (REQUIRES BUFFER!!!!)
        {
            KRATOS_ERROR << "Missing variable " << rNodalLengthVariable;
        }
        
        KRATOS_CATCH(""); 
    }

    /// Destructor.
    virtual ~ALMVariablesCalculationProcess()
    {
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{
    
    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }
    
    ///@}
    ///@name Operations
    ///@{
    
    virtual void Execute() override
    {
        KRATOS_TRY;
        
        /* We compute the penalty factor */
        
        // We initialize the mean values
        double MeanYoungModulusSlave  = 0.0;
        double MeanNodalHSlave        = 0.0;
        double MeanYoungModulusMaster = 0.0;
        double MeanNodalHMaster       = 0.0;
        
        // We initialize the total areas and volumes
        double TotalVolumeSlave  = 0.0;
        double TotalAreaSlave    = 0.0;
        double TotalVolumeMaster = 0.0;
        double TotalAreaMaster   = 0.0;
        
        // Now we iterate over the conditions to calculate the nodal area
        ConditionsArrayType& pConditions = mrThisModelPart.Conditions();
        auto numConditions = pConditions.end() - pConditions.begin();
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numConditions; i++) 
        {
            auto itCond = pConditions.begin() + i;
            
            // We get the condition geometry
            GeometryType& rThisGeometry = itCond->GetGeometry();
            const unsigned NumNodesGeometry = rThisGeometry.size();
            
            // We get the values from the element
            Element::Pointer pElem = itCond->GetValue(ELEMENT_POINTER);
            
            const double YoungModulus = pElem->GetProperties()[YOUNG_MODULUS];
            const double ElementVolume = pElem->GetGeometry().Area();
            
            // We get the values from the condition
            const double ConditionArea = rThisGeometry.Area();
            const double NodalConditionArea = ConditionArea/NumNodesGeometry;
            
            if (itCond->Is(SLAVE) == true)
            {
                #pragma omp atomic
                TotalVolumeSlave += ElementVolume;
                #pragma omp atomic
                TotalAreaSlave += ConditionArea;
                #pragma omp atomic
                MeanYoungModulusSlave += YoungModulus * ElementVolume;
                
                for (unsigned int iNode = 0; iNode < NumNodesGeometry; iNode++)
                {
                    #pragma omp atomic
                    MeanNodalHSlave += rThisGeometry[iNode].FastGetSolutionStepValue(NODAL_H) * NodalConditionArea;
                }
            }
            
            if (itCond->Is(MASTER) == true)
            {

                #pragma omp atomic
                TotalVolumeMaster += ElementVolume;
                #pragma omp atomic
                TotalAreaMaster += ConditionArea;
                #pragma omp atomic
                MeanYoungModulusMaster += YoungModulus * ElementVolume;
                
                for (unsigned int iNode = 0; iNode < NumNodesGeometry; iNode++)
                {
                    #pragma omp atomic
                    MeanNodalHMaster += rThisGeometry[iNode].FastGetSolutionStepValue(NODAL_H) * NodalConditionArea;
                }
            }
        }
        
        // Now we divide between the total areas and volumes
        MeanNodalHSlave /= (TotalAreaSlave + 1.0e-12); 
        MeanYoungModulusSlave /= (TotalVolumeSlave + 1.0e-12);
        
        MeanNodalHMaster /= (TotalAreaMaster  + 1.0e-12); 
        MeanYoungModulusMaster /= (TotalVolumeMaster + 1.0e-12);
        
        // Finally we compute the penalty factor
        const double PenaltyParameterSlave  = mFactorStiffness * MeanYoungModulusSlave/MeanNodalHSlave;
        const double ScaleFactorSlave    = mFactorStiffness * MeanYoungModulusSlave/MeanNodalHSlave;
        const double PenaltyParameterMaster = mFactorStiffness * MeanYoungModulusMaster/MeanNodalHMaster;
        const double ScaleFactorMaster   = mFactorStiffness * MeanYoungModulusMaster/MeanNodalHMaster; 
        
        mrThisModelPart.GetProcessInfo()[PENALTY_PARAMETER] = (PenaltyParameterSlave > PenaltyParameterMaster) ? PenaltyParameterSlave : PenaltyParameterMaster; // NOTE: > or <? , we are supposed to take the largest of the values (more stiff)
        mrThisModelPart.GetProcessInfo()[SCALE_FACTOR]   = (ScaleFactorSlave > ScaleFactorMaster) ? ScaleFactorSlave : ScaleFactorMaster;
        
        KRATOS_CATCH("");
    }
    
    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        return "ALMVariablesCalculationProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ALMVariablesCalculationProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{
    
    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    ModelPart& mrThisModelPart;
    Variable<double>& mrNodalLengthVariable;
    const double mFactorStiffness;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ALMVariablesCalculationProcess& operator=(ALMVariablesCalculationProcess const& rOther);

    /// Copy constructor.
    //ALMVariablesCalculationProcess(ALMVariablesCalculationProcess const& rOther);


    ///@}

}; // Class ALMVariablesCalculationProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   ALMVariablesCalculationProcess& rThis);
// 
// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const ALMVariablesCalculationProcess& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);
// 
//     return rOStream;
// }

}
#endif /* KRATOS_ALM_VARIABLES_CALCULATION_PROCESS defined */
