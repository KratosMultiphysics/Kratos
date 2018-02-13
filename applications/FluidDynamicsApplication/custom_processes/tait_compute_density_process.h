//
//  Main authors:    Riccardo Rossi
//
//

#ifndef KRATOS_TAIT_COMPUTE_DENSITY_PROCESS_PROCESS_H
#define KRATOS_TAIT_COMPUTE_DENSITY_PROCESS_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"
#include "custom_utilities/tait_equation_utility.h"

// Application includes
namespace Kratos
{

class TaitComputeDensityProcess : public Process
{
public:

    /// Pointer definition of TaitComputeDensityProcess
    KRATOS_CLASS_POINTER_DEFINITION(TaitComputeDensityProcess);

    typedef Node<3>                     NodeType;
    typedef Geometry<NodeType>      GeometryType;

    /// Constructor.
    TaitComputeDensityProcess(ModelPart& rModelPart,
                    const int PropertyId): Process(), mrModelPart(rModelPart), mPropertyId(PropertyId)
    {
    }

    /// Destructor.
    ~TaitComputeDensityProcess() {}


    void Execute() override
    {
        KRATOS_TRY;
        
        
        const auto& MaterialProperties = mrModelPart.GetProperties(mPropertyId);
        const Vector& bm  = MaterialProperties[TAIT_PARAMETERS_MOLTEN_STATE];
        const Vector& bs  = MaterialProperties[TAIT_PARAMETERS_SOLID_STATE];
        
#pragma omp parallel for
        for(int i=0; i<mrModelPart.Nodes().size(); ++i)
        {
            auto it = mrModelPart.NodesBegin() + i;
            
            const double T = it->FastGetSolutionStepValue(TEMPERATURE) + 273.15; //transofrming to absolute temperature
            const double P = it->FastGetSolutionStepValue(PRESSURE) + 101325.0; //transforming to absolute pressure
            
            const double density = TaitEquationUtility::CalculateRho(bm,bs,T,P);
            
            it->FastGetSolutionStepValue(DENSITY) = density;
        }

        KRATOS_CATCH("");
    }



    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "TaitComputeDensityProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "TaitComputeDensityProcess";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}


protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ModelPart&                                 mrModelPart;
    unsigned int mPropertyId;


private:



}; // Class TaitComputeDensityProcess

};  // namespace Kratos.

#endif // KRATOS_TAIT_COMPUTE_DENSITY_PROCESS_PROCESS_H
