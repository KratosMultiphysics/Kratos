//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//                    
//

#if !defined(KRATOS_SET_INITIAL_STATE_H_INCLUDED )
#define  KRATOS_SET_INITIAL_STATE_H_INCLUDED



// System includes



// External includes


// Project includes

#include "processes/process.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/initial_state.h"




namespace Kratos
{

///@name Kratos Classes
///@{

/// The SetInitialStateProcess.
/** This Operation is a derived class from the process.h
 *  
*/

class SetInitialStateProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(SetInitialStateProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SetInitialStateProcess(ModelPart& rModelPart,
                           const Vector& rInitialStrain,
                           const Vector& rInitialStress,
                           const Matrix& rInitialF,
                           int InitialImposingType, 
                           ) :
        mrModelPart(model_part), mInitialStrain(rInitialStrain), 
        mInitialStress(rInitialStress), mInitialF(rInitialF), 
        mInitialImposingType(InitialImposingType)
    {
    }

    /// Destructor.
    ~SetInitialStateProcess() override {}


    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{


    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
        const auto it_elem_begin = mrModelPart.ElementsBegin();
        const auto& integration_points = it_elem_begin->GetGeometry().IntegrationPoints(it_elem_begin->GetIntegrationMethod());
        
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {

            auto it_elem = it_elem_begin + i;
            InitialState initial_state = InitialState(mInitialStrain, mInitialStress, mInitialF);

            for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number) {
                it_elem->mConstitutiveLawVector[point_number]->SetInitialState(initial_state);
            }
        }
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

private:
    ///@name Static Member Variables
    ///@{


    ModelPart& mrModelPart;

    bool mInitialStrain;
    bool mInitialStress;
    bool mInitialF;

    int mInitialImposingType = 0;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    SetInitialStateProcess& operator=(Process const& rOther);

    ///@}

}; // Class Process

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


}  // namespace Kratos.

#endif //  defined 


