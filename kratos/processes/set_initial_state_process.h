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
                           const Matrix& rInitialF) :
        mrModelPart(model_part), mInitialStrain(rInitialStrain), 
        mInitialStress(rInitialStress), mInitialF(rInitialF)
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

            Vector aux_initial_strain = mInitialStrain;
            Vector aux_initial_stress = mInitialStrain;
            Vector aux_initial_F      = mInitialF;

            // If the values are set element-wise have priority
            if (it_elem->GetGeometry()->Has(INITIAL_STRAIN_VECTOR)) {
                noalias(aux_initial_strain) = it_elem->GetGeometry()->GetValue(INITIAL_STRAIN_VECTOR)
            }
            if (it_elem->GetGeometry()->Has(INITIAL_STRESS_VECTOR)) {
                noalias(aux_initial_stress) = it_elem->GetGeometry()->GetValue(INITIAL_STRESS_VECTOR)
            }
            if (it_elem->GetGeometry()->Has(INITIAL_DEFORMATION_GRADIENT_MATRIX)) {
                noalias(aux_initial_F) = it_elem->GetGeometry()->GetValue(INITIAL_DEFORMATION_GRADIENT_MATRIX)
            }
            
            InitialState initial_state = InitialState(aux_initial_strain, aux_initial_strain, aux_initial_F);

            // Assign the values to the GP
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


