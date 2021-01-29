//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:         BSD License 
//                   Kratos default license: kratos/license.txt
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
#include "includes/mat_variables.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// The SetInitialStateProcess.
/** This Operation is a derived class from the process.h
 *  
*/
template<SizeType TDim>
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
    SetInitialStateProcess(ModelPart& rModelPart) :
        mrModelPart(rModelPart)
    {
        const SizeType voigt_size = (TDim == 3) ? 6 : 3;

        mInitialStrain.resize(voigt_size, false);
        mInitialStress.resize(voigt_size, false);
        mInitialF.resize(TDim, TDim, false);

        noalias(mInitialStrain) = ZeroVector(voigt_size);
        noalias(mInitialStress) = ZeroVector(voigt_size);
        noalias(mInitialF)      = ZeroMatrix(TDim, TDim);
    }

    /// Full constructor.
    SetInitialStateProcess(ModelPart& rModelPart,
                           const Vector& rInitialStrain,
                           const Vector& rInitialStress,
                           const Matrix& rInitialF) :
        mrModelPart(rModelPart), mInitialStrain(rInitialStrain), 
        mInitialStress(rInitialStress), mInitialF(rInitialF)
    {
    }

    /// Constructor with imposed vector.
    SetInitialStateProcess(ModelPart& rModelPart,
                           const Vector& rInitialStateVector, 
                           const int InitialStateType) :
        mrModelPart(rModelPart)
    {
        const SizeType voigt_size = (TDim == 3) ? 6 : 3;

        mInitialStrain.resize(voigt_size, false);
        mInitialStress.resize(voigt_size, false);
        mInitialF.resize(TDim, TDim, false);

        if (InitialStateType == 0) {
            noalias(mInitialStrain) = rInitialStateVector;
            noalias(mInitialStress) = ZeroVector(voigt_size);
        } else if (InitialStateType == 1) {
            noalias(mInitialStrain) = ZeroVector(voigt_size);
            noalias(mInitialStress) = rInitialStateVector;
        } else {
            noalias(mInitialStrain) = ZeroVector(voigt_size);
            noalias(mInitialStress) = ZeroVector(voigt_size);
        }
        noalias(mInitialF) = ZeroMatrix(TDim, TDim);
    }

    /// Constructor with imposed F.
    SetInitialStateProcess(ModelPart& rModelPart,
                           const Matrix& rInitialStateF) :
        mrModelPart(rModelPart)
    {
        const SizeType voigt_size = (TDim == 3) ? 6 : 3;

        mInitialStrain.resize(voigt_size, false);
        mInitialStress.resize(voigt_size, false);
        mInitialF.resize(TDim, TDim, false);

        noalias(mInitialStrain) = ZeroVector(voigt_size);
        noalias(mInitialStress) = ZeroVector(voigt_size);
        noalias(mInitialF)      = rInitialStateF;
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
        const auto& r_integration_points = it_elem_begin->GetGeometry().IntegrationPoints(it_elem_begin->GetIntegrationMethod());

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
            auto it_elem = it_elem_begin + i;

            Vector aux_initial_strain = mInitialStrain;
            Vector aux_initial_stress = mInitialStress;
            Matrix aux_initial_F      = mInitialF;

            // If the values are set element-wise have priority
            if (it_elem->GetGeometry().Has(INITIAL_STRAIN_VECTOR)) {
                noalias(aux_initial_strain) = (it_elem->GetGeometry()).GetValue(INITIAL_STRAIN_VECTOR);
            }
            if (it_elem->GetGeometry().Has(INITIAL_STRESS_VECTOR)) {
                noalias(aux_initial_stress) = (it_elem->GetGeometry()).GetValue(INITIAL_STRESS_VECTOR);
            }
            if (it_elem->GetGeometry().Has(INITIAL_DEFORMATION_GRADIENT_MATRIX)) {
                noalias(aux_initial_F) = (it_elem->GetGeometry()).GetValue(INITIAL_DEFORMATION_GRADIENT_MATRIX);
            }
            InitialState initial_state = InitialState(aux_initial_strain, aux_initial_stress, aux_initial_F);

            // Assign the values to the GP of the element
            std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector;
            it_elem->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_law_vector, mrModelPart.GetProcessInfo());
            for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
                constitutive_law_vector[point_number]->SetInitialState(initial_state);
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

    Vector mInitialStrain;
    Vector mInitialStress;
    Matrix mInitialF;

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


