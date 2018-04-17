// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder 
//   

#ifndef ADJOINT_STRAIN_ENERGY_RESPONSE_FUNCTION_H
#define ADJOINT_STRAIN_ENERGY_RESPONSE_FUNCTION_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "adjoint_structural_response_function.h"

// ==============================================================================

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

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.

 */

//template<class TDenseSpace>

class AdjointStrainEnergyResponseFunction : public AdjointStructuralResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    typedef AdjointStructuralResponseFunction BaseType;
    typedef array_1d<double, 3> array_3d;



    /// Pointer definition of AdjointStrainEnergyResponseFunction
    KRATOS_CLASS_POINTER_DEFINITION(AdjointStrainEnergyResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AdjointStrainEnergyResponseFunction(ModelPart& model_part, Parameters& responseSettings)
    : AdjointStructuralResponseFunction(model_part, responseSettings)
    {

        // Initialize member variables to NULL
        mCurrentResponseValue = 0.0;

    }

    /// Destructor.
    virtual ~AdjointStrainEnergyResponseFunction()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override
    {
        KRATOS_TRY;

        BaseType::Initialize();

        // It is necessary to initialize the elements since no adjoint problem is solved for this response type.
        // For this response type the elements are only created!
        ModelPart& r_model_part = this->GetModelPart();
    #pragma omp parallel
        {
            ModelPart::ElementIterator elements_begin;
            ModelPart::ElementIterator elements_end;
            OpenMPUtils::PartitionedIterators(r_model_part.Elements(), elements_begin, elements_end);
            for (auto it = elements_begin; it != elements_end; ++it)
                it->Initialize();
        }
        // TODO: Check if initialization is also necessary for conditions!

        KRATOS_CATCH("");
    }

    // ==============================================================================
    double CalculateValue(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        ModelPart& r_model_part = rModelPart; 
        ProcessInfo &r_current_process_info = r_model_part.GetProcessInfo();
        mCurrentResponseValue = 0.0;

        // Check if there are at the time of calling adjoint or primal elements
        KRATOS_ERROR_IF( r_current_process_info[IS_ADJOINT] )
             << "Calculate value for strain energy response is not availible when using adjoint elements" << std::endl;
            
        // Sum all elemental strain energy values calculated as: W_e = u_e^T K_e u_e
        for (auto& elem_i : r_model_part.Elements())
        {
            Matrix LHS;
            Vector RHS;
            Vector disp;

            // Get state solution relevant for energy calculation
            elem_i.GetValuesVector(disp,0);

            elem_i.CalculateLocalSystem(LHS, RHS, r_current_process_info);

            // Compute strain energy
            mCurrentResponseValue += 0.5 * inner_prod(disp, prod(LHS,disp));
         }

        return mCurrentResponseValue;

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
    virtual std::string Info() const
    {
        return "AdjointStrainEnergyResponseFunction";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "AdjointStrainEnergyResponseFunction";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}



    // =============================================================================
    void UpdateSensitivities() override
    {
        KRATOS_TRY;

        BaseType::UpdateSensitivities();

        KRATOS_CATCH("");
    }


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

    // ==============================================================================
    void CalculateSensitivityGradient(Element& rAdjointElem,
                                              const Variable<array_1d<double,3>>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo) override
    {
          KRATOS_TRY

          if (rResponseGradient.size() != rDerivativesMatrix.size1())
              rResponseGradient.resize(rDerivativesMatrix.size1(), false);
        rResponseGradient.clear();

        // There will be a mistake, if body forces are considered. Because the elements are responsible for the body forces!

          KRATOS_CATCH("")
    }

    // ==============================================================================
    void CalculateSensitivityGradient(Element& rAdjointElem,
                                              const Variable<double>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo) override
    {
          KRATOS_TRY

        if (rResponseGradient.size() != rDerivativesMatrix.size1())
              rResponseGradient.resize(rDerivativesMatrix.size1(), false);
        rResponseGradient.clear();

        // There will be a mistake, if body forces are considered. Because the elements are responsible for the body forces!

        KRATOS_CATCH("")
    }

    // ==============================================================================
    void CalculateSensitivityGradient(Condition& rAdjointCondition,
                                              const Variable<array_1d<double,3>>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo) override
    {
        KRATOS_TRY;

        Vector adjoint_variables;

        rAdjointCondition.GetValuesVector(adjoint_variables); // = 0.5*u

        KRATOS_ERROR_IF(adjoint_variables.size() != rDerivativesMatrix.size2())
            << "Size of adjoint vector does not fit to the size of the pseudo load!" << std::endl;

        if (rResponseGradient.size() != rDerivativesMatrix.size2())
            rResponseGradient.resize(adjoint_variables.size(), false);

        noalias(rResponseGradient) = prod(rDerivativesMatrix, adjoint_variables);

        KRATOS_CATCH("");
    }

    // ==============================================================================
    void CalculateSensitivityGradient(Condition& rAdjointCondition,
                                              const Variable<double>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo) override
    {
        KRATOS_TRY;

        Vector adjoint_variables;

        rAdjointCondition.GetValuesVector(adjoint_variables); // = 0.5*u

        KRATOS_ERROR_IF(adjoint_variables.size() != rDerivativesMatrix.size2())
             << "Size of adjoint vector does not fit to the size of the pseudo load!" << std::endl;

        if (rResponseGradient.size() != rDerivativesMatrix.size2())
            rResponseGradient.resize(adjoint_variables.size(), false);

        noalias(rResponseGradient) = prod(rDerivativesMatrix, adjoint_variables);

        KRATOS_CATCH("");
    }

    // ==============================================================================

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

    double mCurrentResponseValue;

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
    //      AdjointStrainEnergyResponseFunction& operator=(AdjointStrainEnergyResponseFunction const& rOther);

    /// Copy constructor.
    //      AdjointStrainEnergyResponseFunction(AdjointStrainEnergyResponseFunction const& rOther);

    ///@}

}; // Class AdjointStrainEnergyResponseFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // ADJOINT_STRAIN_ENERGY_RESPONSE_FUNCTION_H
