// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:   Alejandro Cornejo & Lucia Barbu
//

#if !defined(KRATOS_TANGENT_OPERATOR_CALCULATOR_PROCESS_H_INCLUDED )
#define  KRATOS_TANGENT_OPERATOR_CALCULATOR_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/constitutive_law.h"

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

/**
 * @class TangentOperatorCalculatorProcess
 * @ingroup StructuralMechanicsApplication
 * @brief An algorithm that derives numerically the constitutive tangent tensor at one GP
 * @details The procedure is defined in the PAPER "Caracterización de la delaminación en materiales
 compuestos mediante la teoría de mezclas serie/paralelo" X. Martinez, S. Oller y E. Barbero.
 * @author Alejandro Cornejo & Lucia Barbu
 */

template <class TConstitutiveLawType>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) TangentOperatorCalculatorProcess
    : public Process
{
public:

    /// Pointer definition of TangentOperatorCalculatorProcess
    KRATOS_CLASS_POINTER_DEFINITION(TangentOperatorCalculatorProcess);

    /// Constructor
    TangentOperatorCalculatorProcess(
        ConstitutiveLaw::Parameters& rValues
    )
    {
        mValues = rValues;
    }

    /// Destructor.
    virtual ~TangentOperatorCalculatorProcess() {}

    void Execute()
    {
        Matrix C;
        TConstitutiveLawType::CalculateElasticMatrix(C, mMaterialProperties);

        const Vector& StrainVectorGP = mValues.GetStrainVector();
        const Vector& StressVectorGP = mValues.GetStressVector();

        Matrix& TangentTensor = mValues.GetConstitutiveMatrix();
        TangentTensor.clear();

        const int NumComp = StrainVectorGP.size();

        // Loop over components of the strain
        for (int Component = 0; Component < NumComp; Component++)
        {
            ConstitutiveLaw::Parameters PerturbedValues = mValues;
            Vector& PerturbedStrain = PerturbedValues.GetStrainVector();
            

            double Perturbation;
            this->CalculatePerturbation(PerturbedStrain, Component, Perturbation);
            this->PerturbateStrainVector(PerturbedStrain, StrainVectorGP, Perturbation, Component);
            this->IntegratePerturbedStrain(PerturbedValues);

            Vector& PerturbedIntegratedStress = PerturbedValues.GetStressVector(); // now integrated
            const Vector& DeltaStress = PerturbedIntegratedStress - StressVectorGP; 
            this->AssignComponentsToTangentTensor(TangentTensor, DeltaStress, Perturbation, Component);
        }
    }

    void CalculatePerturbation(
        const Vector& StrainVector, 
        const int Component,
        double& Perturbation
    )
    {
        double Pert1, Pert2;
        if (StrainVector[Component] != 0.0)
        {
            Pert1 = 1e-5*StrainVector[Component];
        }
        else
        {
            double MinStrainComp;
            this->GetMinAbsValue(StrainVector, MinStrainComp);
            Pert1 = 1e-5*MinStrainComp;
        }
        double MaxStrainComp;
        this->GetMaxAbsValue(StrainVector, MaxStrainComp);
        Pert2 = 1e-10*MaxStrainComp;

        Perturbation = std::max(Pert1, Pert2);
    }

    void PerturbateStrainVector(
        Vector& PerturbedStrainVector, 
        const Vector& StrainVectorGP,
        const double Perturbation,
        const int Component
    )
    {
        PerturbedStrainVector = StrainVectorGP;
        PerturbedStrainVector[Component] += Perturbation;
    }

    void IntegratePerturbedStrain(ConstitutiveLaw::Parameters& rValues)
    {
        TConstitutiveLawType::CalculateMaterialResponseCauchy(rValues);
    }

    void GetMaxAbsValue(
        const Vector& ArrayValues, 
        double& MaxValue
    )
    {
        const int Dim = ArrayValues.size();
        double aux = std::abs(ArrayValues[0]);

        for (int i = 1; i < Dim; i++)
        {
            if (std::abs(ArrayValues[i]) > aux) aux = std::abs(ArrayValues[i]);
        }
    }

    void GetMinAbsValue(
        const Vector& ArrayValues, 
        double& MaxValue
    )
    {
        const int Dim = ArrayValues.size();
        double aux = std::abs(ArrayValues[0]);

        for (int i = 1; i < Dim; i++)
        {
            if (std::abs(ArrayValues[i]) < aux) aux = std::abs(ArrayValues[i]);
        }
    }

    void AssignComponentsToTangentTensor(
        Matrix& TangentTensor, 
        const Vector& DeltaStress,
        const double Perturbation,
        const int Component
    )
    {
        const int Dim = DeltaStress.size();
        for (int row = 0; row < Dim; row++)
        {
            TangentTensor(row, Component) = DeltaStress[row] / Perturbation;
        }
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
    // Matrix& mTangentTensor; 
    // const Vector mStressVectorGP;
    // const Vector mStrainVectorGP;
    // const Properties& mMaterialProperties;
    ConstitutiveLaw::Parameters& mValues;

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
    TangentOperatorCalculatorProcess& operator=(TangentOperatorCalculatorProcess const& rOther);

};
}// namespace Kratos.
#endif // KRATOS_TANGENT_OPERATOR_CALCULATOR_PROCESS_H_INCLUDED  defined