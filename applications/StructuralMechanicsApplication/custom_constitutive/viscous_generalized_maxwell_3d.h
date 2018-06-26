// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo & Lucia Barbu 
//

#if !defined (KRATOS_VISCOUS_GENERALIZED_MAXWELL_H_INCLUDED)
#define  KRATOS_VISCOUS_GENERALIZED_MAXWELL_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/properties.h"
#include "utilities/math_utils.h"

#include "includes/constitutive_law.h"
#include "structural_mechanics_application_variables.h"

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
 * @class GenericConstitutiveLawIntegrator
 * @ingroup StructuralMechanicsApplication
 * @brief
 * @details
 * @author Alejandro Cornejo & Lucia Barbu
 */

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ViscousGeneralizedMaxwell3D
    : public ConstitutiveLaw
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of GenericYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(ViscousGeneralizedMaxwell3D);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Default constructor.
    */
    ViscousGeneralizedMaxwell3D()
    {
    }

    /**
    * Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override
    {
        ViscousGeneralizedMaxwell3D::Pointer p_clone
            (new ViscousGeneralizedMaxwell3D(*this));
        return p_clone;
    }

    /**
    * Copy constructor.
    */
    ViscousGeneralizedMaxwell3D (const ViscousGeneralizedMaxwell3D& rOther)
    : ConstitutiveLaw(rOther)
    {
    }
    /**
    * Destructor.
    */
    ~ViscousGeneralizedMaxwell3D() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override
    {
        return 3;
    };

    /**
     * @brief Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return 6;
    };

    int GetVoigtSize() {return 6;}
    int GetWorkingSpaceDimension() {return 3;}

    Vector GetPreviousStressVector() { return mPrevStressVector; }
    void SetPreviousStressVector(Vector toStress) { mPrevStressVector = toStress; }
    Vector GetNonConvPreviousStressVector() { return mNonConvPrevStressVector; }
    void SetNonConvPreviousStressVector(Vector toStress) { mNonConvPrevStressVector = toStress; }

    Vector GetPreviousStrainVector() { return mPrevStrainVector; }
    void SetPreviousStrainVector(Vector toStrain) { mPrevStrainVector = toStrain; }
    Vector GetNonConvPreviousStrainVector() { return mNonConvPrevStrainVector; }
    void SetNonConvPreviousStrainVector(Vector toStrain) { mNonConvPrevStrainVector = toStrain; }


    void CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override
    {
        this->CalculateMaterialResponseCauchy(rValues);
    }
    void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override
    {
        this->CalculateMaterialResponseCauchy(rValues);
    }
    void CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override
    {
        this->CalculateMaterialResponseCauchy(rValues);
    }

    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override
    {
        // Integrate Stress Damage
        const Properties& rMaterialProperties = rValues.GetMaterialProperties();
        const int VoigtSize = this->GetVoigtSize();
        Vector& IntegratedStressVector = rValues.GetStressVector(); // To be updated
        const Vector& StrainVector = rValues.GetStrainVector();
        Matrix& TangentTensor = rValues.GetConstitutiveMatrix(); // todo modify after integration
        const ProcessInfo& ProcessInfo = rValues.GetProcessInfo();
        const double TimeStep = ProcessInfo[DELTA_TIME];

        const double Kvisco    = rMaterialProperties[VISCOUS_PARAMETER]; //  C1/Cinf
        const double DelayTime = rMaterialProperties[DELAY_TIME];

        // Elastic Matrix
        Matrix C;
        this->CalculateElasticMatrix(C, rMaterialProperties);

        const Vector& PreviousStrain  = this->GetPreviousStrainVector();
        const Vector& PreviousStress  = this->GetPreviousStressVector();
        const Vector& StrainIncrement = StrainVector - PreviousStrain;

        const double coef = Kvisco * TimeStep / ((1 + Kvisco)*2.0*DelayTime);
        const Vector& Aux = -(StrainVector - StrainIncrement)*std::exp(-TimeStep/DelayTime)*(1 + coef)
                             + StrainVector*(1 - coef);
        

        noalias(IntegratedStressVector) = PreviousStress*std::exp(-TimeStep/DelayTime) + prod(C, Aux);

        this->SetNonConvPreviousStressVector(IntegratedStressVector);
        this->SetNonConvPreviousStrainVector(StrainVector);

    } // End CalculateMaterialResponseCauchy

    void FinalizeSolutionStep(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo
    ) override
    {
        // Update the required vectors
        this->SetPreviousStrainVector(this->GetNonConvPreviousStrainVector());
        this->SetPreviousStressVector(this->GetNonConvPreviousStressVector());
    }

    void CalculateElasticMatrix(
        Matrix &rElasticityTensor,
        const Properties &rMaterialProperties
    )
    {
        const double E = rMaterialProperties[YOUNG_MODULUS];
        const double poisson_ratio = rMaterialProperties[POISSON_RATIO];
        const double lambda =
            E * poisson_ratio / ((1. + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
        const double mu = E / (2.0 + 2.0 * poisson_ratio);

        if (rElasticityTensor.size1() != 6 || rElasticityTensor.size2() != 6)
            rElasticityTensor.resize(6, 6, false);
        rElasticityTensor.clear();

        rElasticityTensor(0, 0) = lambda + 2.0 * mu;
        rElasticityTensor(0, 1) = lambda;
        rElasticityTensor(0, 2) = lambda;
        rElasticityTensor(1, 0) = lambda;
        rElasticityTensor(1, 1) = lambda + 2.0 * mu;
        rElasticityTensor(1, 2) = lambda;
        rElasticityTensor(2, 0) = lambda;
        rElasticityTensor(2, 1) = lambda;
        rElasticityTensor(2, 2) = lambda + 2.0 * mu;
        rElasticityTensor(3, 3) = mu;
        rElasticityTensor(4, 4) = mu;
        rElasticityTensor(5, 5) = mu;
    }

    void FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
    {
    }
    void FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
    {
    }
    void FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
    {
    }
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
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

    // Converged values
    Vector mPrevStressVector = ZeroVector(6);
    Vector mPrevStrainVector = ZeroVector(6);
    
    // Non Converged values
    Vector mNonConvPrevStressVector = ZeroVector(6);
    Vector mNonConvPrevStrainVector = ZeroVector(6);


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

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
        rSerializer.save("PrevStressVector", mPrevStressVector);
        rSerializer.save("PrevStrainVector", mPrevStrainVector);
        rSerializer.save("NonConvPrevStressVector", mNonConvPrevStressVector);
        rSerializer.save("NonConvPrevStrainVector", mNonConvPrevStrainVector);

    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw)
        rSerializer.load("PrevStressVector", mPrevStressVector);
        rSerializer.load("PrevStrainVector", mPrevStrainVector);
        rSerializer.load("NonConvPrevStressVector", mNonConvPrevStressVector);
        rSerializer.load("NonConvPrevStrainVector", mNonConvPrevStrainVector);
    }

    ///@}

}; // Class GenericYieldSurface

} // namespace kratos
#endif
