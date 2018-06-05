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

#if !defined (KRATOS_VISCOUS_GENERALIZED_KELVIN_H_INCLUDED)
#define  KRATOS_VISCOUS_GENERALIZED_KELVIN_H_INCLUDED

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

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ViscousGeneralizedKelvin3D
    : public ConstitutiveLaw
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of GenericYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(ViscousGeneralizedKelvin3D);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Default constructor.
    */
    ViscousGeneralizedKelvin3D()
    {
    }

    /**
    * Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override
    {
        ViscousGeneralizedKelvin3D::Pointer p_clone
            (new ViscousGeneralizedKelvin3D(*this));
        return p_clone;
    }

    /**
    * Copy constructor.
    */
    ViscousGeneralizedKelvin3D (const ViscousGeneralizedKelvin3D& rOther)
    : ConstitutiveLaw(rOther)
    {
    }
    /**
    * Destructor.
    */
    ~ViscousGeneralizedKelvin3D() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    int GetVoigtSize(){return 6;}
    int GetWorkingSpaceDimension() {return 3;}

    Vector GetPreviousStressVector() { return mPrevStressVector; }
    void SetPreviousStressVector(Vector toStress) { mPrevStressVector = toStress; }
    Vector GetNonConvPreviousStressVector() { return mNonConvPrevStressVector; }
    void SetNonConvPreviousStressVector(Vector toStress) { mNonConvPrevStressVector = toStress; }

    Vector GetPreviousInelasticStrainVector() { return mPrevInelasticStrainVector; }
    void SetPreviousInelasticStrainVector(Vector toStrain) { mPrevInelasticStrainVector = toStrain; }
    Vector GetNonConvPreviousInelasticStrainVector() { return mNonConvPrevInelasticStrainVector; }
    void SetNonConvPreviousInelasticStrainVector(Vector toStrain) { mNonConvPrevInelasticStrainVector = toStrain; }


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

        const double Kvisco    = rMaterialProperties[VISCOUS_PARAMETER]; // C1/Cinf
        const double DelayTime = rMaterialProperties[DELAY_TIME];

        // Elastic Matrix
        Matrix C;
        this->CalculateElasticMatrix(C, rMaterialProperties);

        Vector InelasticStrainVector = this->GetPreviousInelasticStrainVector();
        const Vector& PreviousStress  = this->GetPreviousStressVector();

        const int NumberOfSubIncrements = 10;
        const double dt = TimeStep / NumberOfSubIncrements;

        Vector AuxStressVector;
        noalias(AuxStressVector) = PreviousStress;
        BoundedVector<double, 6> Aux = ZeroVector(6);

        for (int i = 0; i < NumberOfSubIncrements; i++)
        {
            noalias(Aux) = std::exp(-dt/DelayTime)*prod(C, AuxStressVector) / DelayTime;
            InelasticStrainVector = std::exp(-dt/DelayTime)*InelasticStrainVector + Aux;
            noalias(AuxStressVector) = prod(C, StrainVector - InelasticStrainVector);
        }

        noalias(IntegratedStressVector) = AuxStressVector;
        
        this->SetNonConvPreviousStressVector(IntegratedStressVector);
        this->SetNonConvPreviousInelasticStrainVector(InelasticStrainVector);

    } // End CalculateMaterialResponseCauchy

    // void CalculateTangentTensor(Matrix& C) 
    // {

    // }

    void FinalizeSolutionStep(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo
    ) override
    {
        // Update the required vectors
        this->SetPreviousInelasticStrainVector(this->GetNonConvPreviousInelasticStrainVector());
        this->SetPreviousStressVector(this->GetNonConvPreviousStressVector());
    }

    void CalculateElasticMatrix(Matrix &rElasticityTensor,
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
    Vector mPrevInelasticStrainVector = ZeroVector(6);
    
    // Non Converged values
    Vector mNonConvPrevStressVector = ZeroVector(6);
    Vector mNonConvPrevInelasticStrainVector = ZeroVector(6);


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
        rSerializer.save("PrevInelasticStrainVector", mPrevInelasticStrainVector);
        rSerializer.save("NonConvPrevStressVector", mNonConvPrevStressVector);
        rSerializer.save("NonConvPrevInelasticStrainVector", mNonConvPrevInelasticStrainVector);

    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw)
        rSerializer.load("PrevStressVector", mPrevStressVector);
        rSerializer.load("PrevInelasticStrainVector", mPrevInelasticStrainVector);
        rSerializer.load("NonConvPrevStressVector", mNonConvPrevStressVector);
        rSerializer.load("NonConvPrevInelasticStrainVector", mNonConvPrevInelasticStrainVector);
    }

    ///@}

}; // Class GenericYieldSurface

} // namespace kratos
#endif
