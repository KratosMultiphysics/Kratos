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

//**************
// Project includes
#include "includes/process_info.h"
#include "testing/testing.h"

// Application includes

// Integrator
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_damage.h"
// Yield surfaces
#include "custom_constitutive/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/yield_surfaces/rankine_yield_surface.h"
#include "custom_constitutive/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_constitutive/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_constitutive/yield_surfaces/tresca_yield_surface.h"
// Plastic potentials
#include "custom_constitutive/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/tresca_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/drucker_prager_plastic_potential.h"
// Constitutive law
#include "custom_constitutive/generic_small_strain_isotropic_damage_3d.h"
#include "includes/model_part.h"
#include "geometries/tetrahedra_3d_4.h"
//*****************

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

		this->test();

        // Elastic Matrix
        Matrix C, InvC;
        double detC = 0.0;
        this->CalculateElasticMatrix(C, rMaterialProperties);
        MathUtils<double>::InvertMatrix( C, InvC, detC ); 

        Vector InelasticStrainVector  = this->GetPreviousInelasticStrainVector();
        const Vector& PreviousStress  = this->GetPreviousStressVector();

        const int NumberOfSubIncrements = 10;
        const double dt = TimeStep / NumberOfSubIncrements;

        Vector AuxStressVector;
        AuxStressVector = PreviousStress;
        Vector Aux = ZeroVector(6);

        Vector ElasticStrain;
        for (int i = 0; i < NumberOfSubIncrements; i++) {
			Aux = std::exp(-dt/DelayTime) * prod(InvC, AuxStressVector) / DelayTime;
            InelasticStrainVector = std::exp(-dt/DelayTime)*InelasticStrainVector + Aux;
            ElasticStrain = StrainVector - InelasticStrainVector;
            noalias(AuxStressVector) = prod(C, ElasticStrain);
        }

        noalias(IntegratedStressVector) = AuxStressVector;
        noalias(TangentTensor) = C;
        
        this->SetNonConvPreviousStressVector(IntegratedStressVector);
        this->SetNonConvPreviousInelasticStrainVector(InelasticStrainVector);

    } // End CalculateMaterialResponseCauchy

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

	void test()
	{
		typedef Node < 3 > NodeType;

		typedef GenericSmallStrainIsotropicDamage3D
			<GenericConstitutiveLawIntegratorDamage
				<ModifiedMohrCoulombYieldSurface
					<ModifiedMohrCoulombPlasticPotential>>> MC;

		typedef GenericSmallStrainIsotropicDamage3D
			<GenericConstitutiveLawIntegratorDamage
				<VonMisesYieldSurface
					<ModifiedMohrCoulombPlasticPotential>>> VM;

		typedef GenericSmallStrainIsotropicDamage3D
			<GenericConstitutiveLawIntegratorDamage
				<DruckerPragerYieldSurface
					<ModifiedMohrCoulombPlasticPotential>>> DP;

		typedef GenericSmallStrainIsotropicDamage3D
			<GenericConstitutiveLawIntegratorDamage
				<TrescaYieldSurface
					<ModifiedMohrCoulombPlasticPotential>>> T;

		typedef GenericSmallStrainIsotropicDamage3D
			<GenericConstitutiveLawIntegratorDamage
				<SimoJuYieldSurface
					<ModifiedMohrCoulombPlasticPotential>>> SJ;

		typedef GenericSmallStrainIsotropicDamage3D
			<GenericConstitutiveLawIntegratorDamage
				<RankineYieldSurface
					<ModifiedMohrCoulombPlasticPotential>>> R;

		ConstitutiveLaw::Parameters rValues;
		Properties rMaterialProperties;
		Vector rStressVector, rStrainVector;

		ModelPart& TestMdpa = ModelPart();

		NodeType::Pointer Node1 = TestMdpa.CreateNewNode(1, 0.0, 0.0, 0.0);
		NodeType::Pointer Node2 = TestMdpa.CreateNewNode(2, 1.0, 0.0, 0.0);
		NodeType::Pointer Node3 = TestMdpa.CreateNewNode(3, 0.0, 1.0, 0.0);
		NodeType::Pointer Node4 = TestMdpa.CreateNewNode(4, 0.0, 0.0, 1.0);

		Tetrahedra3D4<NodeType> Geom = Tetrahedra3D4<NodeType>(Node1, Node2, Node3, Node4);

		rStressVector = ZeroVector(6);
		rStressVector[0] = 5.40984e+06;
		rStressVector[1] = 5.40984e+06;
		rStressVector[2] = 1.91803e+07;
		rStressVector[3] = 0.0;
		rStressVector[4] = 0.0;
		rStressVector[5] = 1.45804e-10;

		rStrainVector = ZeroVector(6);
		rStrainVector[0] = 0.0;
		rStrainVector[1] = 0.0;
		rStrainVector[2] = 8.0e-5;
		rStrainVector[3] = 0.0;
		rStrainVector[4] = 0.0;
		rStrainVector[5] = 1.6941e-21;

		rMaterialProperties.SetValue(YOUNG_MODULUS, 210e9);
		rMaterialProperties.SetValue(POISSON_RATIO, 0.22);
		rMaterialProperties.SetValue(YIELD_STRESS_COMPRESSION, 3.0e6);
		rMaterialProperties.SetValue(YIELD_STRESS_TENSION, 3.0e6);
		rMaterialProperties.SetValue(FRICTION_ANGLE, 32.0);
		rMaterialProperties.SetValue(DILATANCY_ANGLE, 16.0);
		rMaterialProperties.SetValue(FRACTURE_ENERGY, 1e5);
		rMaterialProperties.SetValue(SOFTENING_TYPE, 1);

		rValues.SetElementGeometry(Geom);
		rValues.SetMaterialProperties(rMaterialProperties);
		rValues.SetStrainVector(rStrainVector);
		rValues.SetStressVector(rStressVector);

		// Create the CL's
		MC MohrCoulombCL = MC();
		VM VonMisesCL = VM();
		DP DruckerPragerCL = DP();
		T TrescaCL = T();
		R RankineCL = R();
		SJ SimoJuCL = SJ();

		std::vector<double> MCres, VMres, DPres, Tres, Rres, SJres;
		MCres = { -9.07262e+06, -9.07262e+06, -1.18548e+07, 0.0, 0.0, -2.94576e-11 };
		VMres = { -9.09508e+06, -9.09508e+06, -1.18098e+07, 0.0, 0.0, -2.87441e-11 };
		DPres = { -5.40984e+06, -5.40984e+06, -1.91803e+07, 0.0, 0.0, -1.45804e-10 };
		Tres = { -9.09508e+06, -9.09508e+06, -1.18098e+07, 0.0, 0.0, -2.87441e-11 };

		Vector TestMC, TestVM, TestDP, TestT, TestR, TestSJ;
		MohrCoulombCL.CalculateMaterialResponseCauchy(rValues);
		TestMC = rValues.GetStressVector();
		KRATOS_WATCH(TestMC)

		VonMisesCL.CalculateMaterialResponseCauchy(rValues);
		TestVM = rValues.GetStressVector();
		KRATOS_WATCH(TestVM)

		DruckerPragerCL.CalculateMaterialResponseCauchy(rValues);
		TestDP = rValues.GetStressVector();
		KRATOS_WATCH(TestDP)

		TrescaCL.CalculateMaterialResponseCauchy(rValues);
		TestT = rValues.GetStressVector();
		KRATOS_WATCH(TestT)

		RankineCL.CalculateMaterialResponseCauchy(rValues);
		TestR = rValues.GetStressVector();
		KRATOS_WATCH(TestR)

        rMaterialProperties.SetValue(FRACTURE_ENERGY, 1.0e15);
		rValues.SetMaterialProperties(rMaterialProperties);
		SimoJuCL.CalculateMaterialResponseCauchy(rValues);
		TestSJ = rValues.GetStressVector();
		KRATOS_WATCH(TestSJ)
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
