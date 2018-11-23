// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Sergio Jiménez/Alejandro Cornejo/Lucia Barbu
//  Collaborator:    
//

#if !defined(KRATOS_GENERIC_SMALL_STRAIN_HIGH_CYCLE_FATIGUE_LAW_H_INCLUDED)
#define KRATOS_GENERIC_SMALL_STRAIN_HIGH_CYCLE_FATIGUE_LAW_H_INCLUDED
// System includes

// External includes

// Project includes
#include "custom_constitutive/generic_small_strain_isotropic_damage.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

// The size type definition
typedef std::size_t SizeType;
    
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
 * @class GenericSmallStrainHighCycleFatigueLaw
 * @ingroup StructuralMechanicsApplication
 * @brief This class is the base class which define all the constitutive laws for damage in small deformation
 * @details This class considers a constitutive law integrator as an intermediate utility to compute the damage
 * @tparam TConstLawIntegratorType The constitutive law integrator considered
 * @author Sergio Jiménez/Alejandro Cornejo/Lucia Barbu
 */
template <class TConstLawIntegratorType>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) GenericSmallStrainHighCycleFatigueLaw
    : public GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of GenericYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(GenericSmallStrainHighCycleFatigueLaw);

    /// The node definition
    typedef Node<3> NodeType;
    
    /// The geometry definition
    typedef Geometry<NodeType> GeometryType;
    
    /// Definition of the machine precision tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Default constructor.
    */
    GenericSmallStrainHighCycleFatigueLaw()
    {
    }

    /**
    * Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>>(*this);
    }

    /**
    * Copy constructor.
    */
    GenericSmallStrainHighCycleFatigueLaw(const GenericSmallStrainHighCycleFatigueLaw &rOther)
        : GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>(rOther),
          mFatigueReductionFactor(rOther.mFatigueReductionFactor)
    {
		//GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::GenericSmallStrainIsotropicDamage(rOther);
    }

    /**
    * Destructor.
    */
    ~GenericSmallStrainHighCycleFatigueLaw() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{


void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues);


void CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues);


void CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues);


void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues);

void FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType &rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo
    );

bool Has(const Variable<double>& rThisVariable);

void SetValue(const Variable<double>& rThisVariable, const double& rValue,
              const ProcessInfo& rCurrentProcessInfo);

double& GetValue(const Variable<double>& rThisVariable, double& rValue);

Matrix& CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    );

Vector& CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    );

double& CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    );

Matrix& GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    );

Vector& GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    );
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
    double GetFatigueReductionFactor() {return mFatigueReductionFactor;}
    void SetFatigueReductionFactor(const double toFred) {mFatigueReductionFactor = toFred;}

    Vector GetPreviousStresses() {return mPreviousStresses;}
    void SetPreviousStresses(const Vector toPreviousStresses) {mPreviousStresses = toPreviousStresses;}

    double GetMaxStress() {return mMaxStress;}
    void SetMaxStress(const double toMaxStress) {mMaxStress = toMaxStress;}

    double GetMinStress() {return mMinStress;}
    void SetMinStress(const double toMinStress) {mMinStress = toMinStress;}

    unsigned int GetNumberOfCycles() {return mNumberOfCycles;}
    void AddCycle() {mNumberOfCycles++;}
    void SetNumberOfCycles(const unsigned int toCycles) {mNumberOfCycles = toCycles;}

    double GetReversionFactor() {return mReversionFactor;}
    void SetReversionFactor(const double toReversionFactor) {mReversionFactor = toReversionFactor;}

    double GetFatigueReductionParameter() {return mFatigueReductionParameter;}
    void SetFatigueReductionParameter(const double toFatigueReductionParameter) {mFatigueReductionParameter = toFatigueReductionParameter;}

    Vector GetStressVector() {return mStressVector;}
    void SetStressVector(const Vector toStressVector) {mStressVector = toStressVector;}

    void ResetCycleCounter(){mHasCountedCycle = false;}
    void SetCycleCounter(const bool tocycle){mHasCountedCycle = tocycle;}
    bool GetCycleCounter() {return mHasCountedCycle;}

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
	double mFatigueReductionFactor = 1.0;
    Vector mPreviousStresses = ZeroVector(2); // [S_t-2, S_t-1]
    double mMaxStress = 0.0;
    double mMinStress = 0.0;
    unsigned int mNumberOfCycles = 0;
    double mReversionFactor = 0.0;  // = mMinStress/mMaxStress
    double mFatigueReductionParameter = 0.0; // B0
    Vector mStressVector = ZeroVector(VoigtSize);
    bool mHasCountedCycle = false;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


    ///@}

}; // Class GenericYieldSurface

} // namespace Kratos
#endif