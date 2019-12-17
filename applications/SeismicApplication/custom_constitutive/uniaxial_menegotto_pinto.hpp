//  KRATOS  ____       _               _
//         / ___|  ___(_)___ _ __ ___ (_) ___
//         \___ \ / _ \ / __| '_ ` _ \| |/ __|
//          ___) |  __/ \__ \ | | | | | | (__
//         |____/ \___|_|___/_| |_| |_|_|\___|
//
//  License:     BSD License
//  license:     structural_mechanics_application/license.txt
//
//  Main authors: Mahmoud Zidan
//

#if !defined(KRATOS_UNIAXIAL_MENGOTTO_PINTO_H_INCLUDED)
#define  KRATOS_UNIAXIAL_MENGOTTO_PINTO_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/constitutive_law.h"

namespace Kratos
{

///@}
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
///@name  Kratos Classes
///@{

/**
 * @class UniaxialMenegottoPintoMaterialLaw
 *
 * @brief A constitutive model for the steel uniaxial fibers of the beam-column element.
 * @details The constitutive law is a Menegotto-Pinto material law
 *          [M. Menegotto and P. E. Pinto. Method of Analysis for Cyclic Loaded R. C. Plane Frame
             Including Changes in Geometry and Non-Elastic Behaviour of Elements under Combined
             Normal Force and Bending, volume 11, pages 15–22. 1973.]
 *
 * @author Mahmoud Zidan
 */

class KRATOS_API(SEISMIC_APPLICATION) UniaxialMenegottoPintoMaterialLaw
    : public ConstitutiveLaw
{

public:

    ///@name Type Definitions
    ///@{

    typedef ConstitutiveLaw    BaseType;
    typedef ProcessInfo ProcessInfoType;
    typedef std::size_t        SizeType;

    ///@}
    ///@name Pointer Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(UniaxialMenegottoPintoMaterialLaw);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Default constructor
     */
    UniaxialMenegottoPintoMaterialLaw() = default;

    /**
     * @return pointer to the constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor
     */
    UniaxialMenegottoPintoMaterialLaw(const UniaxialMenegottoPintoMaterialLaw& rOther) = default;

    /**
     * Destructor
     */
    ~UniaxialMenegottoPintoMaterialLaw() override = default;

    ///@}
    ///@name Operators
    ///@{

    /**
     * Assignment operator
     */
    UniaxialMenegottoPintoMaterialLaw& operator=(const UniaxialMenegottoPintoMaterialLaw& rOther) = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * @return Voigt tensor size
     */
    SizeType GetStrainSize() override
    {
        return 1;
    }

    void InitializeMaterialResponsePK2 (Parameters& rValues) override;
    void CalculateMaterialResponsePK2 (Parameters& rValues) override;
    void FinalizeMaterialResponsePK2 (Parameters& rValues) override;

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
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

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

    // == history variables

    double mTangentModulus = 0.0;    // tangent modulus

    unsigned int mLoadingIndex = 0;  // loading index (0->initial 1->increasing 2->decreasing 3->unchanged)
    double mStrain0 = 0.0;           // strain 0
    double mStress0 = 0.0;           // stress 0
    double mStrainR = 0.0;           // strain reversal
    double mStressR = 0.0;           // stress reversal
    double mStrainPlastic = 0.0;     // plastic strain
    double mStrainMax = 0.0;         // max strain
    double mStrainMin = 0.0;         // min strain

    unsigned int mConvergedLoadingIndex = 0;
    double mConvergedStrain0 = 0.0;
    double mConvergedStress0 = 0.0;
    double mConvergedStrainR = 0.0;
    double mConvergedStressR = 0.0;
    double mConvergedStrainPlastic = 0.0;
    double mConvergedStrainMax = 0.0;
    double mConvergedStrainMin = 0.0;
    double mConvergedStrain = 0.0;
    double mConvergedStress = 0.0;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateStressResponsePK2 (Parameters& rValues, Vector& rStrainVector, Vector& rStressVector);
    void CalculateConstitutiveMatrix (Parameters& rValues, Matrix& rConstitutiveMatrix);

    ///@}
    ///@name Serialization
    ///@{

    friend Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

};  // class UniaxialMenegottoPintoMaterialLaw

/// output stream
inline std::ostream & operator <<(std::ostream& rOStream, const UniaxialMenegottoPintoMaterialLaw& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos

#endif  // #if !defined(KRATOS_UNIAXIAL_MENGOTTO_PINTO_H_INCLUDED)