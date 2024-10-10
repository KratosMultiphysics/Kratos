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

#if !defined(KRATOS_UNIAXIAL_KENT_PARK_H_INCLUDED)
#define  KRATOS_UNIAXIAL_KENT_PARK_H_INCLUDED

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
 * @class UniaxialKentParkMaterialLaw
 *
 * @brief A constitutive model for the concrete uniaxial fibers of the beam-column element.
 * @details The constitutive law is a Kent-Park material law
 *          [D. C. Kent. Inelastic Behavior of Reinforced Concrete Members with Cyclic Loading.
 *           PhD thesis, University of Canterbury, 1969.]
 *
 * @author Mahmoud Zidan
 */

class KRATOS_API(SEISMIC_APPLICATION) UniaxialKentParkMaterialLaw
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

    KRATOS_CLASS_POINTER_DEFINITION(UniaxialKentParkMaterialLaw);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Default constructor
     */
    UniaxialKentParkMaterialLaw() = default;

    /**
     * @return pointer to the constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor
     */
    UniaxialKentParkMaterialLaw(const UniaxialKentParkMaterialLaw& rOther) = default;

    /**
     * Destructor
     */
    ~UniaxialKentParkMaterialLaw() override = default;

    ///@}
    ///@name Operators
    ///@{

    /**
     * Assignment Operator
     */
    UniaxialKentParkMaterialLaw& operator=(const UniaxialKentParkMaterialLaw& rOther) = default;

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

    // material history variables
    double mStrainR        = 0.0;  // strain reversal
    double mStrainP        = 0.0;  // plastic strain
    double mUnloadSlope    = 0.0;  // unload slope
    double mTangentModulus = 0.0;  // tangent modulus

    // material converged history variables
    double mConvergedStress         = 0.0;
    double mConvergedStrain         = 0.0;
    double mConvergedStrainR        = 0.0;
    double mConvergedStrainP        = 0.0;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateStressResponsePK2 (Parameters& rValues, Vector& rStrainVector, Vector& rStressVector);
    void CalculateConstitutiveMatrixPK2 (Parameters& rValues, Matrix& rConstitutiveMatrix);
    void Reload(Parameters& rValues, Vector& rStrainVector, Vector& rStressVector);
    void Envelope(Parameters& rValues, Vector& rStrainVector, Vector& rStressVector);
    void Unload(Parameters& rValues, Vector& rStrainVector, Vector& rStressVector);

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

};  // class UniaxialKentParkMaterialLaw

/// output stream
inline std::ostream & operator <<(std::ostream& rOStream, const UniaxialKentParkMaterialLaw& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos

#endif  // #if !defined(KRATOS_UNIAXIAL_KENT_PARK_H_INCLUDED)