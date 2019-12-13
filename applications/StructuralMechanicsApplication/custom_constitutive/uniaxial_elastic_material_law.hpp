// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//  license:     structural_mechanics_application/license.txt
//
//  Main authors: Mahmoud Zidan
//

#if !defined(KRATOS_UNIAXIAL_ELASTIC_MATERIAL_LAW_H_INCLUDED )
#define  KRATOS_UNIAXIAL_ELASTIC_MATERIAL_LAW_H_INCLUDED

// System includes


// External includes


// Project includes
// #include "includes/define.h"
// #include "includes/variables.h"
// #include "structural_mechanics_application_variables.h"
// #include "custom_constitutive/uniaxial_fiber_beam_column_material_law.hpp"
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
 * @class UniaxialElasticMaterialLaw
 *
 * @brief A constitutive model for the steel uniaxial fibers of the beam-column element.
 * @details The constitutive law is a Menegotto-Pinto material law
 *
 * @author Mahmoud Zidan
 */

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) UniaxialElasticMaterialLaw
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

    KRATOS_CLASS_POINTER_DEFINITION(UniaxialElasticMaterialLaw);

    ///@}
    ///@name Life Cycle
    ///@{

    UniaxialElasticMaterialLaw();
    ConstitutiveLaw::Pointer Clone() const override;
    UniaxialElasticMaterialLaw(const UniaxialElasticMaterialLaw& rOther);
    ~UniaxialElasticMaterialLaw() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return 1;
    }

    // /**
    //  * This function provides the place to perform checks on the completeness of the input.
    //  * It is designed to be called only once (or anyway, not often) typically at the beginning
    //  * of the calculations, so to verify that nothing is missing from the input
    //  * or that no common error is found.
    //  * @param rMaterialProperties: The properties of the material
    //  * @param rElementGeometry: The geometry of the element
    //  * @param rCurrentProcessInfo: The current process info instance
    //  */
    // int Check(
    //     const Properties& rMaterialProperties,
    //     const GeometryType& rElementGeometry,
    //     const ProcessInfo& rCurrentProcessInfo
    // ) override;

    // array_1d<double, 3 > & GetValue(const Variable<array_1d<double, 3 > >& rThisVariable,
    //     array_1d<double, 3 > & rValue) override;

    // double& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
    //     const Variable<double>& rThisVariable,double& rValue) override;

    // Vector& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
    //     const Variable<Vector>& rThisVariable,
    //     Vector& rValue) override;

    // array_1d<double, 3 > & CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
    //     const Variable<array_1d<double, 3 > >& rVariable,
    //     array_1d<double, 3 > & rValue) override;

    // /**
    //  * @brief Returns the value of a specified variable (double)
    //  * @param rThisVariable the variable to be returned
    //  * @param rValue a reference to the returned value
    //  * @return rValue output: the value of the specified variable
    //  */
    // double& GetValue(const Variable<double>& rThisVariable, double& rValue) override;

    void CalculateMaterialResponsePK2(Parameters& rValues) override;

    void FinalizeMaterialResponsePK2(Parameters& rValues) override;

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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

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

};  // class UniaxialElasticMaterialLaw

/// output stream
inline std::ostream & operator <<(std::ostream& rOStream, const UniaxialElasticMaterialLaw& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos

#endif