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

#if !defined(KRATOS_UNIAXIAL_FIBER_BEAM_COLUMN_STEEL_MATERIAL_LAW_H_INCLUDED )
#define  KRATOS_UNIAXIAL_FIBER_BEAM_COLUMN_STEEL_MATERIAL_LAW_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "structural_mechanics_application_variables.h"
#include "custom_constitutive/uniaxial_fiber_beam_column_material_law.hpp"

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
 * @class UniaxialFiberBeamColumnSteelMaterialLaw
 *
 * @brief A constitutive model for the steel uniaxial fibers of the beam-column element.
 * @details The constitutive law is a Menegotto-Pinto material law
 *
 * @author Mahmoud Zidan
 */

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) UniaxialFiberBeamColumnSteelMaterialLaw
    : public UniaxialFiberBeamColumnMaterialLaw
{

public:

    ///@name Type Definitions
    ///@{

    typedef UniaxialFiberBeamColumnMaterialLaw          BaseType;
    typedef BaseType::PropertiesType              PropertiesType;
    typedef std::size_t                                 SizeType;

    ///@}
    ///@name Pointer Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(UniaxialFiberBeamColumnSteelMaterialLaw);

    ///@}
    ///@name Life Cycle
    ///@{

    UniaxialFiberBeamColumnSteelMaterialLaw();
    UniaxialFiberBeamColumnSteelMaterialLaw(PropertiesType::Pointer pProperties);

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override;
    void CalculateMaterialResponse() override;
    void FinalizeMaterialResponse() override;

    ///@}
    ///@name Access
    ///@{

    void SetStrain(const double Strain) override { mStrain = Strain; }

    double& GetStress() override         { return mStress; }
    double& GetTangentModulus() override { return mTangentModulus; }

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

    double mStress;
    double mStrain;
    double mTangentModulus;

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

    unsigned int mLoadingIndex = 0;
    double mStrain0 = 0.0;
    double mStress0 = 0.0;
    double mStrainR = 0.0;
    double mStressR = 0.0;
    double mStrainPlastic = 0.0;
    double mStrainMax = 0.0;
    double mStrainMin = 0.0;

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

};  // class UniaxialFiberBeamColumnSteelMaterialLaw

/// output stream
inline std::ostream & operator <<(std::ostream& rOStream, const UniaxialFiberBeamColumnSteelMaterialLaw& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos

#endif