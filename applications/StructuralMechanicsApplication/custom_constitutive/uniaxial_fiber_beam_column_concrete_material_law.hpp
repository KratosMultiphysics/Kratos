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

#if !defined(KRATOS_UNIAXIAL_FIBER_BEAM_COLUMN_CONCRETE_MATERIAL_LAW_H_INCLUDED )
#define  KRATOS_UNIAXIAL_FIBER_BEAM_COLUMN_CONCRETE_MATERIAL_LAW_H_INCLUDED

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
 * @class FiberBeamColumnElement3D2N
 *
 * @brief A 3D-2node fiber beam-column element for reinforced concrete modeling
 *
 * @author Mahmoud Zidan
 */

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) UniaxialFiberBeamColumnConcreteMaterialLaw
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

    KRATOS_CLASS_POINTER_DEFINITION(UniaxialFiberBeamColumnConcreteMaterialLaw); //FIXME: This gives an error

    ///@}
    ///@name Life Cycle
    ///@{

    UniaxialFiberBeamColumnConcreteMaterialLaw();
    UniaxialFiberBeamColumnConcreteMaterialLaw(PropertiesType::Pointer pProperties);

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

};  // class UniaxialFiberBeamColumnConcreteMaterialLaw

/// output stream
inline std::ostream & operator <<(std::ostream& rOStream, const UniaxialFiberBeamColumnConcreteMaterialLaw& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos

#endif