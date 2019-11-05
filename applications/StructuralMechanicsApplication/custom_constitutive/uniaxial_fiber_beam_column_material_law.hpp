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

#if !defined(KRATOS_UNIAXIAL_FIBER_BEAM_COLUMN_MATERIAL_LAW_H_INCLUDED )
#define  KRATOS_UNIAXIAL_FIBER_BEAM_COLUMN_MATERIAL_LAW_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/properties.h"
#include "includes/define.h"
#include "includes/variables.h"

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

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) UniaxialFiberBeamColumnMaterialLaw
    // : public ConstitutiveLaw
{

public:

    ///@name Type Definitions
    ///@{


    typedef Properties PropertiesType;

    ///@}
    ///@name Pointer Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(UniaxialFiberBeamColumnMaterialLaw); //FIXME: This gives an error

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    UniaxialFiberBeamColumnMaterialLaw();
    UniaxialFiberBeamColumnMaterialLaw( PropertiesType::Pointer pProperties );

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    virtual void Initialize();
    virtual void CalculateMaterialResponse();
    virtual void FinalizeMaterialResponse();

    ///@}
    ///@name Access
    ///@{

    virtual double& GetTangentModulus();
    virtual double& GetStress();
    virtual void SetStrain(const double Strain);

    PropertiesType& GetProperties() { return *mpProperties; }
    PropertiesType const& GetProperties() const { return *mpProperties; }
    PropertiesType::Pointer pGetProperties() { return mpProperties; }
    const PropertiesType::Pointer pGetProperties() const { return mpProperties; }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;

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

    PropertiesType::Pointer mpProperties = nullptr;

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

    virtual void save(Serializer& rSerializer) const;
    virtual void load(Serializer& rSerializer);

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

};  // class UniaxialFiberBeamColumnMaterialLaw

/// output stream
inline std::ostream & operator <<(std::ostream& rOStream, const UniaxialFiberBeamColumnMaterialLaw& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos

#endif