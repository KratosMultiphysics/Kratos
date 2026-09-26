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
//    Co authors: Long Chen
//

#if !defined(KRATOS_FIBER_BEAM_COLUMN_UNIAXIAL_FIBER_H_INCLUDED )
#define  KRATOS_FIBER_BEAM_COLUMN_UNIAXIAL_FIBER_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/element.h"

#include "custom_classes/fiber_beam_column_component.h"

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
 * @class FiberBeamColumnUniaxialFiber
 *
 * @brief A 3D unaxial fiber for the beam-column element for reinforced concrete modeling
 * @details it is stored in the section, in a container of fibers
 *
 * @author Mahmoud Zidan
 */

class KRATOS_API(SEISMIC_APPLICATION) FiberBeamColumnUniaxialFiber : public FiberBeamColumnComponent
{

public:

    ///@name Type Definitions
    ///@{

    typedef FiberBeamColumnComponent       BaseType;
    typedef BaseType::GeometryType     GeometryType;
    typedef BaseType::PropertiesType PropertiesType;
    typedef BaseType::IndexType           IndexType;
    typedef BaseType::SizeType             SizeType;
    typedef BaseType::MatrixType         MatrixType;
    typedef BaseType::VectorType         VectorType;

    ///@}
    ///@name Pointer Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(FiberBeamColumnUniaxialFiber);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Default Constructor
     */
    FiberBeamColumnUniaxialFiber(IndexType NewId = 0);

    /**
     * Constructor
     * @param NewId: id of the fiber
     * @param Y: y coordinate
     * @param Z: z coordinate
     * @param pMaterial: pointer to the constitutive law
     * @param pProperties: pointer to the properties
     */
    FiberBeamColumnUniaxialFiber(
        IndexType NewId,
        double Y,
        double Z,
        double Area,
        ConstitutiveLaw::Pointer pMaterial,
        PropertiesType::Pointer pProperties
    );

    /**
     * Destructor
     */
    ~FiberBeamColumnUniaxialFiber() = default;

    /**
     * Copy Constructor
     */
    FiberBeamColumnUniaxialFiber(const FiberBeamColumnUniaxialFiber& rOther) = default;

    ///@}
    ///@name Operators
    ///@{

    /**
     * Assignement operator
     */
    FiberBeamColumnUniaxialFiber& operator=(const FiberBeamColumnUniaxialFiber& rOther) = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * Initializes the material law
     */
    void Initialize() override;

    /**
     * This method updates the strain in the fiber, calculates the material response,
     * and stores the stress value.
     * @param rSectionDeformationIncrements: the deformation increments of the section
     */
    bool StateDetermination(const Vector& rSectionDeformationIncrements) override;

    /**
     * returns the global stiffness matrix of the fiber.
     * @param rGlobalStiffnessMatrix
     */
    void GetGlobalStiffnessMatrix(Matrix& rGlobalStiffnessMatrix) override;

    /**
     * returns the global internal forces of the fiber.
     * @param rGlobalStiffnessMatrix
     */
    void GetGlobalInternalForces(Vector& rGlobalFiberInternalForces) override;

    /**
     * Called when the non linear iterations have converged
     */
    void FinalizeSolutionStep() override;

    ///@}
    ///@name Access
    ///@{

    double const& GetStress() const {return mStress[0];}
    double const& GetStrain() const {return mStrain[0];}

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;

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

    double mArea = 0;    // area of the fiber
    Vector mTransformationVector = ZeroVector(3);   // coordinates [-y z 1]
    ConstitutiveLaw::Pointer mpMaterial = nullptr;  // pointer to the constitutive law
    PropertiesType::Pointer mpProperties = nullptr; // pointer to the properties
    Vector mStrain = ZeroVector(1);  // strain in the fiber
    Vector mStress = ZeroVector(1);  // stress in the fiber

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
    void save(Serializer& rSerializer) const;
    void load(Serializer& rSerializer);

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

};  // class FiberBeamColumnUniaxialFiber

/// output stream
inline std::ostream & operator <<(std::ostream& rOStream, const FiberBeamColumnUniaxialFiber& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos

#endif