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

#if !defined(KRATOS_FIBER_BEAM_COLUMN_UNIAXIAL_FIBER_H_INCLUDED )
#define  KRATOS_FIBER_BEAM_COLUMN_UNIAXIAL_FIBER_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/element.h"

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

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) FiberBeamColumnUniaxialFiber
{

public:

    ///@name Type Definitions
    ///@{

    typedef Element::GeometryType      GeometryType;
    typedef Element::NodesArrayType  NodesArrayType;
    typedef Element::PropertiesType  PropertiesType;
    typedef Element::IndexType            IndexType;
    typedef Element::SizeType              SizeType;
    typedef Element::MatrixType          MatrixType;
    typedef Element::VectorType          VectorType;

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
    void Initialize();

    /**
     * This method updates the strain in the fiber, calculates the material response,
     * and stores the stress value.
     * @param rSectionDeformationIncrements: the deformation increments of the section
     */
    void StateDetermination(const Vector& rSectionDeformationIncrements);

    /**
     * returns the global stiffness matrix of the fiber.
     * @param rGlobalStiffnessMatrix
     */
    void CreateGlobalFiberStiffnessMatrix(Matrix& rGlobalStiffnessMatrix);

    /**
     * returns the global internal forces of the fiber.
     * @param rGlobalStiffnessMatrix
     */
    void CreateGlobalFiberInternalForces(Vector& rGlobalFiberInternalForces);

    /**
     * Called when the non linear iterations have converged
     */
    void FinalizeSolutionStep();

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

    IndexType mId;       // id of the fiber
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