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

#if !defined(KRATOS_FIBER_BEAM_COLUMN_SECTION_H_INCLUDED )
#define  KRATOS_FIBER_BEAM_COLUMN_SECTION_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/element.h"

#include "custom_classes/fiber_beam_column_component.h"
#include "custom_classes/fiber_beam_column_uniaxial_fiber.h"

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
 * @class FiberBeamColumnSection
 *
 * @brief A 3D section for fiber beam-column element
 * @details
 *
 * @author Mahmoud Zidan
 */

class KRATOS_API(SEISMIC_APPLICATION) FiberBeamColumnSection : public FiberBeamColumnComponent
{

public:

    ///@name Type Definitions
    ///@{

    typedef FiberBeamColumnComponent  BaseType;
    typedef BaseType::IndexType      IndexType;
    typedef BaseType::SizeType        SizeType;
    typedef BaseType::MatrixType    MatrixType;
    typedef BaseType::VectorType    VectorType;

    // by default: integration point is 3-dimensional
    typedef IntegrationPoint<3> IntegrationPointType;

    ///@}
    ///@name Pointer Definitions
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Default Constructor
     */
    FiberBeamColumnSection(IndexType NewId = 0);

    /**
     * Constructor
     * @param NewId: id of the section
     * @param integrationPoint: 1D integrationPoint that holds the position and the weight
     */
    FiberBeamColumnSection(IndexType NewId, IntegrationPointType integrationPoint);

    /**
     * Copy Constructor
     */
    FiberBeamColumnSection(FiberBeamColumnSection const& rOther) = default;

    /**
     * Destructor
     */
    ~FiberBeamColumnSection() = default;

    ///@}
    ///@name Operators
    ///@{

    /**
     * Assignment Operator
     */
    FiberBeamColumnSection& operator=(FiberBeamColumnSection const& rOther) = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * Initializes all the fibers and initializes the local flexibility matrix.
     */
    void Initialize() override;

    /**
     * This method is to get the global flexibility matrix
     * @param rGlobalFlexibilityMatrix: The global flexibility matrix
     */
    void GetGlobalFlexibilityMatrix(Matrix& rGlobalFlexibilityMatrix) override;

    /**
     * This method is to get the global deformation residuals
     * @param rGlobalResiduals: The global deformation residuals
     */
    void GetGlobalDeformationResiduals(Vector& rGlobalResiduals);

    /**
     * This method updates the forces, deformations (and fiber strains).
     * Then, updates the local flexibililty matrix, calculates the unbalance forces,
     * and checks for convergence.
     * @param rElementForceIncrements: the force increments of the element
     * @param Tolerance: tolerance for element equilibrium
     * @return True if norm of the unbalance force is below the tolerance
     */
    bool StateDetermination(const Vector& rElementForceIncrements, double Tolerance) override;

    /**
     * This method resets the deformation residuals to zero.
     */
    void ResetResidual();

    /**
     * Called when the non linear iterations have converged.
     */
    void FinalizeSolutionStep() override;

    ///@}
    ///@name Access
    ///@{

    /// number of fibers
    SizeType Size() { return mFibers.size(); }

    /// numerical integration weight
    double GetWeight() { return mWeight; }

    /// set the fibers of the section
    void SetFibers (std::vector<FiberBeamColumnUniaxialFiber> Fibers) { mFibers = std::move(Fibers); }

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

    double mPosition;    // position of the integration point
    double mWeight;      // weight of the integration point
    std::vector<FiberBeamColumnUniaxialFiber> mFibers;  // the vector the section's fibers

    Matrix mLocalFlexibilityMatrix = ZeroMatrix(3, 3);  // local flexibility matrix
    Vector mForces = ZeroVector(3);                     // forces calculated from the element
    Vector mUnbalanceForces = ZeroVector(3);            // unbalance forces of the section
    Vector mDeformationResiduals = ZeroVector(3);       // deformation residuals

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * This method updates the local flexibility matrix
     * of the section and stores it.
     */
    void UpdateLocalFlexibilityMatrix();

    /**
     * This method calculates the B matrix of the section.
     * @param rBMatrix B matrix
     */
    void CalculateBMatrix(Matrix& rBMatrix);

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

};  // class FiberBeamColumnSection

/// output stream
inline std::ostream & operator <<(std::ostream& rOStream, const FiberBeamColumnSection& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

}  // namespace Kratos


#endif
