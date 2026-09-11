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

#if !defined(KRATOS_FIBER_BEAM_COLUMN_COMPONENT_H_INCLUDED )
#define  KRATOS_FIBER_BEAM_COLUMN_COMPONENT_H_INCLUDED

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
 * @class FiberBeamColumnComponent
 *
 * @brief A base class for fiber beamColumn element
 *
 * @author Mahmoud Zidan
 */

class KRATOS_API(SEISMIC_APPLICATION) FiberBeamColumnComponent
{

public:

    ///@name Type Definitions
    ///@{

    typedef Element::GeometryType      GeometryType;
    typedef Element::PropertiesType  PropertiesType;
    typedef Element::IndexType            IndexType;
    typedef Element::SizeType              SizeType;
    typedef Element::MatrixType          MatrixType;
    typedef Element::VectorType          VectorType;

    ///@}
    ///@name Pointer Definitions
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Default Constructor
     */
    FiberBeamColumnComponent(IndexType NewId = 0);

    /**
     * Copy Constructor
     */
    FiberBeamColumnComponent(FiberBeamColumnComponent const& rOther) = default;

    /**
     * Destructor
     */
    ~FiberBeamColumnComponent() = default;

    ///@}
    ///@name Operators
    ///@{

    /**
     * Assignment Operator
     */
    FiberBeamColumnComponent& operator=(FiberBeamColumnComponent const& rOther) = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * Initialize.
     */
    virtual void Initialize();

    /**
     * This method is to get the global flexibility matrix
     * @param rGlobalFlexibilityMatrix: The global flexibility matrix
     */
    virtual void GetGlobalFlexibilityMatrix(Matrix& rGlobalFlexibilityMatrix);

    /**
     * This method is to get the global stiffness matrix
     * @param rGlobalStiffnessMatrix: The global stiffness matrix
     */
    virtual void GetGlobalStiffnessMatrix(Matrix& rGlobalStiffnessMatrix);

    /**
     * This method is to get the global internal forces
     * @param rGlobalForces: The global internal forces
     */
    virtual void GetGlobalInternalForces(Vector& rGlobalForces);

    /**
     * state determinatopm
     * @param rIncrements: force or deformation increments
     * @return True if in equilibrium
     */
    virtual bool StateDetermination(const Vector& rIncrements);

    /**
     * state determinatopm
     * @param rIncrements: force or deformation increments
     * @param Tolerance: tolerance for equilibrium
     * @return True if in equilibrium
     */
    virtual bool StateDetermination(const Vector& rIncrements, double Tolerance);

    /**
     * Called when the non linear iterations have converged.
     */
    virtual void FinalizeSolutionStep();

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

    IndexType const& Id() const { return mId; }

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

    IndexType mId;  // id of the component

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

};  // class FiberBeamColumnComponent

/// output stream
inline std::ostream & operator <<(std::ostream& rOStream, const FiberBeamColumnComponent& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

}  // namespace Kratos


#endif
