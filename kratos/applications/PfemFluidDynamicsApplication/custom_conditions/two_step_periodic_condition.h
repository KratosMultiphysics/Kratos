
#ifndef KRATOS_TWO_STEP_PERIODIC_CONDITION_H
#define	KRATOS_TWO_STEP_PERIODIC_CONDITION_H

// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "geometries/geometry.h"
#include "includes/properties.h"
#include "includes/process_info.h"
#include "utilities/indexed_object.h"
#include "includes/condition.h"
#include "includes/serializer.h"

// Application includes
/* #include "pfem_fluid_dynamics_application.h" */

namespace Kratos
{
///@addtogroup PfemFluidDynamicsApplication
///@{

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
///@name Kratos Classes
///@{

/// Condition used to assign periodic boundary conditions when using a fractional step fluid solver.
/* /\** */
/*  * A conformant mapping between the two related boundaries is assumed, so this condition */
/*  * will link two nodes, one on each boundary, which are considered to be "images" of one */
/*  * another on both periodic sides.\n */
/*  * Note that the periodic boundary condition is enforced by the builder and solver, this Condition */
/*  * is simply used to provide it the required information to set it up.\n */
/*  * @note TwoStepPeriodicCondition is designed designed to work with a builder and solver */
/*  * that recognizes the presence of periodic conditions, as ResidualBasedBlockBuilderAndSolverPeriodic */
/*  * for shared memory runs or TrilinosResidualBasedBuilderAndSolverMLPeriodic (which is specific */
/*  * for monolithic incompressible flow problems) for an MPI implementation. */
/*  * @see TwoStepPeriodicConditionUtilities,ResidualBasedBlockBuilderAndSolverPeriodic,TrilinosResidualBasedBuilderAndSolverMLPeriodic */
/*  *\/ */
template< unsigned int TDim >
class TwoStepPeriodicCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TwoStepPeriodicCondition
    KRATOS_CLASS_POINTER_DEFINITION(TwoStepPeriodicCondition);

    typedef IndexedObject IndexedObjectType;

    typedef Condition BaseType;

    typedef Node<3> NodeType;

    typedef Properties PropertiesType;

    typedef Geometry<NodeType> GeometryType;

    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    typedef VectorMap<IndexType, DataValueContainer> SolutionStepsConditionalDataContainerType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    /** @param NewId Index number of the new condition (optional)
     */
    TwoStepPeriodicCondition(IndexType NewId = 0);

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    TwoStepPeriodicCondition(IndexType NewId,
                      const NodesArrayType& ThisNodes);

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    TwoStepPeriodicCondition(IndexType NewId,
                      GeometryType::Pointer pGeometry);

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    TwoStepPeriodicCondition(IndexType NewId,
                      GeometryType::Pointer pGeometry,
                      PropertiesType::Pointer pProperties);

    /// Copy constructor.
    TwoStepPeriodicCondition(TwoStepPeriodicCondition const& rOther);


    /// Destructor.
    virtual ~TwoStepPeriodicCondition();


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    TwoStepPeriodicCondition & operator=(TwoStepPeriodicCondition const& rOther);

    ///@}
    ///@name Operations
    ///@{

    /// Create a new TwoStepPeriodicCondition instance
    Condition::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              PropertiesType::Pointer pProperties) const;

    /// Check input to ensure that it makes sense.
    virtual int Check(const ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo);

   virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo);

   virtual void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo);

    /// Provides the global indices for each one of this element's local rows
    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo ProcessInfo instance (unused)
     */
    virtual void EquationIdVector(EquationIdVectorType& rResult,
                                  ProcessInfo& rCurrentProcessInfo);

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo ProcessInfo instance (unused)
     */
    virtual void GetDofList(DofsVectorType& ElementalDofList,
                            ProcessInfo& CurrentProcessInfo);

    /// Returns the values of the unknowns for each node
    virtual void GetValuesVector(Vector& Values, int Step = 0);

    ///@}
    ///@name Conditional Data
    ///@{


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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "TwoStepPeriodicCondition #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "TwoStepPeriodicCondition #" << Id();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        Condition::PrintData(rOStream);
    }


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

    void GetVelocityDofList(DofsVectorType& rElementalDofList,
            ProcessInfo& rCurrentProcessInfo);

    void GetPressureDofList(DofsVectorType& rElementalDofList,
            ProcessInfo& rCurrentProcessInfo);



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
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


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

}; // Class TwoStepPeriodicCondition

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim >
inline std::istream & operator >>(std::istream& rIStream,
                                  TwoStepPeriodicCondition<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim >
inline std::ostream & operator <<(std::ostream& rOStream,
                                  const TwoStepPeriodicCondition<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@}

} // namespace Kratos.


#endif	/* KRATOS_TWO_STEP_PERIODIC_CONDITION_H */

