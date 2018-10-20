//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#ifndef KRATOS_PERIODIC_CONDITION_H
#define	KRATOS_PERIODIC_CONDITION_H

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
#include "includes/variables.h"
#include "containers/periodic_variables_container.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
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

/// Condition used to assign periodic boundary conditions.
/**
 * A conformant mapping between the two related boundaries is assumed, so this condition
 * will link two nodes, one on each boundary, which are considered to be "images" of one
 * another on both periodic sides.\n
 * Note that the periodic boundary condition is enforced by the builder and solver, this Condition
 * is simply used to provide it the required information to set it up.\n
 * This condtion has to be carefully set up in order to work properly, see PeriodicConditionUtilities,
 * which provides some tools to do so.
 * @note PeriodicCondition is designed designed to work with a builder and solver
 * that recognizes the presence of periodic conditions, as ResidualBasedBlockBuilderAndSolverPeriodic
 * for shared memory runs or TrilinosResidualBasedBuilderAndSolverMLPeriodic (which is specific
 * for monolithic incompressible flow problems) for an MPI implementation.
 * @see PeriodicConditionUtilities,ResidualBasedBlockBuilderAndSolverPeriodic,TrilinosResidualBasedBuilderAndSolverMLPeriodic
 */
class KRATOS_API(KRATOS_CORE) PeriodicCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PeriodicCondition
    KRATOS_CLASS_POINTER_DEFINITION(PeriodicCondition);

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
    PeriodicCondition(IndexType NewId = 0);

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    PeriodicCondition(IndexType NewId,
                      const NodesArrayType& ThisNodes);

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    PeriodicCondition(IndexType NewId,
                      GeometryType::Pointer pGeometry);

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    PeriodicCondition(IndexType NewId,
                      GeometryType::Pointer pGeometry,
                      PropertiesType::Pointer pProperties);

    /// Copy constructor.
    PeriodicCondition(PeriodicCondition const& rOther);


    /// Destructor.
    ~PeriodicCondition() override;


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    PeriodicCondition & operator=(PeriodicCondition const& rOther);

    ///@}
    ///@name Operations
    ///@{

    /// Create a new PeriodicCondition instance
    Condition::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              PropertiesType::Pointer pProperties) const override;

    /// Check input to ensure that it makes sense.
    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    /// Returns a matrix of penalty terms for the periodic variables.
    /**
     * The weight of the penalty terms is given by the member variable mWeight,
     * set using SetValueOnIntegrationPoints. The periodic variables are read from
     * the value of PERIODIC_VARIABLES stored in rCurrentProcessInfo.
     * @param rLeftHandSideMatrix Local left hand side matrix (output)
     * @param rRightHandSideVector Local right hand side vector (output)
     * @param rCurrentProcessInfo ProcessInfo instance (unused)
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;

    /// Returns a matrix of penalty terms for the periodic variables.
    /**
     * @param rLeftHandSideMatrix Local left hand side matrix (output)
     * @param rCurrentProcessInfo ProcessInfo instance (unused)
     */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               const ProcessInfo& rCurrentProcessInfo) override;

    /// Returns RHS values for the penalized dofs.
    /**
     * @param rRightHandSideVector Local right hand side vector (output)
     * @param rCurrentProcessInfo ProcessInfo instance (unused)
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                const ProcessInfo& rCurrentProcessInfo) override;

    /// Provides the global indices for each one of this element's local rows
    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo ProcessInfo instance (unused)
     */
    void EquationIdVector(EquationIdVectorType& rResult,
                          const ProcessInfo& rCurrentProcessInfo) const override;

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo ProcessInfo instance (unused)
     */
    void GetDofList(DofsVectorType& ElementalDofList,
                    const ProcessInfo& CurrentProcessInfo) override;

    /// Returns the values of the unknowns for each node
    void GetValuesVector(Vector& Values, int Step = 0) const override;

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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "PeriodicCondition #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PeriodicCondition #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;


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

}; // Class PeriodicCondition

///@}

template class KRATOS_API(KRATOS_CORE) KratosComponents<PeriodicCondition >;

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream & operator >>(std::istream& rIStream,
                                  PeriodicCondition& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream & operator <<(std::ostream& rOStream,
                                  const PeriodicCondition& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@}

} // namespace Kratos.


#endif	/* KRATOS_FLUID_PERTIODIC_CONDITION_H */


