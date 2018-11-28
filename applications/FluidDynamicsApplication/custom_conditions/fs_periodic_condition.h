/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

#ifndef KRATOS_FS_PERIODIC_CONDITION_H
#define	KRATOS_FS_PERIODIC_CONDITION_H

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
#include "fluid_dynamics_application_variables.h"

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

/// Condition used to assign periodic boundary conditions when using a fractional step fluid solver.
/**
 * A conformant mapping between the two related boundaries is assumed, so this condition
 * will link two nodes, one on each boundary, which are considered to be "images" of one
 * another on both periodic sides.\n
 * Note that the periodic boundary condition is enforced by the builder and solver, this Condition
 * is simply used to provide it the required information to set it up.\n
 * @note FSPeriodicCondition is designed designed to work with a builder and solver
 * that recognizes the presence of periodic conditions, as ResidualBasedBlockBuilderAndSolverPeriodic
 * for shared memory runs or TrilinosResidualBasedBuilderAndSolverMLPeriodic (which is specific
 * for monolithic incompressible flow problems) for an MPI implementation.
 * @see FSPeriodicConditionUtilities,ResidualBasedBlockBuilderAndSolverPeriodic,TrilinosResidualBasedBuilderAndSolverMLPeriodic
 */
template< unsigned int TDim >
class FSPeriodicCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FSPeriodicCondition
    KRATOS_CLASS_POINTER_DEFINITION(FSPeriodicCondition);

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

    ;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    /** @param NewId Index number of the new condition (optional)
     */
    FSPeriodicCondition(IndexType NewId = 0);

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    FSPeriodicCondition(IndexType NewId,
                      const NodesArrayType& ThisNodes);

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    FSPeriodicCondition(IndexType NewId,
                      GeometryType::Pointer pGeometry);

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    FSPeriodicCondition(IndexType NewId,
                      GeometryType::Pointer pGeometry,
                      PropertiesType::Pointer pProperties);

    /// Copy constructor.
    FSPeriodicCondition(FSPeriodicCondition const& rOther);


    /// Destructor.
    ~FSPeriodicCondition() override;


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    FSPeriodicCondition & operator=(FSPeriodicCondition const& rOther);

    ///@}
    ///@name Operations
    ///@{

    /// Create a new FSPeriodicCondition instance
    Condition::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              PropertiesType::Pointer pProperties) const override;

    /// Create a new FSPeriodicCondition instance
    Condition::Pointer Create(
		IndexType NewId,
		GeometryType::Pointer pGeom,
		PropertiesType::Pointer pProperties) const override;
    

    /// Check input to ensure that it makes sense.
    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo) override;

   void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) override;

   void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo) override;

    /// Provides the global indices for each one of this element's local rows
    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo ProcessInfo instance (unused)
     */
    void EquationIdVector(EquationIdVectorType& rResult,
                                  ProcessInfo& rCurrentProcessInfo) override;

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo ProcessInfo instance (unused)
     */
    void GetDofList(DofsVectorType& ElementalDofList,
                            ProcessInfo& CurrentProcessInfo) override;

    /// Returns the values of the unknowns for each node
    void GetValuesVector(Vector& Values, int Step = 0) override;

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
        buffer << "FSPeriodicCondition #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FSPeriodicCondition #" << Id();
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

}; // Class FSPeriodicCondition

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim >
inline std::istream & operator >>(std::istream& rIStream,
                                  FSPeriodicCondition<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim >
inline std::ostream & operator <<(std::ostream& rOStream,
                                  const FSPeriodicCondition<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@}

} // namespace Kratos.


#endif	/* KRATOS_FS_PERIODIC_CONDITION_H */

