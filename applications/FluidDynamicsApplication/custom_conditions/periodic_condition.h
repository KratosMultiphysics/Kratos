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
#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/periodic_variables_container.h"

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

    /// Condition used to assign periodic boundary conditions
    /**
     * This condition imposes periodic boundary conditions in a weak sense.
     * A conformant mapping between the two related boundaries is assumed, so this condition
     * will link two nodes, one on each boundary.
     */
    class PeriodicCondition : public Condition
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
        virtual ~PeriodicCondition();


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
                                  PropertiesType::Pointer pProperties) const;

        /// Check input to ensure that it makes sense.
        virtual int Check(const ProcessInfo& rCurrentProcessInfo);

        /// Returns a matrix of penalty terms for the periodic variables.
        /**
         * The weight of the penalty terms is given by the member variable mWeight,
         * set using SetValueOnIntegrationPoints. The periodic variables are read from
         * the value of PERIODIC_VARIABLES stored in rCurrentProcessInfo.
         * @param rLeftHandSideMatrix Local left hand side matrix (output)
         * @param rRightHandSideVector Local right hand side vector (output)
         * @param rCurrentProcessInfo ProcessInfo instance (unused)
         */
        virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          ProcessInfo& rCurrentProcessInfo);

        /// Returns a matrix of penalty terms for the periodic variables.
        /**
         * @param rLeftHandSideMatrix Local left hand side matrix (output)
         * @param rCurrentProcessInfo ProcessInfo instance (unused)
         */
        virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                           ProcessInfo& rCurrentProcessInfo);

        /// Returns RHS values for the penalized dofs.
        /**
         * @param rRightHandSideVector Local right hand side vector (output)
         * @param rCurrentProcessInfo ProcessInfo instance (unused)
         */
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
            buffer << "PeriodicCondition #" << Id();
            return buffer.str();
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "PeriodicCondition #" << Id();
        }

//        /// Print object's data.
//        virtual void PrintData(std::ostream& rOStream) const
//        {
//            BaseType::PrintData(rOStream);
//        }


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

    }; // Class PeriodicCondition

    ///@}

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

