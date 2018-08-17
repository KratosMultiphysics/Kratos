//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

// This is an adapted version of the boundary condition for the stabilized ASGS/OSS element described in
// J. Cotela "Applications of Turbulence Modeling in Civil Engineering" PhD thesis, 2016

#ifndef KRATOS_SKSY_FLUID_CONDITION_H
#define KRATOS_SKSY_FLUID_CONDITION_H

// System includes
#include <string>
#include <iostream>

#include "includes/kratos_flags.h"

// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/process_info.h"


namespace Kratos
{
    ///@addtogroup TurbulentFlowApplication
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

    /// Boundary terms for DSS type elements.
    template< unsigned int TDim, unsigned int TNumNodes = TDim >
    class SkSyFluidCondition : public Condition
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of WallCondition
        KRATOS_CLASS_POINTER_DEFINITION(SkSyFluidCondition);

        typedef Node < 3 > NodeType;

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

        /// Type for shape function values container
        typedef Kratos::Vector ShapeFunctionsType;

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        /** Admits an Id as a parameter.
          @param NewId Index for the new condition
          */
        SkSyFluidCondition(IndexType NewId = 0);

        /// Constructor using an array of nodes
        /**
         @param NewId Index of the new condition
         @param ThisNodes An array containing the nodes of the new condition
         */
        SkSyFluidCondition(IndexType NewId,
                const NodesArrayType& ThisNodes);

        /// Constructor using Geometry
        /**
         @param NewId Index of the new condition
         @param pGeometry Pointer to a geometry object
         */
        SkSyFluidCondition(IndexType NewId,
                GeometryType::Pointer pGeometry);

        /// Constructor using Properties
        /**
         @param NewId Index of the new element
         @param pGeometry Pointer to a geometry object
         @param pProperties Pointer to the element's properties
         */
        SkSyFluidCondition(IndexType NewId,
                GeometryType::Pointer pGeometry,
                PropertiesType::Pointer pProperties);

        /// Copy constructor.
        SkSyFluidCondition(SkSyFluidCondition const& rOther);

        /// Destructor.
        virtual ~SkSyFluidCondition();


        ///@}
        ///@name Operators
        ///@{

        /// Copy constructor
        SkSyFluidCondition & operator=(SkSyFluidCondition const& rOther);

        ///@}
        ///@name Operations
        ///@{

        /// Create a new SkSyFluidCondition object.
        /**
          @param NewId Index of the new condition
          @param ThisNodes An array containing the nodes of the new condition
          @param pProperties Pointer to the element's properties
          */
        virtual Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;


        virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                           ProcessInfo& rCurrentProcessInfo) override;



        /// Calculate wall stress term for all nodes with IS_STRUCTURE != 0.0
        /**
          @param rDampMatrix Left-hand side matrix
          @param rRightHandSideVector Right-hand side vector
          @param rCurrentProcessInfo ProcessInfo instance (unused)
          */
        virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          ProcessInfo& rCurrentProcessInfo) override;


        virtual void CalculateLocalVelocityContribution(MatrixType& rDampMatrix,
                                                        VectorType& rRightHandSideVector,
                                                        ProcessInfo& rCurrentProcessInfo) override;


        /// Check that all data required by this condition is available and reasonable
        virtual int Check(const ProcessInfo& rCurrentProcessInfo) override;


        /// Provides the global indices for each one of this element's local rows.
        /** This determines the elemental equation ID vector for all elemental DOFs
         * @param rResult A vector containing the global Id of each row
         * @param rCurrentProcessInfo the current process info object (unused)
         */
        virtual void EquationIdVector(EquationIdVectorType& rResult,
                                      ProcessInfo& rCurrentProcessInfo) override;


        /// Returns a list of the element's Dofs
        /**
         * @param ElementalDofList the list of DOFs
         * @param rCurrentProcessInfo the current process info instance
         */
        virtual void GetDofList(DofsVectorType& ConditionDofList,
                                ProcessInfo& CurrentProcessInfo) override;


        /// Returns VELOCITY_X, VELOCITY_Y, (VELOCITY_Z,) PRESSURE for each node
        /**
         * @param Values Vector of nodal unknowns
         * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
         */
        virtual void GetValuesVector(Vector& Values,
                                     int Step = 0) override;


        virtual void GetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                                                 std::vector<array_1d<double, 3 > >& rValues,
                                                 const ProcessInfo& rCurrentProcessInfo) override;



        virtual void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                                                 std::vector<double>& rValues,
                                                 const ProcessInfo& rCurrentProcessInfo) override;


        virtual void GetValueOnIntegrationPoints(const Variable<array_1d<double, 6 > >& rVariable,
                                                 std::vector<array_1d<double, 6 > >& rValues,
                                                 const ProcessInfo& rCurrentProcessInfo) override;

        virtual void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                 std::vector<Vector>& rValues,
                                                 const ProcessInfo& rCurrentProcessInfo) override;


        virtual void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                 std::vector<Matrix>& rValues,
                                                 const ProcessInfo& rCurrentProcessInfo) override;


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
        virtual std::string Info() const override;

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const override;

        /// Print object's data.
        virtual void PrintData(std::ostream& rOStream) const override;


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


        void CalculateNormal(array_1d<double,3>& An );

        double CalculateJacobian( double Area ) const;


        /// Apply boundary terms to allow imposing a pressure (normal stress), with a correction to prevent inflow.
        /** This correction should prevent numerical problems arising from inflow in outflow areas, typically due to vortices.
         *  exiting the domain.
         * @param rLocalMatrix Local LHS matrix
         * @param rLocalVector Local RHS vector
         */
        void AddBoundaryTerms(MatrixType& rLocalMatrix,
                              VectorType& rLocalVector,
                              const ProcessInfo& rCurrentProcessInfo);


        /**
         * @brief EvaluateInPoint Interpolate nodal data inside the element.
         * Evaluate a nodal variable in the point where the form functions take the
         * values given by rShapeFunc and write the result to rResult.
         * This is an auxiliary function used to compute values in integration points.
         * @param rResult The variable where the value will be added to
         * @param rVar The nodal variable to be read
         * @param rShapeFunc The values of the form functions in the point
         */
        template< class TVariableType >
        void EvaluateInPoint(TVariableType& rResult,
                             const Kratos::Variable<TVariableType>& Var,
                             const ShapeFunctionsType& rShapeFunc) const
        {
            const GeometryType& rGeom = this->GetGeometry();

            rResult = rShapeFunc[0] * rGeom[0].FastGetSolutionStepValue(Var);

            for(unsigned int i = 1; i < TNumNodes; i++)
            {
                rResult += rShapeFunc[i] * rGeom[i].FastGetSolutionStepValue(Var);
            }
        }


        void ResolvedConvectiveVelocity(array_1d<double,3>& rConvVel,
                                        const ShapeFunctionsType& rN);


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

        virtual void save(Serializer& rSerializer) const override;

        virtual void load(Serializer& rSerializer) override;

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

    }; // Class WallCondition


    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function
    template< unsigned int TDim, unsigned int TNumNodes >
    inline std::istream& operator >> (std::istream& rIStream,
                                      SkSyFluidCondition<TDim,TNumNodes>& rThis)
    {
        return rIStream;
    }

    /// output stream function
    template< unsigned int TDim, unsigned int TNumNodes >
    inline std::ostream& operator << (std::ostream& rOStream,
                                      const SkSyFluidCondition<TDim,TNumNodes>& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }

    ///@}

    ///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_SKSY_FLUID_CONDITION_H
