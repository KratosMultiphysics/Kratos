/*
==============================================================================
Kratos Fluid Dynamics Application
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

#ifndef KRATOS_MONOLITHIC_WALL_CONDITION_H
#define KRATOS_MONOLITHIC_WALL_CONDITION_H

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/process_info.h"

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

    /// Implements a wall condition for the monolithic formulation.
    /** It is intended to be used in combination with ASGS and VMS elements or their derived classes.
      @see ASGS2D,ASGS3D,VMS
     */
    template< unsigned int TDim, unsigned int TNumNodes = TDim >
    class MonolithicWallCondition : public Condition
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of MonolithicWallCondition
        KRATOS_CLASS_POINTER_DEFINITION(MonolithicWallCondition);

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

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        /** Admits an Id as a parameter.
          @param NewId Index for the new condition
          */
        MonolithicWallCondition(IndexType NewId = 0):
            Condition(NewId)
        {
        }

        /// Constructor using an array of nodes
        /**
         @param NewId Index of the new condition
         @param ThisNodes An array containing the nodes of the new condition
         */
        MonolithicWallCondition(IndexType NewId,
                                const NodesArrayType& ThisNodes):
            Condition(NewId,ThisNodes)
        {
        }

        /// Constructor using Geometry
        /**
         @param NewId Index of the new condition
         @param pGeometry Pointer to a geometry object
         */
        MonolithicWallCondition(IndexType NewId,
                                GeometryType::Pointer pGeometry):
            Condition(NewId,pGeometry)
        {
        }

        /// Constructor using Properties
        /**
         @param NewId Index of the new element
         @param pGeometry Pointer to a geometry object
         @param pProperties Pointer to the element's properties
         */
        MonolithicWallCondition(IndexType NewId,
                                GeometryType::Pointer pGeometry,
                                PropertiesType::Pointer pProperties):
            Condition(NewId,pGeometry,pProperties)
        {
        }

        /// Copy constructor.
        MonolithicWallCondition(MonolithicWallCondition const& rOther):
            Condition(rOther)
        {
        }

        /// Destructor.
        virtual ~MonolithicWallCondition(){}


        ///@}
        ///@name Operators
        ///@{

        /// Copy constructor
        MonolithicWallCondition & operator=(MonolithicWallCondition const& rOther)
        {
            Condition::operator=(rOther);

            return *this;
        }

        ///@}
        ///@name Operations
        ///@{

        /// Create a new MonolithicWallCondition object.
        /**
          @param NewId Index of the new condition
          @param ThisNodes An array containing the nodes of the new condition
          @param pProperties Pointer to the element's properties
          */
        virtual Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
        {
            return Condition::Pointer(new MonolithicWallCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
        }


        /// Return local contributions of the correct size, filled with zeros (for compatibility with time schemes).
        /** The actual local contributions are computed in the Damping functions
          @see CalculateLocalVelocityContribution
          */
        virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          ProcessInfo& rCurrentProcessInfo)
        {
            const SizeType BlockSize = TDim + 1;
            const SizeType LocalSize = BlockSize * TNumNodes;

            if (rLeftHandSideMatrix.size1() != LocalSize)
                rLeftHandSideMatrix.resize(LocalSize,LocalSize);

            if (rRightHandSideVector.size() != LocalSize)
                rRightHandSideVector.resize(LocalSize);

            noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);
            noalias(rRightHandSideVector) = ZeroVector(LocalSize);
        }

        /// Return a matrix of the correct size, filled with zeros (for compatibility with time schemes).
        /** The actual local contributions are computed in the Damping functions
          @see DampMatrix
          */
        virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                           ProcessInfo& rCurrentProcessInfo)
        {
            const SizeType BlockSize = TDim + 1;
            const SizeType LocalSize = BlockSize * TNumNodes;

            if (rLeftHandSideMatrix.size1() != LocalSize)
                rLeftHandSideMatrix.resize(LocalSize,LocalSize);

            noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);
        }

        /// Return local right hand side of the correct size, filled with zeros (for compatibility with time schemes).
        /** The actual local contributions are computed in the Damping functions
          @see CalculateLocalVelocityContribution
          */
        virtual void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                            ProcessInfo& rCurrentProcessInfo)
        {
            const SizeType BlockSize = TDim + 1;
            const SizeType LocalSize = BlockSize * TNumNodes;

            if (rRightHandSideVector.size() != LocalSize)
                rRightHandSideVector.resize(LocalSize);

            noalias(rRightHandSideVector) = ZeroVector(LocalSize);
        }



        virtual void DampMatrix(MatrixType& rDampMatrix,
                                ProcessInfo& rCurrentProcessInfo)
        {
            VectorType RHS;
            this->CalculateLocalVelocityContribution(rDampMatrix,RHS,rCurrentProcessInfo);
        }



        virtual void CalculateLocalVelocityContribution(MatrixType& rDampMatrix,
                                                        VectorType& rRightHandSideVector,
                                                        ProcessInfo& rCurrentProcessInfo)
        {
            // Initialize local contributions
            const SizeType LocalSize = (TDim + 1) * TNumNodes;

            if (rDampMatrix.size1() != LocalSize)
                rDampMatrix.resize(LocalSize,LocalSize);
            if (rRightHandSideVector.size() != LocalSize)
                rRightHandSideVector.resize(LocalSize);

            noalias(rDampMatrix) = ZeroMatrix(LocalSize,LocalSize);
            noalias(rRightHandSideVector) = ZeroVector(LocalSize);

            this->ApplyWallLaw(rDampMatrix,rRightHandSideVector);
        }



        virtual int Check(const ProcessInfo& rCurrentProcessInfo)
        {
            KRATOS_TRY;

            int Check = Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area > 0

            if (Check != 0)
            {
                return Check;
            }
            else
            {
                // Add condition-specific checks here
                // we use: VELOCITY, MESH_VELOCITY, PRESSURE, IS_STRUCTURE (non-historic!) Y_WALL (in properties)
                return Check;
            }

            KRATOS_CATCH("");
        }


        /// Provides the global indices for each one of this element's local rows.
        /** This determines the elemental equation ID vector for all elemental DOFs
         * @param rResult A vector containing the global Id of each row
         * @param rCurrentProcessInfo the current process info object (unused)
         */
        virtual void EquationIdVector(EquationIdVectorType& rResult,
                                      ProcessInfo& rCurrentProcessInfo);


        /// Returns a list of the element's Dofs
        /**
         * @param ElementalDofList the list of DOFs
         * @param rCurrentProcessInfo the current process info instance
         */
        virtual void GetDofList(DofsVectorType& ConditionDofList,
                                ProcessInfo& CurrentProcessInfo);


        /// Returns VELOCITY_X, VELOCITY_Y, (VELOCITY_Z,) PRESSURE for each node
        /**
         * @param Values Vector of nodal unknowns
         * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
         */
        virtual void GetFirstDerivativesVector(Vector& Values,
                                               int Step = 0)
        {
            const SizeType LocalSize = (TDim + 1) * TNumNodes;
            unsigned int LocalIndex = 0;

            if (Values.size() != LocalSize)
                Values.resize(LocalSize, false);

            for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
            {
                array_1d<double,3>& rVelocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, Step);
                for (unsigned int d = 0; d < TDim; ++d)
                    Values[LocalIndex++] = rVelocity[d];
                Values[LocalIndex++] = this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE, Step);
            }
        }


        /// Returns ACCELERATION_X, ACCELERATION_Y, (ACCELERATION_Z,) 0 for each node
        /**
         * @param Values Vector of nodal second derivatives
         * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
         */
        virtual void GetSecondDerivativesVector(Vector& Values,
                                                int Step = 0)
        {
            const SizeType LocalSize = (TDim + 1) * TNumNodes;
            unsigned int LocalIndex = 0;

            if (Values.size() != LocalSize)
                Values.resize(LocalSize, false);

            for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
            {
                array_1d<double,3>& rVelocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(ACCELERATION, Step);
                for (unsigned int d = 0; d < TDim; ++d)
                    Values[LocalIndex++] = rVelocity[d];
                Values[LocalIndex++] = 0.0; // No value on pressure positions
            }
        }


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
            buffer << "MonolithicWallCondition" << TDim << "D";
            return buffer.str();
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "MonolithicWallCondition";}

        /// Print object's data.
        virtual void PrintData(std::ostream& rOStream) const {}


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

        void ApplyWallLaw(MatrixType& rLocalMatrix,
                          VectorType& rLocalVector)
        {
            GeometryType& rGeometry = this->GetGeometry();
            const size_t BlockSize = TDim + 1;

            double area = 0.3333333333333 * rGeometry.Area();

            for(size_t itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
            {
                if( rGeometry[itNode].GetValue(IS_STRUCTURE) == 1.0 )
                {
                    array_1d<double,3> Vel = rGeometry[itNode].FastGetSolutionStepValue(VELOCITY);
                    const array_1d<double,3>& VelMesh = rGeometry[itNode].FastGetSolutionStepValue(MESH_VELOCITY);
                    Vel -= VelMesh;

                    const double Ikappa = 1.0/0.41; // inverse of Von Karman's kappa
                    const double B = 5.2;
                    const double limit_yplus = 10.9931899; // limit between linear and log regions

                    const double y = rGeometry[itNode].GetValue(Y_WALL); // wall distance to use in stress calculation

                    const double rho = rGeometry[itNode].FastGetSolutionStepValue(DENSITY);
                    const double nu = rGeometry[itNode].FastGetSolutionStepValue(VISCOSITY);

                    double wall_vel = 0.0;
                    for (size_t d = 0; d < 3; d++)
                    {
                        wall_vel += Vel[d]*Vel[d];
                    }
                    wall_vel = sqrt(wall_vel);

                    if(wall_vel <= 1e-10)
                        wall_vel = 1e-10;

                    if (wall_vel > 1e-12) // do not bother if velocity is zero
                    {

                        // linear region
                        double utau = sqrt(wall_vel * nu / y);
                        double yplus = y * utau / nu;

                        // log region
                        if (yplus > limit_yplus)
                        {

                            // wall_vel / utau = 1/kappa * log(yplus) + B
                            // this requires solving a nonlinear problem:
                            // f(utau) = utau*(1/kappa * log(y*utau/nu) + B) - wall_vel = 0
                            // note that f'(utau) = 1/kappa * log(y*utau/nu) + B + 1/kappa

                            unsigned int iter = 0;
                            double dx = 1e10;
                            const double tol = 1e-6;
                            double uplus = Ikappa * log(yplus) + B;

                            while(iter < 100 && fabs(dx) > tol * utau)
                            {
                                // Newton-Raphson iteration
                                double f = utau * uplus - wall_vel;
                                double df = uplus + Ikappa;
                                dx = f/df;

                                // Update variables
                                utau -= dx;
                                yplus = y * utau / nu;
                                uplus = Ikappa * log(yplus) + B;
                                ++iter;
                            }
                            if (iter == 100)
                            {
                                std::cout << "Warning: wall condition Newton-Raphson did not converge. Residual is " << dx << std::endl;
                            }
                        }

                        for (size_t d = 0; d < 3; d++)
                        {
                            size_t k = itNode*BlockSize+d;
                            rLocalVector[k] -= Vel[d] * area * utau * utau * rho / wall_vel;
                            rLocalMatrix(k,k) += area * utau * utau * rho / wall_vel;
                        }
                    }
                }
            }
        }


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

        virtual void save(Serializer& rSerializer) const
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
        }

        virtual void load(Serializer& rSerializer)
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
        }

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

    }; // Class MonolithicWallCondition


    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function
    template< unsigned int TDim, unsigned int TNumNodes >
    inline std::istream& operator >> (std::istream& rIStream,
                                      MonolithicWallCondition<TDim,TNumNodes>& rThis)
    {
        return rIStream;
    }

    /// output stream function
    template< unsigned int TDim, unsigned int TNumNodes >
    inline std::ostream& operator << (std::ostream& rOStream,
                                      const MonolithicWallCondition<TDim,TNumNodes>& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }

    ///@}

    ///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_MONOLITHIC_WALL_CONDITION_H
