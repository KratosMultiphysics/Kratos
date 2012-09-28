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

#ifndef KRATOS_WALL_CONDITION_H
#define KRATOS_WALL_CONDITION_H

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
    /**
      It is intended to be used in combination with ASGS and VMS elements or their derived classes
      and the ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent time scheme, which supports
      slip conditions.
      This condition will add a wall stress term to all nodes identified with IS_STRUCTURE!=0.0 (in the
      non-historic database, that is, assigned using Node.SetValue()). This stress term is determined
      according to the wall distance provided as Y_WALL.
      @see ASGS2D,ASGS3D,VMS,ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent
     */
    template< unsigned int TDim, unsigned int TNumNodes = TDim >
    class WallCondition : public Condition
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of WallCondition
        KRATOS_CLASS_POINTER_DEFINITION(WallCondition);

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
        WallCondition(IndexType NewId = 0):
            Condition(NewId)
        {
        }

        /// Constructor using an array of nodes
        /**
         @param NewId Index of the new condition
         @param ThisNodes An array containing the nodes of the new condition
         */
        WallCondition(IndexType NewId,
                                const NodesArrayType& ThisNodes):
            Condition(NewId,ThisNodes)
        {
        }

        /// Constructor using Geometry
        /**
         @param NewId Index of the new condition
         @param pGeometry Pointer to a geometry object
         */
        WallCondition(IndexType NewId,
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
        WallCondition(IndexType NewId,
                                GeometryType::Pointer pGeometry,
                                PropertiesType::Pointer pProperties):
            Condition(NewId,pGeometry,pProperties)
        {
        }

        /// Copy constructor.
        WallCondition(WallCondition const& rOther):
            Condition(rOther)
        {
        }

        /// Destructor.
        virtual ~WallCondition(){}


        ///@}
        ///@name Operators
        ///@{

        /// Copy constructor
        WallCondition & operator=(WallCondition const& rOther)
        {
            Condition::operator=(rOther);

            return *this;
        }

        ///@}
        ///@name Operations
        ///@{

        /// Create a new WallCondition object.
        /**
          @param NewId Index of the new condition
          @param ThisNodes An array containing the nodes of the new condition
          @param pProperties Pointer to the element's properties
          */
        virtual Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
        {
            return Condition::Pointer(new WallCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
        }


        virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                           ProcessInfo& rCurrentProcessInfo)
        {
            VectorType RHS;
            this->CalculateLocalSystem(rLeftHandSideMatrix,RHS,rCurrentProcessInfo);
        }



        /// Calculate wall stress term for all nodes with IS_STRUCTURE != 0.0
        /**
          @param rDampMatrix Left-hand side matrix
          @param rRightHandSideVector Right-hand side vector
          @param rCurrentProcessInfo ProcessInfo instance (unused)
          */
        virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          ProcessInfo& rCurrentProcessInfo)
        {
            if (rCurrentProcessInfo[FRACTIONAL_STEP] == 1)
            {
                // Initialize local contributions
                const SizeType LocalSize = TDim * TNumNodes;

                if (rLeftHandSideMatrix.size1() != LocalSize)
                    rLeftHandSideMatrix.resize(LocalSize,LocalSize);
                if (rRightHandSideVector.size() != LocalSize)
                    rRightHandSideVector.resize(LocalSize);

                noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);
                noalias(rRightHandSideVector) = ZeroVector(LocalSize);

                this->ApplyInflowCondition(rLeftHandSideMatrix,rRightHandSideVector);

                this->ApplyWallLaw(rLeftHandSideMatrix,rRightHandSideVector);
            }
            else
            {
                if (rLeftHandSideMatrix.size1() != 0)
                    rLeftHandSideMatrix.resize(0,0,false);

                if (rRightHandSideVector.size() != 0)
                    rRightHandSideVector.resize(0,false);
            }
        }


        /// Check that all data required by this condition is available and reasonable
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
                // Check that all required variables have been registered
                if(VELOCITY.Key() == 0)
                    KRATOS_ERROR(std::invalid_argument,"VELOCITY Key is 0. Check if the application was correctly registered.","");
                if(MESH_VELOCITY.Key() == 0)
                    KRATOS_ERROR(std::invalid_argument,"MESH_VELOCITY Key is 0. Check if the application was correctly registered.","");
                if(NORMAL.Key() == 0)
                    KRATOS_ERROR(std::invalid_argument,"NORMAL Key is 0. Check if the application was correctly registered.","")
                if(IS_STRUCTURE.Key() == 0)
                    KRATOS_ERROR(std::invalid_argument,"IS_STRUCTURE Key is 0. Check if the application was correctly registered.","");
                if(Y_WALL.Key() == 0)
                    KRATOS_ERROR(std::invalid_argument,"Y_WALL Key is 0. Check if the application was correctly registered.","")

                // Checks on nodes

                // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
                for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
                {
                    if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
                        KRATOS_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
                    if(this->GetGeometry()[i].SolutionStepsDataHas(MESH_VELOCITY) == false)
                        KRATOS_ERROR(std::invalid_argument,"missing MESH_VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
                    if(this->GetGeometry()[i].SolutionStepsDataHas(NORMAL) == false)
                        KRATOS_ERROR(std::invalid_argument,"missing NORMAL variable on solution step data for node ",this->GetGeometry()[i].Id());
                    if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
                       this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
                       this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
                        KRATOS_ERROR(std::invalid_argument,"missing VELOCITY component degree of freedom on node ",this->GetGeometry()[i].Id());
                }

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


        /// Returns VELOCITY_X, VELOCITY_Y, (VELOCITY_Z,) for each node
        /**
         * @param Values Vector of nodal unknowns
         * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
         */
        virtual void GetValuesVector(Vector& Values,
                                     int Step = 0)
        {
            const SizeType LocalSize = TDim * TNumNodes;
            unsigned int LocalIndex = 0;

            if (Values.size() != LocalSize)
                Values.resize(LocalSize, false);

            for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
            {
                array_1d<double,3>& rVelocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, Step);
                for (unsigned int d = 0; d < TDim; ++d)
                    Values[LocalIndex++] = rVelocity[d];
            }
        }


//        /// Returns ACCELERATION_X, ACCELERATION_Y, (ACCELERATION_Z,) 0 for each node
//        /**
//         * @param Values Vector of nodal second derivatives
//         * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
//         */
//        virtual void GetSecondDerivativesVector(Vector& Values,
//                                                int Step = 0)
//        {
//            const SizeType LocalSize = (TDim + 1) * TNumNodes;
//            unsigned int LocalIndex = 0;

//            if (Values.size() != LocalSize)
//                Values.resize(LocalSize, false);

//            for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
//            {
//                array_1d<double,3>& rVelocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(ACCELERATION, Step);
//                for (unsigned int d = 0; d < TDim; ++d)
//                    Values[LocalIndex++] = rVelocity[d];
//            }
//        }


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
            buffer << "WallCondition" << TDim << "D";
            return buffer.str();
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "WallCondition";}

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

        /// Commpute the wall stress and add corresponding terms to the system contributions.
        /**
          @param rLocalMatrix Local system matrix
          @param rLocalVector Local right hand side
          */
        void ApplyWallLaw(MatrixType& rLocalMatrix,
                          VectorType& rLocalVector)
        {
            GeometryType& rGeometry = this->GetGeometry();
            const unsigned int BlockSize = TDim;
            const double NodalFactor = 1.0 / double(TDim);

            double area = NodalFactor * rGeometry.DomainSize();
            // DomainSize() is the way to ask the geometry's length/area/volume (whatever is relevant for its dimension) without asking for the number of spatial dimensions first

            for(unsigned int itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
            {
                const NodeType& rConstNode = rGeometry[itNode];
                const double y = rConstNode.GetValue(Y_WALL); // wall distance to use in stress calculation
                if( y > 0.0 && rConstNode.GetValue(IS_STRUCTURE) != 0.0 )
                {
                    array_1d<double,3> Vel = rGeometry[itNode].FastGetSolutionStepValue(VELOCITY);
                    const array_1d<double,3>& VelMesh = rGeometry[itNode].FastGetSolutionStepValue(MESH_VELOCITY);
                    Vel -= VelMesh;

                    const double Ikappa = 1.0/0.41; // inverse of Von Karman's kappa
                    const double B = 5.2;
                    const double limit_yplus = 10.9931899; // limit between linear and log regions

                    const double rho = rGeometry[itNode].FastGetSolutionStepValue(DENSITY);
                    const double nu = rGeometry[itNode].FastGetSolutionStepValue(VISCOSITY);

                    double wall_vel = 0.0;
                    for (size_t d = 0; d < TDim; d++)
                    {
                        wall_vel += Vel[d]*Vel[d];
                    }
                    wall_vel = sqrt(wall_vel);

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
                                std::cout << "WARNING: wall condition Newton-Raphson did not converge. Residual is " << dx << std::endl;
                            }
                        }

                        const double Tmp = area * utau * utau * rho / wall_vel;
                        for (unsigned int d = 0; d < TDim; d++)
                        {
                            unsigned int k = itNode*BlockSize+d;
                            rLocalVector[k] -= Vel[d] * Tmp;
                            rLocalMatrix(k,k) += Tmp;
                        }
                    }
                }
            }
        }

        /// Apply condition to prevent numerical problems due to flow into the domain in unexpected places.
        /** This condition prevents problems arising from inflow in outflow areas, typically due to vortices
         *  exiting the domain.
         * @param rLocalMatrix Local LHS matrix
         * @param rLocalVector Local RHS vector
         */
        void ApplyInflowCondition(MatrixType& rLocalMatrix,
                                  VectorType& rLocalVector)
        {
            const double N = 1.0 / static_cast<double>(TNumNodes);
            const double Weight = N * N * this->GetGeometry().DomainSize();

            unsigned int FirstRow = 0;
            const unsigned int LocalSize = TNumNodes;// + 1;

            for( unsigned int j = 0; j < TNumNodes; j++ )
            {
                const NodeType& rConstNode = this->GetGeometry()[j];
                if ( rConstNode.GetValue(IS_STRUCTURE) == 0.0 )
                {
                    const array_1d<double,3>& rNormal = this->GetGeometry()[j].FastGetSolutionStepValue(NORMAL);
                    const array_1d<double,3>& rVel = this->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY);
                    double Proj = rNormal[0]*rVel[0] + rNormal[1]*rVel[1] + rNormal[2]*rVel[2];

                    if( Proj  < 0 )
                    {
                        double Tij = Proj * Weight;
                        for (unsigned int d = 0; d < TDim;d++)
                            rLocalMatrix(FirstRow+d,FirstRow+d) -= Tij;
                    }
                }

                FirstRow += LocalSize;
            }

            VectorType Values = ZeroVector(rLocalVector.size());
            this->GetValuesVector(Values);
            noalias(rLocalVector) -= prod(rLocalMatrix,Values);
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
                                      WallCondition<TDim,TNumNodes>& rThis)
    {
        return rIStream;
    }

    /// output stream function
    template< unsigned int TDim, unsigned int TNumNodes >
    inline std::ostream& operator << (std::ostream& rOStream,
                                      const WallCondition<TDim,TNumNodes>& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }

    ///@}

    ///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_WALL_CONDITION_H
