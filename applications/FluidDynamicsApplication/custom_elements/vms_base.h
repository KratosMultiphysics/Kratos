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

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: jcotela $
//   Date:                $Date: 2010-10-09 10:34:00 $
//   Revision:            $Revision: 0.1 $
//
//


#if !defined(KRATOS_VMS_BASE_H_INCLUDED )
#define  KRATOS_VMS_BASE_H_INCLUDED

// Choose alternate Tau formula by defining this flag
//#define KRATOS_VMS_ALT_TAU

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/element.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{
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

    /// A stabilized element for the incompressible Navier-Stokes equations
    /**
     * This class implements a common framework for elements implementing a stabilized
     * Variational Multiscale formulation. The the subscales can be modeled
     * using either Algebraic Subgird Scales (ASGS) or Orthogonal Subscales (OSS).
     * In the case of OSS, the projection terms are treated explicitly (computed
     * using the results of the previous iteration) and the subscales are not
     * tracked in time. The choice of subscale model is made based on the Process Info
     * variable OSS_SWITCH (OSS if 1.0, ASGS otherwise)
     * Note that this class does NOT contain a complete implementation, but contains
     * common functions for the 2D and 3D implementations (VMS2D,VMS3D)
     * @see VMS2D,VMS3D
     */
    template< unsigned int TDim,
              unsigned int TNumNodes = TDim + 1,
              unsigned int TBlockSize = TDim + 1,
              unsigned int TLocalSize = TBlockSize * TNumNodes
            >
    class VMSBase : public Element
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of VMSBase
        KRATOS_CLASS_POINTER_DEFINITION(VMSBase);

        ///base type: an IndexedObject that automatically has a unique number
        typedef IndexedObject BaseType;

        ///definition of node type (default is: Node<3>)
        typedef Node < 3 > NodeType;

        /**
         * Properties are used to store any parameters
         * related to the constitutive law
         */
        typedef Properties PropertiesType;

        ///definition of the geometry type with given NodeType
        typedef Geometry<NodeType> GeometryType;

        ///definition of nodes container type, redefined from GeometryType
        typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

        typedef Vector VectorType;

        typedef Matrix MatrixType;

        typedef std::size_t IndexType;

        typedef std::size_t SizeType;

        typedef std::vector<std::size_t> EquationIdVectorType;

        typedef std::vector< Dof<double>::Pointer > DofsVectorType;

        typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

        typedef VectorMap<IndexType, DataValueContainer> SolutionStepsElementalDataContainerType;

        ///@}
        ///@name Life Cycle
        ///@{

        ///Constructors.

        /**
         * Default constuctor.
         * @param NewId Index number of the new element (optional)
         */
        VMSBase(IndexType NewId = 0) :
            Element(NewId)
        {}

        /**
         * Constructor using an array of nodes.
         * @param NewId Index of the new element
         * @param ThisNodes An array containing the nodes of the new element
         */
        VMSBase(IndexType NewId, const NodesArrayType& ThisNodes) :
            Element(NewId, ThisNodes)
        {}

        /**
         * Constructor using a geometry object.
         * @param NewId Index of the new element
         * @param pGeometry Pointer to a geometry object
         */
        VMSBase(IndexType NewId, GeometryType::Pointer pGeometry) :
            Element(NewId, pGeometry)
        {}

        /**
         * Constuctor using geometry and properties.
         * @param NewId Index of the new element
         * @param pGeometry Pointer to a geometry object
         * @param pProperties Pointer to the element's properties
         */
        VMSBase(IndexType NewId, GeometryType::Pointer pGeometry,
                 PropertiesType::Pointer pProperties) :
            Element(NewId, pGeometry, pProperties)
        {}

        /// Destructor.

        virtual ~VMSBase()
        {}


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        /// Create a new element of this type
        /**
         * Returns a pointer to a new VMSBase element, created using given input
         * @param NewId: the ID of the new element
         * @param ThisNodes: the nodes of the new element
         * @param pProperties: the properties assigned to the new element
         * @return a Pointer to the new element
         */
        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                PropertiesType::Pointer pProperties) const
        {
            KRATOS_TRY
//            return Element::Pointer(new VMSBase(NewId, GetGeometry().Create(ThisNodes), pProperties));
            KRATOS_ERROR(std::logic_error,"VMSBase::Create failed. Cannot create an instance of an abstract class. Please use VMS2D or VMS3D elements instead","")
            KRATOS_CATCH("");
        }

                /// Provides local contributions from body forces and OSS projection terms
        /**
         * This is called during the assembly process and provides the terms of the
         * system that are either constant or computed explicitly (from the 'old'
         * iteration variables). In this case this means the body force terms and the
         * OSS projections, that are treated explicitly.
         * @param rLeftHandSideMatrix: the elemental left hand side matrix. Not used here, required for compatibility purposes only.
         * @param rRightHandSideVector: the elemental right hand side
         * @param rCurrentProcessInfo: the current process info
         */
        virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          ProcessInfo& rCurrentProcessInfo)
        {
            // Check sizes and initialize matrix
            if (rLeftHandSideMatrix.size1() != TLocalSize)
                rLeftHandSideMatrix.resize(TLocalSize, TLocalSize, false);

            noalias(rLeftHandSideMatrix) = ZeroMatrix(TLocalSize, TLocalSize);

            // Calculate RHS
            this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
        }

        /// Returns a zero matrix of appropiate size (provided for compatibility with scheme)
        virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                           ProcessInfo& rCurrentProcessInfo)
        {
            if (rLeftHandSideMatrix.size1() != TLocalSize)
                rLeftHandSideMatrix.resize(TLocalSize, TLocalSize);

            noalias(rLeftHandSideMatrix) = ZeroMatrix(TLocalSize, TLocalSize);
        }

        /// Provides local contributions from body forces and projections to the RHS
        /**
         * This is called during the assembly process and provides the RHS terms of the
         * system that are either constant or computed explicitly (from the 'old'
         * iteration variables). In this case this means the body force terms and the
         * OSS projections, that are treated explicitly.
         * @param rRightHandSideVector: the elemental right hand side
         * @param rCurrentProcessInfo: the current process info
         */
        virtual void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                            ProcessInfo& rCurrentProcessInfo)
        {
            // Check sizes and initialize
            if (rRightHandSideVector.size() != TLocalSize)
                rRightHandSideVector.resize(TLocalSize, false);

            noalias(rRightHandSideVector) = ZeroVector(TLocalSize);

            // Calculate this element's geometric parameters
            double Area;
            array_1d<double, TNumNodes> N;
            boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            // Calculate this element's fluid properties
            double Density, KinViscosity;
            this->GetPointContribution(Density, DENSITY, N);
            this->GetPointContribution(KinViscosity, VISCOSITY, N);

            // Calculate Momentum RHS contribution
            this->AddMomentumRHS(rRightHandSideVector, Density, N, Area);

            // For OSS: Add projection of residuals to RHS
            if (rCurrentProcessInfo[OSS_SWITCH] == 1)
            {
                array_1d<double, 3> AdvVel;
                this->GetAdvectiveVel(AdvVel, N);

                // Calculate stabilization parameters
                double TauOne, TauTwo;
                this->CalculateTau(TauOne, TauTwo, AdvVel, Area, KinViscosity, rCurrentProcessInfo);

                this->AddProjectionToRHS(rRightHandSideVector, TauOne, TauTwo, N, DN_DX, Area);
            }
        }

        /// Computes local contributions to the mass matrix
        /**
         * Provides the local contributions to the mass matrix, which is defined here
         * as the matrix associated to velocity derivatives. Note that the mass
         * matrix implemented here is lumped.
         * @param rMassMatrix: the elemental mass matrix
         * @param rCurrentProcessInfo: the current process info instance
         */
        virtual void MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
        {
            // Resize and set to zero
            if (rMassMatrix.size1() != TLocalSize)
                rMassMatrix.resize(TLocalSize, TLocalSize, false);

            rMassMatrix = ZeroMatrix(TLocalSize, TLocalSize);

            // Get the element's geometric parameters
            double Area;
            array_1d<double, TNumNodes> N;
            boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            // Calculate this element's fluid properties
            double Density;
            this->GetPointContribution(Density, DENSITY, N);

            // Add 'classical' mass matrix (lumped)
            double Coeff = Density * Area / TNumNodes; ///@TODO:Optimize
            this->CalculateLumpedMassMatrix(rMassMatrix, Coeff);

            // For ASGS: add dynamic stabilization terms
            // in OSS we treat them explicitly, see CalculateLocalVelocityContribution
            if (rCurrentProcessInfo[OSS_SWITCH] != 1)
            {
                double KinViscosity;
                this->GetPointContribution(KinViscosity, VISCOSITY, N);

                // Get Advective velocity
                array_1d<double, 3> AdvVel;
                this->GetAdvectiveVel(AdvVel, N);

                // Calculate stabilization parameters
                double TauOne, TauTwo;
                this->CalculateTau(TauOne, TauTwo, AdvVel, Area, KinViscosity, rCurrentProcessInfo);

                // Add dynamic stabilization terms ( all terms involving a delta(u) )
                this->AddMassStabTerms<MatrixType> (rMassMatrix, Density, TauOne, N, DN_DX, Area);
            }
        }

        /// Computes the local contribution associated to 'new' velocity and pressure values
        /**
         * Provides local contributions to the system associated to the velocity and
         * pressure terms (convection, diffusion, pressure gradient/velocity divergence
         * and stabilization).
         * @param rDampMatrix: the velocity-proportional "damping" matrix
         * @param rRightHandSideVector: the elemental right hand side vector
         * @param rCurrentProcessInfo: the current process info instance
         */
        virtual void CalculateLocalVelocityContribution(MatrixType& rDampMatrix,
                                                        VectorType& rRightHandSideVector,
                                                        ProcessInfo& rCurrentProcessInfo)
        {
            // Resize and set to zero the matrix
            // Note that we don't clean the RHS because it will already contain body force (and stabilization) contributions
            if (rDampMatrix.size1() != TLocalSize)
                rDampMatrix.resize(TLocalSize, TLocalSize, false);

            noalias(rDampMatrix) = ZeroMatrix(TLocalSize, TLocalSize);

            // Get this element's geometric properties
            double Area;
            array_1d<double, TNumNodes> N;
            boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            // Calculate this element's fluid properties
            double Density, KinViscosity;
            this->GetPointContribution(Density, DENSITY, N);
            this->GetPointContribution(KinViscosity, VISCOSITY, N);

            // Get Advective velocity
            array_1d<double, 3> AdvVel;
            this->GetAdvectiveVel(AdvVel, N);

            // Calculate stabilization parameters
            double TauOne, TauTwo;
            this->CalculateTau(TauOne, TauTwo, AdvVel, Area, KinViscosity, rCurrentProcessInfo);

            this->AddIntegrationPointVelocityContribution(rDampMatrix, rRightHandSideVector, Density, KinViscosity, TauOne, TauTwo, N, DN_DX, Area);

            // Now calculate an additional contribution to the residual: r -= rDampMatrix * (u,p)
            VectorType U = ZeroVector(TLocalSize);
            int LocalIndex = 0;

            for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
            {
                array_1d< double,3 >& rVel = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
                for (unsigned int d = 0; d < TDim; ++d) // Velocity Dofs
                {
                    U[LocalIndex] = rVel[d];
                    ++LocalIndex;
                }
                U[LocalIndex] = this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE); // Pressure Dof
                ++LocalIndex;
            }

            noalias(rRightHandSideVector) -= prod(rDampMatrix, U);

            // The following is only needed for OSS
            if (rCurrentProcessInfo[OSS_SWITCH] == 1)
            {
                LocalIndex = 0;
                boost::numeric::ublas::bounded_matrix<double,TLocalSize,TLocalSize> MassStabilization = ZeroMatrix(TLocalSize, TLocalSize);

                this->AddMassStabTerms<boost::numeric::ublas::bounded_matrix<double,TLocalSize,TLocalSize> > (MassStabilization, Density, TauOne, N, DN_DX, Area);

                for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
                {
                    const array_1d<double, 3 > & rAcceleration = this->GetGeometry()[iNode].FastGetSolutionStepValue(ACCELERATION);
                    for (unsigned int d = 0; d < TDim; ++d) // Acceleration Dofs
                    {
                        U[LocalIndex] = rAcceleration[d];
                        ++LocalIndex;
                    }
                    U[LocalIndex] = 0.0; // Pressure Dof
                    ++LocalIndex;
                }

                noalias(rRightHandSideVector) -= prod(MassStabilization,U);
            }
        }

        /// Implementation of Calculate to compute the local OSS projections
        /**
         * This function computes the OSS projection terms from last iteration's
         * pressure and velocity values. The function signature is inherited from
         * element.h to mantain compatibility.
         * @param rVariable: Unused vector variable reference
         * @param Output: Unused output array
         * @param rCurrentProcessInfo: Process info instance
         */
        virtual void Calculate(const Variable<array_1d<double, 3 > >& rVariable,
                               array_1d<double, 3 > & rOutput,
                               const ProcessInfo& rCurrentProcessInfo)
        {
            // Get the element's geometric parameters
            double Area;
            array_1d<double, TNumNodes> N;
            boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            // Calculate this element's fluid properties
            double Density, Viscosity;
            this->GetPointContribution(Density, DENSITY, N);
            this->GetPointContribution(Viscosity, VISCOSITY, N);

            // Get Advective velocity
            array_1d<double, 3> AdvVel;
            this->GetAdvectiveVel(AdvVel, N);

            // Calculate stabilization parameters
            double TauOne, TauTwo;
            this->CalculateTau(TauOne, TauTwo, AdvVel, Area, Viscosity, rCurrentProcessInfo);

            // Output containers
            array_1d< double,3 > ElementalMomRes(3,0.0);
            double ElementalMassRes(0);

            this->AddProjectionResidualContribution(AdvVel, Density, TauOne, TauTwo, ElementalMomRes,ElementalMassRes, rCurrentProcessInfo, N, DN_DX, Area);

            if( rCurrentProcessInfo[OSS_SWITCH] == 1)
            {
                // Carefully write results to nodal variables, to avoid parallelism problems
                for (unsigned int i = 0; i < TNumNodes; ++i)
                {
                    ///@TODO: Test using atomic
                    this->GetGeometry()[i].SetLock(); // So it is safe to write in the node in OpenMP
                    array_1d< double, 3 > & rAdvProj = this->GetGeometry()[i].FastGetSolutionStepValue(ADVPROJ);
                    for (unsigned int d = 0; d < TDim; ++d)
                        rAdvProj[d] += N[i] * ElementalMomRes[d];

                    this->GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ) += N[i] * ElementalMassRes;
                    this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += Area * N[i];
                    this->GetGeometry()[i].UnSetLock(); // Free the node for other threads
                }
            }

            /// Return output
            rOutput = ElementalMomRes;
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
            buffer << "VMSBase #" << Id();
            return buffer.str();
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "VMSBase #" << Id();
        }

//        /// Print object's data.
//        virtual void PrintData(std::ostream& rOStream) const;

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

        /// Calculate Stabilization parameters
        /**
         * Calculates both tau parameters based on a given advective velocity
         * @param TauOne: First stabilization parameter (momentum equation)
         * @param TauTwo: Second stabilization parameter (mass equation)
         * @param rAdvVel: advection velocity
         * @param Area: Elemental area
         * @param KinViscosity: Elemental kinematic viscosity (nu)
         * @param rCurrentProcessInfo: Process info instance
         */

        #ifndef KRATOS_VMS_ALT_TAU
        virtual void CalculateTau(double& TauOne,
                                  double& TauTwo,
                                  const array_1d< double, 3 > & rAdvVel,
                                  const double Area,
                                  const double KinViscosity,
                                  const ProcessInfo& rCurrentProcessInfo)
        {
           // Compute mean advective velocity norm
            double AdvVelNorm = 0.0;
            for (unsigned int d = 0; d < TDim; ++d)
                AdvVelNorm += rAdvVel[d] * rAdvVel[d];

            AdvVelNorm = sqrt(AdvVelNorm);

            const double Element_Size = this->ElementSize(Area);

            TauOne = 1.0 / (rCurrentProcessInfo[DYNAMIC_TAU] / rCurrentProcessInfo[DELTA_TIME] + 4.0 * KinViscosity / (Element_Size * Element_Size) + 2.0 * AdvVelNorm / Element_Size);
            TauTwo = KinViscosity + 0.5 * Element_Size * AdvVelNorm;
        }

        #else

        virtual void CalculateTau(double& TauOne,
                                  double& TauTwo,
                                  const array_1d< double, 3 > & rAdvVel,
                                  const double Area,
                                  const double KinViscosity,
                                  const ProcessInfo& rCurrentProcessInfo)
        {
            // Compute mean advective velocity norm
            double AdvVelNorm = 0.0;
            for (unsigned int d = 0; d < TDim; ++d)
                AdvVelNorm += rAdvVel[d] * rAdvVel[d];

            AdvVelNorm = sqrt(AdvVelNorm);

            double Element_Size = 0;

            if(AdvVelNorm > 1e-16)
            {
                double tArea;
                array_1d<double, TNumNodes> N;
                boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
                GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, tArea);

                for(unsigned int j = 0; j < TNumNodes; ++j)
                {
                    double Temp = 0;
                    for(unsigned int i = 0; i < TDim; ++i)
                    {
                        Temp += rAdvVel[i]*DN_DX(j,i);
                    }
                    Element_Size += fabs(Temp);
                }

                Element_Size = 2.0*AdvVelNorm / Element_Size;
            }
            else
            {
                Element_Size = this->ElementSize(Area);
            }
//            const double Element_Size = this->ElementSize(Area);

            TauOne = 1.0 / (rCurrentProcessInfo[DYNAMIC_TAU] / rCurrentProcessInfo[DELTA_TIME] + 4.0 * KinViscosity / (Element_Size * Element_Size) + 2.0 * AdvVelNorm / Element_Size);
            TauTwo = KinViscosity + 0.5 * Element_Size * AdvVelNorm;
        }
        #endif

        /// Add the momentum equation contribution to the RHS (body forces)
        void AddMomentumRHS(VectorType& F,
                            const double Density,
                            const array_1d<double, TNumNodes>& rShapeFunc,
                            const double Weight)
        {
            double Coef = Density * Weight;

            array_1d<double, 3> BodyForce(3,0.0);
            this->AddPointContribution(BodyForce, BODY_FORCE, rShapeFunc, Coef);

            // Add the results to the velocity components (Local Dofs are vx, vy, [vz,] p for each node)
            int LocalIndex = 0;

            for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
            {
                for (unsigned int d = 0; d < TDim; ++d)
                {
                    F[LocalIndex++] += rShapeFunc[iNode] * BodyForce[d];
                }
                ++LocalIndex; // Skip pressure Dof
            }
        }

        /// Add OSS projection terms to the RHS
        void AddProjectionToRHS(VectorType& RHS,
                                const double TauOne,
                                const double TauTwo,
                                const array_1d<double, TNumNodes>& rShapeFunc,
                                const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim>& rShapeDeriv,
                                const double Weight)
        {
            array_1d<double, TNumNodes> AGradN;
            this->GetConvectionOperator(AGradN, rShapeFunc, rShapeDeriv); // Get a * grad(Ni)

            // Add to the ouptut vector
            unsigned int FirstRow(0); // Position of the term we want to add to
            double Const1, Const2; // Partial results that remain constant for a given j

            for (unsigned int j = 0; j < TNumNodes; ++j) // loop over nodes (for components of the residual)
            {
                // Compute the 'constant' part of the (nodal) projection terms
                // See this->Calculate() for the computation of the projections
                const array_1d< double, 3 > & rMomProj = this->GetGeometry()[j].FastGetSolutionStepValue(ADVPROJ);

                Const1 = Weight * TauOne * rShapeFunc[j];
                Const2 = Weight * TauTwo * rShapeFunc[j] * this->GetGeometry()[j].FastGetSolutionStepValue(DIVPROJ);

                // Reset row reference
                FirstRow = 0;

                for (unsigned int i = 0; i < TNumNodes; ++i) // loop over nodes (for components of L*(Vh) )
                {
                    for (unsigned int d = 0; d < TDim; ++d)
                    {
                        RHS[FirstRow + d] += Const1 * AGradN[i] * rMomProj[d] + Const2 * rShapeDeriv(i, d); // TauOne * ( a * Grad(v) ) * MomProjection + TauTwo * Div(v) * MassProjection
                        RHS[FirstRow + TDim] += Const1 * rShapeDeriv(i, d) * rMomProj[d]; // TauOne * Grad(q) * MomProjection
                    }
                    // Update row reference
                    FirstRow += TBlockSize;
                }
            }
        }

        /// Add lumped mass matrix
        /**
         * Adds the lumped mass matrix to an elemental LHS matrix. Note that the time factor
         * (typically 1/(k*Dt) ) is added by the scheme outside the element.
         * @param rLHSMatrix: The local matrix where the result will be added
         * @param Mass: The weight assigned to each node (typically Density * Area / NumNodes or Density*Volume / NumNodes)
         */
        void CalculateLumpedMassMatrix(MatrixType& rLHSMatrix,
                                       const double Mass)
        {
            unsigned int DofIndex = 0;
            for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
            {
                for (unsigned int d = 0; d < TDim; ++d)
                {
                    rLHSMatrix(DofIndex, DofIndex) += Mass;
                    ++DofIndex;
                }
                ++DofIndex; // Skip pressure Dof
            }
        }


        /// Add mass-like stabilization terms to LHS
        template < class TMatrixType >
        void AddMassStabTerms(TMatrixType& rLHSMatrix,
                              const double Density, const double TauOne,
                              const array_1d<double, TNumNodes>& rShapeFunc,
                              const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim>& rShapeDeriv,
                              const double Weight)
        {
            double Coef = Density * Weight * TauOne;
            unsigned int FirstRow(0), FirstCol(0);
            double K; // Temporary results

            // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
            array_1d<double, TNumNodes> AGradN;
            this->GetConvectionOperator(AGradN, rShapeFunc, rShapeDeriv); // Get a * grad(Ni)

            // Note: Dof order is (vx,vy,[vz,]p) for each node
            for (unsigned int i = 0; i < TNumNodes; ++i)
            {
                // Loop over columns
                for (unsigned int j = 0; j < TNumNodes; ++j)
                {
                    // Delta(u) * TauOne * [ AdvVel * Grad(v) ] in velocity block
                    K = Coef * AGradN[i] * rShapeFunc[j];

                    for (unsigned int d = 0; d < TDim; ++d) // iterate over dimensions for velocity Dofs in this node combination
                    {
                        rLHSMatrix(FirstRow + d, FirstCol + d) += K;
                        // Delta(u) * TauOne * Grad(q) in q * Div(u) block
                        rLHSMatrix(FirstRow + TDim, FirstCol + d) += Coef * rShapeDeriv(i, d) * rShapeFunc[j];
                    }
                    // Update column index
                    FirstCol += TBlockSize;
                }
                // Update matrix indices
                FirstRow += TBlockSize;
                FirstCol = 0;
            }
        }

        /// Add a the contribution from a single integration point to the velocity contribution
        void AddIntegrationPointVelocityContribution(MatrixType& rDampMatrix,
                                                     VectorType& rDampRHS,
                                                     const double Density,
                                                     const double KinViscosity,
                                                     const double TauOne,
                                                     const double TauTwo,
                                                     const array_1d< double, TNumNodes >& rShapeFunc,
                                                     const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim >& rShapeDeriv,
                                                     const double Weight)
        {
            // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
            array_1d<double, TNumNodes> AGradN;
            this->GetConvectionOperator(AGradN, rShapeFunc, rShapeDeriv); // Get a * grad(Ni)

            // Build the local matrix and RHS
            unsigned int FirstRow(0), FirstCol(0); // position of the first term of the local matrix that corresponds to each node combination
            double K, G, PDivV, L; // Temporary results
            double Coef = Density * Weight;

            // Note that we iterate first over columns, then over rows to read the Body Force only once per node
            for (unsigned int j = 0; j < TNumNodes; ++j) // iterate over colums
            {
                // Get Body Force
                const array_1d<double, 3 > & rBodyForce = this->GetGeometry()[j].FastGetSolutionStepValue(BODY_FORCE);

                for (unsigned int i = 0; i < TNumNodes; ++i) // iterate over rows
                {
                    // Calculate the part of the contributions that is constant for each node combination

                    // Velocity block
                    K = rShapeFunc[i] * AGradN[j]; // Convective term: v * ( a * Grad(u) )
                    K += TauOne * AGradN[i] * AGradN[j]; // Stabilization: (a * Grad(v)) * TauOne * (a * Grad(u))

                    // q-p stabilization block (reset result)
                    L = 0;

                    for (unsigned int m = 0; m < TDim; ++m) // iterate over v components (vx,vy[,vz])
                    {
                        // Velocity block
                        K += KinViscosity * rShapeDeriv(i, m) * rShapeDeriv(j, m); // Diffusive term: Viscosity * Grad(v) * Grad(u)
                        // Note that we are usig kinematic viscosity, as we will multiply it by density later

                        // v * Grad(p) block
                        G = TauOne * AGradN[i] * rShapeDeriv(j, m); // Stabilization: (a * Grad(v)) * TauOne * Grad(p)
                        PDivV = rShapeDeriv(i, m) * rShapeFunc[j]; // Div(v) * p

                        // Write v * Grad(p) component
                        rDampMatrix(FirstRow + m, FirstCol + TDim) = Weight * (G - PDivV);
                        // Use symmetry to write the q * rho * Div(u) component
                        rDampMatrix(FirstCol + TDim, FirstRow + m) = Coef * (G + PDivV);

                        // q-p stabilization block
                        L += rShapeDeriv(i, m) * rShapeDeriv(j, m); // Stabilization: Grad(q) * TauOne * Grad(p)

                        for (unsigned int n = 0; n < TDim; ++n) // iterate over u components (ux,uy[,uz])
                        {
                            // Velocity block
                            rDampMatrix(FirstRow + m, FirstCol + n) = Coef * TauTwo * rShapeDeriv(i, m) * rShapeDeriv(j, n); // Stabilization: Div(v) * TauTwo * Div(u)
                        }

                    }

                    // Write remaining terms to velocity block
                    K *= Coef; // Weight by nodal area and density
                    for (unsigned int d = 0; d < TDim; ++d)
                        rDampMatrix(FirstRow + d, FirstCol + d) += K;

                    // Write q-p stabilization block
                    rDampMatrix(FirstRow + TDim, FirstCol + TDim) = Weight * TauOne * L;

                    // Operate on RHS
                    L = 0; // We reuse one of the temporary variables for the pressure RHS

                    for (unsigned int d = 0; d < TDim; ++d)
                    {
                        rDampRHS[FirstRow + d] += Coef * TauOne * AGradN[i] * rShapeFunc[j] * rBodyForce[d]; // ( a * Grad(v) ) * TauOne * (Density * BodyForce)
                        L += rShapeDeriv(i, d) * rShapeFunc[j] * rBodyForce[d];
                    }
                    rDampRHS[FirstRow + TDim] += Coef * TauOne * L; // Grad(q) * TauOne * (Density * BodyForce)

                    // Update reference row index for next iteration
                    FirstRow += TBlockSize;
                }

                // Update reference indices
                FirstRow = 0;
                FirstCol += TBlockSize;
            }
        }

        /// Assemble the contribution from an integration point to the element's residual
        void AddProjectionResidualContribution(const array_1d< double, 3 > & rAdvVel,
                                               const double Density,
                                               const double TauOne,
                                               const double TauTwo,
                                               array_1d< double,3 >& rElementalMomRes,
                                               double& rElementalMassRes,
                                               const ProcessInfo& rCurrentProcessInfo,
                                               const array_1d< double, TNumNodes >& rShapeFunc,
                                               const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim >& rShapeDeriv,
                                               const double Weight)
        {
            // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
            array_1d<double, TNumNodes> AGradN;
            this->GetConvectionOperator(AGradN, rShapeFunc, rShapeDeriv); // Get a * grad(Ni)

            const double WeightedMass = Weight * Density;

            ///@TODO: Check that this is really safe for multiple Gauss points. Suggestion: return the elemental residual, write to nodes in caller function.
            // Compute contribution to Kij * Uj, with Kij = Ni * Residual(Nj); Uj = (v,p)Node_j (column vector)
            for (unsigned int i = 0; i < TNumNodes; ++i) // Iterate over element nodes
            {

                // Variable references
                const array_1d< double, 3 > & rVelocity = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
                const array_1d< double, 3 > & rAcceleration = this->GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION);
                const array_1d< double, 3 > & rBodyForce = this->GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
                const double& rPressure = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);

                // Compute this node's contribution to the residual (evaluated at inegration point)
                // Pressure contribution to the momentum residual is done separately because it is not computed from array_1d variables
                for (unsigned int d = 0; d < TDim; ++d)
                {
//                    ElementalMomRes[d] += Weight * ( Density * ( /*rShapeFunc[i] * (rAcceleration[d]-rBodyForce[d])*/ + AGradN[i] * rVelocity[d] ) + rShapeDeriv(i,d) * rPressure);
                    ///@TODO: Body force or not Body force?
                    rElementalMomRes[d] += Weight * (Density * (rShapeFunc[i] * (rAcceleration[d] - rBodyForce[d]) + AGradN[i] * rVelocity[d]) + rShapeDeriv(i, d) * rPressure);
                    rElementalMassRes += WeightedMass * rShapeDeriv(i, d) * rVelocity[d];
                }
            }
        }

        /// Write the advective velocity evaluated at this point to an array
        /**
         * Writes the value of the advective velocity evaluated at a point inside
         * the element to an array_1d
         * @param rAdvVel: Output array
         * @param rShapeFunc: Shape functions evaluated at the point of interest
         * @param Step: The time Step (Defaults to 0 = Current)
         */
        virtual void GetAdvectiveVel(array_1d< double, 3 > & rAdvVel,
                             const array_1d< double, TNumNodes >& rShapeFunc,
                             const std::size_t Step = 0)
        {
            // Compute the weighted value of the advective velocity in the (Gauss) Point
            rAdvVel = rShapeFunc[0] * (this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,Step) - this->GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY,Step));
            for (unsigned int iNode = 1; iNode < TNumNodes; ++iNode)
                rAdvVel += rShapeFunc[iNode] * (this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY,Step) - this->GetGeometry()[iNode].FastGetSolutionStepValue(MESH_VELOCITY,Step));
        }

        /// Write the convective operator evaluated at this point (for each nodal funciton) to an array
        /**
         * Evaluate the convective operator for each node's shape function at an arbitrary point
         * @param rResult: Output vector
         * @param rShapeFunc: Shape functions evaluated at the integration point
         * @param rShapeDeriv: Derivatives of shape functions evaluated at the integration point
         */
        void GetConvectionOperator(array_1d< double, TNumNodes >& rResult,
                                   const array_1d< double, TNumNodes>& rShapeFunc,
                                   const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim >& rShapeDeriv)
        {
            // Initialize result
            this->SetToZero<TNumNodes > (rResult);

            // Evaluate the convective velocity at the integration point
            array_1d< double, 3> AGauss;
            this->GetAdvectiveVel(AGauss, rShapeFunc);

            // Evaluate (and weight) the a * Grad(Ni) operator in the integration point, for each node i
            for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode) // Loop over nodes
                for (unsigned int d = 0; d < TDim; ++d) // loop over components
                    rResult[iNode] += AGauss[d] * rShapeDeriv(iNode, d);
        }

        /// Add the weighted value of a variable at a point inside the element to a double
        /**
         * Evaluate a scalar variable in the point where the form functions take the
         * values given by rShapeFunc and add the result, weighted by Weight, to
         * rResult. This is an auxiliary function used to compute values in integration
         * points.
         * @param rResult: The double where the value will be added to
         * @param rVariable: The nodal variable to be read
         * @param rShapeFunc: The values of the form functions in the point
         * @param Step: The time Step (Defaults to 0 = Current)
         * @param Weight: The variable will be weighted by this value before it is added to rResult
         */
        void AddPointContribution(double& rResult,
                                  const Variable< double >& rVariable,
                                  const array_1d< double, TNumNodes >& rShapeFunc,
                                  const std::size_t Step = 0,
                                  const double Weight = 1.0)
        {
            // Compute the weighted value of the nodal variable in the (Gauss) Point
            for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
                rResult += Weight * rShapeFunc[iNode] * this->GetGeometry()[iNode].FastGetSolutionStepValue(rVariable,Step);
        }

        /// Write the value of a variable at a point inside the element to a double
        /**
         * Evaluate a scalar variable in the point where the form functions take the
         * values given by rShapeFunc and write the result to rResult.
         * This is an auxiliary function used to compute values in integration points.
         * @param rResult: The double where the value will be added to
         * @param rVariable: The nodal variable to be read
         * @param rShapeFunc: The values of the form functions in the point
         * @param Step: The time Step (Defaults to 0 = Current)
         */
        void GetPointContribution(double& rResult,
                                  const Variable< double >& rVariable,
                                  const array_1d< double, TNumNodes >& rShapeFunc,
                                  const std::size_t Step = 0)
        {
            // Compute the weighted value of the nodal variable in the (Gauss) Point
            rResult = rShapeFunc[0] * this->GetGeometry()[0].FastGetSolutionStepValue(rVariable,Step);
            for (unsigned int iNode = 1; iNode < TNumNodes; ++iNode)
                rResult += rShapeFunc[iNode] * this->GetGeometry()[iNode].FastGetSolutionStepValue(rVariable,Step);
        }

        /// Add the weighted value of a variable at a point inside the element to a vector
        /**
         * Evaluate a vector variable in the point where the form functions take the
         * values given by rShapeFunc and add the result, weighted by Weight, to
         * rResult. This is an auxiliary function used to compute values in integration
         * points.
         * @param rResult: The vector where the value will be added to
         * @param rVariable: The nodal variable to be read
         * @param rShapeFunc: The values of the form functions in the point
         * @param Weight: The variable will be weighted by this value before it is added to rResult
         */
        inline void AddPointContribution(array_1d< double, 3 > & rResult,
                                         const Variable< array_1d< double, 3 > >& rVariable,
                                         const array_1d< double, TNumNodes>& rShapeFunc,
                                         const double Weight = 1.0)
        {
            // Compute the weighted value of the nodal variable in the (Gauss) Point
            for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
                rResult += Weight * rShapeFunc[iNode] * this->GetGeometry()[iNode].FastGetSolutionStepValue(rVariable);
        }

        /// Write the value of a variable at a point inside the element to a double
        /**
         * Evaluate a scalar variable in the point where the form functions take the
         * values given by rShapeFunc and write the result to rResult.
         * This is an auxiliary function used to compute values in integration points.
         * @param rResult: The double where the value will be added to
         * @param rVariable: The nodal variable to be read
         * @param rShapeFunc: The values of the form functions in the point
         */
        void GetPointContribution(array_1d< double, 3 > & rResult,
                                  const Variable< array_1d< double, 3 > >& rVariable,
                                  const array_1d< double, TNumNodes >& rShapeFunc)
        {
            // Compute the weighted value of the nodal variable in the (Gauss) Point
            rResult = rShapeFunc[0] * this->GetGeometry()[0].FastGetSolutionStepValue(rVariable);
            for (unsigned int iNode = 1; iNode < TNumNodes; ++iNode)
                rResult += rShapeFunc[iNode] * this->GetGeometry()[iNode].FastGetSolutionStepValue(rVariable);
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
        ///@name Private Operators
        ///@{


        ///@}
        ///@name Private Operations
        ///@{

        /// Set an array_1d to zero

        template< unsigned int TSize >
        inline void SetToZero(array_1d< double, TSize >& rArray)
        {
            for (typename array_1d< double, TSize>::iterator itArray = rArray.begin(); itArray != rArray.end(); ++itArray)
                *itArray = 0.0;
        }

        /// Return an estimate for the element size h, used to calculate the stabilization parameters
        /**
         * Estimate the element size from its area or volume, required to calculate stabilization parameters.
         * Note that this is a pure virtual function, as its implementation is different for 2D or 3D elements.
         * @see VMS2D, VMS3D for actual implementation
         * @param Volume (in 3D) or Area (in 2D) of the element
         * @return Element size h
         */
        virtual double ElementSize(const double) = 0;

        ///@}
        ///@name Private  Access
        ///@{


        ///@}
        ///@name Private Inquiry
        ///@{


        ///@}
        ///@name Un accessible methods
        ///@{

        /// Assignment operator.
        VMSBase & operator=(VMSBase const& rOther);

        /// Copy constructor.
        VMSBase(VMSBase const& rOther);

        ///@}

    }; // Class VMSBase

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function
    template< unsigned int TDim,
              unsigned int TNumNodes,
              unsigned int TBlockSize,
              unsigned int TLocalSize
            >
    inline std::istream & operator >>(std::istream& rIStream,
                                      VMSBase<TDim,TNumNodes,TBlockSize,TLocalSize>& rThis)
    {
        return rIStream;
    }

    /// output stream function
    template< unsigned int TDim,
              unsigned int TNumNodes,
              unsigned int TBlockSize,
              unsigned int TLocalSize
            >
    inline std::ostream & operator <<(std::ostream& rOStream,
                                      const VMSBase<TDim,TNumNodes,TBlockSize,TLocalSize>& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}


} // namespace Kratos.

//#undef KRATOS_VMS_ALT_TAU

#endif // KRATOS_VMS_BASE_H_INCLUDED  defined


