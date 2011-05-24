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
//   Date:                $Date: 2011-04-29 12:30:00 $
//   Revision:            $Revision: 0.2 $
//
//


#if !defined(KRATOS_DYNAMIC_VMS_H_INCLUDED )
#define  KRATOS_DYNAMIC_VMS_H_INCLUDED

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/element.h"
#include "vms.h"
#include "includes/serializer.h"
#include "utilities/geometry_utilities.h"

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

    /// A stabilized element for the incompressible Navier-Stokes equations.
    /**
     * @todo Rewrite this documentation
     * This class implements a stabilized formulation based on the
     * Variational Multiscale framework. The the subscales can be modeled
     * using either Algebraic Subgird Scales (ASGS) or Orthogonal Subscales (OSS).
     * In the case of OSS, the projection terms are treated explicitly (computed
     * using the results of the previous iteration) and the subscales are not
     * tracked in time. The choice of subscale model is made based on the ProcessInfo
     * variable OSS_SWITCH (OSS if 1, ASGS otherwise).
     * This class implements both the 2D and 3D versions of the element.
     *
     * The ASGS implementation follows Ramon Codina, A stabilized finite element
     * method for generalized stationary incompressible flows, Computer Methods in
     * Applied Mechanics and Engineering. Vol. 190 (2001), 2681-2706.
     *
     * The OSS implementation corresponds to the case identified as explicit, quasi-
     * static orthogonal subscales in Ramon Codina, Stabilized finite element approximation
     * of transient incompressible flows using orthogonal subscales, Computer Methods
     * in Applied Mechanics and Engineering. Vol. 191 (2002), 4295-4321.
     *
     * In addition to the stabilization, this element implements the Smagorinsky
     * model of turbulence. This turbulent term is only activated if the elemental
     * value C_SMAGORINSKY is set to something other than zero.
     *
     * This class requires at least the following variables:\n
     * On each Node, as solution step variables VELOCITY, PRESSURE, ACCELERATION, MESH_VELOCITY.\n
     * On ProcessInfo OSS_SWITCH, DELTA_TIME.\n
     * If OSS is used, the nodes also require NODAL_AREA, ADVPROJ and DIVPROJ as solution step variables.\n
     * If Smagorinsky is used, C_SMAGORINSKY has to be defined on the elements.\n
     * Error estimation stores ERROR_RATIO on the elements.\n
     * Some additional variables can be used to print results on the element: TAUONE, TAUTWO, MU, VORTICITY.
     *
     * @note Unlike VMS, this class does not use the DYNAMIC_TAU ProcessInfo value, which is always
     * assumed 1.0
     *
     * @see ResidualBasedEliminationBuilderAndSolver compatible monolithic solution strategy.
     * @see PressureSplittingBuilderAndSolver compatible segregated solution strategy.
     * @see TrilinosPressureSplittingBuilderAndSolver compatible mpi strategy.
     * @see DynamicSmagorinskyUtils to set the Smagorinsky parameter dynamically.
     * @see ResidualBasedPredictorCorrectorVelocityBossakScheme time scheme that can use
     * OSS stabilization.
     */
    template< unsigned int TDim,
              unsigned int TNumNodes = TDim + 1 >
    class DynamicVMS : public VMS<TDim,TNumNodes>
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of DynamicVMS
        KRATOS_CLASS_POINTER_DEFINITION(DynamicVMS);

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

        typedef VMS<TDim,TNumNodes> BaseElementType;

        ///@}
        ///@name Life Cycle
        ///@{

        //Constructors.

        /// Default constuctor.
        /**
         * @param NewId Index number of the new element (optional)
         */
        DynamicVMS(IndexType NewId = 0) :
            BaseElementType(NewId),
            mSubscaleVel(3,0.0),
            mOldSubscaleVel(3,0.0),
            mIterCount(0)
        {}

        ///Constructor using an array of nodes.
        /**
         * @param NewId Index of the new element
         * @param ThisNodes An array containing the nodes of the new element
         */
        DynamicVMS(IndexType NewId, const NodesArrayType& ThisNodes) :
            BaseElementType(NewId, ThisNodes),
            mSubscaleVel(3,0.0),
            mOldSubscaleVel(3,0.0),
            mIterCount(0)
        {}

        /// Constructor using a geometry object.
        /**
         * @param NewId Index of the new element
         * @param pGeometry Pointer to a geometry object
         */
        DynamicVMS(IndexType NewId, GeometryType::Pointer pGeometry) :
            BaseElementType(NewId, pGeometry),
            mSubscaleVel(3,0.0),
            mOldSubscaleVel(3,0.0),
            mIterCount(0)
        {}

        /// Constuctor using geometry and properties.
        /**
         * @param NewId Index of the new element
         * @param pGeometry Pointer to a geometry object
         * @param pProperties Pointer to the element's properties
         */
        DynamicVMS(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
            BaseElementType(NewId, pGeometry, pProperties),
            mSubscaleVel(3,0.0),
            mOldSubscaleVel(3,0.0),
            mIterCount(0)
        {}

        /// Destructor.
        virtual ~DynamicVMS()
        {}


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        /// Create a new element of this type
        /**
         * Returns a pointer to a new DynamicVMS element, created using given input
         * @param NewId: the ID of the new element
         * @param ThisNodes: the nodes of the new element
         * @param pProperties: the properties assigned to the new element
         * @return a Pointer to the new element
         */
        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                PropertiesType::Pointer pProperties) const
        {
            return Element::Pointer(new DynamicVMS(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
        }

        /// Update the stored subscale values
        /***/
        virtual void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
        {
            mOldSubscaleVel = mSubscaleVel;
        }

        /// Calculate a new value for the velocity subscale
        /**
         * @param rCurrentProcessInfo ProcessInfo instance containig the time step as DELTA_TIME
         */
        virtual void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
        {
            // Calculate this element's geometric parameters
            GeometryType& rGeom = this->GetGeometry();
            double Area;
            array_1d<double, TNumNodes> N;
            boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
            GeometryUtils::CalculateGeometryData(rGeom, DN_DX, N, Area);

            // Calculate this element's fluid properties
            double Density, KinViscosity;
            this->EvaluateInPoint(Density, DENSITY, N);
            this->EvaluateInPoint(KinViscosity, VISCOSITY, N);

            double Viscosity;
            this->GetEffectiveViscosity(Density,KinViscosity, N, DN_DX, Viscosity, rCurrentProcessInfo);

            // Get Advective velocity
            array_1d<double,3> AdvVel;
            BaseElementType::GetAdvectiveVel(AdvVel, N); // We use the base function here, which does not add the subscale

            // Evaluate the part of the residual that depends on v_h (and will remain constant in this function)
            array_1d<double,3> MomRes(3,0.0);
            if (rCurrentProcessInfo[OSS_SWITCH] == 1)
                // OSS version, includes projection term
                this->OSSMomResidual(AdvVel,Density,MomRes,N,DN_DX,Area);
            else
                // ASGS version, includes dynamic term (which is orthogonal to subscale space for OSS)
                this->ASGSMomResidual(AdvVel,Density,MomRes,N,DN_DX,Area);

            // Some constant values for iteration
            const double DeltaTime = rCurrentProcessInfo[DELTA_TIME];
            double InvStep = 1.0 / DeltaTime;
            const double ElemSize = this->ElementSize(Area);

            MatrixType Gradient = ZeroMatrix(TDim,TDim); // Velocity gradient

            array_1d<double,3> rVelocity;

            for (unsigned int k = 0; k < TNumNodes; ++k)
            {
                rVelocity = rGeom[k].FastGetSolutionStepValue(VELOCITY) - rGeom[k].FastGetSolutionStepValue(MESH_VELOCITY);
                for (unsigned int i = 0; i < TDim; ++i)
                    for (unsigned int j = 0; j < TDim; ++j)
                        Gradient(i,j) += DN_DX(k,j)*rVelocity[i];
            }

            Gradient *= Density;

            // Iteration variables
            VectorType Subscale = ZeroVector(TDim); ///@todo Ensure I'm using dense matrix types here.
            VectorType Vel(TDim);
            for (unsigned int d = 0; d < TDim; ++d)
                Vel[d] = AdvVel[d] + mSubscaleVel[d]; // Here I'm taking the subscale from last iteration as initial guess

            double InvTau;
            double VelNorm;

            double SubscaleNorm = 1e10;
            double SubscaleError = 1e10;
            mIterCount = 0;

            MatrixType J(TDim,TDim);
            VectorType r(TDim);
            VectorType dx = ZeroVector(TDim);

            while ( SubscaleError > mSubscaleTol && SubscaleNorm > mSubscaleTol && mIterCount <= 10 ) ///@todo add some absolute criteria too
            {
                // Update advection velocity
                noalias(Vel) += dx;

                // Update velocity norm
                VelNorm = Vel[0] * Vel[0];
                for (unsigned int d = 1; d < TDim; ++d)
                    VelNorm += Vel[d] * Vel[d];
                VelNorm = sqrt(VelNorm);

                // Update Tau
                InvTau = this->InverseTau(Viscosity,VelNorm,ElemSize,DeltaTime); // InvTau includes the 1/Dt term
                InvTau *= Density;

                // Build system using uBLAS types
                noalias(J) = Gradient;
                for (unsigned int d = 0; d < TDim; ++d)
                {
                    // Build LHS
                    J(d,d) += InvTau;
                    // Build RHS
                    r[d] = MomRes[d] + Density * mOldSubscaleVel[d] * InvStep; ///@todo This is constant and can be moved out of the solution loop
                }
                // Finish RHS: r -= J * Subscale
                noalias(r) -= prod(J,Subscale);

                // Solve system
                this->DenseSystemSolve(J,dx,r);

                // Update solution
                noalias(Subscale) += dx;

                // Update Error
                SubscaleError = dx[0]*dx[0];
                SubscaleNorm = Subscale[0]*Subscale[0];
                for(unsigned int d = 1; d < TDim; ++d)
                {
                    SubscaleError += dx[d]*dx[d];
                    SubscaleNorm += Subscale[0]*Subscale[0];
                }
                SubscaleError /= SubscaleNorm;
                
                // Iteration counter
                ++mIterCount;
            }

            // Store the converged subscale
            for(unsigned int d = 0; d < TDim; ++d)
                mSubscaleVel[d] = Subscale[d];
        }


        /// Provides local contributions from body forces and projections to the RHS
        /**
         * This is called during the assembly process and provides the RHS terms of the
         * system that are either constant or treated explicitly (from the 'old'
         * iteration variables).
         * @param rRightHandSideVector Will be filled with the elemental right hand side
         * @param rCurrentProcessInfo ProcessInfo instance from the ModelPart. It is
         * expected to contain values for OSS_SWITCH and DELTA_TIME
         */
        virtual void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                            ProcessInfo& rCurrentProcessInfo)
        {
            const unsigned int LocalSize = (TDim + 1) * TNumNodes;

            // Check sizes and initialize
            if (rRightHandSideVector.size() != LocalSize)
                rRightHandSideVector.resize(LocalSize, false);

            noalias(rRightHandSideVector) = ZeroVector(LocalSize);

            // Calculate this element's geometric parameters
            double Area;
            array_1d<double, TNumNodes> N;
            boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            // Calculate this element's fluid properties
            double Density;
            this->EvaluateInPoint(Density, DENSITY, N);

            // Calculate Momentum RHS contribution
            this->AddMomentumRHS(rRightHandSideVector, Density, N, Area);

            if (rCurrentProcessInfo[OSS_SWITCH] == 1)
            {
                // OSS: Add projection of residuals to RHS
                array_1d<double, 3 > AdvVel;
                this->GetAdvectiveVel(AdvVel, N);

                double KinViscosity;
                this->EvaluateInPoint(KinViscosity, VISCOSITY, N);

                double Viscosity;
                this->GetEffectiveViscosity(Density,KinViscosity, N, DN_DX, Viscosity, rCurrentProcessInfo);

                // Calculate stabilization parameters
                double TauOne, TauTwo, MassFactor;
                this->CalculateTau(TauOne,TauTwo,MassFactor,AdvVel,Area,Viscosity,rCurrentProcessInfo);

                this->AddProjectionToRHS(rRightHandSideVector,AdvVel,TauOne,TauTwo,MassFactor,N,DN_DX,Area);
            }
            else
            {
                // ASGS: Include velocity subscale derivative term
                this->AddDynamicSubscaleTerm(rRightHandSideVector,Density,rCurrentProcessInfo[DELTA_TIME],N,Area);
            }
        }


        /// Computes the local contribution associated to 'new' velocity and pressure values
        /**
         * Provides local contributions to the system associated to the velocity and
         * pressure terms (convection, diffusion, pressure gradient/velocity divergence
         * and stabilization).
         * @param rDampMatrix Will be filled with the velocity-proportional "damping" matrix
         * @param rRightHandSideVector the elemental right hand side vector
         * @param rCurrentProcessInfo the current process info instance
         */
        virtual void CalculateLocalVelocityContribution(MatrixType& rDampMatrix,
                                                        VectorType& rRightHandSideVector,
                                                        ProcessInfo& rCurrentProcessInfo)
        {
            const unsigned int LocalSize = (TDim + 1) * TNumNodes;

            // Resize and set to zero the matrix
            // Note that we don't clean the RHS because it will already contain body force (and stabilization) contributions
            if (rDampMatrix.size1() != LocalSize)
                rDampMatrix.resize(LocalSize, LocalSize, false);

            noalias(rDampMatrix) = ZeroMatrix(LocalSize, LocalSize);

            // Get this element's geometric properties
            double Area;
            array_1d<double, TNumNodes> N;
            boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            // Calculate this element's fluid properties
            double Density, KinViscosity;
            this->EvaluateInPoint(Density, DENSITY, N);
            this->EvaluateInPoint(KinViscosity, VISCOSITY, N);

            double Viscosity;
            this->GetEffectiveViscosity(Density,KinViscosity, N, DN_DX, Viscosity, rCurrentProcessInfo);

            // Get Advective velocity
            array_1d<double, 3 > AdvVel;
            this->GetAdvectiveVel(AdvVel, N);

            // Calculate stabilization parameters
            double TauOne, TauTwo, MassFactor;
            this->CalculateTau(TauOne, TauTwo, MassFactor, AdvVel, Area, Viscosity, rCurrentProcessInfo);

            this->AddIntegrationPointVelocityContribution(rDampMatrix,rRightHandSideVector,Density,Viscosity,AdvVel,TauOne,TauTwo,rCurrentProcessInfo[DELTA_TIME],MassFactor,N,DN_DX,Area);

            // Now calculate an additional contribution to the residual: r -= rDampMatrix * (u,p)
            VectorType U = ZeroVector(LocalSize);
            int LocalIndex = 0;

            for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
            {
                array_1d< double, 3 > & rVel = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
                for (unsigned int d = 0; d < TDim; ++d) // Velocity Dofs
                {
                    U[LocalIndex] = rVel[d];
                    ++LocalIndex;
                }
                U[LocalIndex] = this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE); // Pressure Dof
                ++LocalIndex;
            }

            noalias(rRightHandSideVector) -= prod(rDampMatrix, U);
        }


        /// Implementation of Calculate to compute an error estimate.
        /**
         * If rVariable == ERROR_RATIO, this function will provide an a posteriori
         * estimate of the norm of the subscale velocity, calculated as TauOne*||MomentumResidual||.
         * Note that the residual of the momentum equation is evaluated at the element center
         * and that the result has units of velocity (L/T).
         * The error estimate both saved as the elemental ERROR_RATIO variable and returned as rOutput.
         * @param rVariable Use ERROR_RATIO
         * @param rOutput Returns the error estimate
         * @param rCurrentProcessInfo Process info instance (will be checked for OSS_SWITCH)
         * @see MarkForRefinement for a use of the error ratio
         */
        virtual void Calculate(const Variable<double>& rVariable,
                               double& rOutput,
                               const ProcessInfo& rCurrentProcessInfo)
        {
            if (rVariable == ERROR_RATIO)
            {
                rOutput = 0.0;

                for (unsigned int i = 0; i < TDim; ++i)
                {
                    rOutput += mSubscaleVel[i] * mSubscaleVel[i];
//                    UNorm += UGauss[i] * UGauss[i];
                }
                rOutput = sqrt(rOutput); // / UNorm);
                this->SetValue(ERROR_RATIO, rOutput);
            }
        }

//        /// Implementation of Calculate to compute the local OSS projections.
//        /**
//         * If rVariable == ADVPROJ, This function computes the OSS projection
//         * terms using pressure and velocity values from the previous iteration. The
//         * projections are then added to the nodal variables ADVPROJ (Momentum residual)
//         * and DIVPROJ (Mass continuity residual). It is assumed that the scheme will
//         * divide the result by the assembled NODAL_AREA, which is equivalent to a
//         * nodal interpolation using a lumped mass matrix.
//         * @param rVariable Use ADVPROJ
//         * @param Output Will be overwritten with the elemental momentum error
//         * @param rCurrentProcessInfo Process info instance (unused)
//         */
//        virtual void Calculate(const Variable<array_1d<double, 3 > >& rVariable,
//                               array_1d<double, 3 > & rOutput,
//                               const ProcessInfo& rCurrentProcessInfo)
//        {
//            if (rVariable == ADVPROJ) // Compute residual projections for OSS
//            {
//                // Get the element's geometric parameters
//                double Area;
//                array_1d<double, TNumNodes> N;
//                boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
//                GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);
//
//                // Calculate this element's fluid properties
//                double Density, KinViscosity;
//                this->EvaluateInPoint(Density, DENSITY, N);
//                this->EvaluateInPoint(KinViscosity, VISCOSITY, N);
//
//                double Viscosity;
//                this->GetEffectiveViscosity(Density,KinViscosity, N, DN_DX, Viscosity, rCurrentProcessInfo);
//
//                // Get Advective velocity
//                array_1d<double, 3 > AdvVel;
//                this->GetAdvectiveVel(AdvVel, N);
//
//                // Output containers
//                array_1d< double, 3 > ElementalMomRes(3, 0.0);
//                double ElementalMassRes(0);
//
//                this->AddProjectionResidualContribution(AdvVel, Density, ElementalMomRes, ElementalMassRes, N, DN_DX, Area);
//
//                if (rCurrentProcessInfo[OSS_SWITCH] == 1)
//                {
//                    // Carefully write results to nodal variables, to avoid parallelism problems
//                    for (unsigned int i = 0; i < TNumNodes; ++i)
//                    {
//                        this->GetGeometry()[i].SetLock(); // So it is safe to write in the node in OpenMP
//                        array_1d< double, 3 > & rAdvProj = this->GetGeometry()[i].FastGetSolutionStepValue(ADVPROJ);
//                        for (unsigned int d = 0; d < TDim; ++d)
//                            rAdvProj[d] += N[i] * ElementalMomRes[d];
//
//                        this->GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ) += N[i] * ElementalMassRes;
//                        this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += Area * N[i];
//                        this->GetGeometry()[i].UnSetLock(); // Free the node for other threads
//                    }
//                }
//
//                /// Return output
//                rOutput = ElementalMomRes;
//            }
//        }


        /// Obtain an array_1d<double,3> elemental variable, evaluated on gauss points.
        /**
         * If the variable is VORTICITY, computes the vorticity (rotational of the velocity)
         * based on the current velocity values. Otherwise, it assumes that the input
         * variable is an elemental value and retrieves it. Implemented for a
         * single gauss point only.
         * @param rVariable Kratos vector variable to get
         * @param Output Will be filled with the values of the variable on integrartion points
         * @param rCurrentProcessInfo Process info instance
         */
        virtual void GetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                                                 std::vector<array_1d<double, 3 > >& rOutput,
                                                 const ProcessInfo& rCurrentProcessInfo);

        /// Obtain a double elemental variable, evaluated on gauss points.
        /**
         * If the variable is TAUONE or TAUTWO, calculates the corresponding stabilization
         * parameter for the element, based on rCurrentProcessInfo's DELTA_TIME and
         * DYNAMIC_TAU. If the variable is MU, calculates the effective viscosity at the
         * element center due to Smagorinsky (in 'dynamic' units). Otherwise, it assumes
         * that the input variable is an elemental value and retrieves it.
         * Implemented for a single gauss point only.
         * @param rVariable Kratos vector variable to compute
         * @param Output Will be filled with the values of the variable on integrartion points
         * @param rCurrentProcessInfo Process info instance
         */
        virtual void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                                                 std::vector<double>& rValues,
                                                 const ProcessInfo& rCurrentProcessInfo)
        {
            if (rVariable == NL_ITERATION_NUMBER)
            {
                rValues.resize(1, false);
                rValues[0] = mIterCount;
            }
            else
            {
                BaseElementType::GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
            }
        }

        /// Empty implementation of unused CalculateOnIntegrationPoints overloads to avoid compilation warning
        virtual void GetValueOnIntegrationPoints(const Variable<array_1d<double, 6 > >& rVariable,
                                                 std::vector<array_1d<double, 6 > >& rValues,
                                                 const ProcessInfo& rCurrentProcessInfo)
        {}

        /// Empty implementation of unused CalculateOnIntegrationPoints overloads to avoid compilation warning
        virtual void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                 std::vector<Vector>& rValues,
                                                 const ProcessInfo& rCurrentProcessInfo)
        {}

        /// Empty implementation of unused CalculateOnIntegrationPoints overloads to avoid compilation warning
        virtual void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                 std::vector<Matrix>& rValues,
                                                 const ProcessInfo& rCurrentProcessInfo)
        {}

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
            buffer << "DynamicVMS #" << Element::Id();
            return buffer.str();
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "DynamicVMS" << TDim << "D";
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

        /// Calculate Stabilization parameters.
        /**
         * Calculates both tau parameters based on a given advective velocity.
         * @param TauOne First stabilization parameter (momentum equation)
         * @param TauTwo Second stabilization parameter (mass equation)
         * @param rAdvVel Advective (linearized) velocity
         * @param Area Elemental area
         * @param KinViscosity Elemental kinematic viscosity (nu)
         * @param rCurrentProcessInfo Process info instance
         * (containing the time step as DELTA_TIME)
         */
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

            TauOne = 1.0 / ( 1.0 / rCurrentProcessInfo[DELTA_TIME] + 4.0 * KinViscosity / (Element_Size * Element_Size) + 2.0 * AdvVelNorm / Element_Size);
            TauTwo = KinViscosity + 0.5 * Element_Size * AdvVelNorm;
        }

        /// Calculate Stabilization parameters.
        /**
         * Calculates both tau parameters and their ratio (as MassResidualFactor)
         * based on a given advective velocity.
         * @param TauOne First stabilization parameter (momentum equation)
         * @param TauTwo Second stabilization parameter (mass equation)
         * @param MassResidualFactor Factor that multiplies the dynamic part of
         * the pressure subscale equation, equivalent to StaticTauOne * TauTwo / Dt
         * (Where StaticTauOne does not include the 1/Dt term)
         * @param rAdvVel Advective (linearized) velocity
         * @param Area Elemental area
         * @param KinViscosity Elemental kinematic viscosity (nu)
         * @param rCurrentProcessInfo Process info instance
         * (containing the time step as DELTA_TIME)
         * @see DynamicVMS::CalculateRightHandSide
         */
        virtual void CalculateTau(double& TauOne,
                                  double& TauTwo,
                                  double& MassResidualFactor,
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

            TauOne = 1.0 / ( 1.0 / rCurrentProcessInfo[DELTA_TIME] + 4.0 * KinViscosity / (Element_Size * Element_Size) + 2.0 * AdvVelNorm / Element_Size);
            TauTwo = KinViscosity + 0.5 * Element_Size * AdvVelNorm;
            MassResidualFactor = Element_Size * Element_Size / ( 4.0 * rCurrentProcessInfo[DELTA_TIME] );
        }


        /// Add OSS projection terms to the RHS
        void AddProjectionToRHS(VectorType& RHS,
                                const array_1d<double, 3 > & rAdvVel,
                                const double TauOne,
                                const double TauTwo,
                                const double MassFactor,
                                const array_1d<double, TNumNodes>& rShapeFunc,
                                const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim>& rShapeDeriv,
                                const double Weight)
        {
            const unsigned int BlockSize = TDim + 1;

            array_1d<double, TNumNodes> AGradN;
            this->GetConvectionOperator(AGradN, rAdvVel, rShapeDeriv); // Get a * grad(Ni)

            // Add to the ouptut vector
            unsigned int FirstRow(0); // Position of the term we want to add to
            double Const1, Const2, Const3; // Partial results that remain constant for a given j

            for (unsigned int j = 0; j < TNumNodes; ++j) // loop over nodes (for components of the residual)
            {
                // Compute the 'constant' part of the (nodal) projection terms
                // See this->Calculate() for the computation of the projections
                const array_1d< double, 3 > & rMomProj = this->GetGeometry()[j].FastGetSolutionStepValue(ADVPROJ);

                Const1 = Weight * TauOne * rShapeFunc[j];
                Const2 = Weight * (TauTwo + MassFactor) * rShapeFunc[j] * this->GetGeometry()[j].FastGetSolutionStepValue(DIVPROJ);
                Const3 = Weight * MassFactor * rShapeFunc[j] * this->GetGeometry()[j].FastGetSolutionStepValue(DIVPROJ,1);

                // Reset row reference
                FirstRow = 0;

                for (unsigned int i = 0; i < TNumNodes; ++i) // loop over nodes (for components of L*(Vh) )
                {
                    for (unsigned int d = 0; d < TDim; ++d)
                    {
                        // TauOne * ( a * Grad(v) ) * MomProjection + TauTwo * Div(v) * MassProjection
                        RHS[FirstRow + d] += Const1 * AGradN[i] * rMomProj[d] + (Const2 - Const3) * rShapeDeriv(i, d);
                        // TauOne * Grad(q) * MomProjection
                        RHS[FirstRow + TDim] += Const1 * rShapeDeriv(i, d) * rMomProj[d];
                    }
                    // Update row reference
                    FirstRow += BlockSize;
                }
            }
        }

        /// Add the stabilization term coming from the derivative of the subscale velocity
        /**
         * Note that this term is zero (and therefore not computed) in OSS, as subscales
         * are orthogonal to test functions.
         * @param rRightHandSideVector Vector storing RHS for the elemental system
         * @param Density Fluid density on integration point
         * @param TimeStep
         * @param rShapeFunc Shape functions evaluated on integration point
         * @param Weight Integration weight (as a fraction of total area or volume)
         */
        virtual void AddDynamicSubscaleTerm(VectorType& rRightHandSideVector,
                                            const double Density,
                                            const double TimeStep,
                                            const array_1d<double,TNumNodes>& rShapeFunc,
                                            const double Weight)
        {
            const double Coef = Density * Weight / TimeStep;

            // Add the results to the velocity components (Local Dofs are vx, vy, [vz,] p for each node)
            int LocalIndex = 0;

            for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
            {
                for (unsigned int d = 0; d < TDim; ++d)
                {
                    // RHS -= Density * v * (SubscaleVel - OldSubscaleVel) / Dt
                    rRightHandSideVector[LocalIndex++] += Coef * rShapeFunc[iNode] * (mOldSubscaleVel[d] - mSubscaleVel[d]);
                }
                ++LocalIndex; // Skip pressure Dof
            }
        }


        /// Add a the contribution from a single integration point to the velocity contribution
        void AddIntegrationPointVelocityContribution(MatrixType& rDampMatrix,
                                                     VectorType& rDampRHS,
                                                     const double Density,
                                                     const double Viscosity,
                                                     const array_1d< double, 3 > & rAdvVel,
                                                     const double TauOne,
                                                     const double TauTwo,
                                                     const double MassFactor,
                                                     const double TimeStep,
                                                     const array_1d< double, TNumNodes >& rShapeFunc,
                                                     const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim >& rShapeDeriv,
                                                     const double Weight)
        {
            const unsigned int BlockSize = TDim + 1;

            // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
            array_1d<double, TNumNodes> AGradN;
            this->GetConvectionOperator(AGradN, rAdvVel, rShapeDeriv); // Get a * grad(Ni)

            // Build the local matrix and RHS
            unsigned int FirstRow(0), FirstCol(0); // position of the first term of the local matrix that corresponds to each node combination
            double K, G, PDivV, L, DivDiv; // Temporary results
            const double MassStab = TauTwo + MassFactor; // "Dynamic TauTwo" including extra TauOne*TauTwo/Dt term coming from dynamic tracking of subscales
            const double InvTimeStep = 1.0 / TimeStep; // Velocity subscale approximated as TauOne*(MomResProj + Density * OldSubscale/DeltaTime)
            double Coef = Density * Weight;

            // Note that we iterate first over columns, then over rows to read the Body Force only once per node
            for (unsigned int j = 0; j < TNumNodes; ++j) // iterate over colums
            {
                // Get Body Force
                const array_1d<double, 3 > & rBodyForce = this->GetGeometry()[j].FastGetSolutionStepValue(BODY_FORCE);
                const array_1d<double, 3 > & rOldVel = this->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY,1);

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
//                        K += Viscosity * rShapeDeriv(i, m) * rShapeDeriv(j, m); // Diffusive term: Viscosity * Grad(v) * Grad(u)
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
                            DivDiv = rShapeDeriv(i, m) * rShapeDeriv(j, n);
                            // Velocity block
                            rDampMatrix(FirstRow + m, FirstCol + n) = Coef * MassStab * DivDiv; // Stabilization: Div(v) * TauTwo * ( Div(u) + StaticTauOne * Div(u)/Dt
                            rDampRHS(FirstRow+m) += Coef * MassFactor * DivDiv * rOldVel[n];// Stab RHS: Div(v) * TauTwo * StaticTauOne * Div(u_old)/Dt
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
                        // ( a * Grad(v) ) * TauOne * Density * (BodyForce - OldSubscaleVel/Dt)
                        rDampRHS[FirstRow + d] += Coef * TauOne * AGradN[i] * (rShapeFunc[j] * rBodyForce[d] - mOldSubscaleVel[d] * InvTimeStep);
                        L += rShapeDeriv(i, d) * (rShapeFunc[j] * rBodyForce[d] - mOldSubscaleVel[d] * InvTimeStep);
                    }
                    rDampRHS[FirstRow + TDim] += Coef * TauOne * L; // Grad(q) * TauOne * Density * (BodyForce - OldSubscaleVel/Dt)

                    // Update reference row index for next iteration
                    FirstRow += BlockSize;
                }

                // Update reference indices
                FirstRow = 0;
                FirstCol += BlockSize;
            }

//            this->AddBTransCB(rDampMatrix,rShapeDeriv,Viscosity*Coef);
            this->AddViscousTerm(rDampMatrix,rShapeDeriv,Viscosity*Coef);
        }


        /// Write the advective velocity evaluated at this point to an array
        /**
         * Writes the value of the advective velocity evaluated at a point inside
         * the element to an array_1d. Note that the advective velocity includes
         * the effect of the subscale velocity.
         * @param rAdvVel: Output array
         * @param rShapeFunc: Shape functions evaluated at the point of interest
         * @param Step: The time Step (Defaults to 0 = Current)
         */
        virtual void GetAdvectiveVel(array_1d< double,3> & rAdvVel,
                                     const array_1d< double, TNumNodes >& rShapeFunc,
                                     const std::size_t Step = 0)
        {
            // Compute the weighted value of the advective velocity in the (Gauss) Point
            rAdvVel = mSubscaleVel;
            for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
                rAdvVel += rShapeFunc[iNode] * (this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, Step) - this->GetGeometry()[iNode].FastGetSolutionStepValue(MESH_VELOCITY, Step));
        }


        /// Efficient evaluation of 1/TauOne to use in the subscale velocity calculation
        /**
         *
         * @param Viscosity Effective kinematic viscosity
         * @param VelNorm Advection velocity norm (including subscale)
         * @param ElemSize Element size (h)
         * @param DeltaTime Time step
         * @return 1/TauOne
         * @see DynamicVMS::CalculateTau, DynamicVMS::InitializeNonLinearIteration
         */
        double InverseTau(const double Viscosity,
                          const double VelNorm,
                          const double ElemSize,
                          const double DeltaTime)
        {
            return 1.0 / DeltaTime + 4.0 * Viscosity / (ElemSize*ElemSize) + 2.0 * VelNorm / ElemSize;
        }

        /// Solve a linear system of TDim linear equations by computing the algebraic inverse.
        /**
         * This function is a wrapper for calls to MathUtils::InvertMatrix.
         * Matrix and vector sizes are not checked.
         * @param rA System matrix of size TDim x TDim
         * @param rx vector of unknowns (size TDim)
         * @param rb right hand side vector (size TDim)
         */
        void DenseSystemSolve(const MatrixType& rA, VectorType& rx, const VectorType& rb);


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
        
        /// Number of Dofs on each node
        static const unsigned int msLocalSize;
        
        /// Total number of DOFs on the element
        static const unsigned int msBlockSize;

        /// Tolerance for subscale iterations
        static const double mSubscaleTol;

        ///@}
        ///@name Member Variables
        ///@{

        /// Subscale velocity evaluated on the integration point
        array_1d<double,3> mSubscaleVel;

        /// Subscale velocity on integration point obtained on last time step
        array_1d<double,3> mOldSubscaleVel;

        /// Iteration count for the non-linear velocity subscale loop
        unsigned int mIterCount;


        ///@}
        ///@name Serialization
        ///@{

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const;

        virtual void load(Serializer& rSerializer)
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseElementType);
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

        /// Assignment operator.
        DynamicVMS & operator=(DynamicVMS const& rOther);

        /// Copy constructor.
        DynamicVMS(DynamicVMS const& rOther);

        ///@}

    }; // Class DynamicVMS

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function
    template< unsigned int TDim,
              unsigned int TNumNodes >
    inline std::istream& operator >>(std::istream& rIStream,
                                     DynamicVMS<TDim, TNumNodes>& rThis)
    {
        return rIStream;
    }

    /// output stream function
    template< unsigned int TDim,
              unsigned int TNumNodes >
    inline std::ostream& operator <<(std::ostream& rOStream,
                                     const DynamicVMS<TDim, TNumNodes>& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}

    ///@} // Fluid Dynamics Application group
    
    template< unsigned int TDim, unsigned int TNumNodes >
    const unsigned int DynamicVMS<TDim,TNumNodes>::msBlockSize = TDim + 1; // TDim velocity dofs + 1 pressure dof

    template< unsigned int TDim, unsigned int TNumNodes >
    const unsigned int DynamicVMS<TDim,TNumNodes>::msLocalSize = TNumNodes * (TDim + 1); // NumNodes * LocalSize

    template< unsigned int TDim, unsigned int TNumNodes >
    const double DynamicVMS<TDim,TNumNodes>::mSubscaleTol = 1e-6;

} // namespace Kratos.

#endif // KRATOS_DYNAMIC_VMS__H_INCLUDED  defined
