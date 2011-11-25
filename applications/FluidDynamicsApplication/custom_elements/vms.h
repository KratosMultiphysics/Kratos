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


#if !defined(KRATOS_VMS_H_INCLUDED )
#define  KRATOS_VMS_H_INCLUDED

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/element.h"
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
     * On each Node, as solution step variables VELOCITY, PRESSURE, ACCELERATION, MESH_VELOCITY, DENSITY, VISCOSITY.\n
     * On ProcessInfo OSS_SWITCH, DYNAMIC_TAU, DELTA_TIME.\n
     * If OSS is used, the nodes also require NODAL_AREA, ADVPROJ and DIVPROJ as solution step variables.\n
     * If Smagorinsky is used, C_SMAGORINSKY has to be defined on the elements.\n
     * Error estimation stores ERROR_RATIO on the elements.\n
     * Some additional variables can be used to print results on the element: SUBSCALE_VELOCITY, SUBSCALE_PRESSURE, TAUONE, TAUTWO, MU, VORTICITY.
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
    class VMS : public Element
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of VMS
        KRATOS_CLASS_POINTER_DEFINITION(VMS);

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

        //Constructors.

        /// Default constuctor.
        /**
         * @param NewId Index number of the new element (optional)
         */
        VMS(IndexType NewId = 0) :
            Element(NewId)
        {}

        ///Constructor using an array of nodes.
        /**
         * @param NewId Index of the new element
         * @param ThisNodes An array containing the nodes of the new element
         */
        VMS(IndexType NewId, const NodesArrayType& ThisNodes) :
            Element(NewId, ThisNodes)
        {}

        /// Constructor using a geometry object.
        /**
         * @param NewId Index of the new element
         * @param pGeometry Pointer to a geometry object
         */
        VMS(IndexType NewId, GeometryType::Pointer pGeometry) :
            Element(NewId, pGeometry)
        {}

        /// Constuctor using geometry and properties.
        /**
         * @param NewId Index of the new element
         * @param pGeometry Pointer to a geometry object
         * @param pProperties Pointer to the element's properties
         */
        VMS(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
            Element(NewId, pGeometry, pProperties)
        {}

        /// Destructor.
        virtual ~VMS()
        {}


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        /// Create a new element of this type
        /**
         * Returns a pointer to a new VMS element, created using given input
         * @param NewId: the ID of the new element
         * @param ThisNodes: the nodes of the new element
         * @param pProperties: the properties assigned to the new element
         * @return a Pointer to the new element
         */
        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                PropertiesType::Pointer pProperties) const
        {
            return Element::Pointer(new VMS(NewId, GetGeometry().Create(ThisNodes), pProperties));
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
            const unsigned int LocalSize = (TDim + 1) * TNumNodes;

            // Check sizes and initialize matrix
            if (rLeftHandSideMatrix.size1() != LocalSize)
                rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

            noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);

            // Calculate RHS
            this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
        }

        /// Returns a zero matrix of appropiate size (provided for compatibility with scheme)
        /**
         * @param rLeftHandSideMatrix Local matrix, will be filled with zeros
         * @param rCurrentProcessInfo Process info instance
         */
        virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                           ProcessInfo& rCurrentProcessInfo)
        {
            const unsigned int LocalSize = (TDim + 1) * TNumNodes;

            if (rLeftHandSideMatrix.size1() != LocalSize)
                rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

            noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
        }

        /// Provides local contributions from body forces and projections to the RHS
        /**
         * This is called during the assembly process and provides the RHS terms of the
         * system that are either constant or computed explicitly (from the 'old'
         * iteration variables). In this case this means the body force terms and the
         * OSS projections, that are treated explicitly.
         * @param rRightHandSideVector Will be filled with the elemental right hand side
         * @param rCurrentProcessInfo ProcessInfo instance from the ModelPart. It is
         * expected to contain values for OSS_SWITCH, DYNAMIC_TAU and DELTA_TIME
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

            // For OSS: Add projection of residuals to RHS
            if (rCurrentProcessInfo[OSS_SWITCH] == 1)
            {
                array_1d<double, 3 > AdvVel;
                this->GetAdvectiveVel(AdvVel, N);

                double KinViscosity;
                this->EvaluateInPoint(KinViscosity, VISCOSITY, N);

                double Viscosity;
                this->GetEffectiveViscosity(Density,KinViscosity, N, DN_DX, Viscosity, rCurrentProcessInfo);

                // Calculate stabilization parameters
                double TauOne, TauTwo;
                this->CalculateTau(TauOne, TauTwo, AdvVel, Area, Viscosity, rCurrentProcessInfo);

                this->AddProjectionToRHS(rRightHandSideVector, AdvVel, TauOne, TauTwo, N, DN_DX, Area,rCurrentProcessInfo[DELTA_TIME]);
            }
            else if (this->GetValue(TRACK_SUBSCALES) == 1)// Experimental: Dynamic tracking of ASGS subscales (see Codina 2002 Stabilized finite element ... using orthogonal subscales)
            {
                /* We want to evaluate v * d(u_subscale)/dt. This term is zero in OSS, due the orthogonality of the two terms. in ASGS, we approximate it as
                 * d(u_s)/dt = u_s(last iteration) - u_s(previous time step) / DeltaTime where u_s(last iteration) = TauOne(MomResidual(last_iteration) + u_s(previous step)/DeltaTime)
                 */
                array_1d<double, 3 > AdvVel;
                this->GetAdvectiveVel(AdvVel, N);

                double KinViscosity;
                this->EvaluateInPoint(KinViscosity, VISCOSITY, N);

                double Viscosity;
                this->GetEffectiveViscosity(Density,KinViscosity, N, DN_DX, Viscosity, rCurrentProcessInfo);

                double StaticTauOne;
                this->CalculateStaticTau(StaticTauOne,AdvVel,Area,Viscosity);

                double TauOne,TauTwo;
                this->CalculateTau(TauOne,TauTwo,AdvVel,Area,Viscosity,rCurrentProcessInfo);

                array_1d<double,3> ElemMomRes;
                this->ASGSMomResidual(AdvVel,Density,ElemMomRes,N,DN_DX,1.0);

                const array_1d<double,3>& rOldSubscale = this->GetValue(SUBSCALE_VELOCITY);
                const double DeltaTime = rCurrentProcessInfo[DELTA_TIME];

//                array_1d<double,3> Subscale = TauOne * (ElemMomRes + Density * rOldSubscale / DeltaTime);
                const double C1 = 1.0 - TauOne/StaticTauOne;
                const double C2 = TauOne / (StaticTauOne*DeltaTime);

                unsigned int LocalIndex = 0;
                for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
                {
                    for(unsigned int d = 0; d < TDim; d++)
                        rRightHandSideVector[LocalIndex++] -= N[iNode] * Area * ( C1 * ElemMomRes[d] - C2 * rOldSubscale[d] );
//                        rRightHandSideVector[LocalIndex++] -= N[iNode] * Area * ( Subscale[d] - rOldSubscale[d] ) / DeltaTime;
                    LocalIndex++; // Pressure Dof
                }
            }
        }

        /// Computes local contributions to the mass matrix
        /**
         * Provides the local contributions to the mass matrix, which is defined here
         * as the matrix associated to velocity derivatives. Note that the mass
         * matrix implemented here is lumped.
         * @param rMassMatrix Will be filled with the elemental mass matrix
         * @param rCurrentProcessInfo the current process info instance
         */
        virtual void MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
        {
            const unsigned int LocalSize = (TDim + 1) * TNumNodes;

            // Resize and set to zero
            if (rMassMatrix.size1() != LocalSize)
                rMassMatrix.resize(LocalSize, LocalSize, false);

            rMassMatrix = ZeroMatrix(LocalSize, LocalSize);

            // Get the element's geometric parameters
            double Area;
            array_1d<double, TNumNodes> N;
            boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            // Calculate this element's fluid properties
            double Density;
            this->EvaluateInPoint(Density, DENSITY, N);

            // Add 'classical' mass matrix (lumped)
            double Coeff = Density * Area / TNumNodes; //Optimize!
            this->CalculateLumpedMassMatrix(rMassMatrix, Coeff);

            /* For ASGS: add dynamic stabilization terms.
             These terms are not used in OSS, as they belong to the finite element
             space and cancel out with their projections.
             */
            if (rCurrentProcessInfo[OSS_SWITCH] != 1)
            {
                double KinViscosity;
                this->EvaluateInPoint(KinViscosity, VISCOSITY, N);

                double Viscosity;
                this->GetEffectiveViscosity(Density,KinViscosity, N, DN_DX, Viscosity, rCurrentProcessInfo);

                // Get Advective velocity
                array_1d<double, 3 > AdvVel;
                this->GetAdvectiveVel(AdvVel, N);

                // Calculate stabilization parameters
                double TauOne, TauTwo;
                this->CalculateTau(TauOne, TauTwo, AdvVel, Area, Viscosity, rCurrentProcessInfo);

                // Add dynamic stabilization terms ( all terms involving a delta(u) )
                this->AddMassStabTerms<MatrixType > (rMassMatrix, Density, AdvVel, TauOne, N, DN_DX, Area);
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
            double TauOne, TauTwo;
            this->CalculateTau(TauOne, TauTwo, AdvVel, Area, Viscosity, rCurrentProcessInfo);

            if(this->GetValue(TRACK_SUBSCALES)==1)
            {
                const double DeltaTime = rCurrentProcessInfo[DELTA_TIME];
                this->AddIntegrationPointVelocityContribution(rDampMatrix, rRightHandSideVector, Density, Viscosity, AdvVel, TauOne, TauTwo, N, DN_DX, Area,DeltaTime);
            }
            else
            {
                this->AddIntegrationPointVelocityContribution(rDampMatrix, rRightHandSideVector, Density, Viscosity, AdvVel, TauOne, TauTwo, N, DN_DX, Area);
            }

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

        virtual void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
        {
            if(this->GetValue(TRACK_SUBSCALES) == 1)
            {
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
                double StaticTauOne;
                this->CalculateStaticTau(StaticTauOne,AdvVel,Area,Viscosity);

                array_1d<double,3> ElementalMomRes(3,0.0);

                if ( rCurrentProcessInfo[OSS_SWITCH] != 1 ) // ASGS
                {
                    this->ASGSMomResidual(AdvVel,Density,ElementalMomRes,N,DN_DX,1.0);
                }
                else // OSS
                {
                    this->OSSMomResidual(AdvVel,Density,ElementalMomRes,N,DN_DX,1.0);;
                }

                // Update subscale term
                const double DeltaTime = rCurrentProcessInfo.GetValue(DELTA_TIME);
                array_1d<double,3>& OldSubscaleVel = this->GetValue(SUBSCALE_VELOCITY);
                array_1d<double,3> Tmp = ( 1.0/( 1.0/DeltaTime + 1.0/StaticTauOne)) * (Density*OldSubscaleVel/DeltaTime + ElementalMomRes);
                OldSubscaleVel = Tmp;
            }
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
                // Get the element's geometric parameters
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

                // Output container
                array_1d< double, 3 > ElementalMomRes(3, 0.0);

                // Calculate stabilization parameter. Note that to estimate the subscale velocity, the dynamic coefficient in TauOne is assumed zero.
                double TauOne;
                this->CalculateStaticTau(TauOne, AdvVel, Area, Viscosity);

                if ( rCurrentProcessInfo[OSS_SWITCH] != 1 ) // ASGS
                {
                    this->ASGSMomResidual(AdvVel,Density,ElementalMomRes,N,DN_DX,1.0);
                    ElementalMomRes *= TauOne;
                }
                else // OSS
                {
                    this->OSSMomResidual(AdvVel,Density,ElementalMomRes,N,DN_DX,1.0);;
                    ElementalMomRes *= TauOne;
                }

                // Error estimation ( ||U'|| / ||Uh_gauss|| ), taking ||U'|| = TauOne ||MomRes||
                double ErrorRatio(0.0);//, UNorm(0.0);
//                array_1d< double, 3 > UGauss(3, 0.0);
//                this->AddPointContribution(UGauss, VELOCITY, N);

                for (unsigned int i = 0; i < TDim; ++i)
                {
                    ErrorRatio += ElementalMomRes[i] * ElementalMomRes[i];
//                    UNorm += UGauss[i] * UGauss[i];
                }
                ErrorRatio = sqrt(ErrorRatio); // / UNorm);
                ErrorRatio /= Density;
                this->SetValue(ERROR_RATIO, ErrorRatio);
                rOutput = ErrorRatio;
            }
        }

        /// Implementation of Calculate to compute the local OSS projections.
        /**
         * If rVariable == ADVPROJ, This function computes the OSS projection
         * terms using pressure and velocity values from the previous iteration. The
         * projections are then added to the nodal variables ADVPROJ (Momentum residual)
         * and DIVPROJ (Mass continuity residual). It is assumed that the scheme will
         * divide the result by the assembled NODAL_AREA, which is equivalent to a
         * nodal interpolation using a lumped mass matrix.
         * @param rVariable Use ADVPROJ
         * @param Output Will be overwritten with the elemental momentum error
         * @param rCurrentProcessInfo Process info instance (unused)
         */
        virtual void Calculate(const Variable<array_1d<double, 3 > >& rVariable,
                               array_1d<double, 3 > & rOutput,
                               const ProcessInfo& rCurrentProcessInfo)
        {
            if (rVariable == ADVPROJ) // Compute residual projections for OSS
            {
                // Get the element's geometric parameters
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

                // Output containers
                array_1d< double, 3 > ElementalMomRes(3, 0.0);
                double ElementalMassRes(0);

                this->AddProjectionResidualContribution(AdvVel, Density, ElementalMomRes, ElementalMassRes, N, DN_DX, Area);

                if (rCurrentProcessInfo[OSS_SWITCH] == 1)
                {
                    // Carefully write results to nodal variables, to avoid parallelism problems
                    for (unsigned int i = 0; i < TNumNodes; ++i)
                    {
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
            else if (rVariable == SUBSCALE_VELOCITY)
            {
                // Get the element's geometric parameters
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

                // Output containers
                array_1d< double, 3 > ElementalMomRes(3,0.0);
                double ElementalMassRes(0.0);

                this->AddProjectionResidualContribution(AdvVel, Density, ElementalMomRes, ElementalMassRes, N, DN_DX, Area);

                if (rCurrentProcessInfo[OSS_SWITCH] == 1)
                {
                    /* Projections of the elemental residual are computed with
                     * Newton-Raphson iterations of type M(lumped) dx = ElemRes - M(consistent) * x
                     */
                    const double Weight = ConsitentMassCoef(Area); // Consistent mass matrix is Weigth * ( Ones(TNumNodes,TNumNodes) + Identity(TNumNodes,TNumNodes) )
                    // Carefully write results to nodal variables, to avoid parallelism problems
                    for (unsigned int i = 0; i < TNumNodes; ++i)
                    {
                        this->GetGeometry()[i].SetLock(); // So it is safe to write in the node in OpenMP

                        // Add elemental residual to RHS
                        array_1d< double, 3 > & rMomRHS = this->GetGeometry()[i].GetValue(ADVPROJ);
                        double& rMassRHS = this->GetGeometry()[i].GetValue(DIVPROJ);
                        for (unsigned int d = 0; d < TDim; ++d)
                            rMomRHS[d] += N[i] * ElementalMomRes[d];

                        rMassRHS += N[i] * ElementalMassRes;

                        // Write nodal area
                        this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += Area * N[i];

                        // Substract M(consistent)*x(i-1) from RHS
                        for(unsigned int j = 0; j < TNumNodes; ++j) // RHS -= Weigth * Ones(TNumNodes,TNumNodes) * x(i-1)
                        {
                            for(unsigned int d = 0; d < TDim; ++d)
                                rMomRHS[d] -= Weight * this->GetGeometry()[j].FastGetSolutionStepValue(ADVPROJ)[d];
                            rMassRHS -= Weight * this->GetGeometry()[j].FastGetSolutionStepValue(DIVPROJ);
                        }
                        for(unsigned int d = 0; d < TDim; ++d) // RHS -= Weigth * Identity(TNumNodes,TNumNodes) * x(i-1)
                            rMomRHS[d] -= Weight * this->GetGeometry()[i].FastGetSolutionStepValue(ADVPROJ)[d];
                        rMassRHS -= Weight * this->GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ);

                        this->GetGeometry()[i].UnSetLock(); // Free the node for other threads
                    }
                }

                /// Return output
                rOutput = ElementalMomRes;
            }
        }

        // The following methods have different implementations depending on TDim
        /// Provides the global indices for each one of this element's local rows
        /**
         * this determines the elemental equation ID vector for all elemental
         * DOFs
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
        virtual void GetDofList(DofsVectorType& rElementalDofList,
                                ProcessInfo& rCurrentProcessInfo);

        /// Returns VELOCITY_X, VELOCITY_Y, (VELOCITY_Z,) PRESSURE for each node
        /**
         * @param Values Vector of nodal unknowns
         * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
         */
        virtual void GetFirstDerivativesVector(Vector& Values, int Step = 0);

        /// Returns ACCELERATION_X, ACCELERATION_Y, (ACCELERATION_Z,) 0 for each node
        /**
         * @param Values Vector of nodal second derivatives
         * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
         */
        virtual void GetSecondDerivativesVector(Vector& Values, int Step = 0);

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
            if (rVariable == TAUONE || rVariable == TAUTWO || rVariable == MU)
            {
                double TauOne, TauTwo;
                double Area;
                array_1d<double, TNumNodes> N;
                boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
                GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

                array_1d<double, 3 > AdvVel;
                this->GetAdvectiveVel(AdvVel, N);

                double Density,KinViscosity;
                this->EvaluateInPoint(Density, DENSITY, N);
                this->EvaluateInPoint(KinViscosity, VISCOSITY, N);

                double Viscosity;
                this->GetEffectiveViscosity(Density,KinViscosity, N, DN_DX, Viscosity, rCurrentProcessInfo);

                this->CalculateTau(TauOne, TauTwo, AdvVel, Area, Viscosity, rCurrentProcessInfo);

                rValues.resize(1, false);
                if (rVariable == TAUONE)
                {
                    rValues[0] = TauOne;
                }
                else if (rVariable == TAUTWO)
                {
                    rValues[0] = TauTwo;
                }
                else if (rVariable == MU)
                {
                    rValues[0] = Density * Viscosity;
                }
            }
            else if(rVariable == SUBSCALE_PRESSURE)
            {
                double TauOne, TauTwo;
                double Area;
                array_1d<double, TNumNodes> N;
                boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
                GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

                array_1d<double, 3 > AdvVel;
                this->GetAdvectiveVel(AdvVel, N);

                double Density,KinViscosity;
                this->EvaluateInPoint(Density, DENSITY, N);
                this->EvaluateInPoint(KinViscosity, VISCOSITY, N);

                double Viscosity;
                this->GetEffectiveViscosity(Density,KinViscosity, N, DN_DX, Viscosity, rCurrentProcessInfo);

                this->CalculateTau(TauOne, TauTwo, AdvVel, Area, Viscosity, rCurrentProcessInfo);

                double DivU = 0.0;
                for(unsigned int i=0; i < TNumNodes;i++)
                {
                    for(unsigned int d = 0; d < TDim; d++)
                        DivU -= DN_DX(i,d) * this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[d];
                }

                rValues.resize(1, false);
                rValues[0] = TauTwo * DivU;// *Density?? decide on criteria and use the same for SUBSCALE_VELOCITY
                if(rCurrentProcessInfo[OSS_SWITCH]==1)
                {
                    double Proj = 0.0;
                    for(unsigned int i=0; i < TNumNodes;i++)
                    {
                        Proj += N[i]*this->GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ);
                    }
                    rValues[0] -= TauTwo*Proj;
                }
            }
            else // Default behaviour (returns elemental data)
            {
                rValues.resize(1, false);
                /*
                 The cast is done to avoid modification of the element's data. Data modification
                 would happen if rVariable is not stored now (would initialize a pointer to &rVariable
                 with associated value of 0.0). This is catastrophic if the variable referenced
                 goes out of scope.
                 */
                const VMS<TDim, TNumNodes>* const_this = static_cast<const VMS<TDim, TNumNodes>*> (this);
                rValues[0] = const_this->GetValue(rVariable);
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
        ///@name Elemental Data
        ///@{

        /// Checks the input and that all required Kratos variables have been registered.
        /**
         * This function provides the place to perform checks on the completeness of the input.
         * It is designed to be called only once (or anyway, not often) typically at the beginning
         * of the calculations, so to verify that nothing is missing from the input
         * or that no common error is found.
         * @param rCurrentProcessInfo The ProcessInfo of the ModelPart that contains this element.
         * @return 0 if no errors were found.
         */
        virtual int Check(const ProcessInfo& rCurrentProcessInfo)
        {
            KRATOS_TRY

            // Perform basic element checks
            int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
            if(ErrorCode != 0) return ErrorCode;

            // Check that all required variables have been registered
            if(VELOCITY.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"VELOCITY Key is 0. Check if the application was correctly registered.","")
            if(MESH_VELOCITY.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"MESH_VELOCITY Key is 0. Check if the application was correctly registered.","")
            if(ACCELERATION.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"ACCELERATION Key is 0. Check if the application was correctly registered.","")
            if(PRESSURE.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"PRESSURE Key is 0. Check if the application was correctly registered.","")
            if(DENSITY.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"DENSITY Key is 0. Check if the application was correctly registered.","")
            if(VISCOSITY.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"VISCOSITY Key is 0. Check if the application was correctly registered.","")
            if(OSS_SWITCH.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"OSS_SWITCH Key is 0. Check if the application was correctly registered.","")
            if(DYNAMIC_TAU.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"DYNAMIC_TAU Key is 0. Check if the application was correctly registered.","")
            if(DELTA_TIME.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"DELTA_TIME Key is 0. Check if the application was correctly registered.","")
            if(ADVPROJ.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"ADVPROJ Key is 0. Check if the application was correctly registered.","")
            if(DIVPROJ.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"DIVPROJ Key is 0. Check if the application was correctly registered.","")
            if(NODAL_AREA.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"NODAL_AREA Key is 0. Check if the application was correctly registered.","")
            if(C_SMAGORINSKY.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"C_SMAGORINSKY Key is 0. Check if the application was correctly registered.","")
            if(ERROR_RATIO.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"ERROR_RATIO Key is 0. Check if the application was correctly registered.","")
            // Additional variables, only required to print results:
            // SUBSCALE_VELOCITY, SUBSCALE_PRESSURE, TAUONE, TAUTWO, MU, VORTICITY.

            // Checks on nodes

            // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
            for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
            {
                if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
                    KRATOS_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
                if(this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
                    KRATOS_ERROR(std::invalid_argument,"missing PRESSURE variable on solution step data for node ",this->GetGeometry()[i].Id());
                if(this->GetGeometry()[i].SolutionStepsDataHas(MESH_VELOCITY) == false)
                    KRATOS_ERROR(std::invalid_argument,"missing MESH_VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
                if(this->GetGeometry()[i].SolutionStepsDataHas(ACCELERATION) == false)
                    KRATOS_ERROR(std::invalid_argument,"missing ACCELERATION variable on solution step data for node ",this->GetGeometry()[i].Id());
                if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
                   this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
                   this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
                    KRATOS_ERROR(std::invalid_argument,"missing VELOCITY component degree of freedom on node ",this->GetGeometry()[i].Id());
                if(this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
                    KRATOS_ERROR(std::invalid_argument,"missing PRESSURE component degree of freedom on node ",this->GetGeometry()[i].Id());
            }
            // Not checking OSS related variables NODAL_AREA, ADVPROJ, DIVPROJ, which are only required as SolutionStepData if OSS_SWITCH == 1

            // If this is a 2D problem, check that nodes are in XY plane
            if (this->GetGeometry().WorkingSpaceDimension() == 2)
            {
                for (unsigned int i=0; i<this->GetGeometry().size(); ++i)
                {
                    if (this->GetGeometry()[i].Z() != 0.0)
                        KRATOS_ERROR(std::invalid_argument,"Node with non-zero Z coordinate found. Id: ",this->GetGeometry()[i].Id());
                }
            }

            return 0;

            KRATOS_CATCH("");
        }

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
            buffer << "VMS #" << Id();
            return buffer.str();
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "VMS" << TDim << "D";
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
         * Takes time step and dynamic coefficient from given ProcessInfo instance.
         * ProcessInfo variables DELTA_TIME and DYNAMIC_TAU will be used.
         * @param TauOne First stabilization parameter (momentum equation)
         * @param TauTwo Second stabilization parameter (mass equation)
         * @param rAdvVel advection velocity
         * @param Area Elemental area
         * @param KinViscosity Elemental kinematic viscosity (nu)
         * @param rCurrentProcessInfo Process info instance
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

            TauOne = 1.0 / (rCurrentProcessInfo[DYNAMIC_TAU] / rCurrentProcessInfo[DELTA_TIME] + 4.0 * KinViscosity / (Element_Size * Element_Size) + 2.0 * AdvVelNorm / Element_Size);
            TauTwo = KinViscosity + 0.5 * Element_Size * AdvVelNorm;
        }

        /// Calculate momentum stabilization parameter (without time term).
        /**
         * Calculates the momentum tau parameter based on a given advective velocity.
         * The dynamic term is not taken into account. This function
         * is intended for error estimation only. In other cases use CalculateTau
         * @param TauOne First stabilization parameter (momentum equation)
         * @param rAdvVel advection velocity
         * @param Area Elemental area
         * @param KinViscosity Elemental kinematic viscosity (nu)
         */
        virtual void CalculateStaticTau(double& TauOne,
                                        const array_1d< double, 3 > & rAdvVel,
                                        const double Area,
                                        const double KinViscosity)
        {
            // Compute mean advective velocity norm
            double AdvVelNorm = 0.0;
            for (unsigned int d = 0; d < TDim; ++d)
                AdvVelNorm += rAdvVel[d] * rAdvVel[d];

            AdvVelNorm = sqrt(AdvVelNorm);

            const double Element_Size = this->ElementSize(Area);

            TauOne = 1.0 / (4.0 * KinViscosity / (Element_Size * Element_Size) + 2.0 * AdvVelNorm / Element_Size);
        }

        /// Add the momentum equation contribution to the RHS (body forces)
        void AddMomentumRHS(VectorType& F,
                            const double Density,
                            const array_1d<double, TNumNodes>& rShapeFunc,
                            const double Weight)
        {
            double Coef = Density * Weight;

            array_1d<double, 3 > BodyForce(3, 0.0);
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
                                const array_1d<double, 3 > & rAdvVel,
                                const double TauOne,
                                const double TauTwo,
                                const array_1d<double, TNumNodes>& rShapeFunc,
                                const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim>& rShapeDeriv,
                                const double Weight,
                                const double DeltaTime = 1.0)
        {

            // TO BE REMOVED ---------------------------------------------------------------------------------------------------
            // Experimental: Dynamic tracking of subscales (see Codina 2002 Stabilized finite element ... using orthogonal subscales)
            const array_1d<double,3>& OldSubscale = this->GetValue(SUBSCALE_VELOCITY);
            //------------------------------------------------------------------------------------------------------------------

            const unsigned int BlockSize = TDim + 1;

            array_1d<double, TNumNodes> AGradN;
            this->GetConvectionOperator(AGradN, rAdvVel, rShapeDeriv); // Get a * grad(Ni)

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
                        RHS[FirstRow + d] -= Const1 * AGradN[i] * rMomProj[d] + Const2 * rShapeDeriv(i, d); // TauOne * ( a * Grad(v) ) * MomProjection + TauTwo * Div(v) * MassProjection
                        RHS[FirstRow + TDim] -= Const1 * rShapeDeriv(i, d) * rMomProj[d]; // TauOne * Grad(q) * MomProjection
                        if(this->GetValue(TRACK_SUBSCALES)==1)
                        {
                            RHS[FirstRow + d] += Const1 * AGradN[i] * TauOne * OldSubscale[d]/DeltaTime;
                            RHS[FirstRow + TDim] += Const1 * rShapeDeriv(i, d) * TauOne * OldSubscale[d]/DeltaTime;
                        }
                    }
                    // Update row reference
                    FirstRow += BlockSize;
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


        /// Add mass-like stabilization terms to LHS.
        /**
         * This function is only used in ASGS. For OSS, we avoid computing these
         * terms, as they shoud cancel out with the dynamic part of the projection
         * (which is not computed either)
         * @param rLHSMatrix Left hand side of the velocity-pressure system
         * @param Density Density on integration point
         * @param rAdvVel Advective velocity on integration point
         * @param TauOne Stabilization parameter for momentum equation
         * @param rShapeFunc Shape funcitions evaluated on integration point
         * @param rShapeDeriv Shape function derivatives evaluated on integration point
         * @param Weight Area (or volume) times integration point weight
         */
        template < class TMatrixType >
        void AddMassStabTerms(TMatrixType& rLHSMatrix,
                              const double Density,
                              const array_1d<double, 3 > & rAdvVel,
                              const double TauOne,
                              const array_1d<double, TNumNodes>& rShapeFunc,
                              const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim>& rShapeDeriv,
                              const double Weight)
        {
            const unsigned int BlockSize = TDim + 1;

            double Coef = Density * Weight * TauOne;
            unsigned int FirstRow(0), FirstCol(0);
            double K; // Temporary results

            // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
            array_1d<double, TNumNodes> AGradN;
            this->GetConvectionOperator(AGradN, rAdvVel, rShapeDeriv); // Get a * grad(Ni)

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
                    FirstCol += BlockSize;
                }
                // Update matrix indices
                FirstRow += BlockSize;
                FirstCol = 0;
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
//                        K += Viscosity * rShapeDeriv(i, m) * rShapeDeriv(j, m); // Diffusive term: Viscosity * Grad(v) * Grad(u)
                        // Note that we are usig kinematic viscosity, as we will multiply it by density later

                        // v * Grad(p) block
                        G = TauOne * AGradN[i] * rShapeDeriv(j, m); // Stabilization: (a * Grad(v)) * TauOne * Grad(p)
                        PDivV = rShapeDeriv(i, m) * rShapeFunc[j]; // Div(v) * p

                        // Write v * Grad(p) component
                        rDampMatrix(FirstRow + m, FirstCol + TDim) += Weight * (G - PDivV);
                        // Use symmetry to write the q * rho * Div(u) component
                        rDampMatrix(FirstCol + TDim, FirstRow + m) += Coef * (G + PDivV);

                        // q-p stabilization block
                        L += rShapeDeriv(i, m) * rShapeDeriv(j, m); // Stabilization: Grad(q) * TauOne * Grad(p)

                        for (unsigned int n = 0; n < TDim; ++n) // iterate over u components (ux,uy[,uz])
                        {
                            // Velocity block
                            rDampMatrix(FirstRow + m, FirstCol + n) += Coef * TauTwo * rShapeDeriv(i, m) * rShapeDeriv(j, n); // Stabilization: Div(v) * TauTwo * Div(u)
                        }

                    }

                    // Write remaining terms to velocity block
                    K *= Coef; // Weight by nodal area and density
                    for (unsigned int d = 0; d < TDim; ++d)
                        rDampMatrix(FirstRow + d, FirstCol + d) += K;

                    // Write q-p stabilization block
                    rDampMatrix(FirstRow + TDim, FirstCol + TDim) += Weight * TauOne * L;

                    // Operate on RHS
                    L = 0; // We reuse one of the temporary variables for the pressure RHS

                    for (unsigned int d = 0; d < TDim; ++d)
                    {
                        rDampRHS[FirstRow + d] += Coef * TauOne * AGradN[i] * rShapeFunc[j] * rBodyForce[d]; // ( a * Grad(v) ) * TauOne * (Density * BodyForce)
                        L += rShapeDeriv(i, d) * rShapeFunc[j] * rBodyForce[d];
                    }
                    rDampRHS[FirstRow + TDim] += Coef * TauOne * L; // Grad(q) * TauOne * (Density * BodyForce)

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

        /// Add a the contribution from a single integration point to the velocity contribution
        void AddIntegrationPointVelocityContribution(MatrixType& rDampMatrix,
                                                     VectorType& rDampRHS,
                                                     const double Density,
                                                     const double Viscosity,
                                                     const array_1d< double, 3 > & rAdvVel,
                                                     const double TauOne,
                                                     const double TauTwo,
                                                     const array_1d< double, TNumNodes >& rShapeFunc,
                                                     const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim >& rShapeDeriv,
                                                     const double Weight,
                                                     const double DeltaTime)
        {
            const unsigned int BlockSize = TDim + 1;

            const double TauCoef = TauOne*TauTwo / DeltaTime;

            // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
            array_1d<double, TNumNodes> AGradN;
            this->GetConvectionOperator(AGradN, rAdvVel, rShapeDeriv); // Get a * grad(Ni)

            // Build the local matrix and RHS
            unsigned int FirstRow(0), FirstCol(0); // position of the first term of the local matrix that corresponds to each node combination
            double K, G, PDivV, L; // Temporary results
            double Coef = Density * Weight;

            // Note that we iterate first over columns, then over rows to read the Body Force only once per node
            for (unsigned int j = 0; j < TNumNodes; ++j) // iterate over colums
            {
                // Get Body Force
                const array_1d<double, 3 > & rBodyForce = this->GetGeometry()[j].FastGetSolutionStepValue(BODY_FORCE);

                const array_1d<double,3>& OldVelocity = this->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY,1);

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
                        rDampMatrix(FirstRow + m, FirstCol + TDim) += Weight * (G - PDivV);
                        // Use symmetry to write the q * rho * Div(u) component
                        rDampMatrix(FirstCol + TDim, FirstRow + m) += Coef * (G + PDivV);

                        // q-p stabilization block
                        L += rShapeDeriv(i, m) * rShapeDeriv(j, m); // Stabilization: Grad(q) * TauOne * Grad(p)

                        for (unsigned int n = 0; n < TDim; ++n) // iterate over u components (ux,uy[,uz])
                        {
                            // Velocity block
                            rDampMatrix(FirstRow + m, FirstCol + n) += Coef * (TauTwo + TauCoef) * rShapeDeriv(i, m) * rShapeDeriv(j, n); // Stabilization: Div(v) * TauTwo *( 1+TauOne/Dt) * Div(u)
                            rDampRHS[FirstRow + m] -= Coef * TauCoef * rShapeDeriv(i, m) * rShapeDeriv(j, n) * OldVelocity[n]; // Stabilization: Div(v) * TauTwo*TauOne/Dt * Div(u_old)
                        }

                    }

                    // Write remaining terms to velocity block
                    K *= Coef; // Weight by nodal area and density
                    for (unsigned int d = 0; d < TDim; ++d)
                        rDampMatrix(FirstRow + d, FirstCol + d) += K;

                    // Write q-p stabilization block
                    rDampMatrix(FirstRow + TDim, FirstCol + TDim) += Weight * TauOne * L;

                    // Operate on RHS
                    L = 0; // We reuse one of the temporary variables for the pressure RHS

                    for (unsigned int d = 0; d < TDim; ++d)
                    {
                        rDampRHS[FirstRow + d] += Coef * TauOne * AGradN[i] * rShapeFunc[j] * rBodyForce[d]; // ( a * Grad(v) ) * TauOne * (Density * BodyForce)
                        L += rShapeDeriv(i, d) * rShapeFunc[j] * rBodyForce[d];
                    }
                    rDampRHS[FirstRow + TDim] += Coef * TauOne * L; // Grad(q) * TauOne * (Density * BodyForce)

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

        /// Assemble the contribution from an integration point to the element's residual.
        /** Note that the dynamic term is not included in the momentum equation.
         *  If OSS_SWITCH = 1, we don't take into account the 'dynamic' stabilization
         *  terms, as it they belong to the finite element space.
         */
        void AddProjectionResidualContribution(const array_1d< double, 3 > & rAdvVel,
                                               const double Density,
                                               array_1d< double, 3 > & rElementalMomRes,
                                               double& rElementalMassRes,
                                               const array_1d< double, TNumNodes >& rShapeFunc,
                                               const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim >& rShapeDeriv,
                                               const double Weight)
        {
            // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
            array_1d<double, TNumNodes> AGradN;
            this->GetConvectionOperator(AGradN, rAdvVel, rShapeDeriv); // Get a * grad(Ni)

            const double WeightedMass = Weight * Density;

            // Compute contribution to Kij * Uj, with Kij = Ni * Residual(Nj); Uj = (v,p)Node_j (column vector)
            for (unsigned int i = 0; i < TNumNodes; ++i) // Iterate over element nodes
            {

                // Variable references
                const array_1d< double, 3 > & rVelocity = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
                const array_1d< double, 3 > & rBodyForce = this->GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
                const double& rPressure = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);

                // Compute this node's contribution to the residual (evaluated at inegration point)
                for (unsigned int d = 0; d < TDim; ++d)
                {
                    rElementalMomRes[d] += Weight * (Density * (rShapeFunc[i] * rBodyForce[d] - AGradN[i] * rVelocity[d]) - rShapeDeriv(i, d) * rPressure);
                    rElementalMassRes -= WeightedMass * rShapeDeriv(i, d) * rVelocity[d];
                }
            }
        }

        /// Assemble the contribution from an integration point to the element's residual.
        /**
         * ASGS version. Note that rElementalMomRes should be initialized before calling this.
         * @param rAdvVel Convection velocity (not including subscale)
         * @param Density Fluid density evaluated at integration point
         * @param rElementalMomRes Result
         * @param rShapeFunc Shape functions evaluated at integration point
         * @param rShapeDeriv Shape function derivatives evaluated at integration point
         * @param Weight Integration point weight (as a fraction of area or volume)
         */
        void ASGSMomResidual(const array_1d< double, 3 > & rAdvVel,
                             const double Density,
                             array_1d< double, 3 > & rElementalMomRes,
                             const array_1d< double, TNumNodes >& rShapeFunc,
                             const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim >& rShapeDeriv,
                             const double Weight)
        {
            // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
            array_1d<double, TNumNodes> AGradN;
            this->GetConvectionOperator(AGradN, rAdvVel, rShapeDeriv); // Get a * grad(Ni)

            // Compute contribution to Kij * Uj, with Kij = Ni * Residual(Nj); Uj = (v,p)Node_j (column vector)
            for (unsigned int i = 0; i < TNumNodes; ++i) // Iterate over element nodes
            {

                // Variable references
                const array_1d< double, 3 > & rVelocity = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
                const array_1d< double, 3 > & rAcceleration = this->GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION);
                const array_1d< double, 3 > & rBodyForce = this->GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
                const double& rPressure = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);

                // Compute this node's contribution to the residual (evaluated at inegration point)
                for (unsigned int d = 0; d < TDim; ++d)
                {
                    rElementalMomRes[d] += Weight * (Density * (rShapeFunc[i] * (rBodyForce[d] - rAcceleration[d]) - AGradN[i] * rVelocity[d]) - rShapeDeriv(i, d) * rPressure);
                }
            }
        }

        /// Assemble the contribution from an integration point to the element's residual.
        /**
         * OSS version. Note that rElementalMomRes should be initialized before calling this.
         * @param rAdvVel Convection velocity (not including subscale)
         * @param Density Fluid density evaluated at integration point
         * @param rElementalMomRes Result
         * @param rShapeFunc Shape functions evaluated at integration point
         * @param rShapeDeriv Shape function derivatives evaluated at integration point
         * @param Weight Integration point weight (as a fraction of area or volume)
         */
        void OSSMomResidual(const array_1d< double, 3 > & rAdvVel,
                            const double Density,
                            array_1d< double, 3 > & rElementalMomRes,
                            const array_1d< double, TNumNodes >& rShapeFunc,
                            const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim >& rShapeDeriv,
                            const double Weight)
        {
            // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
            array_1d<double, TNumNodes> AGradN;
            this->GetConvectionOperator(AGradN, rAdvVel, rShapeDeriv); // Get a * grad(Ni)

            // Compute contribution to Kij * Uj, with Kij = Ni * Residual(Nj); Uj = (v,p)Node_j (column vector)
            for (unsigned int i = 0; i < TNumNodes; ++i) // Iterate over element nodes
            {

                // Variable references
                const array_1d< double, 3 > & rVelocity = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
                const array_1d< double, 3 > & rBodyForce = this->GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
                const array_1d< double, 3 > & rProjection = this->GetGeometry()[i].FastGetSolutionStepValue(ADVPROJ);
                const double& rPressure = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);

                // Compute this node's contribution to the residual (evaluated at inegration point)
                for (unsigned int d = 0; d < TDim; ++d)
                {
                    rElementalMomRes[d] += Weight * (Density * (rShapeFunc[i] * rBodyForce[d] - AGradN[i] * rVelocity[d]) - rShapeDeriv(i, d) * rPressure);
                    rElementalMomRes[d] -= Weight * rShapeFunc[i] * rProjection[d];
                }
            }
        }


        /// Add the contribution from Smagorinsky model to the (kinematic) viscosity
        /**
         * @param Denstity Density evaluated at the integration point (unused)
         * @param MolecularViscosity Viscosity of the fluid, in kinematic units (m2/s)
         * @param rShapeFunc Elemental shape functions, evaluated on the integration point
         * @param rShapeDeriv Shape function derivatives, evaluated on the integration point
         * @param TotalViscosity Effective viscosity (output)
         * @param rCurrentProcessInfo ProcessInfo instance (unused)
         */
        virtual void GetEffectiveViscosity(const double Density,
                                           const double MolecularViscosity,
                                           const array_1d<double, TNumNodes>& rShapeFunc,
                                           const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim >& rShapeDeriv,
                                           double& TotalViscosity,
                                           const ProcessInfo& rCurrentProcessInfo)
        {
            const double C = this->GetValue(C_SMAGORINSKY);

            TotalViscosity = MolecularViscosity;
            if (C != 0.0 )
            {
                // The filter width in Smagorinsky is typically the element size h. We will store the square of h, as the final formula involves the squared filter width
                const double FilterWidth = this->FilterWidth(rShapeDeriv);

                const double NormS = this->SymmetricGradientNorm(rShapeDeriv);

                // Total Viscosity
                TotalViscosity += 2.0 * C * C * FilterWidth * NormS;
            }
        }

        /// Write the advective velocity evaluated at this point to an array
        /**
         * Writes the value of the advective velocity evaluated at a point inside
         * the element to an array_1d
         * @param rAdvVel: Output array
         * @param rShapeFunc: Shape functions evaluated at the point of interest
         */
        virtual void GetAdvectiveVel(array_1d< double, 3 > & rAdvVel,
                                     const array_1d< double, TNumNodes >& rShapeFunc)
        {
            // Compute the weighted value of the advective velocity in the (Gauss) Point
            rAdvVel = rShapeFunc[0] * (this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY) - this->GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY));
            for (unsigned int iNode = 1; iNode < TNumNodes; ++iNode)
                rAdvVel += rShapeFunc[iNode] * (this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY) - this->GetGeometry()[iNode].FastGetSolutionStepValue(MESH_VELOCITY));
        }

        /// Write the advective velocity evaluated at this point to an array
        /**
         * Writes the value of the advective velocity evaluated at a point inside
         * the element to an array_1d
         * @param rAdvVel: Output array
         * @param rShapeFunc: Shape functions evaluated at the point of interest
         * @param Step: The time Step
         */
        virtual void GetAdvectiveVel(array_1d< double, 3 > & rAdvVel,
                                     const array_1d< double, TNumNodes >& rShapeFunc,
                                     const std::size_t Step)
        {
            // Compute the weighted value of the advective velocity in the (Gauss) Point
            rAdvVel = rShapeFunc[0] * (this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY, Step) - this->GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY, Step));
            for (unsigned int iNode = 1; iNode < TNumNodes; ++iNode)
                rAdvVel += rShapeFunc[iNode] * (this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, Step) - this->GetGeometry()[iNode].FastGetSolutionStepValue(MESH_VELOCITY, Step));
        }

        /// Write the convective operator evaluated at this point (for each nodal funciton) to an array
        /**
         * Evaluate the convective operator for each node's shape function at an arbitrary point
         * @param rResult: Output vector
         * @param rVelocity: Velocity evaluated at the integration point
         * @param rShapeDeriv: Derivatives of shape functions evaluated at the integration point
         * @see GetAdvectiveVel provides rVelocity
         */
        void GetConvectionOperator(array_1d< double, TNumNodes >& rResult,
                                   const array_1d< double, 3 > & rVelocity,
                                   const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim >& rShapeDeriv)
        {
            // Evaluate (and weight) the a * Grad(Ni) operator in the integration point, for each node i
            for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode) // Loop over nodes
            {
                // Initialize result
                rResult[iNode] = rVelocity[0] * rShapeDeriv(iNode, 0);
                for (unsigned int d = 1; d < TDim; ++d) // loop over components
                    rResult[iNode] += rVelocity[d] * rShapeDeriv(iNode, d);
            }
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
                rResult += Weight * rShapeFunc[iNode] * this->GetGeometry()[iNode].FastGetSolutionStepValue(rVariable, Step);
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
        void EvaluateInPoint(double& rResult,
                             const Variable< double >& rVariable,
                             const array_1d< double, TNumNodes >& rShapeFunc,
                             const std::size_t Step = 0)
        {
            // Compute the weighted value of the nodal variable in the (Gauss) Point
            rResult = rShapeFunc[0] * this->GetGeometry()[0].FastGetSolutionStepValue(rVariable, Step);
            for (unsigned int iNode = 1; iNode < TNumNodes; ++iNode)
                rResult += rShapeFunc[iNode] * this->GetGeometry()[iNode].FastGetSolutionStepValue(rVariable, Step);
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
        void EvaluateInPoint(array_1d< double, 3 > & rResult,
                             const Variable< array_1d< double, 3 > >& rVariable,
                             const array_1d< double, TNumNodes >& rShapeFunc)
        {
            // Compute the weighted value of the nodal variable in the (Gauss) Point
            rResult = rShapeFunc[0] * this->GetGeometry()[0].FastGetSolutionStepValue(rVariable);
            for (unsigned int iNode = 1; iNode < TNumNodes; ++iNode)
                rResult += rShapeFunc[iNode] * this->GetGeometry()[iNode].FastGetSolutionStepValue(rVariable);
        }

        /// Return an estimate for the element size h, used to calculate the stabilization parameters
        /**
         * Estimate the element size from its area or volume, required to calculate stabilization parameters.
         * Note that its implementation is different for 2D or 3D elements.
         * @see VMS2D, VMS3D for actual implementation
         * @param Volume (in 3D) or Area (in 2D) of the element
         * @return Element size h
         */
        double ElementSize(const double);

        /// Return the filter width for the smagorinsky model (Delta squared)
        double FilterWidth();

        /// Return the filter width for the smagorinsky model (Delta squared)
        double FilterWidth(const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim >& DN_DX);

        /// Compute the norm of the symmetric gradient of velocity
        /**
         * @param rShapeDeriv derivatives of the shape functions
         * @return Norm of the symmetric gradient
         */
        double SymmetricGradientNorm(const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim >& rShapeDeriv)
        {
            const unsigned int GradientSize = (TDim*(TDim+1))/2; // Number of different terms in the symmetric gradient matrix
            array_1d<double,GradientSize> GradientVector( GradientSize, 0.0 );
            unsigned int Index;

            // Compute Symmetric Grad(u). Note that only the lower half of the matrix is calculated
            for (unsigned int k = 0; k < TNumNodes; ++k)
            {
                const array_1d< double, 3 > & rNodeVel = this->GetGeometry()[k].FastGetSolutionStepValue(VELOCITY);
                Index = 0;
                for (unsigned int i = 0; i < TDim; ++i)
                {
                    for (unsigned int j = 0; j < i; ++j) // Off-diagonal
                        GradientVector[Index++] += 0.5 * (rShapeDeriv(k, j) * rNodeVel[i] + rShapeDeriv(k, i) * rNodeVel[j]);
                    GradientVector[Index++] += rShapeDeriv(k, i) * rNodeVel[i]; // Diagonal
                }
            }

            // Norm[ Symmetric Grad(u) ] = ( 2 * Sij * Sij )^(1/2)
            Index = 0;
            double NormS(0.0);
            for (unsigned int i = 0; i < TDim; ++i)
            {
                for (unsigned int j = 0; j < i; ++j)
                {
                    NormS += 2.0 * GradientVector[Index] * GradientVector[Index]; // Using symmetry, lower half terms of the matrix are added twice
                    ++Index;
                }
                NormS += GradientVector[Index] * GradientVector[Index]; // Diagonal terms
                ++Index; // Diagonal terms
            }

            NormS = sqrt( 2.0 * NormS );
            return NormS;
        }

        /// Adds the contribution of the viscous term to the momentum equation.
        /**
         * The viscous term is written in stress-divergence (Cauchy) form.
         * @param rDampMatrix Elemental Damping matrix
         * @param rShapeDeriv Elemental shape function derivatives
         * @param Weight Effective viscosity, in dynamic units, weighted by the integration point area
         */
        virtual void AddViscousTerm(MatrixType& rDampMatrix,
                                    const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim >& rShapeDeriv,
                                    const double Weight);

        /// Adds the contribution of the viscous term to the momentum equation (alternate).
        /**
         * This function is an alternate implementation of VMS::AddViscousTerm.
         * This version works with ublas matrices, using the relationship between stress and
         * rate of strain given by VMS::CalculateC. It is currently unused (as VMS::AddViscousTerm
         * is a more efficient implementation of the Cauchy equation) but it is left here so derived
         * classes can use it to implement other constitutive equations.
         * @param rDampMatrix Elemental Damping matrix
         * @param rShapeDeriv Elemental shape function derivatives
         * @param Weight Effective viscosity, in dynamic units, weighted by the integration point area
         */
        void AddBTransCB(MatrixType& rDampMatrix,
                         const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim >& rShapeDeriv,
                         const double Weight)
        {
            boost::numeric::ublas::bounded_matrix<double, (TDim * TNumNodes)/2, TDim*TNumNodes > B;
            boost::numeric::ublas::bounded_matrix<double, (TDim * TNumNodes)/2, (TDim*TNumNodes)/2 > C;
            this->CalculateB(B,rShapeDeriv);
            this->CalculateC(C,Weight);

            const unsigned int BlockSize = TDim + 1;
            const unsigned int StrainSize = (TDim*TNumNodes)/2;

            vector<unsigned int> aux(TDim*TNumNodes);
            for(unsigned int i=0; i<TNumNodes; i++)
            {
                int base_index = TDim*i;
                int aux_index = BlockSize*i;
                for(unsigned int j=0; j<TDim; j++)
                {
                    aux[base_index+j] = aux_index+j;
                }
            }

            for(unsigned int k=0; k< StrainSize; k++)
            {
                for(unsigned int l=0; l< StrainSize; l++)
                {
                    const double Ckl = C(k,l);
                    for (unsigned int i = 0; i < TDim*TNumNodes; ++i) // iterate over v components (vx,vy[,vz])
                    {
                        const double Bki=B(k,i);
                        for (unsigned int j = 0; j < TDim*TNumNodes; ++j) // iterate over u components (ux,uy[,uz])
                        {
                            rDampMatrix(aux[i],aux[j]) += Bki*Ckl*B(l,j);
                        }

                    }

                }
            }
        }

        /// Calculate the strain rate matrix
        /**
         * Unused, left to support derived classes. @see VMS::AddBTransCB
         * @param rB Strain rate matrix
         * @param rShapeDeriv Nodal shape funcion derivatives
         */
        void CalculateB( boost::numeric::ublas::bounded_matrix<double, (TDim * TNumNodes) / 2, TDim * TNumNodes >& rB,
                         const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim >& rShapeDeriv);

        /// Calculate a matrix that provides the stress given the strain rate
        /**
         * Unused, left to support derived classes. @see VMS::AddBTransCB.
         * Note that only non-zero terms are written, so the output matrix should be
         * initialized before calling this.
         * @param rC Matrix representation of the stress tensor (output)
         * @param Viscosity Effective viscosity, in dynamic units, weighted by the integration point area
         */
        virtual void CalculateC( boost::numeric::ublas::bounded_matrix<double, (TDim * TNumNodes) / 2, (TDim * TNumNodes) / 2 >& rC,
                                 const double Viscosity);

        double ConsitentMassCoef(const double Area);

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
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
        }

        virtual void load(Serializer& rSerializer)
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
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
        VMS & operator=(VMS const& rOther);

        /// Copy constructor.
        VMS(VMS const& rOther);

        ///@}

    }; // Class VMS

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
                                     VMS<TDim, TNumNodes>& rThis)
    {
        return rIStream;
    }

    /// output stream function
    template< unsigned int TDim,
              unsigned int TNumNodes >
    inline std::ostream& operator <<(std::ostream& rOStream,
                                     const VMS<TDim, TNumNodes>& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}

    ///@} // Fluid Dynamics Application group

} // namespace Kratos.

#endif // KRATOS_VMS_H_INCLUDED  defined
