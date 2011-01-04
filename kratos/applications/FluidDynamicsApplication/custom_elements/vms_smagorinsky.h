/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis

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
//   Date:                $Date: 2010-12-09 17:01:00 $
//   Revision:            $Revision: 0.1 $
//
//

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "custom_elements/vms_base.h"
#include "fluid_dynamics_application_variables.h"

#ifndef KRATOS_VMS_SMAGORINSKY_H_DEFINED
#define	KRATOS_VMS_SMAGORINSKY_H_DEFINED

namespace Kratos
{
    /// A stabilized finite element formulation for Incompressible Navier-Stokes with Smagorinsky eddy viscosity
    /**
     */
    template< unsigned int TDim,
              unsigned int TNumNodes = TDim +1
            >
    class VMSSmagorinsky : public VMSBase< TDim, TNumNodes, TDim+1, TNumNodes*(TDim+1)>
    {
    public:
        /// Pointer definition of VMSSmagorinsky
        typedef VMSSmagorinsky<TDim, TNumNodes> VMSSmagorinskyType;
        KRATOS_CLASS_POINTER_DEFINITION( VMSSmagorinskyType );

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

        ///Type definition for integration methods
        typedef GeometryData::IntegrationMethod IntegrationMethod;

        typedef VMSBase< TDim, TNumNodes, TDim+1, TNumNodes*(TDim+1)> BaseElement;

        ///@}
        ///@name Life Cycle
        ///@{

        ///Constructors.

        /**
         * Default constuctor.
         * @param NewId Index number of the new element (optional)
         */
        VMSSmagorinsky(IndexType NewId = 0) :
            VMSBase<TDim,TNumNodes>(NewId)
        {}

        /**
         * Constructor using an array of nodes.
         * @param NewId Index of the new element
         * @param ThisNodes An array containing the nodes of the new element
         */
        VMSSmagorinsky(IndexType NewId, const NodesArrayType& ThisNodes) :
            VMSBase<TDim,TNumNodes>(NewId, ThisNodes)
        {}

        /**
         * Constructor using a geometry object.
         * @param NewId Index of the new element
         * @param pGeometry Pointer to a geometry object
         */
        VMSSmagorinsky(IndexType NewId, GeometryType::Pointer pGeometry) :
            VMSBase<TDim,TNumNodes>(NewId, pGeometry)
        {}

        /**
         * Constuctor using geometry and properties.
         * @param NewId Index of the new element
         * @param pGeometry Pointer to a geometry object
         * @param pProperties Pointer to the element's properties
         */
        VMSSmagorinsky(IndexType NewId, GeometryType::Pointer pGeometry,
                 PropertiesType::Pointer pProperties) :
            VMSBase<TDim,TNumNodes>(NewId, pGeometry, pProperties)
        {}

        /// Destructor.
        virtual ~VMSSmagorinsky()
        {}


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        /// Create a new element of this type
        /**
         * Returns a pointer to a new VMS2D element, created using given input
         * @param NewId: the ID of the new element
         * @param ThisNodes: the nodes of the new element
         * @param pProperties: the properties assigned to the new element
         * @return a Pointer to the new element
         */
        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                PropertiesType::Pointer pProperties) const
        {
            KRATOS_TRY
            KRATOS_ERROR(std::logic_error,"VMSSmagorinsky::Create failed. Cannot create an instance of an abstract class. Please use VMS2DSmagorinsky or VMS3DSmagorinsky elements instead","")
            KRATOS_CATCH("");
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
            const unsigned int LocalSize = TNumNodes * (TDim + 1);

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
            this->GetPointContribution(Density, DENSITY, N);
            this->GetPointContribution(KinViscosity, VISCOSITY, N);

            // Get Advective velocity
            array_1d<double, 3> AdvVel;
            this->GetAdvectiveVel(AdvVel, N);

            // Calculate stabilization parameters
            double TauOne, TauTwo;
            this->CalculateTau(TauOne, TauTwo, AdvVel, Area, KinViscosity, rCurrentProcessInfo);

            // Get Smagorinsky coefficient
            const double Cs = rCurrentProcessInfo[C_SMAGORINSKY];

            this->AddIntegrationPointVelocityContribution(rDampMatrix, rRightHandSideVector, Density, KinViscosity, TauOne, TauTwo, Cs, N, DN_DX, Area);

            // Now calculate an additional contribution to the residual: r -= rDampMatrix * (u,p)
            VectorType U = ZeroVector(LocalSize);
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
                boost::numeric::ublas::bounded_matrix<double,LocalSize,LocalSize> MassStabilization = ZeroMatrix(LocalSize, LocalSize);

                this->AddMassStabTerms(MassStabilization, Density, TauOne, N, DN_DX, Area);

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
            VMSBase<TDim,TNumNodes>::Calculate(rVariable,rOutput,rCurrentProcessInfo);
//            // Get the element's geometric parameters
//            double Area;
//            array_1d<double, TNumNodes> N;
//            boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
//            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);
//
//            // Calculate this element's fluid properties
//            double Density, Viscosity;
//            this->GetPointContribution(Density, DENSITY, N);
//            this->GetPointContribution(Viscosity, VISCOSITY, N);
//
//            // Get Advective velocity
//            array_1d<double, 3> AdvVel;
//            this->GetAdvectiveVel(AdvVel, N);
//
//            // Calculate stabilization parameters
//            double TauOne, TauTwo;
//            this->CalculateTau(TauOne, TauTwo, AdvVel, Area, Viscosity, rCurrentProcessInfo);
//
//            // Output containers
//            array_1d< double,3 > ElementalMomRes(3,0.0);
//            double ElementalMassRes(0);
//
//            this->AddProjectionResidualContribution(AdvVel, Density, TauOne, TauTwo, ElementalMomRes,ElementalMassRes, rCurrentProcessInfo, N, DN_DX, Area);
////            if (rVariable == COARSE_VELOCITY)
////                this->AddCoarseResidualContribution(AdvVel, Density, TauOne, TauTwo, ElementalMomRes,ElementalMassRes, rCurrentProcessInfo, N, DN_DX, Area);
////            else
////                this->AddResidualContribution(AdvVel, Density, TauOne, TauTwo, ElementalMomRes,ElementalMassRes, rCurrentProcessInfo, N, DN_DX, Area);
//
//            // Checking on variable to ensure that we only write on nodal variables when callling from scheme (residual projection)
//            if (rVariable == ADVPROJ) //if( rCurrentProcessInfo[OSS_SWITCH] == 1)
//            {
//                // Carefully write results to nodal variables, to avoid parallelism problems
//                for (unsigned int i = 0; i < TNumNodes; ++i)
//                {
//                    ///@TODO: Test using atomic
//                    this->GetGeometry()[i].SetLock(); // So it is safe to write in the node in OpenMP
//                    array_1d< double, 3 > & rAdvProj = this->GetGeometry()[i].FastGetSolutionStepValue(ADVPROJ);
//                    for (unsigned int d = 0; d < TDim; ++d)
//                        rAdvProj[d] += N[i] * ElementalMomRes[d];
//
//                    this->GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ) += N[i] * ElementalMassRes;
//                    this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += Area * N[i];
//                    this->GetGeometry()[i].UnSetLock(); // Free the node for other threads
//                }
//            }
//
//            /// Return output
//            rOutput = ElementalMomRes;
        }

        virtual void Calculate( const Variable<Vector>& rVariable,
                                Vector& rOutput,
                                const ProcessInfo& rCurrentProcessInfo)
        {
            rOutput.resize(2,false);
            rOutput = ZeroVector(2);

            boost::numeric::ublas::bounded_matrix<double, TDim, TDim> GradUc = ZeroMatrix(TDim,TDim); // Grad(u) in coarse mesh
            boost::numeric::ublas::bounded_matrix<double, TDim, TDim> GradUf = ZeroMatrix(TDim,TDim); // Grad(u) in fine mesh

            // Get the element's geometric parameters
            double Area;
            array_1d<double, TNumNodes> N;
            boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            // Calculate residual on both meshes
            array_1d< double,3 > FineRes(3,0.0), CoarseRes(3,0.0);
            this->Calculate(VELOCITY,FineRes,rCurrentProcessInfo);
            this->Calculate(COARSE_VELOCITY,CoarseRes,rCurrentProcessInfo);

            // Compute Grad(u) on both meshes and u on coarse mesh
            array_1d< double,3 > ElemCoarseVel(3,0.0);
            for (unsigned int k = 0; k < TNumNodes; ++k)
            {
                const array_1d< double,3 >& rFineVel = this->GetGeometry()[k].FastGetSolutionStepValue(VELOCITY);
                const array_1d< double,3 >& rCoarseVel = this->GetGeometry()[k].GetValue(COARSE_VELOCITY);
                for (unsigned int i = 0; i < TDim; ++i)
                {
                    ElemCoarseVel[i] += rCoarseVel[i] * N[k];
                    for (unsigned int j = 0; j < TDim; ++j)
                    {
                        GradUc(i,j) += DN_DX(k,j) * rCoarseVel[i];
                        GradUf(i,j) += DN_DX(k,j) * rFineVel[i];
                    }
                }
            }

            // Compute Norm[ Grad(u) ] on both meshes
            double NormSc(0.0),NormSf(0.0);
            for (unsigned int i = 0; i < TDim; ++i)
                for (unsigned int j = 0; j < TDim; ++j)
                {
                    NormSc += 2.0 * pow(GradUc(i,j),2);
                    NormSf += 2.0 * pow(GradUf(i,j),2);
                }

            NormSc = sqrt( NormSc );
            NormSf = sqrt( NormSf );

            // Calculate square of Filter Width for both meshes (h,H)
            double sqh,sqH;
            if( this->GetGeometry().WorkingSpaceDimension() == 2)
            {
                sqh = Area;
                sqH = 16 * Area; // Area_coarse = 4 * Area_fine, by construction
            }
            else // Here the variable Area contains the elemental volume
            {
                sqh = pow(Area,1.0/3.0);
                sqh = pow(sqh,2.0);
                sqH = pow(8.0 * Area,1.0/3.0); // Vol_coarse = 8 * Vol_fine
                sqH = pow(sqH,2.0);
            }

            double SmaC(0.0),SmaF(0.0);

            // Calculate local contributions to Variational Germano identity
            for(unsigned int i = 0; i < TDim; ++i)
            {
                rOutput[0] += ElemCoarseVel[i] * (CoarseRes[i] - FineRes[i]);
                for(unsigned int j = 0; j < TDim; ++j)
                {
                    SmaC += 2.0 * GradUc(i,j) * GradUc(i,j);
                    SmaF += 2.0 * fabs( GradUc(i,j) * GradUf(i,j) );
                }
            }
            SmaC = sqH * NormSc * sqrt( SmaC );
            SmaF = sqh * NormSf * sqrt( SmaF );

            rOutput[1] = SmaC - SmaF;

//            //loop on coarse mesh elements
//            this->CalculateRHS();
//            double output = 0.0;

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
            buffer << "VMSSmagorinsky #" << this->Id();
            return buffer.str();
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "VMSSmagorinsky #" << this->Id();
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

        /// Add a the contribution from a single integration point to the velocity contribution
        void AddIntegrationPointVelocityContribution(MatrixType& rDampMatrix,
                                                     VectorType& rDampRHS,
                                                     const double Density,
                                                     const double KinViscosity,
                                                     const double TauOne,
                                                     const double TauTwo,
                                                     const double Cs,
                                                     const array_1d< double, TNumNodes >& rShapeFunc,
                                                     const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim >& rShapeDeriv,
                                                     const double Weight)
        {
            // Dofs per node
            const unsigned int BlockSize = TDim + 1;

            // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
            array_1d<double, TNumNodes> AGradN;
            this->GetConvectionOperator(AGradN, rShapeFunc, rShapeDeriv); // Get a * grad(Ni)

            // Build the local matrix and RHS
            unsigned int FirstRow(0), FirstCol(0); // position of the first term of the local matrix that corresponds to each node combination
            double K, G, PDivV, L; // Temporary results
            double Coef = Density * Weight;

            // Get total viscosity (as given by Smagorinsky)
            double Viscosity;
            this->GetEffectiveViscosity(KinViscosity,Cs,rShapeDeriv,Viscosity); // Smagorinsky viscosity (in "Kinematic" units, m^2/s)

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
                        K += Viscosity * rShapeDeriv(i, m) * rShapeDeriv(j, m); // Diffusive term: Viscosity * Grad(v) * Grad(u)
                        // Effective Viscosity, as given by Smagorinsky

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
                    FirstRow += BlockSize;
                }

                // Update reference indices
                FirstRow = 0;
                FirstCol += BlockSize;
            }
        }

//        /// Assemble the contribution from an integration point to the element's residual
//        void AddResidualContribution(const array_1d< double, 3 > & rAdvVel,
//                                     const double Density,
//                                     const double TauOne,
//                                     const double TauTwo,
//                                     array_1d< double,3 >& rElementalMomRes,
//                                     double& rElementalMassRes,
//                                     const ProcessInfo& rCurrentProcessInfo,
//                                     const array_1d< double, TNumNodes >& rShapeFunc,
//                                     const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim >& rShapeDeriv,
//                                     const double Weight)
//        {
//            // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
//            array_1d<double, TNumNodes> AGradN;
//            this->GetConvectionOperator(AGradN, rShapeFunc, rShapeDeriv); // Get a * grad(Ni)
//
//            const double WeightedMass = Weight * Density;
//
//            ///@TODO: Check that this is really safe for multiple Gauss points. Suggestion: return the elemental residual, write to nodes in caller function.
//            // Compute contribution to Kij * Uj, with Kij = Ni * Residual(Nj); Uj = (v,p)Node_j (column vector)
//            for (unsigned int i = 0; i < TNumNodes; ++i) // Iterate over element nodes
//            {
//
//                // Variable references
//                const array_1d< double, 3 > & rVelocity = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
//                const double& rPressure = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);
//
//                // Compute this node's contribution to the residual (evaluated at inegration point)
//                // Pressure contribution to the momentum residual is done separately because it is not computed from array_1d variables
//                for (unsigned int d = 0; d < TDim; ++d)
//                {
////                    ElementalMomRes[d] += Weight * ( Density * ( /*rShapeFunc[i] * (rAcceleration[d]-rBodyForce[d])*/ + AGradN[i] * rVelocity[d] ) + rShapeDeriv(i,d) * rPressure);
//                    ///@TODO: Body force or not Body force?
//                    rElementalMomRes[d] += Weight * (Density * (AGradN[i] * rVelocity[d]) + rShapeDeriv(i, d) * rPressure);
//                    rElementalMassRes += WeightedMass * rShapeDeriv(i, d) * rVelocity[d];
//                }
//            }
//        }

        /// Assemble the contribution from an integration point to the element's residual
        void AddCoarseResidualContribution(const array_1d< double, 3 > & rAdvVel,
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
                const array_1d< double, 3 > & rVelocity = this->GetGeometry()[i].GetValue(COARSE_VELOCITY);
                const double& rPressure = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);

                // Compute this node's contribution to the residual (evaluated at inegration point)
                // Pressure contribution to the momentum residual is done separately because it is not computed from array_1d variables
                for (unsigned int d = 0; d < TDim; ++d)
                {
//                    ElementalMomRes[d] += Weight * ( Density * ( /*rShapeFunc[i] * (rAcceleration[d]-rBodyForce[d])*/ + AGradN[i] * rVelocity[d] ) + rShapeDeriv(i,d) * rPressure);
                    ///@TODO: Body force or not Body force?
                    rElementalMomRes[d] += Weight * (Density * (AGradN[i] * rVelocity[d]) + rShapeDeriv(i, d) * rPressure);
                    rElementalMassRes += WeightedMass * rShapeDeriv(i, d) * rVelocity[d];
                }
            }
        }

        void GetEffectiveViscosity(const double MolecularViscosity,
                                   const double C,
                                   const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim >& rShapeDeriv,
                                   double& TotalViscosity)
        {
            double FilterWidth;
            if(TDim == 2)
            {
                FilterWidth = GeometryUtils::CalculateVolume2D(this->GetGeometry());
            }
            else
            {
                FilterWidth = GeometryUtils::CalculateVolume3D(this->GetGeometry());
                FilterWidth = pow(FilterWidth, 2.0/3.0);
            }

            boost::numeric::ublas::bounded_matrix<double,TDim,TDim> dv_dx = ZeroMatrix(TDim,TDim);

            TotalViscosity = MolecularViscosity;

            // Compute Grad(u)
            for (unsigned int k = 0; k < TNumNodes; ++k)
            {
                const array_1d< double,3 >& rNodeVel = this->GetGeometry()[k].FastGetSolutionStepValue(VELOCITY);
                for (unsigned int i = 0; i < TDim; ++i)
                    for (unsigned int j = 0; j < TDim; ++j)
                        dv_dx(i,j) += rShapeDeriv(k,j) * rNodeVel[i];
            }

            // Norm[ Grad(u) ]
            double NormS(0.0);
            for (unsigned int i = 0; i < TDim; ++i)
                for (unsigned int j = 0; j < TDim; ++j)
                    NormS += dv_dx(i,j) * dv_dx(i,j);

            NormS = sqrt(NormS);

            // Total Viscosity
            TotalViscosity += 2.0 * C * C * FilterWidth * NormS;
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
                rAdvVel += rShapeFunc[iNode] * (this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY,Step) - this->GetGeometry()[iNode].FastGetSolutionStepValue(MESH_VELOCITY,Step) + this->GetGeometry()[iNode].FastGetSolutionStepValue(ADVPROJ,Step));
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

        /// Estimate element size
        virtual double ElementSize(const double Area) = 0;

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
        VMSSmagorinsky & operator=(VMSSmagorinsky const& rOther);

        /// Copy constructor.
        VMSSmagorinsky(VMSSmagorinsky const& rOther);

        ///@}

    }; // Class VMSSmagorinsky

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function
    template< unsigned int TDim, unsigned int TNumNodes>
    inline std::istream & operator >>(std::istream& rIStream,
                                      VMSSmagorinsky<TDim,TNumNodes>& rThis)
    {
        return rIStream;
    }

    /// output stream function
    template< unsigned int TDim, unsigned int TNumNodes>
    inline std::ostream & operator <<(std::ostream& rOStream,
                                      const VMSSmagorinsky<TDim,TNumNodes>& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}
}

#endif	/* KRATOS_VMS_SMAGORINSKY_H_DEFINED */

