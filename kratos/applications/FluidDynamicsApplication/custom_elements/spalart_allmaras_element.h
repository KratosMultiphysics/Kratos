/*
==============================================================================
KratosFluidDynamicsApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

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
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:32 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_SPALART_ALLMARAS_ELEMENT_H_INCLUDED )
#define  KRATOS_SPALART_ALLMARAS_ELEMENT_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/serializer.h"
#include "utilities/geometry_utilities.h"
#include "utilities/math_utils.h"
#include "fluid_dynamics_application_variables.h"

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

    /// An element that implements the spalart Allmaras turubulence model.
    /**
     * @todo add detailed documentation
     */
    template< unsigned int TDim, unsigned int TNumNodes = TDim + 1>
    class SpalartAllmaras : public Element
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Counted pointer of SpalartAllmaras
        KRATOS_CLASS_POINTER_DEFINITION(SpalartAllmaras);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        SpalartAllmaras(IndexType NewId = 0)
        : Element(NewId)
        {
            //DO NOT ADD DOFS HERE!!!
        }

        SpalartAllmaras(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
        {
            //DO NOT ADD DOFS HERE!!!
        }

        SpalartAllmaras(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
        {}

        /// Destructor.
        virtual ~SpalartAllmaras()
        {}


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties)  const
        {
            return Element::Pointer(new SpalartAllmaras<TDim,TNumNodes>(NewId, GetGeometry().Create(ThisNodes), pProperties));
        }


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
            if(VISCOSITY.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"VISCOSITY Key is 0. Check if the application was correctly registered.","")
            if(MOLECULAR_VISCOSITY.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"MOLECULAR_VISCOSITY Key is 0. Check if the application was correctly registered.","")
            if(TURBULENT_VISCOSITY.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"TURBULENT_VISCOSITY Key is 0. Check if the application was correctly registered.","")
            if(TEMP_CONV_PROJ.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"TEMP_CONV_PROJ Key is 0. Check if the application was correctly registered.","")

            // Checks on nodes

            // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
            for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
            {
                if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
                    KRATOS_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
                if(this->GetGeometry()[i].SolutionStepsDataHas(MESH_VELOCITY) == false)
                    KRATOS_ERROR(std::invalid_argument,"missing MESH_VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
                if(this->GetGeometry()[i].SolutionStepsDataHas(VISCOSITY) == false)
                    KRATOS_ERROR(std::invalid_argument,"missing VISCOSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
                if(this->GetGeometry()[i].SolutionStepsDataHas(MOLECULAR_VISCOSITY) == false)
                    KRATOS_ERROR(std::invalid_argument,"missing MOLECULAR_VISCOSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
                if(this->GetGeometry()[i].SolutionStepsDataHas(TURBULENT_VISCOSITY) == false)
                    KRATOS_ERROR(std::invalid_argument,"missing TURBULENT_VISCOSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
                if(this->GetGeometry()[i].SolutionStepsDataHas(TEMP_CONV_PROJ) == false)
                    KRATOS_ERROR(std::invalid_argument,"missing TEMP_CONV_PROJ variable on solution step data for node ",this->GetGeometry()[i].Id());


                if(this->GetGeometry()[i].HasDofFor(TURBULENT_VISCOSITY) == false)
                    KRATOS_ERROR(std::invalid_argument,"missing TURBULENT_VISCOSITY degree of freedom on node ",this->GetGeometry()[i].Id());
            }

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

        /// Evaluate the elemental contribution to the problem for turbulent viscosity.
        /**
         *
         * @param rLeftHandSideMatrix Elemental left hand side matrix
         * @param rRightHandSideVector Elemental right hand side vector
         * @param rCurrentProcessInfo Reference to the ProcessInfo from the ModelPart containg the element
         */
        void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
        {
            KRATOS_TRY;

            const double sigma = 2.0 / 3.0;

            const double lumping_factor = 1.00 / double(TNumNodes);

            if (rLeftHandSideMatrix.size1() != TNumNodes)
                rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);

            if (rRightHandSideVector.size() != TNumNodes)
                rRightHandSideVector.resize(TNumNodes, false);


            boost::numeric::ublas::bounded_matrix<double,TNumNodes,TNumNodes> MassFactors = lumping_factor * IdentityMatrix(TNumNodes,TNumNodes);
            boost::numeric::ublas::bounded_matrix<double,TNumNodes,TDim> DN_DX;
            array_1d<double,TNumNodes> N;
            array_1d<double,TDim> vel_gauss;
            array_1d<double,TNumNodes> temp_vec_np;
            array_1d<double,TNumNodes> u_DN;
            array_1d<double,TDim> grad_g;
            boost::numeric::ublas::bounded_matrix<double,TDim,TDim> First;
            boost::numeric::ublas::bounded_matrix<double,TDim,TDim> Second;
            boost::numeric::ublas::bounded_matrix<double,TDim,TNumNodes> Third;


            //getting data for the given geometry
            double Area;
            GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);

            //calculating viscosity
            double molecular_viscosity = GetGeometry()[0].FastGetSolutionStepValue(MOLECULAR_VISCOSITY);
            double turbulent_viscosity = GetGeometry()[0].FastGetSolutionStepValue(TURBULENT_VISCOSITY);
            double proj = GetGeometry()[0].FastGetSolutionStepValue(TEMP_CONV_PROJ);
            const array_1d<double, 3 > & v = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3 > & w = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);


            for (unsigned int j = 0; j < TDim; j++)
                vel_gauss[j] = v[j] - w[j];

            for (unsigned int i = 1; i < TNumNodes; i++)
            {
                molecular_viscosity += GetGeometry()[i].FastGetSolutionStepValue(MOLECULAR_VISCOSITY);
                turbulent_viscosity += GetGeometry()[i].FastGetSolutionStepValue(TURBULENT_VISCOSITY);
                proj += GetGeometry()[i].FastGetSolutionStepValue(TEMP_CONV_PROJ);

                const array_1d<double, 3 > & v = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
                const array_1d<double, 3 > & w = GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY);
                for (unsigned int j = 0; j < TDim; j++)
                    vel_gauss[j] += v[j] - w[j];
            }

            molecular_viscosity *= lumping_factor;
            turbulent_viscosity *= lumping_factor;
            const double conductivity = (molecular_viscosity + turbulent_viscosity) / sigma;
            proj *= lumping_factor;
            vel_gauss *= lumping_factor;

            const double c1 = 4.00;
            const double c2 = 2.00;
            const double h = this->ElemSize(Area);
            const double norm_u = norm_2(vel_gauss);
            double tau1 = (h * h) / (c1 * conductivity + c2 * norm_u * h);

            const double node_turb_visc = GetGeometry()[0].FastGetSolutionStepValue(TURBULENT_VISCOSITY);
            for (unsigned int d = 0; d < TDim; ++d)
                grad_g[d] = DN_DX(0,d) * node_turb_visc;

            for (unsigned int i = 1; i < TNumNodes; ++i)
            {
                const double node_turb_visc = GetGeometry()[i].FastGetSolutionStepValue(TURBULENT_VISCOSITY);
                for (unsigned int d = 0; d < TDim; ++d)
                {
                    grad_g[d] += DN_DX(i,d) * node_turb_visc;
                }
            }

            bool UseRotationCorrection = false;
            if (rCurrentProcessInfo[TAUONE] == 1.0) UseRotationCorrection = true;

            const double source_term = this->CalculateSourceTerm(DN_DX, N, molecular_viscosity, turbulent_viscosity, UseRotationCorrection); ///@todo: decide on best option and remove last parameter

            double res = (inner_prod(vel_gauss, grad_g));
	    res -= proj;
            double norm_grad = norm_2(grad_g);
            double k_aux = fabs(res) / (norm_grad + 0.00001);
	    k_aux = 0.0;

            noalias(First) = outer_prod(vel_gauss, trans(vel_gauss));
            First /= ((norm_u + 0.0000000001)*(norm_u + 0.0000000001));
            noalias(Second) = IdentityMatrix(TDim,TDim) - First;
            noalias(Third) = prod(Second, trans(DN_DX));

            //getting the BDF2 coefficients (not fixed to allow variable time step)
            //the coefficients INCLUDE the time step
            const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

            //CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
            noalias(u_DN) = prod(DN_DX, vel_gauss);
            noalias(rLeftHandSideMatrix) = outer_prod(N, u_DN);

            //CONVECTION STABILIZING CONTRIBUTION (Suu)
            noalias(rLeftHandSideMatrix) += tau1 * outer_prod(u_DN, u_DN);

            //VISCOUS CONTRIBUTION TO THE STIFFNESS MATRIX
            noalias(rLeftHandSideMatrix) += (conductivity * prod(DN_DX, trans(DN_DX)) + k_aux * h * prod(DN_DX, Third));

            //INERTIA CONTRIBUTION
            noalias(rLeftHandSideMatrix) += BDFcoeffs[0] * MassFactors;

            // RHS = Fext
            noalias(rRightHandSideVector) = source_term*N;

            //RHS += Suy * proj[component]
            noalias(rRightHandSideVector) += (tau1 * proj) * u_DN;

            //adding the inertia ter
            // RHS += M*vhistory
            //calculating the historical velocity
            for (unsigned int iii = 0; iii < TNumNodes; iii++)
                temp_vec_np[iii] = BDFcoeffs[1] * GetGeometry()[iii].FastGetSolutionStepValue(TURBULENT_VISCOSITY, 1);
            for (unsigned int step = 2; step < BDFcoeffs.size(); step++)
            {
                for (unsigned int iii = 0; iii < TNumNodes; iii++)
                    temp_vec_np[iii] += BDFcoeffs[step] * GetGeometry()[iii].FastGetSolutionStepValue(TURBULENT_VISCOSITY, step);
            }
            noalias(rRightHandSideVector) -= prod(MassFactors, temp_vec_np);

            //subtracting the dirichlet term
            // RHS -= LHS*temperatures
            for (unsigned int iii = 0; iii < TNumNodes; iii++)
                temp_vec_np[iii] = GetGeometry()[iii].FastGetSolutionStepValue(TURBULENT_VISCOSITY);

            noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, temp_vec_np);

            rRightHandSideVector *= Area;
            rLeftHandSideMatrix *= Area;

            KRATOS_CATCH("");
        }

        void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
        {
            KRATOS_ERROR(std::logic_error, "SplartAllmaras::CalculateRightHandSide method not implemented", "");
        }

        //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

        ///Provides the global indices for each one of this element's local rows. @see NoNewtonianASGS2D
        void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
        {
            if (rResult.size() != TNumNodes)
                rResult.resize(TNumNodes, false);

            for (unsigned int i = 0; i < TNumNodes; i++)
                rResult[i] = GetGeometry()[i].GetDof(TURBULENT_VISCOSITY).EquationId();
        }


        ///Returns a list of the element's Dofs. @see NoNewtonianASGS2D
        void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
        {
            if (ElementalDofList.size() != TNumNodes)
                ElementalDofList.resize(TNumNodes);

            for (unsigned int i = 0; i < TNumNodes; i++)
                ElementalDofList[i] = GetGeometry()[i].pGetDof(TURBULENT_VISCOSITY);
        }

        /// Calculates the projection term for stabilization
        void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
        {
            KRATOS_TRY
            int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];

            boost::numeric::ublas::bounded_matrix<double,TNumNodes,TDim> DN_DX;
            array_1d<double,TNumNodes> N;
            array_1d<double,TDim> vel_gauss;
            array_1d<double,TNumNodes> temp_vec_np;
            array_1d<double,TNumNodes> u_DN;

            //getting data for the given geometry
            double Area;
            GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);


            if (FractionalStepNumber == 2) //calculation of convective projection
            {
                const double lumping_factor = 1.00 / double(TNumNodes);

                //calculating viscosity
                temp_vec_np[0] = GetGeometry()[0].FastGetSolutionStepValue(TURBULENT_VISCOSITY);
                const array_1d<double, 3 > & v = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
                const array_1d<double, 3 > & w = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);

                for (unsigned int j = 0; j < TDim; j++)
                    vel_gauss[j] = v[j] - w[j];

                for (unsigned int i = 1; i < TNumNodes; i++)
                {
                    temp_vec_np[i] = GetGeometry()[i].FastGetSolutionStepValue(TURBULENT_VISCOSITY);
                    const array_1d<double, 3 > & v = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
                    const array_1d<double, 3 > & w = GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY);
                    for (unsigned int j = 0; j < TDim; j++)
                        vel_gauss[j] += v[j] - w[j];
                }

                vel_gauss *= lumping_factor;

                //calculating convective auxiliary vector
                noalias(u_DN) = prod(DN_DX, vel_gauss);
                double temp_conv = inner_prod(u_DN, temp_vec_np);
                temp_conv *= Area;

                for (unsigned int i = 0; i < TNumNodes; i++)
                {
                    GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += lumping_factor*Area;
                    GetGeometry()[i].FastGetSolutionStepValue(TEMP_CONV_PROJ) += lumping_factor*temp_conv;
                }
            }
            KRATOS_CATCH("");
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

        virtual std::string Info() const {
            std::stringstream buffer;
            buffer << "Spalart-Allmaras Element #" << Id();
            return buffer.str();
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "SpalartAllmaras" << TDim << "D";
        }

        /// Print object's data.
        //      virtual void PrintData(std::ostream& rOStream) const;


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

        double CalculateSourceTerm(const boost::numeric::ublas::bounded_matrix<double,TNumNodes,TDim> & DN_DX,
                                   const array_1d<double,TNumNodes> & N,
                                   const double molecular_viscosity,
                                   const double turbulent_viscosity,
                                   bool UseRotationCorrection)
        {
            const double sigma = 2.0 / 3.0;
            const double cb1 = 0.1355;
            const double cb2 = 0.622;
            const double kappa = 0.41; // Von Karman constant
            const double cw1 = cb1 / (kappa * kappa) + (1.0 + cb2) / sigma;
            const double cw2 = 0.3;
            const double cw3 = 2.0;
            const double cv1 = 7.1;
            //const double ct1 = 1.0; // For ft1 (trip term, not implemented)
            //const double ct2 = 2.0;
            const double ct3 = 1.2; // Original value is 1.1, correction by Spalart
            const double ct4 = 0.5; // Original value is 2.0, correction by Spalart

            const double xi = turbulent_viscosity/molecular_viscosity;
            const double fv1 = xi * xi * xi / (xi * xi * xi + cv1 * cv1 * cv1);
            const double fv2 = 1.0 - xi / ( 1.0 + xi * fv1);

            double S;
            if (UseRotationCorrection)
                S = this->CorrectedRotationRate(DN_DX);
            else
                S = this->AntimetricGradientNorm(DN_DX);

            Element::GeometryType& rGeom = this->GetGeometry();

            double distance = N[0] * rGeom[0].FastGetSolutionStepValue(DISTANCE);
            for (unsigned int i = 1; i < TNumNodes; ++i)
                distance += N[i] * rGeom[i].FastGetSolutionStepValue(DISTANCE);

            double S_hat = S + fv2 * turbulent_viscosity / (distance * distance * kappa * kappa);
            if(S_hat < 1e-9)
                S_hat = 1e-9;

            const double ft2 = ct3 * exp(-ct4*xi*xi);

            double r = turbulent_viscosity / ( S_hat * kappa * kappa * distance * distance);
            if(r>10.0) r=10.0;
//            else if(r<0.0) r = 10.0; // This happens if S_hat < 0.0, which is not allowed, take S_hat = 0; r = min(inf,10)

            const double g = r + cw2 * (pow(r,6) - r);
            const double fw = g * pow( (1.0+pow(cw3,6)) / ( pow(g,6) + pow(cw3,6) ) , 1.0/6.0);

            array_1d<double,TDim> grad_nu(TDim,0.0);
            for (unsigned int i = 0; i < TNumNodes; ++i)
            {
                const double node_turb_visc = rGeom[i].FastGetSolutionStepValue(TURBULENT_VISCOSITY);
                for (unsigned int d = 0; d < TDim; ++d)
                    grad_nu[d] += DN_DX(i,d) * node_turb_visc;
            }

            double norm2_grad_nu = grad_nu[0]*grad_nu[0];
            for (unsigned int d = 1; d < TDim; ++d)
                norm2_grad_nu += grad_nu[d]*grad_nu[d];

            double source_term = cb1 * (1.0 - ft2) * S_hat * turbulent_viscosity;
            source_term += cb2 * norm2_grad_nu / sigma;

            source_term -= (cw1 * fw - ft2 * cb1/(kappa*kappa) ) * pow(turbulent_viscosity/distance,2);

            return source_term;
        }

        double AntimetricGradientNorm(const boost::numeric::ublas::bounded_matrix<double,TNumNodes,TDim> & rShapeDeriv)
        {
            boost::numeric::ublas::bounded_matrix<double,2,2> grad = ZeroMatrix(TDim,TDim);

            // Compute Antimetric Grad(u)
            for (unsigned int k = 0; k < TNumNodes; ++k)
            {
                const array_1d<double, 3> & rNodeVel = this->GetGeometry()[k].FastGetSolutionStepValue(VELOCITY);
                for (unsigned int i = 0; i < TDim; ++i)
                    for (unsigned int j = 0; j < TDim; ++j)
                        grad(i,j) += 0.5 * (rShapeDeriv(k, j) * rNodeVel[i] - rShapeDeriv(k, i) * rNodeVel[j]);
            }

            double NormS(0.0);
            for (unsigned int i = 0; i < TDim; ++i)
                for (unsigned int j = 0; j < TDim; ++j)
                    NormS += grad(i,j) * grad(i,j);

            NormS = sqrt( 2.0 * NormS );
            return NormS;
        }

        double CorrectedRotationRate(const boost::numeric::ublas::bounded_matrix<double,TNumNodes,TDim> & rShapeDeriv)
        {
            const double Crot = 2.0; // As proposed in reference, this value could be adjusted if necessary

            boost::numeric::ublas::bounded_matrix<double,2,2> sGrad = ZeroMatrix(TDim,TDim);
            boost::numeric::ublas::bounded_matrix<double,2,2> aGrad = ZeroMatrix(TDim,TDim);

            // Compute Antimetric Grad(u)
            for (unsigned int k = 0; k < TNumNodes; ++k)
            {
                const array_1d<double, 3> & rNodeVel = this->GetGeometry()[k].FastGetSolutionStepValue(VELOCITY);
                for (unsigned int i = 0; i < TDim; ++i)
                    for (unsigned int j = 0; j < TDim; ++j)
                    {
                        sGrad(i,j) += 0.5 * (rShapeDeriv(k, j) * rNodeVel[i] + rShapeDeriv(k, i) * rNodeVel[j]);
                        aGrad(i,j) += 0.5 * (rShapeDeriv(k, j) * rNodeVel[i] - rShapeDeriv(k, i) * rNodeVel[j]);
                    }
            }

            double NormS(0.0);
            double NormA(0.0);
            for (unsigned int i = 0; i < TDim; ++i)
                for (unsigned int j = 0; j < TDim; ++j)
                {
                    NormS += sGrad(i,j) * sGrad(i,j);
                    NormA += aGrad(i,j) * aGrad(i,j);
                }

            // S = NormA + 2.0 * min(0,NormS - NormA)
            if (NormA < NormS)
            {
                NormS = sqrt( 2.0 * NormS );
                NormA = sqrt( 2.0 * NormA );

                NormS += Crot * (NormS - NormA);
                return NormS;
            }
            else
                return sqrt( 2.0 * NormS );
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

        virtual void save(Serializer& rSerializer) const {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
        }

        virtual void load(Serializer& rSerializer) {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
        }

        ///@}
        ///@name Private Operators
        ///@{

        ///@}
        ///@name Private Operations
        ///@{

        double ElemSize(const double Area);

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
        SpalartAllmaras& operator=(const SpalartAllmaras& rOther);

        /// Copy constructor.
        SpalartAllmaras(const SpalartAllmaras& rOther);


        ///@}

    }; // Class SpalartAllmaras2D


    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    // input stream function
    template< unsigned int TDim, unsigned int TNumNodes >
    inline std::istream& operator >> (std::istream& rIStream,
                                      SpalartAllmaras<TDim,TNumNodes>& rThis)
    {
        return rIStream;
    }

    /// output stream function
    template< unsigned int TDim, unsigned int TNumNodes >
    inline std::ostream& operator << (std::ostream& rOStream,
                                      const SpalartAllmaras<TDim,TNumNodes>& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }

    ///@}

} // namespace Kratos.

#endif // KRATOS_SPALART_ALLMARAS_ELEMENT_H_INCLUDED  defined


