/*
==============================================================================
KratosPhaseFieldApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
hgbk2008@gmail.com
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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
/* *********************************************************
 *
 *   Last Modified by:    $Author: hbui $
 *   Date:                $Date: 4 Dec 2015 $
 *   Revision:            $Revision: 1.0 $
 *
 * ***********************************************************/

// System includes


// External includes

// Project includes
#include "phase_field_kinematic_linear.h"
#include "phase_field_application.h"
#include "custom_utilities/eig/eig3.h"
#include "structural_application/custom_utilities/sd_math_utils.h"


#define CHECK_NAN


namespace Kratos
{

    extern Variable<double> FRACTURE_ENERGY;
    extern Variable<Vector> STRESSES;
    extern Variable<array_1d<double, 3> > PRESCRIBED_DELTA_DISPLACEMENT;

    PhaseFieldKinematicLinear::PhaseFieldKinematicLinear(IndexType NewId, GeometryType::Pointer pGeometry) : Element(NewId, pGeometry)
    {
        mIsInitialized = false;
    }

    PhaseFieldKinematicLinear::PhaseFieldKinematicLinear(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : Element(NewId, pGeometry, pProperties)
    {
        mIsInitialized = false;
        mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
    }

    Element::Pointer PhaseFieldKinematicLinear::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new PhaseFieldKinematicLinear(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    Element::Pointer PhaseFieldKinematicLinear::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new PhaseFieldKinematicLinear(NewId, pGeom, pProperties));
    }

    PhaseFieldKinematicLinear::~PhaseFieldKinematicLinear()
    {
    }

    void PhaseFieldKinematicLinear::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
    {
        DofsVectorType DofList;
        this->GetDofList(DofList, CurrentProcessInfo);
        if(rResult.size() != DofList.size())
            rResult.resize(DofList.size());
        for(unsigned int i = 0; i < DofList.size(); ++i)
            rResult[i] = DofList[i]->EquationId();
    }

    void PhaseFieldKinematicLinear::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        ElementalDofList.resize( 0 );

        if(CurrentProcessInfo[FRACTIONAL_STEP] == 0) // solve for displacements
        {
            for(unsigned int i = 0; i < GetGeometry().size(); ++i)
            {
                ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
                ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
                if(dim == 3)
                    ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
            }
        }
        else if(CurrentProcessInfo[FRACTIONAL_STEP] == 1) // solve for phase field
        {
            for(unsigned int i = 0; i < GetGeometry().size(); ++i)
                ElementalDofList.push_back(GetGeometry()[i].pGetDof(PHASE_FIELD));
        }
    }

    void PhaseFieldKinematicLinear::Initialize()
    {
        if(mIsInitialized == false)
        {
            unsigned int dim = GetGeometry().WorkingSpaceDimension();
            unsigned int strain_size = dim * (dim + 1) / 2;

            // extract the integration points
            const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

            // initialize the reference energy density
            mReferenceEnergyDensity.resize(integration_points.size());
            std::fill(mReferenceEnergyDensity.begin(), mReferenceEnergyDensity.end(), 0.0);

            // initialize the stresses container
            mCurrentStresses.resize(integration_points.size());
            for(unsigned int i = 0; i < integration_points.size(); ++i)
                mCurrentStresses[i].resize(strain_size);

            mIsInitialized = true;
        }
    }

    void PhaseFieldKinematicLinear::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
        if(CurrentProcessInfo[FRACTIONAL_STEP] == 1) // solve for phase field
        {
            /* compute the reference energy density */
            unsigned int number_of_nodes = GetGeometry().size();
            unsigned int dim = GetGeometry().WorkingSpaceDimension();
            unsigned int strain_size = dim * (dim + 1) / 2;
            unsigned int mat_size = dim * number_of_nodes;
            Matrix B(strain_size, mat_size);
            Vector StrainVector(strain_size);
            Matrix DN_DX(number_of_nodes, dim);
            Matrix CurrentDisp(number_of_nodes, dim);
            Matrix InvJ(dim, dim);
            double DetJ;
            std::vector<array_1d<double, 3> > n(3);
            double e[3];
            Matrix e_plus(3, 3);

            double E = GetProperties()[YOUNG_MODULUS];
            double NU = GetProperties()[POISSON_RATIO];
            double lambda = E * NU / (1.0 + NU) / (1.0 - 2.0 * NU);
            double mu = 0.5 * E / (1.0 + NU);

            #ifdef ENABLE_BEZIER_GEOMETRY
            // initialize geometry
            GetGeometry().Initialize(mThisIntegrationMethod);
            #endif

            // extract the integration points
            const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

            // extract the shape function local gradients
            const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

            // compute Jacobian
            GeometryType::JacobiansType J(integration_points.size());
            J = GetGeometry().Jacobian(J, mThisIntegrationMethod);

            // extract current displacements
            for(unsigned int node = 0; node < number_of_nodes; ++node)
                noalias(row(CurrentDisp, node)) = GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT);

            // loop through integration to compute LHS & RHS contribution
            for(unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
            {
                // compute inverse of Jacobian
                MathUtils<double>::InvertMatrix(J[PointNumber], InvJ, DetJ);

                // compute the gradient w.r.t global coordinates
                noalias(DN_DX) = prod(DN_De[PointNumber], InvJ);

                // compute the B Operator
                this->CalculateBoperator(B, DN_DX);

                // compute the strain
                this->CalculateStrain(StrainVector, B, CurrentDisp);

                double trace_strain_plus = 0.0;
                for(unsigned int i = 0; i < dim; ++i)
                    trace_strain_plus += StrainVector(i);
                if(trace_strain_plus < 0.0)
                    trace_strain_plus = 0.0;

                // perform the spectral decomposition
                if(dim == 2)
                    spectral_decomposition(StrainVector(0), StrainVector(1), 0.0, 0.5 * StrainVector(2), 0.0, 0.0, e[0], e[1], e[2], n[0], n[1], n[2]);
                else if(dim == 3)
                    spectral_decomposition(StrainVector(0), StrainVector(1), StrainVector(2), 0.5 * StrainVector(3), 0.5 * StrainVector(4), 0.5 * StrainVector(5), e[0], e[1], e[2], n[0], n[1], n[2]);

                // compute positive part of strain tensor
                noalias(e_plus) = ZeroMatrix(3, 3);
                for(unsigned int i = 0; i < 3; ++i)
                    if(e[i] > 0.0)
                    {
                        for(unsigned int j = 0; j < 3; ++j)
                            for(unsigned int k = 0; k < 3; ++k)
                                e_plus(j, k) += e[i] * n[i](j) * n[i](k);
                    }

                // compute reference energy density at the integration point
                double psi_plus = 0.0;
                psi_plus += 0.5 * lambda * pow(trace_strain_plus, 2);
                for(unsigned int j = 0; j < 3; ++j)
                    for(unsigned int k = 0; k < 3; ++k)
                        psi_plus += mu * pow(e_plus(j, k), 2);

                // update the history field
                if(psi_plus > mReferenceEnergyDensity[PointNumber])
                    mReferenceEnergyDensity[PointNumber] = psi_plus;
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            // clean the geometry internal data
            GetGeometry().Clean();
            #endif
        }
    }

    void PhaseFieldKinematicLinear::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
    {
        if(CurrentProcessInfo[FRACTIONAL_STEP] == 0) // solve for displacement
        {
            // reset all resistant forces at node
            for ( unsigned int i = 0; i < GetGeometry().size(); ++i )
            {
                GetGeometry()[i].GetSolutionStepValue( REACTION_X ) = 0.0;
                GetGeometry()[i].GetSolutionStepValue( REACTION_Y ) = 0.0;
                GetGeometry()[i].GetSolutionStepValue( REACTION_Z ) = 0.0;
            }
        }
        else if(CurrentProcessInfo[FRACTIONAL_STEP] == 1) // solve for phase field
        {
            // reset all dual phase field at node
            for ( unsigned int i = 0; i < GetGeometry().size(); ++i )
                GetGeometry()[i].GetSolutionStepValue( PHASE_FIELD_DUAL_VARIABLE ) = 0.0;
        }
    }

    void PhaseFieldKinematicLinear::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, true, true, true);
    }

    void PhaseFieldKinematicLinear::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        MatrixType temp = Matrix();

        bool need_calculate_stiffness = false;
        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
        {
            if(GetGeometry()[node].IsFixed(DISPLACEMENT_X)
            || GetGeometry()[node].IsFixed(DISPLACEMENT_Y)
            || GetGeometry()[node].IsFixed(DISPLACEMENT_Z))
            {
                need_calculate_stiffness = true;
                break;
            }
        }

        CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, need_calculate_stiffness, true, false);
    }

    void PhaseFieldKinematicLinear::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, bool CalculateStiffnessMatrixFlag, bool CalculateResidualVectorFlag, bool MaterialUpdateFlag)
    {
        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        #ifdef ENABLE_BEZIER_GEOMETRY
        // initialize the geometry
        GetGeometry().Initialize(mThisIntegrationMethod);
        #endif

        // extract the integration points
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

        // extract the shape function values
        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);

        // extract the shape function local gradients
        const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

        // compute Jacobian
        GeometryType::JacobiansType J(integration_points.size());
        J = GetGeometry().Jacobian(J, mThisIntegrationMethod);

        if(rCurrentProcessInfo[FRACTIONAL_STEP] == 0) // solve for displacement
        {
            // declare local variables
            unsigned int strain_size = dim * (dim + 1) / 2;
            unsigned int mat_size = dim * number_of_nodes;
            Matrix B(strain_size, mat_size);
            Matrix TanC(strain_size, strain_size);
            Vector StrainVector(strain_size);
            Vector StressVector(strain_size);
            Matrix DN_DX(number_of_nodes, dim);
            Matrix CurrentDisp(number_of_nodes, dim);
            Vector InternalForces(dim);
            double C;
            Matrix InvJ(dim, dim);
            double DetJ;
            double e[3];
            Matrix e_plus(3, 3);

            // material parameters
            double E = GetProperties()[YOUNG_MODULUS];
            double NU = GetProperties()[POISSON_RATIO];
            double lambda = E * NU / (1.0 + NU) / (1.0 - 2.0 * NU);
            double mu = 0.5 * E / (1.0 + NU);
            double kappa = GetProperties()[KAPPA];

            // reset the LHS
            if(CalculateStiffnessMatrixFlag == true)
            {
                if(rLeftHandSideMatrix.size1() != mat_size || rLeftHandSideMatrix.size2() != mat_size)
                    rLeftHandSideMatrix.resize(mat_size, mat_size);
                noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);
            }

            // reset the RHS
            if(CalculateResidualVectorFlag == true)
            {
                if(rRightHandSideVector.size() != mat_size)
                    rRightHandSideVector.resize(mat_size);
                noalias(rRightHandSideVector) = ZeroVector(mat_size);
            }

            // extract current displacements
            for(unsigned int node = 0; node < number_of_nodes; ++node)
                noalias(row(CurrentDisp, node)) = GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT);

            // loop through integration to compute LHS & RHS contribution
            for(unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
            {
                // compute weight at the integration_points
                double Weight = integration_points[PointNumber].Weight();
                if(dim == 2)
                    Weight *= GetProperties()[THICKNESS];

                // compute inverse of Jacobian
                MathUtils<double>::InvertMatrix(J[PointNumber], InvJ, DetJ);

                // compute the gradient w.r.t global coordinates
                noalias(DN_DX) = prod(DN_De[PointNumber], InvJ);

                // compute the B Operator
                this->CalculateBoperator(B, DN_DX);

                // compute the strain
                this->CalculateStrain(StrainVector, B, CurrentDisp);


                /***********VARIANT 1: classsical ******************/
//                this->CalculateElasticMatrix(TanC, E, NU);
//                noalias(StressVector) = prod(TanC, StrainVector);
                /***************************************************/

                /*******VARIANT 2: decompose strain****************/
//                CLParams Params;
//                Params.PointNumber = PointNumber;
//                Params.lambda = lambda;
//                Params.mu = mu;
//                this->CalculateStress(StrainVector, StressVector, Params);
//                this->CalculateTangent(StrainVector, TanC, Params);
                /***************************************************/

                /*******VARIANT 3: numerical tangent****************/
                CLParams Params;
                Params.PointNumber = PointNumber;
                Params.lambda = lambda;
                Params.mu = mu;
                Params.pertubation = 1.0e-8;
                this->CalculateStress(StrainVector, StressVector, Params);
                this->CalculateNumericalTangent(StrainVector, StressVector, TanC, Params);
                /***************************************************/

                // store the stresses
                noalias(mCurrentStresses[PointNumber]) = StressVector;

                /* compute the contribution of external forces to RHS */
                if(CalculateResidualVectorFlag == true)
                {
                    // body force
                    if(GetProperties().Has(BODY_FORCE))
                    {
                        const Vector& BodyForce = GetProperties()[BODY_FORCE];
                        for(unsigned int i = 0; i < number_of_nodes; ++i)
                            for(unsigned int j = 0; j < dim; ++j)
                                rRightHandSideVector(dim * i + j) += Ncontainer(PointNumber, i) * BodyForce[j] * Weight * DetJ;
                    }

                    // gravity
                    if(GetProperties().Has(GRAVITY) && GetProperties().Has(DENSITY))
                    {
                        double density = GetProperties()[DENSITY];
                        const Vector& Gravity = GetProperties()[GRAVITY];
                        for(unsigned int i = 0; i < number_of_nodes; ++i)
                            for(unsigned int j = 0; j < dim; ++j)
                                rRightHandSideVector(dim * i + j) += Ncontainer(PointNumber, i) * density * Gravity( i ) * Weight * DetJ;
                    }

                    // internal force
                    for(unsigned int i = 0; i < number_of_nodes; ++i)
                    {
                        for( unsigned int j = 0; j < dim; ++j)
                        {
                            InternalForces(j) = 0.0;
                            for(unsigned int k = 0; k < strain_size; ++k)
                                InternalForces(j) += B(k, dim * i + j) * StressVector(k) * Weight * DetJ;

                            rRightHandSideVector(dim * i + j) -= InternalForces(j);
                        }
                    }
                }

                /* compute the contribution to LHS */
                if(CalculateStiffnessMatrixFlag == true)
                    noalias( rLeftHandSideMatrix ) += prod(trans(B), (Weight * DetJ) * Matrix(prod(TanC, B)));
            }

            if(CalculateResidualVectorFlag == true)
            {
                // modify the right hand side to account for prescribed displacement
                // according to the book of Bazant & Jirasek, this scheme is more stable than the current scheme for prescribing displacement.
                // // However, I have to temporarily disable it to keep the consistency.
                for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
                {
                    if(GetGeometry()[node].IsFixed(DISPLACEMENT_X))
                    {
                        double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X);
                        for( unsigned int i = 0; i < mat_size; ++i )
                            rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * dim) * temp;
                    }
                    if(GetGeometry()[node].IsFixed(DISPLACEMENT_Y))
                    {
                        double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y);
                        for( unsigned int i = 0; i < mat_size; ++i )
                            rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * dim + 1) * temp;
                    }
                    if(GetGeometry()[node].IsFixed(DISPLACEMENT_Z) && dim == 3)
                    {
                        double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Z);
                        for( unsigned int i = 0; i < mat_size; ++i )
                            rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * dim + 2) * temp;
                    }
                }
            }

            if(Id() == 1)
            {
                KRATOS_WATCH(rRightHandSideVector)
                KRATOS_WATCH(rLeftHandSideMatrix)
            }
        }
        else if(rCurrentProcessInfo[FRACTIONAL_STEP] == 1) // solve for phase field
        {
            unsigned int mat_size = number_of_nodes;
            Matrix DN_DX(number_of_nodes, dim); // shape functions gradient in global coordinates
            double C; // temporary variable for crack phase field
            Vector D_C(dim); // temporary variable for crack phase field gradient
            Matrix InvJ(dim, dim); // temporary variable of inverse of Jacobian at integration point
            double DetJ; // temporary variable of determinant of Jacobian at integration point

            // fracturing parameters
            double kappa = GetProperties()[KAPPA]; // stabilized parameter
            double l0 = GetProperties()[LENGTH_SCALE]; // length scale of the crack phase field
            double Gc = GetProperties()[FRACTURE_ENERGY]; // fracture energy per unit volume

            // reset the LHS
            if(CalculateStiffnessMatrixFlag == true)
            {
                if(rLeftHandSideMatrix.size1() != mat_size || rLeftHandSideMatrix.size2() != mat_size)
                    rLeftHandSideMatrix.resize(mat_size, mat_size);
                noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);
            }

            // reset the RHS
            if(CalculateResidualVectorFlag == true)
            {
                if(rRightHandSideVector.size() != mat_size)
                    rRightHandSideVector.resize(mat_size);
                noalias(rRightHandSideVector) = ZeroVector(mat_size);
            }

            // loop through integration_points
            for(unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
            {
                // compute weight at the integration_points
                double Weight = integration_points[PointNumber].Weight();
                if(dim == 2)
                    Weight *= GetProperties()[THICKNESS];

                // compute determinant of Jacobian
                MathUtils<double>::InvertMatrix(J[PointNumber], InvJ, DetJ);

                // compute the gradient w.r.t global coordinates
                noalias(DN_DX) = prod(DN_De[PointNumber], InvJ);

                // compute the phase field at integration point
                C = 0.0;
                for(unsigned int i = 0; i < number_of_nodes; ++i)
                    C += Ncontainer(PointNumber, i) * GetGeometry()[i].GetSolutionStepValue(PHASE_FIELD);

                // compute the phase field gradient
                for(unsigned int j = 0; j < dim; ++j)
                {
                    D_C(j) = 0.0;
                    for(unsigned int i = 0; i < number_of_nodes; ++i)
                        D_C(j) += DN_DX(i, j) * GetGeometry()[i].GetSolutionStepValue(PHASE_FIELD);
                }

                if(GetProperties()[PHASE_FIELD_ORDER] == 2)
                {
                    double aux = 4 * l0 * (1 - kappa) * mReferenceEnergyDensity[PointNumber] / Gc + 1.0;

                    /* compute contribution to RHS */
                    if(CalculateResidualVectorFlag == true)
                    {
                        for(unsigned int i = 0; i < number_of_nodes; ++i)
                        {
                            rRightHandSideVector(i) += Ncontainer(PointNumber, i) * Weight * DetJ;
                            rRightHandSideVector(i) -= aux * Ncontainer(PointNumber, i) * C * Weight * DetJ;
                            double DNxDC = 0.0;
                            for(unsigned int j = 0; j < dim; ++j)
                                DNxDC += DN_DX(i, j) * D_C(j);
                            rRightHandSideVector(i) -= 4.0 * pow(l0, 2) * DNxDC * Weight * DetJ;
                        }
                    }

                    /* compute contribution to LHS */
                    if(CalculateStiffnessMatrixFlag == true)
                    {
                        for(unsigned int i = 0; i < number_of_nodes; ++i)
                        {
                            for(unsigned int j = 0; j < number_of_nodes; ++j)
                            {
                                rLeftHandSideMatrix(i, j) += aux * Ncontainer(PointNumber, i) * Ncontainer(PointNumber, j) * Weight * DetJ;
                                double DNixDNj = 0.0;
                                for(unsigned int k = 0; k < dim ; ++k)
                                    DNixDNj += DN_DX(i, k) * DN_DX(j, k);
                                rLeftHandSideMatrix(i, j) += 4.0 * pow(l0, 2) * DNixDNj * Weight * DetJ;
                            }
                        }
                    }
                }
                else if(GetProperties()[PHASE_FIELD_ORDER] == 4)
                {
                    // compute the second derivatives w.r.t physical coordinates
                    std::vector<Vector> D2N_DX2;
                    CalculateSecondDerivatives(D2N_DX2, integration_points[PointNumber]);

                    // compute the phase field Laplace
                    Vector D_C2(dim); // this array contains [D^2C/dX^2 D^2C/dY^2 D^2C/dZ^2]
                    double Delta_C = 0.0;
                    for(unsigned int j = 0; j < dim; ++j)
                    {
                        D_C2(j) = 0.0;
                        for(unsigned int i = 0; i < number_of_nodes; ++i)
                            D_C2(j) += D2N_DX2[i](j) * GetGeometry()[i].GetSolutionStepValue(PHASE_FIELD);
                        Delta_C += D_C2(j);
                    }

                    double aux = 4 * l0 * (1 - kappa) * mReferenceEnergyDensity[PointNumber] / Gc + 1.0;

                    /* compute contribution to RHS */
                    if(CalculateResidualVectorFlag == true)
                    {
                        for(unsigned int i = 0; i < number_of_nodes; ++i)
                        {
                            rRightHandSideVector(i) += Ncontainer(PointNumber, i) * Weight * DetJ;
                            rRightHandSideVector(i) -= aux * Ncontainer(PointNumber, i) * C * Weight * DetJ;

                            double DNxDC = 0.0;
                            for(unsigned int j = 0; j < dim; ++j)
                                DNxDC += DN_DX(i, j) * D_C(j);
                            rRightHandSideVector(i) -= 2.0 * pow(l0, 2) * DNxDC * Weight * DetJ;

//                            double D2NxD2C = 0.0;
//                            for(unsigned int j = 0; j < dim; ++j)
//                                D2NxD2C += D2N_DX2[i](j) * D_C2(j);
//                            rRightHandSideVector(i) -= pow(l0, 4) * D2NxD2C * Weight * DetJ;
                            double Delta_N = 0.0;
                            for(unsigned int j = 0; j < dim; ++j)
                                Delta_N += D2N_DX2[i](j);
                            rRightHandSideVector(i) -= pow(l0, 4) * Delta_N * Delta_C * Weight * DetJ;
                        }
                    }

                    /* compute contribution to LHS */
                    if(CalculateStiffnessMatrixFlag == true)
                    {
                        for(unsigned int i = 0; i < number_of_nodes; ++i)
                        {
                            for(unsigned int j = 0; j < number_of_nodes; ++j)
                            {
                                rLeftHandSideMatrix(i, j) += aux * Ncontainer(PointNumber, i) * Ncontainer(PointNumber, j) * Weight * DetJ;

                                double DNixDNj = 0.0;
                                for(unsigned int k = 0; k < dim ; ++k)
                                    DNixDNj += DN_DX(i, k) * DN_DX(j, k);
                                rLeftHandSideMatrix(i, j) += 2.0 * pow(l0, 2) * DNixDNj * Weight * DetJ;

//                                double D2NixD2Nj = 0.0;
//                                for(unsigned int k = 0; k < dim ; ++k)
//                                    D2NixD2Nj += D2N_DX2[i](k) * D2N_DX2[j](k);
//                                rLeftHandSideMatrix(i, j) += pow(l0, 4) * D2NixD2Nj * Weight * DetJ;
                                double Delta_Ni = 0.0;
                                double Delta_Nj = 0.0;
                                for(unsigned int k = 0; k < dim; ++j)
                                {
                                    Delta_Ni += D2N_DX2[i](k);
                                    Delta_Nj += D2N_DX2[j](k);
                                }
                                rLeftHandSideMatrix(i, j) += pow(l0, 4) * Delta_Ni * Delta_Nj * Weight * DetJ;
                            }
                        }
                    }
                }
                else
                    KRATOS_THROW_ERROR(std::logic_error, "THis PHASE_FIELD_ORDER is not or not yet supported: ", GetProperties()[PHASE_FIELD_ORDER])
            }
//            KRATOS_WATCH(rRightHandSideVector)
//            KRATOS_WATCH(rLeftHandSideMatrix)
        }
        else // solve for nothing (should not happen)
        {
            noalias(rLeftHandSideMatrix) = ZeroMatrix(rLeftHandSideMatrix.size1(), rLeftHandSideMatrix.size2());
            noalias(rRightHandSideVector) = ZeroVector(rRightHandSideVector.size());
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        // clean the geometry internal data
        GetGeometry().Clean();
        #endif
    }

    void PhaseFieldKinematicLinear::CalculateStress(const Vector& StrainVector, Vector& StressVector, CLParams& Params)
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        // perform the spectral decomposition for strain
        Matrix StrainTensor(dim, dim);
        if(dim == 2)
        {
            IsotropicTensorUtility<2>::StrainVectorToTensor(StrainVector, StrainTensor);
            Params.decomp_params = IsotropicTensorUtility<2>::SpectralDecomposition(StrainTensor, Params.e, Params.E);
        }
        else if(dim == 3)
        {
            IsotropicTensorUtility<3>::StrainVectorToTensor(StrainVector, StrainTensor);
            Params.decomp_params = IsotropicTensorUtility<3>::SpectralDecomposition(StrainTensor, Params.e, Params.E);
        }
        if(Params.decomp_params.flag != 0)
        {
            KRATOS_WATCH(StrainVector)
            KRATOS_WATCH(Params.e)
            for(int i = 0; i < Params.E.size(); ++i)
                KRATOS_WATCH(Params.E[i])
            std::stringstream ss;
            ss << "Error running decomposition at element " << Id() << ", integration point " << Params.PointNumber << std::endl;
            ss << "Error code: " << Params.decomp_params.flag;
            KRATOS_THROW_ERROR(std::runtime_error, ss.str(), "")
        }

        // compute the trace of the strain
        double trace_strain = 0.0;
        for(unsigned int i = 0; i < dim; ++i)
            trace_strain += StrainVector(i);

        // compute contribution to stress tensor from trace part of strain
        Matrix StressTensor = ZeroMatrix(dim, dim);
        noalias(StressTensor) += Params.lambda * trace_strain * IdentityMatrix(dim);

        // compute contribution to stress tensor from positive part of strain
        Matrix StrainTensor_decompose = ZeroMatrix(dim, dim);
        if(dim == 2)
        {
            MyIsotropicTensorFunction<2> Phi;
            IsotropicTensorUtility<2>::Value(StrainTensor, Params.e, Params.E, Params.decomp_params, Phi, StrainTensor_decompose);
        }
        else if(dim == 3)
        {
            MyIsotropicTensorFunction<3> Phi;
            IsotropicTensorUtility<3>::Value(StrainTensor, Params.e, Params.E, Params.decomp_params, Phi, StrainTensor_decompose);
        }
        noalias(StressTensor) += 2.0 * Params.mu * StrainTensor_decompose;

        if(dim == 2)
            IsotropicTensorUtility<2>::StressTensorToVector(StressTensor, StressVector);
        else if(dim == 3)
            IsotropicTensorUtility<3>::StressTensorToVector(StressTensor, StressVector);
    }

    void PhaseFieldKinematicLinear::CalculateTangent(const Vector& StrainVector, Matrix& TanC, CLParams& Params)
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int strain_size = dim * (dim + 1) / 2;

        /* compute the tangent */
        // compute tangent from trace part of strain
        noalias(TanC) = ZeroMatrix(strain_size, strain_size);
        for(unsigned int i = 0; i < dim; ++i)
            for(unsigned int j = 0; j < dim; ++j)
                TanC(i, j) += Params.lambda;

        // compute the tangent from decomposed part of strain
        Matrix Ct;
        Matrix StrainTensor(dim, dim);
        if(dim == 2)
        {
            IsotropicTensorUtility<2>::StrainVectorToTensor(StrainVector, StrainTensor);
            MyIsotropicTensorFunction<2> Phi;
            IsotropicTensorUtility<2>::Fourth_Order_Tensor D;
            IsotropicTensorUtility<2>::InitializeFourthOrderTensor(D);
            IsotropicTensorUtility<2>::Derivative(StrainTensor, Params.e, Params.E, Params.decomp_params, Phi, D);
            IsotropicTensorUtility<2>::FourthOrderTensorToMatrix(D, Ct);
        }
        else if(dim == 3)
        {
            IsotropicTensorUtility<3>::StrainVectorToTensor(StrainVector, StrainTensor);
            MyIsotropicTensorFunction<3> Phi;
            IsotropicTensorUtility<3>::Fourth_Order_Tensor D;
            IsotropicTensorUtility<3>::InitializeFourthOrderTensor(D);
            IsotropicTensorUtility<3>::Derivative(StrainTensor, Params.e, Params.E, Params.decomp_params, Phi, D);
            IsotropicTensorUtility<3>::FourthOrderTensorToMatrix(D, Ct);
        }
        noalias(TanC) += 2.0 * Params.mu * Ct;
    }

    void PhaseFieldKinematicLinear::CalculateNumericalTangent(const Vector& StrainVector, const Vector& StressVector, Matrix& TanC, CLParams& Params)
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int strain_size = dim * (dim + 1) / 2;

        /* compute the tangent */
        Vector PertubedStrainVector(strain_size);
        Vector PertubedStressVector(strain_size);
        for(std::size_t col = 0; col < strain_size; ++col)
        {
            // pertube the strain component
            noalias(PertubedStrainVector) = StrainVector;
            PertubedStrainVector(col) += Params.pertubation;

            // compute the resulting pertubed stress
            this->CalculateStress(PertubedStrainVector, PertubedStressVector, Params);

            // compute the associated column of the numerical tangent
            for(std::size_t row = 0; row < strain_size; ++row)
                TanC(row, col) = (PertubedStressVector(row) - StressVector(row)) / Params.pertubation;
        }
    }

    void PhaseFieldKinematicLinear::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
    {
    }

    void PhaseFieldKinematicLinear::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
//        unsigned int dim = GetGeometry().WorkingSpaceDimension();
//        double DetJ;
// 
//        /* compute the total energy density */
//        std::vector<double> EnergyDensity;
//        this->GetValueOnIntegrationPoints(ENERGY_FUNCTIONAL_DENSITY, EnergyDensity, CurrentProcessInfo);
// 
//        // extract the integration_points
//        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);
// 
//        // compute Jacobian
//        GeometryType::JacobiansType J(integration_points.size());
//        J = GetGeometry().Jacobian(J, mThisIntegrationMethod);
// 
//        double TotalEnergy = 0.0;
//        for(int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
//        {
//            DetJ = MathUtils<double>::Det(J[PointNumber]);
// 
//            double Weight = integration_points[PointNumber].Weight();
//            if(dim == 2)
//                Weight *= GetProperties()[THICKNESS];
// 
//            TotalEnergy += EnergyDensity[PointNumber] * Weight * DetJ;
//        }
// 
//        this->SetValue(ENERGY_FUNCTIONAL, TotalEnergy);
    }

    int PhaseFieldKinematicLinear::Check(const Kratos::ProcessInfo& rCurrentProcessInfo)
    {
        if( this->GetProperties().Has( YOUNG_MODULUS ) == false )
        {
            KRATOS_THROW_ERROR( std::logic_error, "YOUNG_MODULUS not provided for property", this->GetProperties().Id() )
        }

        if( this->GetProperties().Has( POISSON_RATIO ) == false )
        {
            KRATOS_THROW_ERROR( std::logic_error, "POISSON_RATIO not provided for property", this->GetProperties().Id() )
        }

        if( this->GetProperties().Has( KAPPA ) == false )
        {
            KRATOS_THROW_ERROR( std::logic_error, "KAPPA (stabilized parameter for phase field calculation) not provided for property", this->GetProperties().Id() )
        }

        if( this->GetProperties().Has( LENGTH_SCALE ) == false )
        {
            KRATOS_THROW_ERROR( std::logic_error, "LENGTH_SCALE not provided for property", this->GetProperties().Id() )
        }

        if( this->GetProperties().Has( FRACTURE_ENERGY ) == false )
        {
            KRATOS_THROW_ERROR( std::logic_error, "FRACTURE_ENERGY not provided for property", this->GetProperties().Id() )
        }

        if( this->GetProperties().Has( PHASE_FIELD_ORDER ) == false )
        {
            KRATOS_THROW_ERROR( std::logic_error, "PHASE_FIELD_ORDER not provided for property", this->GetProperties().Id() )
        }

        return 0;
    }

    ////////////////ACCESS////////////////////////////////

    void PhaseFieldKinematicLinear::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
        if(rVariable == PHASE_FIELD)
        {
            #ifdef ENABLE_BEZIER_GEOMETRY
            // initialize the geometry
            GetGeometry().Initialize(mThisIntegrationMethod);
            #endif

            // extract the integration_points
            const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

            // extract the shape function values
            const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);

            // resize as needed
            if(rValues.size() != integration_points.size())
                rValues.resize(integration_points.size());

            // loop through integration_points to compute the phase field at the integration_points
            for(unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
            {
                rValues[PointNumber] = 0.0;
                for(unsigned int i = 0; i < GetGeometry().size(); ++i)
                    rValues[PointNumber] += Ncontainer(PointNumber, i) * GetGeometry()[i].GetSolutionStepValue(PHASE_FIELD);
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            // clean the geometry internal data
            GetGeometry().Clean();
            #endif
        }

        if(rVariable == REFERENCE_ENERGY_DENSITY) // Note that this is psi_plus_0 term and not the energy functional
        {
            if(rValues.size() != mReferenceEnergyDensity.size())
                rValues.resize(mReferenceEnergyDensity.size());
            std::copy(mReferenceEnergyDensity.begin(), mReferenceEnergyDensity.end(), rValues.begin());
        }

        if(rVariable == ENERGY_FUNCTIONAL_DENSITY)
        {
            unsigned int number_of_nodes = GetGeometry().size();
            unsigned int dim = GetGeometry().WorkingSpaceDimension();
            unsigned int strain_size = dim * (dim + 1) / 2;
            unsigned int mat_size = dim * number_of_nodes;
            Matrix B(strain_size, mat_size);
            Vector StrainVector(strain_size);
            Matrix DN_DX(number_of_nodes, dim);
            Matrix CurrentDisp(number_of_nodes, dim);
            Matrix InvJ(dim, dim);
            double DetJ;
            std::vector<array_1d<double, 3> > n(3);
            double e[3];
            Matrix e_plus(3, 3);
            Matrix e_minus(3, 3);
            double C;
            Vector D_C(3);

            double E = GetProperties()[YOUNG_MODULUS];
            double NU = GetProperties()[POISSON_RATIO];
            double lambda = E * NU / (1.0 + NU) / (1.0 - 2.0 * NU);
            double mu = 0.5 * E / (1.0 + NU);
            double kappa = GetProperties()[KAPPA]; // stabilized parameter
            double l0 = GetProperties()[LENGTH_SCALE]; // length scale of the crack phase field
            double Gc = GetProperties()[FRACTURE_ENERGY]; // fracture energy per unit volume

            #ifdef ENABLE_BEZIER_GEOMETRY
            // initialize geometry
            GetGeometry().Initialize(mThisIntegrationMethod);
            #endif

            // extract the integration points
            const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

            // resize as needed
            if(rValues.size() != integration_points.size())
                rValues.resize(integration_points.size());

            // extract the shape function values
            const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);

            // extract the shape function local gradients
            const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

            // compute Jacobian
            GeometryType::JacobiansType J(integration_points.size());
            J = GetGeometry().Jacobian(J, mThisIntegrationMethod);

            // extract current displacements
            for(unsigned int node = 0; node < number_of_nodes; ++node)
                noalias(row(CurrentDisp, node)) = GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT);

            // extract the current phase field

            // loop through integration
            for(unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
            {
                // compute inverse of Jacobian
                MathUtils<double>::InvertMatrix(J[PointNumber], InvJ, DetJ);

                // compute the gradient w.r.t global coordinates
                noalias(DN_DX) = prod(DN_De[PointNumber], InvJ);

                // compute the B Operator
                this->CalculateBoperator(B, DN_DX);

                // compute the strain
                this->CalculateStrain(StrainVector, B, CurrentDisp);

                double trace_strain_plus = 0.0;
                double trace_strain_minus = 0.0;
                for(unsigned int i = 0; i < dim; ++i)
                    trace_strain_plus += StrainVector(i);
                if(trace_strain_plus < 0.0)
                {
                    trace_strain_minus = trace_strain_plus;
                    trace_strain_plus = 0.0;
                }

                // perform the spectral decomposition
                if(dim == 2)
                    spectral_decomposition(StrainVector(0), StrainVector(1), 0.0, 0.5 * StrainVector(2), 0.0, 0.0, e[0], e[1], e[2], n[0], n[1], n[2]);
                else if(dim == 3)
                    spectral_decomposition(StrainVector(0), StrainVector(1), StrainVector(2), 0.5 * StrainVector(3), 0.5 * StrainVector(4), 0.5 * StrainVector(5), e[0], e[1], e[2], n[0], n[1], n[2]);

                // compute positive part of strain tensor
                noalias(e_plus) = ZeroMatrix(3, 3);
                noalias(e_minus) = ZeroMatrix(3, 3);
                for(unsigned int i = 0; i < 3; ++i)
                {
                    if(e[i] > 0.0)
                    {
                        for(unsigned int j = 0; j < 3; ++j)
                            for(unsigned int k = 0; k < 3; ++k)
                                e_plus(j, k) += e[i] * n[i](j) * n[i](k);
                    }
                    else
                    {
                        for(unsigned int j = 0; j < 3; ++j)
                            for(unsigned int k = 0; k < 3; ++k)
                                e_minus(j, k) += e[i] * n[i](j) * n[i](k);
                    }
                }

                // compute energy density at the integration point
                double psi_plus = 0.0;
                psi_plus += 0.5 * lambda * pow(trace_strain_plus, 2);
                for(unsigned int j = 0; j < 3; ++j)
                    for(unsigned int k = 0; k < 3; ++k)
                        psi_plus += mu * pow(e_plus(j, k), 2);

                double psi_minus = 0.0;
                psi_minus += 0.5 * lambda * pow(trace_strain_minus, 2);
                for(unsigned int j = 0; j < 3; ++j)
                    for(unsigned int k = 0; k < 3; ++k)
                        psi_minus += mu * pow(e_minus(j, k), 2);

                // compute the phase field at integration point
                C = 0.0;
                for(unsigned int i = 0; i < number_of_nodes; ++i)
                    C += Ncontainer(PointNumber, i) * GetGeometry()[i].GetSolutionStepValue(PHASE_FIELD);

                // compute the phase field gradient
                for(unsigned int j = 0; j < dim; ++j)
                {
                    D_C(j) = 0.0;
                    for(unsigned int i = 0; i < number_of_nodes; ++i)
                        D_C(j) += DN_DX(i, j) * GetGeometry()[i].GetSolutionStepValue(PHASE_FIELD);
                }

                // compute the fracture energy density
                double Gamma_C = 0.0;
                if(GetProperties()[PHASE_FIELD_ORDER] == 2)
                {
                    Gamma_C = 0.25 / l0 * ( pow(C - 1, 2) + 4 * pow(l0 * norm_2(D_C), 2));
                }
                else if(GetProperties()[PHASE_FIELD_ORDER] == 4)
                {
                    // compute the second derivatives w.r.t physical coordinates
                    std::vector<Vector> D2N_DX2;
                    CalculateSecondDerivatives(D2N_DX2, integration_points[PointNumber]);

                    // compute the phase field Laplace
                    Vector D_C2(dim); // this array contains [D^2C/dX^2 D^2C/dY^2 D^2C/dZ^2]
                    double Delta_C = 0.0;
                    for(unsigned int j = 0; j < dim; ++j)
                    {
                        D_C2(j) = 0.0;
                        for(unsigned int i = 0; i < number_of_nodes; ++i)
                            D_C2(j) += D2N_DX2[i](j) * GetGeometry()[i].GetSolutionStepValue(PHASE_FIELD);
                        Delta_C += D_C2(j);
                    }

                    Gamma_C = 0.25 / l0 * ( pow(C - 1, 2) + 2 * pow(l0 * norm_2(D_C), 2) + pow(l0, 4) * pow(Delta_C, 2));
                }

                // compute the total energy density
                rValues[PointNumber] = ((1 - kappa) * pow(C, 2) + kappa) * psi_plus + psi_minus + Gc * Gamma_C;
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            // clean the geometry internal data
            GetGeometry().Clean();
            #endif
        }
    }

    void PhaseFieldKinematicLinear::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int strain_size = dim * (dim + 1) / 2;

        // extract the integration points
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

        // resize as needed
        if ( rValues.size() != integration_points.size() )
            rValues.resize( integration_points.size() );

        if(rVariable == STRESSES)
        {
            for ( unsigned int i = 0; i < integration_points.size(); ++i )
            {
                if ( rValues[i].size() != strain_size )
                    rValues[i].resize( strain_size );
                noalias( rValues[i] ) = mCurrentStresses[i];
            }
        }
    }

    void PhaseFieldKinematicLinear::GetValueOnIntegrationPoints(const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
        // extract the integration points
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

        // resize as needed
        if ( rValues.size() != integration_points.size() )
            rValues.resize( integration_points.size() );

        if(rVariable == INTEGRATION_POINT_GLOBAL_COORDINATES)
        {
            #ifdef ENABLE_BEZIER_GEOMETRY
            // initialize the geometry
            GetGeometry().Initialize(mThisIntegrationMethod);
            #endif

            // extract the shape function values
            const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);

            for ( unsigned int i = 0; i < integration_points.size(); ++i )
            {
                noalias(rValues[i]) = ZeroVector(3);

                for( unsigned int j = 0; j < GetGeometry().size(); ++j )
                    rValues[i](0) += Ncontainer(i, j) * GetGeometry()[j].X0();

                for( unsigned int j = 0; j < GetGeometry().size(); ++j )
                    rValues[i](1) += Ncontainer(i, j) * GetGeometry()[j].Y0();

                for( unsigned int j = 0; j < GetGeometry().size(); ++j )
                    rValues[i](2) += Ncontainer(i, j) * GetGeometry()[j].Z0();
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            // clean the geometry internal data
            GetGeometry().Clean();
            #endif
        }
    }

    void PhaseFieldKinematicLinear::CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
        GetValueOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    }

    void PhaseFieldKinematicLinear::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
        GetValueOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    }

    void PhaseFieldKinematicLinear::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
        GetValueOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    }

    ///////////////////ASSIGN/////////////////////////////

    void PhaseFieldKinematicLinear::SetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
        if(rVariable == REFERENCE_ENERGY_DENSITY)
        {
            if(rValues.size() != mReferenceEnergyDensity.size())
                KRATOS_THROW_ERROR(std::logic_error, "The size of given values and number of internal values are incompatible", "")
            std::copy(rValues.begin(), rValues.end(), mReferenceEnergyDensity.begin());
        }
    }

    ///////////////PRIVATE SUBROUTINES////////////////////

    void PhaseFieldKinematicLinear::CalculateBoperator(Matrix& B_Operator, const Matrix& DN_DX)
    {
        unsigned int number_of_nodes = GetGeometry().PointsNumber();
        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int strain_size = dim * (dim + 1) / 2;

        if(B_Operator.size1() != strain_size || B_Operator.size2() != number_of_nodes * dim)
            B_Operator.resize(strain_size, number_of_nodes * dim);
        noalias( B_Operator ) = ZeroMatrix( strain_size, number_of_nodes * dim );

        if(dim == 2)
        {
            for(unsigned int i = 0; i < number_of_nodes; ++i)
            {
                B_Operator( 0, i*2 ) = DN_DX( i, 0 );
                B_Operator( 1, i*2 + 1 ) = DN_DX( i, 1 );
                B_Operator( 2, i*2 ) = DN_DX( i, 1 );
                B_Operator( 2, i*2 + 1 ) = DN_DX( i, 0 );
            }
        }
        else if(dim == 3)
        {
            for(unsigned int i = 0; i < number_of_nodes; ++i)
            {
                B_Operator( 0, i*3 ) = DN_DX( i, 0 );
                B_Operator( 1, i*3 + 1 ) = DN_DX( i, 1 );
                B_Operator( 2, i*3 + 2 ) = DN_DX( i, 2 );
                B_Operator( 3, i*3 ) = DN_DX( i, 1 );
                B_Operator( 3, i*3 + 1 ) = DN_DX( i, 0 );
                B_Operator( 4, i*3 + 1 ) = DN_DX( i, 2 );
                B_Operator( 4, i*3 + 2 ) = DN_DX( i, 1 );
                B_Operator( 5, i*3 ) = DN_DX( i, 2 );
                B_Operator( 5, i*3 + 2 ) = DN_DX( i, 0 );
            }
        }
    }

    void PhaseFieldKinematicLinear::CalculateStrain(Vector& StrainVector, const Matrix& B, const Matrix& Displacements)
    {
        unsigned int Dim = GetGeometry().WorkingSpaceDimension();
        unsigned int strain_size = Dim * (Dim + 1) / 2;

        noalias(StrainVector) = ZeroVector(strain_size);
        for(unsigned int node = 0; node < GetGeometry().size(); ++node)
            for(unsigned int item = 0; item < strain_size; ++item)
                for(unsigned int dim = 0; dim < Dim; ++dim)
                    StrainVector[item] += B(item, Dim * node + dim) * Displacements(node, dim);
    }

    void PhaseFieldKinematicLinear::CalculateElasticMatrix(Matrix& C, double E, double NU)
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int strain_size = dim * (dim + 1) / 2;
        if(C.size1() != strain_size || C.size2() != strain_size)
            C.resize(strain_size, strain_size);

        if(dim == 2) // plane strain case
        {
            double c1 = E * ( 1.00 - NU ) / (( 1.00 + NU ) * ( 1.00 - 2 * NU ) );
            double c2 = E * NU / (( 1.00 + NU ) * ( 1.00 - 2 * NU ) );
            double c3 = 0.5 * E / ( 1 + NU );

            C( 0, 0 ) = c1;
            C( 0, 1 ) = c2;
            C( 0, 2 ) = 0.0;
            C( 1, 0 ) = c2;
            C( 1, 1 ) = c1;
            C( 1, 2 ) = 0.0;
            C( 2, 0 ) = 0.0;
            C( 2, 1 ) = 0.0;
            C( 2, 2 ) = c3;
        }
        else if(dim == 3)
        {
            double c1 = E / (( 1.00 + NU ) * ( 1 - 2 * NU ) );
            double c2 = c1 * ( 1 - NU );
            double c3 = c1 * NU;
            double c4 = c1 * 0.5 * ( 1 - 2 * NU );

            C( 0, 0 ) = c2;
            C( 0, 1 ) = c3;
            C( 0, 2 ) = c3;
            C( 0, 3 ) = 0.0;
            C( 0, 4 ) = 0.0;
            C( 0, 5 ) = 0.0;
            C( 1, 0 ) = c3;
            C( 1, 1 ) = c2;
            C( 1, 2 ) = c3;
            C( 1, 3 ) = 0.0;
            C( 1, 4 ) = 0.0;
            C( 1, 5 ) = 0.0;
            C( 2, 0 ) = c3;
            C( 2, 1 ) = c3;
            C( 2, 2 ) = c2;
            C( 2, 3 ) = 0.0;
            C( 2, 4 ) = 0.0;
            C( 2, 5 ) = 0.0;
            C( 3, 0 ) = 0.0;
            C( 3, 1 ) = 0.0;
            C( 3, 2 ) = 0.0;
            C( 3, 3 ) = c4;
            C( 3, 4 ) = 0.0;
            C( 3, 5 ) = 0.0;
            C( 4, 0 ) = 0.0;
            C( 4, 1 ) = 0.0;
            C( 4, 2 ) = 0.0;
            C( 4, 3 ) = 0.0;
            C( 4, 4 ) = c4;
            C( 4, 5 ) = 0.0;
            C( 5, 0 ) = 0.0;
            C( 5, 1 ) = 0.0;
            C( 5, 2 ) = 0.0;
            C( 5, 3 ) = 0.0;
            C( 5, 4 ) = 0.0;
            C( 5, 5 ) = c4;
        }
    }

    /// REMARKS: the arrangement is like this
    /// in 2D: [d^2/dx^2 d^2/dy^2 d^2/dxdy]
    /// in 3D: [d^2/dx^2 d^2/dy^2 d^2/dz^2 d^2/dxdy d^2/dydz d^2/dxdz]
    void PhaseFieldKinematicLinear::CalculateSecondDerivatives(std::vector<Vector>& rD2N_DX2, const CoordinatesArrayType& rPoint)
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int number_of_nodes = GetGeometry().size();

        #ifdef ENABLE_BEZIER_GEOMETRY
        // initialize the geometry
        GetGeometry().Initialize(mThisIntegrationMethod);
        #endif

        // compute the Jacobian
        Matrix J;
        J = GetGeometry().Jacobian(J, rPoint);

        // compute inverse of Jacobian
        Matrix InvJ;
        double DetJ;
        MathUtils<double>::InvertMatrix(J, InvJ, DetJ);

        // compute the shape function local gradients
        Matrix DN_De;
        DN_De = GetGeometry().ShapeFunctionsLocalGradients(DN_De, rPoint);

        // compute the shape function gradients w.r.t physical coordinates
        Matrix DN_DX(number_of_nodes, dim);
        noalias(DN_DX) = prod(DN_De, InvJ);

        // compute the shape function local second derivatives
        ShapeFunctionsSecondDerivativesType D2N_De2;
        D2N_De2 = GetGeometry().ShapeFunctionsSecondDerivatives(D2N_De2, rPoint);

        // compute the second derivatives of physical coordinates
        Matrix D2X_De2(dim, dim);
        Matrix D2Y_De2(dim, dim);
        Matrix D2Z_De2(dim, dim);
        for(unsigned int i = 0; i < dim; ++i)
            for(unsigned int j = 0; j < dim; ++j)
            {
                D2X_De2(i, j) = 0.0;
                for(unsigned int k = 0; k < number_of_nodes; ++k)
                    D2X_De2(i, j) += GetGeometry()[k].X0() * D2N_De2[k](i, j);

                D2Y_De2(i, j) = 0.0;
                for(unsigned int k = 0; k < number_of_nodes; ++k)
                    D2Y_De2(i, j) += GetGeometry()[k].Y0() * D2N_De2[k](i, j);

                if(dim > 2)
                {
                    D2Z_De2(i, j) = 0.0;
                    for(unsigned int k = 0; k < number_of_nodes; ++k)
                        D2Z_De2(i, j) += GetGeometry()[k].Z0() * D2N_De2[k](i, j);
                }
            }

        // compute the second derivatives in physical space
        unsigned int mat_size = dim * (dim + 1) / 2;
        Matrix D2(mat_size, mat_size);
        if(dim == 2)
        {
            D2(0, 0) = pow(J(0, 0), 2);
            D2(0, 1) = pow(J(0, 1), 2);
            D2(0, 2) = J(0, 0) * J(0, 1);

            D2(1, 0) = pow(J(1, 0), 2);
            D2(1, 1) = pow(J(1, 1), 2);
            D2(1, 2) = J(1, 0) * J(1, 1);

            D2(2, 0) = 2.0 * J(0, 0) * J(1, 0);
            D2(2, 1) = 2.0 * J(0, 1) * J(1, 1);
            D2(2, 2) = J(0, 0) * J(1, 1) + J(1, 0) * J(0, 1);
        }
        else if(dim == 3)
        {
            D2(0, 0) = pow(J(0, 0), 2);
            D2(0, 1) = pow(J(0, 1), 2);
            D2(0, 2) = pow(J(0, 2), 2);
            D2(0, 3) = J(0, 0) * J(0, 1);
            D2(0, 4) = J(0, 1) * J(0, 2);
            D2(0, 5) = J(0, 0) * J(0, 2);

            D2(1, 0) = pow(J(1, 0), 2);
            D2(1, 1) = pow(J(1, 1), 2);
            D2(1, 2) = pow(J(1, 2), 2);
            D2(1, 3) = J(1, 0) * J(1, 1);
            D2(1, 4) = J(1, 1) * J(1, 2);
            D2(1, 5) = J(1, 0) * J(1, 2);

            D2(2, 0) = pow(J(2, 0), 2);
            D2(2, 1) = pow(J(2, 1), 2);
            D2(2, 2) = pow(J(2, 2), 2);
            D2(2, 3) = J(2, 0) * J(2, 1);
            D2(2, 4) = J(2, 1) * J(2, 2);
            D2(2, 5) = J(2, 0) * J(2, 2);

            D2(3, 0) = 2.0 * J(0, 0) * J(1, 0);
            D2(3, 1) = 2.0 * J(0, 1) * J(1, 1);
            D2(3, 2) = 2.0 * J(0, 2) * J(1, 2);
            D2(3, 3) = J(0, 0) * J(1, 1) + J(1, 0) * J(0, 1);
            D2(3, 4) = J(0, 1) * J(1, 2) + J(1, 2) * J(0, 1);
            D2(3, 5) = J(0, 0) * J(1, 2) + J(1, 2) * J(0, 0);

            D2(4, 0) = 2.0 * J(1, 0) * J(2, 0);
            D2(4, 1) = 2.0 * J(1, 1) * J(2, 1);
            D2(4, 2) = 2.0 * J(1, 2) * J(2, 2);
            D2(4, 3) = J(1, 0) * J(2, 1) + J(2, 0) * J(1, 1);
            D2(4, 4) = J(1, 1) * J(2, 2) + J(2, 2) * J(1, 1);
            D2(4, 5) = J(1, 0) * J(2, 2) + J(2, 2) * J(1, 0);

            D2(5, 0) = 2.0 * J(0, 0) * J(2, 0);
            D2(5, 1) = 2.0 * J(0, 1) * J(2, 1);
            D2(5, 2) = 2.0 * J(0, 2) * J(2, 2);
            D2(5, 3) = J(0, 0) * J(2, 1) + J(2, 0) * J(0, 1);
            D2(5, 4) = J(0, 1) * J(2, 2) + J(2, 2) * J(0, 1);
            D2(5, 5) = J(0, 0) * J(2, 2) + J(2, 2) * J(0, 0);
        }

        Matrix InvD2(mat_size, mat_size);
        SD_MathUtils<double>::InvertMatrix(D2, InvD2);

        if(rD2N_DX2.size() != number_of_nodes)
            rD2N_DX2.resize(number_of_nodes);

        Vector R2(mat_size);
        for(unsigned int i = 0; i < number_of_nodes; ++i)
        {
            if(rD2N_DX2[i].size() != mat_size)
                rD2N_DX2[i].resize(mat_size);

            if(dim == 2)
            {
                R2(0) = D2N_De2[i](0, 0) - DN_DX(i, 0) * D2X_De2(0, 0) - DN_DX(i, 1) * D2Y_De2(0, 0);
                R2(1) = D2N_De2[i](1, 1) - DN_DX(i, 0) * D2X_De2(1, 1) - DN_DX(i, 1) * D2Y_De2(1, 1);
                R2(2) = D2N_De2[i](0, 1) - DN_DX(i, 0) * D2X_De2(0, 1) - DN_DX(i, 1) * D2Y_De2(0, 1);
            }
            else if(dim == 3)
            {
                R2(0) = D2N_De2[i](0, 0) - DN_DX(i, 0) * D2X_De2(0, 0) - DN_DX(i, 1) * D2Y_De2(0, 0) - DN_DX(i, 2) * D2Z_De2(0, 0);
                R2(1) = D2N_De2[i](1, 1) - DN_DX(i, 0) * D2X_De2(1, 1) - DN_DX(i, 1) * D2Y_De2(1, 1) - DN_DX(i, 2) * D2Z_De2(1, 1);
                R2(2) = D2N_De2[i](2, 2) - DN_DX(i, 0) * D2X_De2(2, 2) - DN_DX(i, 1) * D2Y_De2(2, 2) - DN_DX(i, 2) * D2Z_De2(2, 2);
                R2(3) = D2N_De2[i](0, 1) - DN_DX(i, 0) * D2X_De2(0, 1) - DN_DX(i, 1) * D2Y_De2(0, 1) - DN_DX(i, 2) * D2Z_De2(0, 1);
                R2(4) = D2N_De2[i](1, 2) - DN_DX(i, 0) * D2X_De2(1, 2) - DN_DX(i, 1) * D2Y_De2(1, 2) - DN_DX(i, 2) * D2Z_De2(1, 2);
                R2(5) = D2N_De2[i](0, 2) - DN_DX(i, 0) * D2X_De2(0, 2) - DN_DX(i, 1) * D2Y_De2(0, 2) - DN_DX(i, 2) * D2Z_De2(0, 2);
            }

            noalias(rD2N_DX2[i]) = prod(R2, InvD2);
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        // clean the geometry internal data
        GetGeometry().Clean();
        #endif
    }

} // Namespace Kratos

#undef ENABLE_DEBUG_CONSTITUTIVE_LAW
#undef CHECK_NAN

