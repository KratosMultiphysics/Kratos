//
//   Project Name:        Kratos
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes
#include "includes/define.h"
#include "custom_elements/monolithic_2fluid_2d_deniz_nitsche.h"
#include "pfem_2_application_variables.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "utilities/enrichment_utilities.h"

#include "processes/process.h"
#include "includes/node.h"
//#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"

#include "utilities/divide_triangle_2d_3.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"

namespace Kratos
{

    //************************************************************************************
    //************************************************************************************
    MonolithicPFEM22DDenizNitsche::MonolithicPFEM22DDenizNitsche(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************
    MonolithicPFEM22DDenizNitsche::MonolithicPFEM22DDenizNitsche(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
    }

    Element::Pointer MonolithicPFEM22DDenizNitsche::Create(IndexType NewId, NodesArrayType const &ThisNodes, PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new MonolithicPFEM22DDenizNitsche(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    MonolithicPFEM22DDenizNitsche::~MonolithicPFEM22DDenizNitsche()
    {
    }

    //************************************************************************************
    //************************************************************************************

    //************************************************************************************
    //************************************************************************************

    void MonolithicPFEM22DDenizNitsche::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector, const ProcessInfo &rCurrentProcessInfo)
    {
        KRATOS_TRY

        const int TDim = 2;

        const SizeType NumNodes = TDim + 1;
        const SizeType LocalSize = (TDim + 1) * (TDim + 1);

        //struct to pass around the data
        element_data data;
        double Volume;
        GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, Volume);

        const Vector &BDFVector = rCurrentProcessInfo[BDF_COEFFICIENTS];
        data.bdf0 = BDFVector[0];
        data.bdf1 = BDFVector[1];
        data.bdf2 = BDFVector[2];

        array_1d<double, NumNodes> distances;
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            const array_1d<double, TDim> &vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, TDim> &body_force = GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
            const array_1d<double, TDim> &body_force_nitsche = GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE_NITSCHE);
            const array_1d<double, TDim> &body_force_n = GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE, 1);
            const array_1d<double, TDim> &body_force_nitsche_n = GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE_NITSCHE, 1);
            const array_1d<double, TDim> &vel_n = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1);
            for (unsigned int k = 0; k < TDim; k++)
            {
                data.v(i, k) = vel[k];
                data.vn(i, k) = vel_n[k];
                data.f(i, k) = body_force[k];
                data.f_n(i, k) = body_force_n[k];
                data.f_nitsche(i, k) = body_force_nitsche[k];
                data.f_nitsche_n(i, k) = body_force_nitsche_n[k];
            }
            data.p[i] = GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);
            distances[i] = GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
        }

        unsigned int npos = 0, nneg = 0;
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            if (distances[i] > 0.0)
                npos++;
            if (distances[i] < 0.0)
                nneg++;
        }

        if (npos == NumNodes) //all AIR
        {
            //           KRATOS_WATCH("air");
            //allocate memory needed
            if (rLeftHandSideMatrix.size1() != LocalSize)
                rLeftHandSideMatrix.resize(LocalSize, LocalSize, false); //false says not to preserve existing storage!!
            if (rRightHandSideVector.size() != LocalSize)
                rRightHandSideVector.resize(LocalSize, false); //false says not to preserve existing storage!!

            bounded_matrix<double, LocalSize, LocalSize> lhs_local = ZeroMatrix(LocalSize, LocalSize);
            array_1d<double, LocalSize> rhs_local = ZeroVector(LocalSize);
            this->ComputeElementAsAIR(lhs_local, rhs_local, rLeftHandSideMatrix, rRightHandSideVector, Volume, data, rCurrentProcessInfo);
        }
        else if (nneg == NumNodes) //all FLUID
                                   //         else  //all FLUID
        {
            //                         KRATOS_WATCH("fluid");
            if (rLeftHandSideMatrix.size1() != LocalSize)
                rLeftHandSideMatrix.resize(LocalSize, LocalSize, false); //false says not to preserve existing storage!!
            if (rRightHandSideVector.size() != LocalSize)
                rRightHandSideVector.resize(LocalSize, false); //false says not to preserve existing storage!!

            bounded_matrix<double, LocalSize, LocalSize> lhs_local = ZeroMatrix(LocalSize, LocalSize);
            array_1d<double, LocalSize> rhs_local = ZeroVector(LocalSize);
            this->ComputeElementAsWATER(lhs_local, rhs_local, rLeftHandSideMatrix, rRightHandSideVector, Volume, data, rCurrentProcessInfo);
        }
        else //element includes both FLUID and AIR
        {
            //           KRATOS_WATCH("mixed");
            if (rLeftHandSideMatrix.size1() != LocalSize * 2)
                rLeftHandSideMatrix.resize(LocalSize * 2, LocalSize * 2, false); //false says not to preserve existing storage!!
            if (rRightHandSideVector.size() != LocalSize * 2)
                rRightHandSideVector.resize(LocalSize * 2, false); //false says not to preserve existing storage!!

            bounded_matrix<double, LocalSize * 2, LocalSize * 2> lhs_local = ZeroMatrix(LocalSize * 2, LocalSize * 2);
            array_1d<double, LocalSize * 2> rhs_local = ZeroVector(LocalSize * 2);
            this->ComputeElementAsMIXED(lhs_local, rhs_local, rLeftHandSideMatrix, rRightHandSideVector, Volume, data, rCurrentProcessInfo, distances);
        }

        KRATOS_CATCH("Error in MonolithicPFEM22DDenizNitsche Element");
    }

    //************************************************************************************
    //************************************************************************************

    void MonolithicPFEM22DDenizNitsche::EquationIdVector(EquationIdVectorType &rResult, const ProcessInfo &CurrentProcessInfo) const
    {

        KRATOS_TRY

        const int TDim = 2;

        const SizeType NumNodes = TDim + 1;
        const SizeType LocalSize = (TDim + 1) * (TDim + 1);
        const GeometryType &rGeom = this->GetGeometry();

        array_1d<double, NumNodes> distances;

        for (unsigned int i = 0; i < NumNodes; i++)
        {
            distances[i] = rGeom[i].FastGetSolutionStepValue(DISTANCE);
        }

        unsigned int npos = 0, nneg = 0;
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            if (distances[i] > 0.0)
                npos++;
            if (distances[i] < 0.0)
                nneg++;
        }

        if (npos == NumNodes or nneg == NumNodes) //non-interface element
        {

            if (rResult.size() != LocalSize)
                rResult.resize(LocalSize, false);

            for (unsigned int i = 0; i < NumNodes; i++)
            {
                rResult[i * (TDim + 1)] = rGeom[i].GetDof(VELOCITY_X).EquationId();
                rResult[i * (TDim + 1) + 1] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
                rResult[i * (TDim + 1) + 2] = rGeom[i].GetDof(PRESSURE).EquationId();
            }
        }
        else //interface element
        {

            if (rResult.size() != LocalSize * 2)
                rResult.resize(LocalSize * 2, false);

            for (unsigned int i = 0; i < NumNodes; i++)
            {
                rResult[i * (TDim + 1)] = rGeom[i].GetDof(VELOCITY_X).EquationId();
                rResult[i * (TDim + 1) + 1] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
                rResult[i * (TDim + 1) + 2] = rGeom[i].GetDof(PRESSURE).EquationId();
            }

            for (unsigned int i = 0; i < NumNodes; i++)
            {
                rResult[i * (TDim + 1) + 9] = rGeom[i].GetDof(VELOCITY_NITSCHE_X).EquationId();
                rResult[i * (TDim + 1) + 10] = rGeom[i].GetDof(VELOCITY_NITSCHE_Y).EquationId();
                rResult[i * (TDim + 1) + 11] = rGeom[i].GetDof(PRESSURE_NITSCHE).EquationId();
            }
        }

        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    void MonolithicPFEM22DDenizNitsche::GetDofList(DofsVectorType &ElementalDofList, const ProcessInfo &CurrentProcessInfo) const
    {

        KRATOS_TRY

        const int TDim = 2;

        const SizeType NumNodes = TDim + 1;
        const SizeType LocalSize = (TDim + 1) * (TDim + 1);
        const GeometryType &rGeom = this->GetGeometry();

        array_1d<double, NumNodes> distances;

        for (unsigned int i = 0; i < NumNodes; i++)
        {
            distances[i] = rGeom[i].FastGetSolutionStepValue(DISTANCE);
        }

        unsigned int npos = 0, nneg = 0;
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            if (distances[i] > 0.0)
                npos++;
            if (distances[i] < 0.0)
                nneg++;
        }

        if (npos == NumNodes or nneg == NumNodes) //non-interface element
        {

            //           KRATOS_WATCH("all air or water");

            if (ElementalDofList.size() != LocalSize)
                ElementalDofList.resize(LocalSize);

            for (unsigned int i = 0; i < NumNodes; i++)
            {
                ElementalDofList[i * (TDim + 1)] = rGeom[i].pGetDof(VELOCITY_X);
                ElementalDofList[i * (TDim + 1) + 1] = rGeom[i].pGetDof(VELOCITY_Y);
                ElementalDofList[i * (TDim + 1) + 2] = rGeom[i].pGetDof(PRESSURE);
            }
        }
        else //interface element
        {

            if (ElementalDofList.size() != LocalSize * 2)
                ElementalDofList.resize(LocalSize * 2);

            for (unsigned int i = 0; i < NumNodes; i++)
            {
                ElementalDofList[i * (TDim + 1)] = rGeom[i].pGetDof(VELOCITY_X);
                ElementalDofList[i * (TDim + 1) + 1] = rGeom[i].pGetDof(VELOCITY_Y);
                ElementalDofList[i * (TDim + 1) + 2] = rGeom[i].pGetDof(PRESSURE);
            }
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                ElementalDofList[i * (TDim + 1) + 9] = rGeom[i].pGetDof(VELOCITY_NITSCHE_X);
                ElementalDofList[i * (TDim + 1) + 10] = rGeom[i].pGetDof(VELOCITY_NITSCHE_Y);
                ElementalDofList[i * (TDim + 1) + 11] = rGeom[i].pGetDof(PRESSURE_NITSCHE);
            }
        }

        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************

    void MonolithicPFEM22DDenizNitsche::FinalizeSolutionStep(const ProcessInfo &rCurrentProcessInfo)
    {
        KRATOS_TRY
        const int TDim = 2;
        const SizeType NumNodes = TDim + 1;
        const SizeType LocalSize = (TDim + 1) * (TDim + 1);

        //struct to pass around the data
        element_data data;
        double Volume;
        GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, Volume);

        array_1d<double, (TDim)> utilde_n_aux = ZeroVector(TDim);  
        array_1d<double, (TDim)> utilde_nn_aux = ZeroVector(TDim);                                 //u tilde at the previous time step
        array_1d<double, (TDim * (TDim + 1))> body_force_vector = ZeroVector(TDim * (TDim + 1)); //body force at tn+1 on the gauss point
        boost::numeric::ublas::bounded_matrix<double, (TDim), (TDim * (TDim + 1))> dummyN = ZeroMatrix((TDim), TDim * (TDim + 1));
        array_1d<double, TDim *(TDim + 1)> v = ZeroVector(TDim * (TDim + 1));   //v vector
        array_1d<double, TDim *(TDim + 1)> v_n = ZeroVector(TDim * (TDim + 1)); //v_n vector
        array_1d<double, (TDim + 1)> p = ZeroVector((TDim + 1)); 
        array_1d<double, TDim *(TDim + 1)> arhs = ZeroVector(TDim * (TDim + 1)); //arhs vector_bdf1*vn+bdf2*vnn
        array_1d<double, TDim> arhs2 = ZeroVector(TDim);
        boost::numeric::ublas::bounded_matrix<double, (TDim + 1), TDim *(TDim + 1)> Kpv_ = ZeroMatrix(TDim + 1, TDim * (TDim + 1));   
        boost::numeric::ublas::bounded_matrix<double, (TDim + 1), (TDim + 1)> Kpp_ = ZeroMatrix(TDim + 1, TDim + 1);   


        array_1d<double, (TDim + 1)> rp_ = ZeroVector((TDim + 1));

        const double mu = GetProperties()[VISCOSITY_WATER];
        const double rho = GetProperties()[DENSITY_WATER];
        const Vector &BDFVector = rCurrentProcessInfo[BDF_COEFFICIENTS];
        double bdf0 = BDFVector[0];
        double bdf1 = BDFVector[1];
        double bdf2 = BDFVector[2];
        double dt = rCurrentProcessInfo[DELTA_TIME];

        //get shape function values
        const bounded_matrix<double, NumNodes, TDim> &DN = data.DN_DX;
        const array_1d<double, NumNodes> &N = data.N;
        for (unsigned int i = 0; i < (TDim + 1); i++)
        {
            dummyN(0, i * TDim) = N(i);
            dummyN(1, i * TDim + 1) = N(i);
            if (TDim == 3)
            {
                dummyN(2, i * TDim + 2) = N(i);
            }
        }

        for (unsigned int i = 0; i < (TDim + 1); i++)
        {
            v[TDim * i] = (GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X));
            v[TDim * i + 1] = (GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y));
            v_n[TDim * i] = (GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X, 1));
            v_n[TDim * i + 1] = (GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y, 1));
            p(i) = (GetGeometry()[i].FastGetSolutionStepValue(PRESSURE));
            arhs[TDim * i] = bdf1 * v_n(TDim * i);
            arhs[TDim * i + 1] = bdf1 * v_n(TDim * i + 1);
            arhs2[0] += N[i] * arhs[TDim * i];
            arhs2[1] += N[i] * arhs[TDim * i + 1];
            // body_force_vector[TDim * i] = f(i, 0);
            // body_force_vector[TDim * i + 1] = f(i, 1);
            // body_force_vector_n[TDim * i] = f_n(i, 0);
            // body_force_vector_n[TDim * i + 1] = f_n(i, 1);
        }

        //Tau for stabilization
        double sqElemSize = 0.0;
        double mElemSize;
        array_1d<double, 3> Edge(3, 0.0);
        Edge = this->GetGeometry()[1].Coordinates() - this->GetGeometry()[0].Coordinates();
        sqElemSize = Edge[0] * Edge[0];
        for (SizeType d = 1; d < TDim; d++)
            sqElemSize += Edge[d] * Edge[d];

        for (SizeType i = 2; i < TDim + 1; i++)
        {
            for (SizeType j = 0; j < i; j++)
            {
                Edge = this->GetGeometry()[i].Coordinates() - this->GetGeometry()[j].Coordinates();
                double Length = Edge[0] * Edge[0];
                for (SizeType d = 1; d < TDim; d++)
                    Length += Edge[d] * Edge[d];
                if (Length < sqElemSize)
                    sqElemSize = Length;
            }
        }

        mElemSize = sqrt(sqElemSize);
        double tau = 0.0;
        if (bdf0 != 0.0)
        {
            tau = 1.0 / ((1.0 * rho / (dt) + 4.0 * mu / (mElemSize * mElemSize)));
        }
        else
        {
            tau = 1.0 / ((4.0 * mu / (mElemSize * mElemSize)));
        }



        //Part-10_
        utilde_n_aux += -tau*rho*bdf0*prod(dummyN,v);
        // KRATOS_WATCH(mElemSize);
        // KRATOS_WATCH(tau);
        // KRATOS_WATCH(rho);
        // KRATOS_WATCH(bdf0);
        // KRATOS_WATCH(dummyN);
        // KRATOS_WATCH(v);
        // KRATOS_WATCH(utilde_n_aux);
        // KRATOS_WATCH(-tau*rho*bdf0*prod(dummyN,v));
        utilde_n_aux += -tau*rho*arhs2;

        

        //Part-07_
        utilde_n_aux += -tau*prod(trans(DN),p);

        //Part-12_
        utilde_nn_aux = this->utilde_n;
        utilde_n_aux += rho * tau * (1.0 / dt) * utilde_nn_aux;

        this->utilde_n = utilde_n_aux;

        // KRATOS_WATCH(utilde_n_aux);



        KRATOS_CATCH("");
    }

    template <class T>
    bool MonolithicPFEM22DDenizNitsche::InvertMatrix(const T &input, T &inverse)
    {
        typedef permutation_matrix<std::size_t> pmatrix;

        // create a working copy of the input
        T A(input);

        // create a permutation matrix for the LU-factorization
        pmatrix pm(A.size1());

        // perform LU-factorization
        int res = lu_factorize(A, pm);
        if (res != 0)
            return false;

        // create identity matrix of "inverse"
        inverse.assign(identity_matrix<double>(A.size1()));

        // backsubstitute to get the inverse
        lu_substitute(A, pm, inverse);

        return true;
    }

    void MonolithicPFEM22DDenizNitsche::ComputeElementAsWATER(bounded_matrix<double, 9, 9> &lhs_local,
                                                              array_1d<double, 9> &rhs_local,
                                                              Matrix &rLeftHandSideMatrix,
                                                              Vector &rRightHandSideVector,
                                                              const double &Volume,
                                                              element_data &data,
                                                              const ProcessInfo &rCurrentProcessInfo)
    {

        const int TDim = 2;
        const SizeType NumNodes = TDim + 1;

        const double fluid_density = GetProperties()[DENSITY_WATER];
        const double fluid_mu = GetProperties()[VISCOSITY_WATER];
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            data.rho[i] = fluid_density;
        }

        const double weight = Volume;
        noalias(rLeftHandSideMatrix) = ZeroMatrix(9, 9);
        noalias(rRightHandSideVector) = ZeroVector(9);

        for (unsigned int igauss = 0; igauss < 1; igauss++)
        {
            //              noalias(data.N) = row(Ncontainer, igauss);

            MonolithicPFEM22DDenizNitsche::ComputeConstitutiveResponse(data, fluid_density, fluid_mu, rCurrentProcessInfo);
            MonolithicPFEM22DDenizNitsche::ComputeGaussPointLHSandRHSContribution(lhs_local, rhs_local, data, rCurrentProcessInfo, weight);
            //here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/NumNodes
            noalias(rLeftHandSideMatrix) += weight * lhs_local;
            noalias(rRightHandSideVector) += weight * rhs_local;
        }

        array_1d<double, TDim *(TDim + 1)> v_aux = ZeroVector(TDim * (TDim + 1));
        array_1d<double, (TDim + 1)> p_aux = ZeroVector(TDim + 1);
        array_1d<double, (TDim + 1) * (TDim + 1)> temp_vector;

        v_aux[0] = (GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X));
        v_aux[1] = (GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y));
        v_aux[2] = (GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_X));
        v_aux[3] = (GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_Y));
        v_aux[4] = (GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_X));
        v_aux[5] = (GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_Y));
        p_aux[0] = (GetGeometry()[0].FastGetSolutionStepValue(PRESSURE));
        p_aux[1] = (GetGeometry()[1].FastGetSolutionStepValue(PRESSURE));
        p_aux[2] = (GetGeometry()[2].FastGetSolutionStepValue(PRESSURE));
        int counter = 0;
        for (unsigned int iii = 0; iii < (TDim + 1); iii++)
        {
            temp_vector[counter++] = v_aux[TDim * iii];
            temp_vector[counter++] = v_aux[TDim * iii + 1];
            temp_vector[counter++] = p_aux[iii];
        }
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, temp_vector);
    }

    void MonolithicPFEM22DDenizNitsche::ComputeElementAsAIR(bounded_matrix<double, 9, 9> &lhs_local,
                                                            array_1d<double, 9> &rhs_local,
                                                            Matrix &rLeftHandSideMatrix,
                                                            Vector &rRightHandSideVector,
                                                            const double &Volume,
                                                            element_data &data,
                                                            const ProcessInfo &rCurrentProcessInfo)
    {

        const int TDim = 2;
        const SizeType NumNodes = TDim + 1;

        const double air_density = GetProperties()[DENSITY_AIR];
        const double air_mu = GetProperties()[VISCOSITY_AIR];
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            data.rho[i] = air_density;
        }

        const double weight = Volume;
        noalias(rLeftHandSideMatrix) = ZeroMatrix(9, 9);
        noalias(rRightHandSideVector) = ZeroVector(9);

        for (unsigned int igauss = 0; igauss < 1; igauss++)
        {
            //              noalias(data.N) = row(Ncontainer, igauss);

            MonolithicPFEM22DDenizNitsche::ComputeConstitutiveResponse(data, air_density, air_mu, rCurrentProcessInfo);
            MonolithicPFEM22DDenizNitsche::ComputeGaussPointLHSandRHSContribution(lhs_local, rhs_local, data, rCurrentProcessInfo, weight);
            //here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/NumNodes
            noalias(rLeftHandSideMatrix) += weight * lhs_local;
            noalias(rRightHandSideVector) += weight * rhs_local;
        }

        array_1d<double, TDim *(TDim + 1)> v_aux = ZeroVector(TDim * (TDim + 1));
        array_1d<double, (TDim + 1)> p_aux = ZeroVector(TDim + 1);
        array_1d<double, (TDim + 1) * (TDim + 1)> temp_vector;

        v_aux[0] = (GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X));
        v_aux[1] = (GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y));
        v_aux[2] = (GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_X));
        v_aux[3] = (GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_Y));
        v_aux[4] = (GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_X));
        v_aux[5] = (GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_Y));
        p_aux[0] = (GetGeometry()[0].FastGetSolutionStepValue(PRESSURE));
        p_aux[1] = (GetGeometry()[1].FastGetSolutionStepValue(PRESSURE));
        p_aux[2] = (GetGeometry()[2].FastGetSolutionStepValue(PRESSURE));
        int counter = 0;
        for (unsigned int iii = 0; iii < (TDim + 1); iii++)
        {
            temp_vector[counter++] = v_aux[TDim * iii];
            temp_vector[counter++] = v_aux[TDim * iii + 1];
            temp_vector[counter++] = p_aux[iii];
        }
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, temp_vector);
    }

    void MonolithicPFEM22DDenizNitsche::ComputeElementAsMIXED(bounded_matrix<double, 18, 18> &lhs_local, array_1d<double, 18> &rhs_local,
                                                              Matrix &rLeftHandSideMatrix,
                                                              Vector &rRightHandSideVector,
                                                              const double &Volume,
                                                              element_data &data,
                                                              const ProcessInfo &rCurrentProcessInfo,
                                                              array_1d<double, 3> &distances)

    {

        const int TDim = 2;
        const SizeType NumNodes = TDim + 1;

        const GeometryType &geom = this->GetGeometry();
        double weight = 0.0;
        noalias(rLeftHandSideMatrix) = ZeroMatrix(18, 18);
        noalias(rRightHandSideVector) = ZeroVector(18);

        // Set the elemental distances vector
        Geometry<Node<3>>::Pointer p_geometry = this->pGetGeometry();

        array_1d<double, 3> distances_vector;
        for (unsigned int i = 0; i < p_geometry->size(); ++i)
        {
            distances_vector(i) = (*p_geometry)[i].FastGetSolutionStepValue(DISTANCE);
        }

        //             KRATOS_WATCH(distances_vector);
        this->SetValue(ELEMENTAL_DISTANCES, distances_vector);

        const Vector &r_elemental_distances = this->GetValue(ELEMENTAL_DISTANCES);

        const double air_density = GetProperties()[DENSITY_AIR];
        const double air_mu = GetProperties()[VISCOSITY_AIR];
        const double water_density = GetProperties()[DENSITY_WATER];
        const double water_mu = GetProperties()[VISCOSITY_WATER];
        ComputeConstitutiveResponse(data, air_density, air_mu, rCurrentProcessInfo);
        const Matrix &C_air = data.C;
        ComputeConstitutiveResponse(data, water_density, water_mu, rCurrentProcessInfo);
        const Matrix &C_water = data.C;

        // Call the modified shape functions calculator
        Triangle2D3ModifiedShapeFunctions triangle_shape_functions(p_geometry, r_elemental_distances);
        Matrix positive_side_sh_func, negative_side_sh_func;
        ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_gradients, negative_side_sh_func_gradients;
        Vector positive_side_weights, negative_side_weights;

        triangle_shape_functions.ComputePositiveSideShapeFunctionsAndGradientsValues(
            positive_side_sh_func,
            positive_side_sh_func_gradients,
            positive_side_weights,
            GeometryData::GI_GAUSS_1);

        triangle_shape_functions.ComputeNegativeSideShapeFunctionsAndGradientsValues(
            negative_side_sh_func,
            negative_side_sh_func_gradients,
            negative_side_weights,
            GeometryData::GI_GAUSS_1);

        // Call the interface modified shape functions calculator
        Matrix positive_interface_side_sh_func, negative_interface_side_sh_func;
        ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_interface_side_sh_func_gradients, negative_interface_side_sh_func_gradients;
        Vector positive_interface_side_weights, negative_interface_side_weights;

        triangle_shape_functions.ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
            positive_interface_side_sh_func,
            positive_interface_side_sh_func_gradients,
            positive_interface_side_weights,
            GeometryData::GI_GAUSS_2);

        triangle_shape_functions.ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
            negative_interface_side_sh_func,
            negative_interface_side_sh_func_gradients,
            negative_interface_side_weights,
            GeometryData::GI_GAUSS_2);

        // Call the interface outwards normal area vector calculator
        std::vector<Vector> positive_side_area_normals, negative_side_area_normals;

        triangle_shape_functions.ComputePositiveSideInterfaceAreaNormals(
            positive_side_area_normals,
            GeometryData::GI_GAUSS_2);

        triangle_shape_functions.ComputeNegativeSideInterfaceAreaNormals(
            negative_side_area_normals,
            GeometryData::GI_GAUSS_2);

        //Calculation of the enlargement matrix T
        boost::numeric::ublas::bounded_matrix<double, (NumNodes * NumNodes * 2), (NumNodes * NumNodes * 2)> Tdublicate_positive = ZeroMatrix((NumNodes * NumNodes * 2), (NumNodes * NumNodes * 2));
        boost::numeric::ublas::bounded_matrix<double, (NumNodes * NumNodes * 2), (NumNodes * NumNodes * 2)> Tdublicate_negative = ZeroMatrix((NumNodes * NumNodes * 2), (NumNodes * NumNodes * 2));
        boost::numeric::ublas::bounded_matrix<double, (NumNodes * 2), (NumNodes * 2)> Tpositive = ZeroMatrix((NumNodes * 2), (NumNodes * 2));
        boost::numeric::ublas::bounded_matrix<double, (NumNodes * 2), (NumNodes * 2)> Tnegative = ZeroMatrix((NumNodes * 2), (NumNodes * 2));
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            if (distances[i] > 0.0)
            {
                Tdublicate_positive(i * NumNodes, i * NumNodes) = 1.0;
                Tdublicate_positive(i * NumNodes + 1, i * NumNodes + 1) = 1.0;
                Tdublicate_positive(i * NumNodes + 2, i * NumNodes + 2) = 1.0;
                Tpositive(i, i) = 1.0;

                Tdublicate_negative(i * NumNodes + 9, i * NumNodes) = 1.0;
                Tdublicate_negative(i * NumNodes + 10, i * NumNodes + 1) = 1.0;
                Tdublicate_negative(i * NumNodes + 11, i * NumNodes + 2) = 1.0;
                Tnegative(i + 3, i) = 1.0;
            }
            if (distances[i] < 0.0)
            {
                Tdublicate_negative(i * NumNodes, i * NumNodes) = 1.0;
                Tdublicate_negative(i * NumNodes + 1, i * NumNodes + 1) = 1.0;
                Tdublicate_negative(i * NumNodes + 2, i * NumNodes + 2) = 1.0;
                Tnegative(i, i) = 1.0;

                Tdublicate_positive(i * NumNodes + 9, i * NumNodes) = 1.0;
                Tdublicate_positive(i * NumNodes + 10, i * NumNodes + 1) = 1.0;
                Tdublicate_positive(i * NumNodes + 11, i * NumNodes + 2) = 1.0;
                Tpositive(i + 3, i) = 1.0;
            }
        }

        //First we start with the positive CalculateLeftandRightHandSide
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            data.rho[i] = air_density;
        }

        //First we calculate the domain terms
        for (unsigned int igauss = 0; igauss < positive_side_sh_func.size1(); igauss++)
        {
            ComputeConstitutiveResponse(data, air_density, air_mu, rCurrentProcessInfo);
            noalias(data.N) = row(positive_side_sh_func, igauss);
            noalias(data.DN_DX) = positive_side_sh_func_gradients[igauss];
            weight = positive_side_weights[igauss];

            MonolithicPFEM22DDenizNitsche::ComputeGaussPointLHSandRHSContribution_NitscheDomainTerms(lhs_local, rhs_local, data, rCurrentProcessInfo, weight, Tdublicate_positive);

            noalias(rLeftHandSideMatrix) += weight * lhs_local;
            noalias(rRightHandSideVector) += weight * rhs_local;
        }

        //Then we calculate the Nitsche boundary terms
        for (unsigned int igauss = 0; igauss < positive_interface_side_sh_func.size1(); igauss++)
        {
            ComputeConstitutiveResponse(data, air_density, air_mu, rCurrentProcessInfo);
            noalias(data.N) = row(positive_interface_side_sh_func, igauss);
            noalias(data.DN_DX) = positive_interface_side_sh_func_gradients[igauss];
            weight = positive_interface_side_weights[igauss];

            array_1d<double, 3> normal = ZeroVector(3);
            normal = positive_side_area_normals[igauss];
            double norm = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
            norm = sqrt(norm);
            normal /= norm;

            MonolithicPFEM22DDenizNitsche::ComputeGaussPointLHSandRHSContribution_NitscheBoundaryTerms(lhs_local, rhs_local, data, rCurrentProcessInfo, weight, normal, Tdublicate_positive, Tdublicate_negative, Tpositive, Tnegative, C_air, C_water);

            noalias(rLeftHandSideMatrix) += weight * lhs_local;
            noalias(rRightHandSideVector) += weight * rhs_local;
        }

        //Second we start with the negative CalculateLeftandRightHandSide
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            data.rho[i] = water_density;
        }

        //First we calculate the domain terms
        for (unsigned int igauss = 0; igauss < negative_side_sh_func.size1(); igauss++)
        {
            ComputeConstitutiveResponse(data, water_density, water_mu, rCurrentProcessInfo);
            noalias(data.N) = row(negative_side_sh_func, igauss);
            noalias(data.DN_DX) = negative_side_sh_func_gradients[igauss];
            weight = negative_side_weights[igauss];

            MonolithicPFEM22DDenizNitsche::ComputeGaussPointLHSandRHSContribution_NitscheDomainTerms(lhs_local, rhs_local, data, rCurrentProcessInfo, weight, Tdublicate_negative);

            noalias(rLeftHandSideMatrix) += weight * lhs_local;
            noalias(rRightHandSideVector) += weight * rhs_local;
        }

        //Then we calculate the Nitsche boundary terms
        for (unsigned int igauss = 0; igauss < negative_interface_side_sh_func.size1(); igauss++)
        {
            ComputeConstitutiveResponse(data, water_density, water_mu, rCurrentProcessInfo);
            noalias(data.N) = row(negative_interface_side_sh_func, igauss);
            noalias(data.DN_DX) = negative_interface_side_sh_func_gradients[igauss];
            weight = negative_interface_side_weights[igauss];

            array_1d<double, 3> normal = ZeroVector(3);
            normal = negative_side_area_normals[igauss];
            double norm = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
            norm = sqrt(norm);
            normal /= norm;

            MonolithicPFEM22DDenizNitsche::ComputeGaussPointLHSandRHSContribution_NitscheBoundaryTerms(lhs_local, rhs_local, data, rCurrentProcessInfo, weight, normal, Tdublicate_negative, Tdublicate_positive, Tnegative, Tpositive, C_water, C_air);

            noalias(rLeftHandSideMatrix) += weight * lhs_local;
            noalias(rRightHandSideVector) += weight * rhs_local;
        }

        array_1d<double, TDim *(TDim + 1)> v_aux = ZeroVector(TDim * (TDim + 1));
        array_1d<double, TDim *(TDim + 1)> v_aux_nitsche = ZeroVector(TDim * (TDim + 1));
        array_1d<double, (TDim + 1)> p_aux = ZeroVector(TDim + 1);
        array_1d<double, (TDim + 1)> p_aux_nitsche = ZeroVector(TDim + 1);
        array_1d<double, ((TDim + 1) * (TDim + 1) * 2)> temp_vector;

        v_aux[0] = (GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X));
        v_aux[1] = (GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y));
        v_aux[2] = (GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_X));
        v_aux[3] = (GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_Y));
        v_aux[4] = (GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_X));
        v_aux[5] = (GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_Y));
        p_aux[0] = (GetGeometry()[0].FastGetSolutionStepValue(PRESSURE));
        p_aux[1] = (GetGeometry()[1].FastGetSolutionStepValue(PRESSURE));
        p_aux[2] = (GetGeometry()[2].FastGetSolutionStepValue(PRESSURE));

        v_aux_nitsche[0] = (GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_NITSCHE_X));
        v_aux_nitsche[1] = (GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_NITSCHE_Y));
        v_aux_nitsche[2] = (GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_NITSCHE_X));
        v_aux_nitsche[3] = (GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_NITSCHE_Y));
        v_aux_nitsche[4] = (GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_NITSCHE_X));
        v_aux_nitsche[5] = (GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_NITSCHE_Y));
        p_aux_nitsche[0] = (GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_NITSCHE));
        p_aux_nitsche[1] = (GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_NITSCHE));
        p_aux_nitsche[2] = (GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_NITSCHE));

        int counter = 0;
        for (unsigned int iii = 0; iii < (TDim + 1); iii++)
        {
            temp_vector[counter++] = v_aux[TDim * iii];
            temp_vector[counter++] = v_aux[TDim * iii + 1];
            temp_vector[counter++] = p_aux[iii];
        }

        for (unsigned int iii = 0; iii < (TDim + 1); iii++)
        {
            temp_vector[counter++] = v_aux_nitsche[TDim * iii];
            temp_vector[counter++] = v_aux_nitsche[TDim * iii + 1];
            temp_vector[counter++] = p_aux_nitsche[iii];
        }
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, temp_vector);
    }

    void MonolithicPFEM22DDenizNitsche::ComputeConstitutiveResponse(element_data &data, const double rho, const double mu, const ProcessInfo &rCurrentProcessInfo) //Deniz-changed nu to mu
    {
        const unsigned int NumNodes = 3;
        const unsigned int dim = 2;
        const unsigned int strain_size = 3;

        if (data.C.size1() != strain_size)
            data.C.resize(strain_size, strain_size, false);

        if (data.stress.size() != strain_size)
            data.stress.resize(strain_size, false);

        const bounded_matrix<double, NumNodes, dim> &v = data.v;
        const bounded_matrix<double, NumNodes, dim> &DN = data.DN_DX;

        //compute strain
        //         Vector strain(strain_size);
        //         strain[0] = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
        //         strain[1] = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
        //         strain[2] = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
        //         strain[3] = DN(0,0)*v(0,1) + DN(0,1)*v(0,0) + DN(1,0)*v(1,1) + DN(1,1)*v(1,0) + DN(2,0)*v(2,1) + DN(2,1)*v(2,0) + DN(3,0)*v(3,1) + DN(3,1)*v(3,0);
        //         strain[4] = DN(0,1)*v(0,2) + DN(0,2)*v(0,1) + DN(1,1)*v(1,2) + DN(1,2)*v(1,1) + DN(2,1)*v(2,2) + DN(2,2)*v(2,1) + DN(3,1)*v(3,2) + DN(3,2)*v(3,1);
        //         strain[5] = DN(0,0)*v(0,2) + DN(0,2)*v(0,0) + DN(1,0)*v(1,2) + DN(1,2)*v(1,0) + DN(2,0)*v(2,2) + DN(2,2)*v(2,0) + DN(3,0)*v(3,2) + DN(3,2)*v(3,0);

        //here we shall call the constitutive law
        data.C.clear();

        //         data.C(0,0) = 2.0*mu;    //Deniz-changed nu to mu
        //         data.C(1,1) = 2.0*mu;    //Deniz-changed nu to mu
        //         data.C(2,2) = 2.0*mu;    //Deniz-changed nu to mu
        //         data.C(3,3) = mu;        //Deniz-changed nu to mu
        //         data.C(4,4) = mu;        //Deniz-changed nu to mu
        //         data.C(5,5) = mu;        //Deniz-changed nu to mu
        //         data.C(0,0) = 4.0/3.0*mu;  data.C(0,1) = -2.0/3.0*mu; data.C(1,0) = -2.0/3.0*mu;
        //         data.C(1,1) = 4.0/3.0*mu;
        //         data.C(2,2) = mu;

        data.C(0, 0) = 2.0 * mu;
        data.C(0, 1) = 0.0;
        data.C(0, 2) = 0.0;
        data.C(1, 1) = 2.0 * mu;
        data.C(1, 0) = 0.0;
        data.C(1, 2) = 0.0;
        data.C(2, 2) = mu;
        data.C(2, 0) = 0.0;
        data.C(2, 1) = 0.0;

        //         const double c2 = mu; //Deniz-changed nu to mu
        //         const double c1 = 2.0*c2;
        //         data.stress[0] =  c1*strain[0];
        //         data.stress[1] =  c1*strain[1];
        //         data.stress[2] =  c1*strain[2];
        //         data.stress[3] =  c2*strain[3];
        //         data.stress[4] =  c2*strain[4];
        //         data.stress[5] =  c2*strain[5];
    }

    void MonolithicPFEM22DDenizNitsche::ComputeGaussPointLHSandRHSContribution(bounded_matrix<double, 9, 9> &lhs, array_1d<double, 9> &rhs, const element_data &data, const ProcessInfo &rCurrentProcessInfo, const double &weight)
    {

        const unsigned int NumNodes = 3;
        const unsigned int TDim = 2;
        lhs = ZeroMatrix(9, 9);
        rhs = ZeroVector(9);

        boost::numeric::ublas::bounded_matrix<double, (TDim + 1) * (TDim) / 2, TDim *(TDim + 1)> B = ZeroMatrix((TDim + 1) * (TDim) / 2, TDim * (TDim + 1));       //B matrix
        boost::numeric::ublas::bounded_matrix<double, (TDim + 1) * (TDim) / 2, TDim *(TDim + 1)> Prod_CB = ZeroMatrix((TDim + 1) * (TDim) / 2, TDim * (TDim + 1)); //Prod_CB matrix
        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1), TDim *(TDim + 1)> Prod_BtCB = ZeroMatrix(TDim * (TDim + 1), TDim * (TDim + 1));            //Prod_BtCB matrix
        boost::numeric::ublas::bounded_matrix<double, TDim, (TDim + 1)> Prod_TrDNtr = ZeroMatrix(TDim, (TDim + 1));
        boost::numeric::ublas::bounded_matrix<double, (TDim), (TDim * (TDim + 1))> dummyN = ZeroMatrix((TDim), TDim * (TDim + 1));

        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1), (TDim + 1)> Kvp_aux = ZeroMatrix(TDim * (TDim + 1), (TDim + 1)); //Kvp_aux matrix
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> Tau = ZeroMatrix(TDim, TDim);                                          //Stabilization Tau vector
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim *(TDim + 1)> Prod_TaudummyN = ZeroMatrix(TDim, TDim * (TDim + 1));      //Prod_TaudummyN
        boost::numeric::ublas::bounded_matrix<double, (TDim + 1), TDim> Prod_DNDX_Tau = ZeroMatrix(TDim + 1, TDim);
        array_1d<double, TDim> tauf = ZeroVector(TDim);                          //An auxilary stabilization vector
        array_1d<double, TDim *(TDim + 1)> arhs = ZeroVector(TDim * (TDim + 1)); //arhs vector_bdf1*vn+bdf2*vnn
        array_1d<double, TDim> arhs2 = ZeroVector(TDim);
        array_1d<double, (TDim)> body_force_on_gauss_point = ZeroVector(TDim);                     //body force at tn+1 on the gauss point
        array_1d<double, (TDim)> body_force_on_gauss_point_n = ZeroVector(TDim);                   //body force at tn on the gauss point
        array_1d<double, (TDim)> a_n_on_gauss_point = ZeroVector(TDim); 
        array_1d<double, (TDim)> utilde_n = ZeroVector(TDim);                                      //u tilde at the previous time step
        array_1d<double, (TDim * (TDim + 1))> body_force_vector = ZeroVector(TDim * (TDim + 1));   //body force at tn+1 on the gauss point
        array_1d<double, (TDim * (TDim + 1))> body_force_vector_n = ZeroVector(TDim * (TDim + 1)); //body force at tn+1 on the gauss point
        array_1d<double, TDim *(TDim + 1)> v_n = ZeroVector(TDim * (TDim + 1));                    //v_n vector
        array_1d<double, (TDim + 1)> p_n = ZeroVector((TDim + 1));                                 //p_n vector
        array_1d<double, TDim *(TDim + 1)> a_n = ZeroVector(TDim * (TDim + 1));                    //a_n vector

        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1), TDim *(TDim + 1)> Kvv = ZeroMatrix(TDim * (TDim + 1), TDim * (TDim + 1));     //Kvv matrix
        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1), TDim *(TDim + 1)> Kvv_aux = ZeroMatrix(TDim * (TDim + 1), TDim * (TDim + 1)); //Kvv_aux matrix
        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1), (TDim + 1)> Kvp = ZeroMatrix(TDim * (TDim + 1), (TDim + 1));                  //Kvp matrix
        boost::numeric::ublas::bounded_matrix<double, (TDim + 1), TDim *(TDim + 1)> Kpv = ZeroMatrix(TDim + 1, TDim * (TDim + 1));                    //Kpv matrix
        boost::numeric::ublas::bounded_matrix<double, (TDim + 1), (TDim + 1)> Kpp = ZeroMatrix(TDim + 1, TDim + 1);                                   //Kpp matrix
        array_1d<double, TDim *(TDim + 1)> rv = ZeroVector(TDim * (TDim + 1));
        array_1d<double, (TDim + 1)> rp = ZeroVector((TDim + 1));
        array_1d<double, (TDim * (TDim + 1))> div_vec = ZeroVector(TDim * (TDim + 1)); //div_vec vector

        const bool momentumstabilize = false;
        const bool dynamicsubgridstabilize = false;
        const bool acceleratedparticles = false;
        //get constitutive matrix
        const Matrix &C = data.C;
        //get density
        const double rho = inner_prod(data.N, data.rho);
        //get bdf coefficients and dt
        const double &bdf0 = data.bdf0;
        const double &bdf1 = data.bdf1;
        const double &bdf2 = data.bdf2;
        double dt = rCurrentProcessInfo[DELTA_TIME];
        const double theta = GetProperties()[TIMEDIS_THETA];

        //      KRATOS_WATCH(dt);
        //      KRATOS_WATCH(bdf0);
        //      KRATOS_WATCH(bdf1);
        //      KRATOS_WATCH(bdf2);

        //get shape function values
        const bounded_matrix<double, NumNodes, TDim> &DN = data.DN_DX;
        const array_1d<double, NumNodes> &N = data.N;
        for (unsigned int i = 0; i < (TDim + 1); i++)
        {
            dummyN(0, i * TDim) = N(i);
            dummyN(1, i * TDim + 1) = N(i);
            if (TDim == 3)
            {
                dummyN(2, i * TDim + 2) = N(i);
            }
        }
        //calculate body force on gauss points for tn+1 and tn
        const bounded_matrix<double, NumNodes, TDim> &f = data.f;
        const bounded_matrix<double, NumNodes, TDim> &f_n = data.f_n;
        body_force_on_gauss_point = prod(trans(f), N);
        body_force_on_gauss_point_n = prod(trans(f_n), N);

        for (unsigned int i = 0; i < (TDim + 1); i++)
        {
            v_n[TDim * i] = (GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X, 1));
            v_n[TDim * i + 1] = (GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y, 1));
            a_n[TDim * i] = (GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION_AUX_X, 1));
            a_n[TDim * i + 1] = (GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION_AUX_Y, 1));
            arhs[TDim * i] = bdf1 * v_n(TDim * i);
            arhs[TDim * i + 1] = bdf1 * v_n(TDim * i + 1);
            arhs2[0] += N[i] * arhs[TDim * i];
            arhs2[1] += N[i] * arhs[TDim * i + 1];
            a_n_on_gauss_point[0] += N[i] * a_n[TDim * i];
            a_n_on_gauss_point[1] += N[i] * a_n[TDim * i + 1];
            body_force_vector[TDim * i] = f(i, 0);
            body_force_vector[TDim * i + 1] = f(i, 1);
            body_force_vector_n[TDim * i] = f_n(i, 0);
            body_force_vector_n[TDim * i + 1] = f_n(i, 1);
            if (TDim == 3)
            {
                v_n[TDim * i + 2] = (GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Z, 1));
                a_n[TDim * i + 3] = (GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION_AUX_Z, 1));
                body_force_vector[TDim * i + 2] = f(i, 2);
                body_force_vector_n[TDim * i + 2] = f_n(i, 2);
                arhs[TDim * i + 2] = bdf1 * v_n(TDim * i + 2);
                arhs2[2] += N[i] * arhs[TDim * i + 2];
            }
            p_n(i) = (GetGeometry()[i].FastGetSolutionStepValue(PRESSURE, 1));
        }

        //      KRATOS_WATCH(v_n);
        //      KRATOS_WATCH(p_n);
        //      KRATOS_WATCH(arhs);
        //      KRATOS_WATCH(arhs2);

        //Tau for stabilization
        double sqElemSize = 0.0;
        double mElemSize;
        array_1d<double, 3> Edge(3, 0.0);
        Edge = this->GetGeometry()[1].Coordinates() - this->GetGeometry()[0].Coordinates();
        sqElemSize = Edge[0] * Edge[0];
        for (SizeType d = 1; d < TDim; d++)
            sqElemSize += Edge[d] * Edge[d];

        for (SizeType i = 2; i < TDim + 1; i++)
        {
            for (SizeType j = 0; j < i; j++)
            {
                Edge = this->GetGeometry()[i].Coordinates() - this->GetGeometry()[j].Coordinates();
                double Length = Edge[0] * Edge[0];
                for (SizeType d = 1; d < TDim; d++)
                    Length += Edge[d] * Edge[d];
                if (Length < sqElemSize)
                    sqElemSize = Length;
            }
        }

        mElemSize = sqrt(sqElemSize);
        double tau = 0.0;
        const double tau2 = C(2, 2);
        if (bdf0 != 0.0)
        {
            tau = 1.0 / ((1.0 * rho / (dt) + 4.0 * C(2, 2) / (mElemSize * mElemSize)));
        }
        else
        {
            tau = 1.0 / ((4.0 * C(2, 2) / (mElemSize * mElemSize)));
        }


        for (unsigned int i = 0; i < (TDim); i++)
        {
            Tau(i, i) = tau;
        }

        //Part-01
        for (unsigned int i = 0; i < (TDim + 1); i++)
        {
            double tmp = (rho * bdf0) * N(i);
            Kvv(i * TDim, i * TDim) += tmp;
            Kvv(i * TDim + 1, i * TDim + 1) += tmp;
        }
        for (unsigned int i = 0; i < (TDim + 1); i++)
        {
            rv(i * TDim) -= (rho * N(i)) * arhs(i * TDim);
            rv(i * TDim + 1) -= (rho * N(i)) * arhs(i * TDim + 1);
        }

        //Part-02_ThetaMethodsIntegrated
        for (unsigned int i = 0; i < (TDim + 1); i++)
        {
            unsigned int index = i * TDim;
            B(0, index + 0) = DN(i, 0);
            B(0, index + 1) = 0;
            B(1, index + 0) = 0;
            B(1, index + 1) = DN(i, 1);
            B(2, index + 0) = DN(i, 1);
            B(2, index + 1) = DN(i, 0);
        }
        noalias(Prod_CB) = prod(C, B);
        noalias(Prod_BtCB) = prod(trans(B), Prod_CB);
        Kvv += theta * Prod_BtCB;
        rv -= (1.0 - theta) * prod(Prod_BtCB, v_n);

        if (momentumstabilize == true)
        {
            //Part-06_ThetaMethodsIntegrated
            rv -= (theta * rho * tau * prod(Prod_BtCB, body_force_vector) + (1.0 - theta) * rho * tau * prod(Prod_BtCB, body_force_vector_n));

            //Part-09_ThetaMethodsIntegrated
            for (int i = 0; i < (TDim + 1); i++)
            {
                for (int k = 0; k < (TDim); k++)
                {
                    div_vec(TDim * i + k) = DN(i, k);
                }
            }
            Kvv_aux += outer_prod(div_vec, div_vec);
            Kvv -= theta * tau2 * Kvv_aux;
            rv += (1.0 - theta) * tau2 * prod(Kvv_aux, v_n);
            //Part-11_ThetaMethodsIntegrated
            Kvv -= theta * tau * rho * bdf0 * Prod_BtCB;
            rv += (1.0 - theta) * tau * rho * prod(Prod_BtCB, arhs2);
        }

        //Part-03_ThetaMethodsIntegrated
        for (unsigned int i = 0; i < (TDim + 1); i++)
        {
            for (unsigned int j = 0; j < (TDim + 1); j++)
            {
                for (unsigned int k = 0; k < TDim; k++)
                {
                    Kvp_aux(i * TDim + k, j) += -DN(i, k) * N(j);
                }
            }
        }
        Kvp += theta * Kvp_aux;
        rv -= (1.0 - theta) * prod(Kvp_aux, p_n);

        //Part-05_ThetaMethodsIntegrated
        for (int i = 0; i < (TDim + 1); i++)
        {
            for (int k = 0; k < TDim; k++)
            {
                rv[i * TDim + k] += theta * rho * N[i] * body_force_on_gauss_point[k];
                rv[i * TDim + k] += (1.0 - theta) * rho * N[i] * body_force_on_gauss_point_n[k];
            }
        }

        //Part-04
        for (int i = 0; i < (TDim + 1); i++)
        {
            for (int j = 0; j < (TDim + 1); j++)
            {
                for (int k = 0; k < TDim; k++)
                {
                    Kpv(i, j * TDim + k) += N(i) * DN(j, k);
                }
            }
        }

        //      //Part-07---
        //      Prod_TrDNtr = prod(Tau, trans(DN));
        //      Kpp += prod(DN, Prod_TrDNtr);

        //Part-07
        Kpp += tau * prod(DN, trans(DN));

        //      //Part-08---
        //      tauf = rho * prod(Tau, body_force_on_gauss_point);
        //      rp  += prod(DN, tauf);

        //Part-08
        rp += rho * tau * prod(DN, body_force_on_gauss_point);

        //      //Part-10---
        //      noalias(Prod_TaudummyN) = prod(Tau, dummyN);
        //      Kpv += bdf0 *  rho  * prod(DN, Prod_TaudummyN);
        //      noalias(Prod_DNDX_Tau) = prod(DN, Tau);
        //      rp += -rho *  prod(Prod_DNDX_Tau, arhs2);

        //Part-10
        Kpv += bdf0 * rho * tau * prod(DN, dummyN);
        rp += -rho * tau * prod(DN, arhs2);

        if (dynamicsubgridstabilize)
        {
            //Part-12
            utilde_n = this->utilde_n;
            rp += rho * tau * (1.0 / dt) * prod(DN, utilde_n);
        }

        if (acceleratedparticles)
        { 
         for (int i = 0; i < (TDim + 1); i++)
         {
            for (int k = 0; k < TDim; k++)
            {
                rv[i * TDim + k] -= N[i] * a_n_on_gauss_point[k];
            }
         }
        }

        

        for (int i = 0; i < (TDim + 1); i++)
        {
            for (int j = 0; j < (TDim + 1); j++)
            {
                for (int k = 0; k < TDim; k++)
                {
                    for (int l = 0; l < TDim; l++)
                    {
                        lhs(i * (TDim + 1) + k, j * (TDim + 1) + l) += Kvv(i * (TDim) + k, j * (TDim) + l);
                    }
                }
            }
        }

        for (int i = 0; i < (TDim + 1); i++)
        {
            for (int j = 0; j < (TDim + 1); j++)
            {
                for (int k = 0; k < TDim; k++)
                {
                    lhs(i * (TDim + 1) + k, j * (TDim + 1) + TDim) += Kvp(i * TDim + k, j);
                }
            }
        }

        for (int i = 0; i < (TDim + 1); i++)
        {
            for (int j = 0; j < (TDim + 1); j++)
            {
                for (int k = 0; k < TDim; k++)
                {
                    lhs(i * (TDim + 1) + TDim, j * (TDim + 1) + k) += Kpv(i, j * TDim + k);
                }
            }
        }

        for (int i = 0; i < (TDim + 1); i++)
        {
            for (int j = 0; j < (TDim + 1); j++)
            {
                lhs(i * (TDim + 1) + TDim, j * (TDim + 1) + TDim) += Kpp(i, j);
            }
        }

        int counter = 0;
        for (int i = 0; i < (TDim + 1); i++)
        {
            for (int j = 0; j < TDim; j++)
            {
                rhs(counter++) += rv(i * TDim + j);
            }
            rhs(counter++) += rp(i);
        }
    }

    void MonolithicPFEM22DDenizNitsche::ComputeGaussPointLHSandRHSContribution_NitscheDomainTerms(bounded_matrix<double, 18, 18> &lhs, array_1d<double, 18> &rhs, const element_data &data, const ProcessInfo &rCurrentProcessInfo, const double &weight, boost::numeric::ublas::bounded_matrix<double, (18), (18)> &Tdublicate)
    {

        const unsigned int NumNodes = 3;
        const unsigned int TDim = 2;
        lhs = ZeroMatrix(18, 18);
        rhs = ZeroVector(18);

        boost::numeric::ublas::bounded_matrix<double, (TDim + 1) * (TDim) / 2, TDim *(TDim + 1)> B = ZeroMatrix((TDim + 1) * (TDim) / 2, TDim * (TDim + 1));       //B matrix
        boost::numeric::ublas::bounded_matrix<double, (TDim + 1) * (TDim) / 2, TDim *(TDim + 1)> Prod_CB = ZeroMatrix((TDim + 1) * (TDim) / 2, TDim * (TDim + 1)); //Prod_CB matrix
        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1), TDim *(TDim + 1)> Prod_BtCB = ZeroMatrix(TDim * (TDim + 1), TDim * (TDim + 1));            //Prod_BtCB matrix
        boost::numeric::ublas::bounded_matrix<double, TDim, (TDim + 1)> Prod_TrDNtr = ZeroMatrix(TDim, (TDim + 1));
        boost::numeric::ublas::bounded_matrix<double, (TDim), (TDim * (TDim + 1))> dummyN = ZeroMatrix((TDim), TDim * (TDim + 1));

        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1), (TDim + 1)> Kvp_aux = ZeroMatrix(TDim * (TDim + 1), (TDim + 1)); //Kvp_aux matrix
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> Tau = ZeroMatrix(TDim, TDim);                                          //Stabilization Tau vector
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim *(TDim + 1)> Prod_TaudummyN = ZeroMatrix(TDim, TDim * (TDim + 1));      //Prod_TaudummyN
        boost::numeric::ublas::bounded_matrix<double, (TDim + 1), TDim> Prod_DNDX_Tau = ZeroMatrix(TDim + 1, TDim);
        array_1d<double, TDim> tauf = ZeroVector(TDim);                          //An auxilary stabilization vector
        array_1d<double, TDim *(TDim + 1)> arhs = ZeroVector(TDim * (TDim + 1)); //arhs vector_bdf1*vn+bdf2*vnn
        array_1d<double, TDim> arhs2 = ZeroVector(TDim);
        array_1d<double, (TDim)> body_force_on_gauss_point = ZeroVector(TDim);                     //body force at tn+1 on the gauss point
        array_1d<double, (TDim)> body_force_on_gauss_point_nitsche = ZeroVector(TDim);             //body force at tn+1 on the gauss point
        array_1d<double, (TDim)> body_force_on_gauss_point_n = ZeroVector(TDim);                   //body force at tn on the gauss point
        array_1d<double, (TDim)> body_force_on_gauss_point_nitsche_n = ZeroVector(TDim);           //body force at tn on the gauss point
        array_1d<double, (TDim * (TDim + 1))> body_force_vector = ZeroVector(TDim * (TDim + 1));   //body force at tn+1 on the gauss point
        array_1d<double, (TDim * (TDim + 1))> body_force_vector_n = ZeroVector(TDim * (TDim + 1)); //body force at tn+1 on the gauss point
        array_1d<double, TDim *(TDim + 1)> v_n = ZeroVector(TDim * (TDim + 1));                    //v_n vector
        array_1d<double, (TDim + 1)> p_n = ZeroVector((TDim + 1));                                 //p_n vector

        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1), TDim *(TDim + 1)> Kvv = ZeroMatrix(TDim * (TDim + 1), TDim * (TDim + 1));     //Kvv matrix
        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1), TDim *(TDim + 1)> Kvv_aux = ZeroMatrix(TDim * (TDim + 1), TDim * (TDim + 1)); //Kvv_aux matrix
        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1), (TDim + 1)> Kvp = ZeroMatrix(TDim * (TDim + 1), (TDim + 1));                  //Kvp matrix
        boost::numeric::ublas::bounded_matrix<double, (TDim + 1), TDim *(TDim + 1)> Kpv = ZeroMatrix(TDim + 1, TDim * (TDim + 1));                    //Kpv matrix
        boost::numeric::ublas::bounded_matrix<double, (TDim + 1), (TDim + 1)> Kpp = ZeroMatrix(TDim + 1, TDim + 1);                                   //Kpp matrix
        array_1d<double, TDim *(TDim + 1)> rv = ZeroVector(TDim * (TDim + 1));
        array_1d<double, (TDim + 1)> rp = ZeroVector((TDim + 1));
        array_1d<double, (TDim * (TDim + 1))> div_vec = ZeroVector(TDim * (TDim + 1)); //div_vec vector

        const double theta = GetProperties()[TIMEDIS_THETA];
        const bool momentumstabilize = false;
        //get constitutive matrix
        const Matrix &C = data.C;
        //get density
        const double rho = inner_prod(data.N, data.rho);
        //get bdf coefficients and dt
        const double &bdf0 = data.bdf0;
        const double &bdf1 = data.bdf1;
        const double &bdf2 = data.bdf2;
        double dt = rCurrentProcessInfo[DELTA_TIME];

        //get shape function values
        const bounded_matrix<double, NumNodes, TDim> &DN = data.DN_DX;
        const array_1d<double, NumNodes> &N = data.N;
        for (unsigned int i = 0; i < (TDim + 1); i++)
        {
            dummyN(0, i * TDim) = N(i);
            dummyN(1, i * TDim + 1) = N(i);
            if (TDim == 3)
                dummyN(2, i * TDim + 2) = N(i);
        }
        //calculate body force on gauss points for tn+1 and tn
        const bounded_matrix<double, NumNodes, TDim> &f = data.f;
        const bounded_matrix<double, NumNodes, TDim> &f_nitsche = data.f_nitsche;
        const bounded_matrix<double, NumNodes, TDim> &f_n = data.f_n;
        const bounded_matrix<double, NumNodes, TDim> &f_nitsche_n = data.f_nitsche_n;
        //      body_force_on_gauss_point = prod(trans(f), N);
        //      body_force_on_gauss_point_n = prod(trans(f_n), N);

        for (unsigned int i = 0; i < (TDim + 1); i++)
        {
            if (Tdublicate(i * (TDim + 1), i * (TDim + 1)) == 1.0)
            {
                v_n[TDim * i] = (GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X, 1));
                v_n[TDim * i + 1] = (GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y, 1));
                p_n(i) = (GetGeometry()[i].FastGetSolutionStepValue(PRESSURE, 1));
                body_force_on_gauss_point(0) += N[i] * (GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE_X));
                body_force_on_gauss_point_n(0) += N[i] * (GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE_X, 1));
                body_force_on_gauss_point(1) += N[i] * (GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE_Y));
                body_force_on_gauss_point_n(1) += N[i] * (GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE_Y, 1));
                body_force_vector[TDim * i] = f(i, 0);
                body_force_vector[TDim * i + 1] = f(i, 1);
                body_force_vector_n[TDim * i] = f_n(i, 0);
                body_force_vector_n[TDim * i + 1] = f_n(i, 1);
            }
            else
            {
                v_n[TDim * i] = (GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_NITSCHE_X, 1));
                v_n[TDim * i + 1] = (GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_NITSCHE_Y, 1));
                p_n(i) = (GetGeometry()[i].FastGetSolutionStepValue(PRESSURE_NITSCHE, 1));
                body_force_on_gauss_point(0) += N[i] * (GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE_NITSCHE_X));
                body_force_on_gauss_point_n(0) += N[i] * (GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE_NITSCHE_X, 1));
                body_force_on_gauss_point(1) += N[i] * (GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE_NITSCHE_Y));
                body_force_on_gauss_point_n(1) += N[i] * (GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE_NITSCHE_Y, 1));
                body_force_vector[TDim * i] = f_nitsche(i, 0);
                body_force_vector[TDim * i + 1] = f_nitsche(i, 1);
                body_force_vector_n[TDim * i] = f_nitsche_n(i, 0);
                body_force_vector_n[TDim * i + 1] = f_nitsche_n(i, 1);
            }
            arhs[TDim * i] = bdf1 * v_n(TDim * i);
            arhs[TDim * i + 1] = bdf1 * v_n(TDim * i + 1);
            arhs2[0] += N[i] * arhs[TDim * i];
            arhs2[1] += N[i] * arhs[TDim * i + 1];
            if (TDim == 3)
            {
                if (Tdublicate(i * (TDim + 1), i * (TDim + 1)) == 1.0)
                {
                    v_n[TDim * i + 2] = (GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Z, 1));
                    body_force_on_gauss_point(2) += N[i] * (GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE_Z));
                    body_force_on_gauss_point_n(2) += N[i] * (GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE_Z, 1));
                    body_force_vector[TDim * i + 2] = f(i, 1);
                    body_force_vector_n[TDim * i + 2] = f_n(i, 1);
                }
                else
                {
                    v_n[TDim * i + 2] = (GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_NITSCHE_Z, 1));
                    body_force_on_gauss_point(2) += N[i] * (GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE_NITSCHE_Z));
                    body_force_on_gauss_point_n(2) += N[i] * (GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE_NITSCHE_Z, 1));
                    body_force_vector[TDim * i + 2] = f_nitsche(i, 2);
                    body_force_vector_n[TDim * i + 2] = f_nitsche_n(i, 2);
                }
                arhs[TDim * i + 2] = bdf1 * v_n(TDim * i + 2);
                arhs2[2] += N[i] * arhs[TDim * i + 2];
            }
        }

        //Tau for stabilization
        double sqElemSize = 0.0;
        double mElemSize;
        array_1d<double, 3> Edge(3, 0.0);
        Edge = this->GetGeometry()[1].Coordinates() - this->GetGeometry()[0].Coordinates();
        sqElemSize = Edge[0] * Edge[0];
        for (SizeType d = 1; d < TDim; d++)
            sqElemSize += Edge[d] * Edge[d];

        for (SizeType i = 2; i < TDim + 1; i++)
            for (SizeType j = 0; j < i; j++)
            {
                Edge = this->GetGeometry()[i].Coordinates() - this->GetGeometry()[j].Coordinates();
                double Length = Edge[0] * Edge[0];
                for (SizeType d = 1; d < TDim; d++)
                    Length += Edge[d] * Edge[d];
                if (Length < sqElemSize)
                    sqElemSize = Length;
            }
        mElemSize = sqrt(sqElemSize);
        double tau = 0.0;
        const double tau2 = C(2, 2);
        if (bdf0 != 0.0)
        {
            tau = 1.0 / ((1.0 * rho / (dt) + 4.0 * C(2, 2) / (mElemSize * mElemSize)));
        }
        else
        {
            tau = 1.0 / ((4.0 * C(2, 2) / (mElemSize * mElemSize)));
        }

        for (unsigned int i = 0; i < (TDim); i++)
        {
            Tau(i, i) = tau;
        }

        //Part-01
        for (unsigned int i = 0; i < (TDim + 1); i++)
        {
            double tmp = (rho * bdf0) * N(i);
            Kvv(i * TDim, i * TDim) += tmp;
            Kvv(i * TDim + 1, i * TDim + 1) += tmp;
        }
        for (unsigned int i = 0; i < (TDim + 1); i++)
        {
            rv(i * TDim) -= (rho * N(i)) * arhs(i * TDim);
            rv(i * TDim + 1) -= (rho * N(i)) * arhs(i * TDim + 1);
        }

        //Part-02_ThetaMethodsIntegrated
        for (unsigned int i = 0; i < (TDim + 1); i++)
        {
            unsigned int index = i * TDim;
            B(0, index + 0) = DN(i, 0);
            B(0, index + 1) = 0;
            B(1, index + 0) = 0;
            B(1, index + 1) = DN(i, 1);
            B(2, index + 0) = DN(i, 1);
            B(2, index + 1) = DN(i, 0);
        }
        noalias(Prod_CB) = prod(C, B);
        noalias(Prod_BtCB) = prod(trans(B), Prod_CB);
        Kvv += theta * Prod_BtCB;
        rv -= (1.0 - theta) * prod(Prod_BtCB, v_n);

        if (momentumstabilize == true)
        {
            //Part-06_ThetaMethodsIntegrated
            rv -= (theta * rho * tau * prod(Prod_BtCB, body_force_vector) + (1.0 - theta) * rho * tau * prod(Prod_BtCB, body_force_vector_n));

            //Part-09_ThetaMethodsIntegrated
            for (int i = 0; i < (TDim + 1); i++)
            {
                for (int k = 0; k < (TDim); k++)
                {
                    div_vec(TDim * i + k) = DN(i, k);
                }
            }
            Kvv_aux += outer_prod(div_vec, div_vec);
            Kvv -= theta * tau2 * Kvv_aux;
            rv += (1.0 - theta) * tau2 * prod(Kvv_aux, v_n);
            //Part-11_ThetaMethodsIntegrated
            Kvv -= theta * tau * rho * bdf0 * Prod_BtCB;
            rv += (1.0 - theta) * tau * rho * prod(Prod_BtCB, arhs2);
        }

        //Part-03_ThetaMethodsIntegrated
        for (unsigned int i = 0; i < (TDim + 1); i++)
        {
            for (unsigned int j = 0; j < (TDim + 1); j++)
            {
                for (unsigned int k = 0; k < TDim; k++)
                {
                    Kvp_aux(i * TDim + k, j) += -DN(i, k) * N(j);
                }
            }
        }
        Kvp += theta * Kvp_aux;
        rv -= (1.0 - theta) * prod(Kvp_aux, p_n);

        //Part-05_ThetaMethodsIntegrated
        for (int i = 0; i < (TDim + 1); i++)
        {
            for (int k = 0; k < TDim; k++)
            {
                rv[i * TDim + k] += theta * rho * N[i] * body_force_on_gauss_point[k];
                rv[i * TDim + k] += (1.0 - theta) * rho * N[i] * body_force_on_gauss_point_n[k];
            }
        }

        //Part-04
        for (int i = 0; i < (TDim + 1); i++)
        {
            for (int j = 0; j < (TDim + 1); j++)
            {
                for (int k = 0; k < TDim; k++)
                {
                    Kpv(i, j * TDim + k) += N(i) * DN(j, k);
                }
            }
        }

        //      //Part-07---
        //      Prod_TrDNtr = prod(Tau, trans(DN));
        //      Kpp += prod(DN, Prod_TrDNtr);

        //Part-07
        Kpp += tau * prod(DN, trans(DN));

        //      //Part-08---
        //      tauf = rho * prod(Tau, body_force_on_gauss_point);
        //      rp  += prod(DN, tauf);

        //Part-08
        rp += rho * tau * prod(DN, body_force_on_gauss_point);

        //      //Part-10---
        //      noalias(Prod_TaudummyN) = prod(Tau, dummyN);
        //      Kpv += bdf0 *  rho  * prod(DN, Prod_TaudummyN);
        //      noalias(Prod_DNDX_Tau) = prod(DN, Tau);
        //      rp += -rho *  prod(Prod_DNDX_Tau, arhs2);

        //Part-10
        Kpv += bdf0 * rho * tau * prod(DN, dummyN);
        rp += -rho * tau * prod(DN, arhs2);

        for (int i = 0; i < (TDim + 1); i++)
        {
            for (int j = 0; j < (TDim + 1); j++)
            {
                for (int k = 0; k < TDim; k++)
                {
                    for (int l = 0; l < TDim; l++)
                    {
                        lhs(i * (TDim + 1) + k, j * (TDim + 1) + l) += Kvv(i * (TDim) + k, j * (TDim) + l);
                    }
                }
            }
        }

        for (int i = 0; i < (TDim + 1); i++)
        {
            for (int j = 0; j < (TDim + 1); j++)
            {
                for (int k = 0; k < TDim; k++)
                {
                    lhs(i * (TDim + 1) + k, j * (TDim + 1) + TDim) += Kvp(i * TDim + k, j);
                }
            }
        }

        for (int i = 0; i < (TDim + 1); i++)
        {
            for (int j = 0; j < (TDim + 1); j++)
            {
                for (int k = 0; k < TDim; k++)
                {
                    lhs(i * (TDim + 1) + TDim, j * (TDim + 1) + k) += Kpv(i, j * TDim + k);
                }
            }
        }

        for (int i = 0; i < (TDim + 1); i++)
        {
            for (int j = 0; j < (TDim + 1); j++)
            {
                lhs(i * (TDim + 1) + TDim, j * (TDim + 1) + TDim) += Kpp(i, j);
            }
        }

        int counter = 0;
        for (int i = 0; i < (TDim + 1); i++)
        {
            for (int j = 0; j < TDim; j++)
            {
                rhs(counter++) += rv(i * TDim + j);
            }
            rhs(counter++) += rp(i);
        }

        lhs = prod(lhs, trans(Tdublicate));
        lhs = prod(Tdublicate, lhs);
        rhs = prod(Tdublicate, rhs);
    }

    void MonolithicPFEM22DDenizNitsche::ComputeGaussPointLHSandRHSContribution_NitscheBoundaryTerms(bounded_matrix<double, 18, 18> &lhs, array_1d<double, 18> &rhs, const element_data &data, const ProcessInfo &rCurrentProcessInfo, const double &weight, array_1d<double, 3> &normal, boost::numeric::ublas::bounded_matrix<double, (18), (18)> &Tdublicate, boost::numeric::ublas::bounded_matrix<double, (18), (18)> &Tdublicate_neighbor, boost::numeric::ublas::bounded_matrix<double, (6), (6)> &Tself, boost::numeric::ublas::bounded_matrix<double, (6), (6)> &Tneighbor, const Matrix &C_self, const Matrix &C_neighbor)
    {

        const unsigned int NumNodes = 3;
        const unsigned int TDim = 2;
        lhs = ZeroMatrix(18, 18);
        rhs = ZeroVector(18);

        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1), (TDim + 1) * (TDim + 1)> Kvup_nitsche = ZeroMatrix(TDim * (TDim + 1), (TDim + 1) * (TDim + 1));
        boost::numeric::ublas::bounded_matrix<double, (TDim + 1) * (TDim + 1), TDim *(TDim + 1)> Kupv_nitsche = ZeroMatrix((TDim + 1) * (TDim + 1), TDim * (TDim + 1));
        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1), TDim *(TDim + 1)> Kvv_nitsche = ZeroMatrix(TDim * (TDim + 1), TDim * (TDim + 1));
        boost::numeric::ublas::bounded_matrix<double, 1, TDim *(TDim + 1)> Kvv_nitsche_Fmatrix = ZeroMatrix(1, TDim * (TDim + 1));
        boost::numeric::ublas::bounded_matrix<double, 1, (TDim + 1) * (TDim + 1)> Kvup_nitsche_Bmatrix = ZeroMatrix(1, (TDim + 1) * (TDim + 1));
        boost::numeric::ublas::bounded_matrix<double, 1, TDim *(TDim + 1)> Kupv_nitsche_Dmatrix = ZeroMatrix(1, TDim * (TDim + 1));
        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1), TDim> Kvu_nitsche_Cmatrix = ZeroMatrix(TDim * (TDim + 1), TDim);
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim *(TDim + 1)> Kvu_nitsche_Dmatrix = ZeroMatrix(TDim, TDim * (TDim + 1));
        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1) * 3, TDim *(TDim + 1) * 3> Kvu_nitsche_large = ZeroMatrix(TDim * (TDim + 1) * 3, TDim * (TDim + 1) * 3);
        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1), 1> Kvp_nitsche_Ematrix = ZeroMatrix(TDim * (TDim + 1), 1);
        boost::numeric::ublas::bounded_matrix<double, 1, TDim + 1> Kvp_nitsche_Fmatrix = ZeroMatrix(1, TDim + 1);
        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1) * 3, TDim *(TDim + 1) * 3> Kvp_nitsche_large = ZeroMatrix(TDim * (TDim + 1) * 3, TDim * (TDim + 1) * 3);
        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1), TDim + 1> Kvp_nitsche = ZeroMatrix(TDim * (TDim + 1), TDim + 1);
        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1), TDim *(TDim + 1)> Kvu_nitsche = ZeroMatrix(TDim * (TDim + 1), TDim * (TDim + 1));

        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1), TDim *(TDim + 1)> Kuv_nitsche = ZeroMatrix(TDim * (TDim + 1), TDim * (TDim + 1));
        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1) * 3, TDim *(TDim + 1) * 3> Kuv_nitsche_self_large = ZeroMatrix(TDim * (TDim + 1) * 3, TDim * (TDim + 1) * 3);
        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1) * 3, TDim *(TDim + 1) * 3> Kuv_nitsche_neighbor_large = ZeroMatrix(TDim * (TDim + 1) * 3, TDim * (TDim + 1) * 3);

        boost::numeric::ublas::bounded_matrix<double, TDim + 1, TDim *(TDim + 1)> Kpv_nitsche = ZeroMatrix(TDim + 1, TDim * (TDim + 1));
        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1) * 3, TDim *(TDim + 1) * 3> Kpv_nitsche_self_large = ZeroMatrix(TDim * (TDim + 1) * 3, TDim * (TDim + 1) * 3);
        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1) * 3, TDim *(TDim + 1) * 3> Kpv_nitsche_neighbor_large = ZeroMatrix(TDim * (TDim + 1) * 3, TDim * (TDim + 1) * 3);
        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1), TDim *(TDim + 1)> Kvu_nitsche_alpha_self = ZeroMatrix(TDim * (TDim + 1), TDim * (TDim + 1));

        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1) * 3, TDim *(TDim + 1) * 3> Kvu_nitsche_alpha_self_large = ZeroMatrix(TDim * (TDim + 1) * 3, TDim * (TDim + 1) * 3);
        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1), TDim *(TDim + 1)> Kvu_nitsche_alpha_neighbor = ZeroMatrix(TDim * (TDim + 1), TDim * (TDim + 1));
        boost::numeric::ublas::bounded_matrix<double, TDim *(TDim + 1) * 3, TDim *(TDim + 1) * 3> Kvu_nitsche_alpha_neighbor_large = ZeroMatrix(TDim * (TDim + 1) * 3, TDim * (TDim + 1) * 3);

        boost::numeric::ublas::bounded_matrix<double, 2, TDim + 1> Kvp_nitsche_Dmatrix = ZeroMatrix(2, TDim + 1);

        array_1d<double, TDim *(TDim + 1)> Kvup_nitsche_Avector = ZeroVector(TDim * (TDim + 1));
        array_1d<double, (TDim + 1) * (TDim + 1)> Kupv_nitsche_Cvector = ZeroVector((TDim + 1) * (TDim + 1));
        array_1d<double, TDim *(TDim + 1)> Kvv_nitsche_Evector = ZeroVector(TDim * (TDim + 1));
        array_1d<double, (TDim + 1) * (TDim + 1)> rgup_nitsche = ZeroVector((TDim + 1) * (TDim + 1));
        array_1d<double, TDim *(TDim + 1)> rvv_nitsche = ZeroVector(TDim * (TDim + 1));

        const double alph = GetProperties()[NITSCHE_ALPHA];
        //     const double delt = GetProperties()[NITSCHE_DELTA];
        const double delt = 1.0;
        const int casenumber = GetProperties()[CASE_NUMBER];
        const double theta = GetProperties()[TIMEDIS_THETA];

        //     array_1d<double,(TDim+1)> Vwall = ZeroVector(TDim+1);
        //     Vwall[0]= 0.0;
        //     Vwall[1]= 0.0;
        //     Vwall[2]= 0.0;
        //     double g = inner_prod(Vwall, normal);
        double g = 0.0;

        //get shape function values
        const bounded_matrix<double, NumNodes, TDim> &DN = data.DN_DX;
        const array_1d<double, NumNodes> &N = data.N;

        array_1d<double, NumNodes * 2> N_small = ZeroVector(NumNodes * 2);
        array_1d<double, NumNodes * 2> N_self = ZeroVector(NumNodes * 2);
        array_1d<double, NumNodes * 2> N_neighbor = ZeroVector(NumNodes * 2);

        for (int i = 0; i < NumNodes; i++)
            N_small[i] = N[i];

        N_self = prod(Tself, N_small);
        N_neighbor = prod(Tneighbor, N_small);
        const double rho = inner_prod(data.N, data.rho);

        //get constitutive matrix
        const Matrix &C = data.C;
        double mu = C(2, 2);
        double mu_self = C_self(2, 2);
        double mu_neighbor = C_neighbor(2, 2);

        double sqElemSize = 0.0;
        double mElemSize;
        array_1d<double, 3> Edge(3, 0.0);
        Edge = this->GetGeometry()[1].Coordinates() - this->GetGeometry()[0].Coordinates();
        sqElemSize = Edge[0] * Edge[0];
        for (SizeType d = 1; d < TDim; d++)
            sqElemSize += Edge[d] * Edge[d];

        for (SizeType i = 2; i < TDim + 1; i++)
            for (SizeType j = 0; j < i; j++)
            {
                Edge = this->GetGeometry()[i].Coordinates() - this->GetGeometry()[j].Coordinates();
                double Length = Edge[0] * Edge[0];
                for (SizeType d = 1; d < TDim; d++)
                    Length += Edge[d] * Edge[d];
                if (Length < sqElemSize)
                    sqElemSize = Length;
            }
        mElemSize = sqrt(sqElemSize);

        //Free-Slip Condition
        if (casenumber == 1)
        {

            //Kvup_nitsche calculation
            Kvup_nitsche_Avector[0] = N[0] * normal[0];
            Kvup_nitsche_Avector[1] = N[0] * normal[1];
            Kvup_nitsche_Avector[2] = N[1] * normal[0];
            Kvup_nitsche_Avector[3] = N[1] * normal[1];
            Kvup_nitsche_Avector[4] = N[2] * normal[0];
            Kvup_nitsche_Avector[5] = N[2] * normal[1];

            Kvup_nitsche_Bmatrix(0, 0) = 2.0 * mu * normal[0] * normal[0] * DN(0, 0) + 2.0 * mu * normal[0] * normal[1] * DN(0, 1);
            Kvup_nitsche_Bmatrix(0, 1) = 2.0 * mu * normal[0] * normal[1] * DN(0, 0) + 2.0 * mu * normal[1] * normal[1] * DN(0, 1);
            Kvup_nitsche_Bmatrix(0, 2) = -N[0];

            Kvup_nitsche_Bmatrix(0, 3) = 2.0 * mu * normal[0] * normal[0] * DN(1, 0) + 2.0 * mu * normal[0] * normal[1] * DN(1, 1);
            Kvup_nitsche_Bmatrix(0, 4) = 2.0 * mu * normal[0] * normal[1] * DN(1, 0) + 2.0 * mu * normal[1] * normal[1] * DN(1, 1);
            Kvup_nitsche_Bmatrix(0, 5) = -N[1];

            Kvup_nitsche_Bmatrix(0, 6) = 2.0 * mu * normal[0] * normal[0] * DN(2, 0) + 2.0 * mu * normal[0] * normal[1] * DN(2, 1);
            Kvup_nitsche_Bmatrix(0, 7) = 2.0 * mu * normal[0] * normal[1] * DN(2, 0) + 2.0 * mu * normal[1] * normal[1] * DN(2, 1);
            Kvup_nitsche_Bmatrix(0, 8) = -N[2];

            for (unsigned int i = 0; i < TDim * (TDim + 1); i++)
            {
                for (unsigned int j = 0; j < (TDim + 1) * (TDim + 1); j++)
                {
                    Kvup_nitsche(i, j) = -Kvup_nitsche_Avector[i] * Kvup_nitsche_Bmatrix(0, j);
                }
            }

            //Kupv_nitsche calculation
            Kupv_nitsche_Cvector[0] = 2.0 * mu * normal[0] * normal[0] * DN(0, 0) + 2.0 * mu * normal[0] * normal[1] * DN(0, 1);
            Kupv_nitsche_Cvector[1] = 2.0 * mu * normal[0] * normal[1] * DN(0, 0) + 2.0 * mu * normal[1] * normal[1] * DN(0, 1);
            Kupv_nitsche_Cvector[2] = -N[0];
            Kupv_nitsche_Cvector[3] = 2.0 * mu * normal[0] * normal[0] * DN(1, 0) + 2.0 * mu * normal[0] * normal[1] * DN(1, 1);
            Kupv_nitsche_Cvector[4] = 2.0 * mu * normal[0] * normal[1] * DN(1, 0) + 2.0 * mu * normal[1] * normal[1] * DN(1, 1);
            Kupv_nitsche_Cvector[5] = -N[1];
            Kupv_nitsche_Cvector[6] = 2.0 * mu * normal[0] * normal[0] * DN(2, 0) + 2.0 * mu * normal[0] * normal[1] * DN(2, 1);
            Kupv_nitsche_Cvector[7] = 2.0 * mu * normal[0] * normal[1] * DN(2, 0) + 2.0 * mu * normal[1] * normal[1] * DN(2, 1);
            Kupv_nitsche_Cvector[8] = -N[2];

            Kupv_nitsche_Dmatrix(0, 0) = N[0] * normal[0];
            Kupv_nitsche_Dmatrix(0, 1) = N[0] * normal[1];

            Kupv_nitsche_Dmatrix(0, 2) = N[1] * normal[0];
            Kupv_nitsche_Dmatrix(0, 3) = N[1] * normal[1];

            Kupv_nitsche_Dmatrix(0, 4) = N[2] * normal[0];
            Kupv_nitsche_Dmatrix(0, 5) = N[2] * normal[1];

            for (unsigned int i = 0; i < (TDim + 1) * (TDim + 1); i++)
            {
                for (unsigned int j = 0; j < TDim * (TDim + 1); j++)
                {
                    Kupv_nitsche(i, j) = -delt * Kupv_nitsche_Cvector[i] * Kupv_nitsche_Dmatrix(0, j);
                }
            }

            //Kvv_nitsche calculation
            Kvv_nitsche_Evector[0] = N[0] * normal[0];
            Kvv_nitsche_Evector[1] = N[0] * normal[1];
            Kvv_nitsche_Evector[2] = N[1] * normal[0];
            Kvv_nitsche_Evector[3] = N[1] * normal[1];
            Kvv_nitsche_Evector[4] = N[2] * normal[0];
            Kvv_nitsche_Evector[5] = N[2] * normal[1];

            Kvv_nitsche_Fmatrix(0, 0) = N[0] * normal[0];
            Kvv_nitsche_Fmatrix(0, 1) = N[0] * normal[1];

            Kvv_nitsche_Fmatrix(0, 2) = N[1] * normal[0];
            Kvv_nitsche_Fmatrix(0, 3) = N[1] * normal[1];

            Kvv_nitsche_Fmatrix(0, 4) = N[2] * normal[0];
            Kvv_nitsche_Fmatrix(0, 5) = N[2] * normal[1];

            for (unsigned int i = 0; i < TDim * (TDim + 1); i++)
            {
                for (unsigned int j = 0; j < TDim * (TDim + 1); j++)
                {
                    Kvv_nitsche(i, j) += (1 / mElemSize) * (1 / alph) * Kvv_nitsche_Evector[i] * Kvv_nitsche_Fmatrix(0, j);
                }
            }

            for (unsigned int i = 0; i < (TDim + 1) * (TDim + 1); i++)
            {
                rgup_nitsche[i] = -delt * g * Kupv_nitsche_Cvector[i];
            }

            for (unsigned int i = 0; i < TDim * (TDim + 1); i++)
            {
                rvv_nitsche[i] = (1 / mElemSize) * (1 / alph) * g * Kvv_nitsche_Evector[i];
            }

            for (int i = 0; i < (TDim + 1); i++)
            {
                for (int j = 0; j < (TDim + 1) * (TDim + 1); j++)
                {
                    for (int k = 0; k < TDim; k++)
                    {
                        lhs(i * (TDim + 1) + k, j) += Kvup_nitsche(i * TDim + k, j);
                    }
                }
            }

            for (int i = 0; i < (TDim + 1) * (TDim + 1); i++)
            {
                for (int j = 0; j < (TDim + 1); j++)
                {
                    for (int k = 0; k < TDim; k++)
                    {
                        lhs(i, j * (TDim + 1) + k) += Kupv_nitsche(i, j * TDim + k);
                    }
                }
            }

            for (int i = 0; i < (TDim + 1); i++)
            {
                for (int j = 0; j < (TDim + 1); j++)
                {
                    for (int k = 0; k < TDim; k++)
                    {
                        for (int l = 0; l < TDim; l++)
                        {
                            lhs(i * (TDim + 1) + k, j * (TDim + 1) + l) += Kvv_nitsche(i * (TDim) + k, j * (TDim) + l);
                        }
                    }
                }
            }

            int counter = 0;
            for (int i = 0; i < (TDim + 1) * (TDim + 1); i++)
            {
                rhs(counter++) += rgup_nitsche[i];
            }

            counter = 0;
            for (int i = 0; i < (TDim + 1); i++)
            {
                for (int j = 0; j < TDim; j++)
                {
                    rhs(counter++) += rvv_nitsche(i * TDim + j);
                }
                rhs(counter++) += 0.0;
            }

            lhs = prod(lhs, trans(Tdublicate));
            lhs = prod(Tdublicate, lhs);
            rhs = prod(Tdublicate, rhs);
        }

        if (casenumber == 3)
        {

            //Kvu_nitsche calculation-start
            //Term-1
            for (unsigned int i = 0; i < (TDim + 1); i++)
            {
                Kvu_nitsche_Cmatrix(i * TDim, 0) = N[i];
                Kvu_nitsche_Cmatrix(i * TDim + 1, 1) = N[i];

                Kvu_nitsche_Dmatrix(0, i * TDim) = 2.0 * mu_neighbor * (DN(i, 0) * normal(0) + 0.5 * DN(i, 1) * normal(1));
                Kvu_nitsche_Dmatrix(0, i * TDim + 1) = 2.0 * mu_neighbor * (0.5 * DN(i, 0) * normal(1));

                Kvu_nitsche_Dmatrix(1, i * TDim) = 2.0 * mu_neighbor * (0.5 * DN(i, 1) * normal(0));
                Kvu_nitsche_Dmatrix(1, i * TDim + 1) = 2.0 * mu_neighbor * (DN(i, 1) * normal(1) + 0.5 * DN(i, 0) * normal(0));
            }

            Kvu_nitsche = -1.0 * prod(Kvu_nitsche_Cmatrix, Kvu_nitsche_Dmatrix);

            for (int i = 0; i < (TDim + 1); i++)
            {
                for (int j = 0; j < (TDim + 1); j++)
                {
                    for (int k = 0; k < TDim; k++)
                    {
                        for (int l = 0; l < TDim; l++)
                        {
                            Kvu_nitsche_large(i * (TDim + 1) + k, j * (TDim + 1) + l) += Kvu_nitsche(i * (TDim) + k, j * (TDim) + l);
                        }
                    }
                }
            }

            Kvu_nitsche_large = prod(Kvu_nitsche_large, trans(Tdublicate_neighbor));
            //     Kvu_nitsche_large = prod(Kvu_nitsche_large,trans(Tdublicate));
            Kvu_nitsche_large = prod(Tdublicate, Kvu_nitsche_large);
            //Kvu_nitsche calculation-end

            //Kvp_nitsche calculation-start
            //Term-2
            //     Kvp_nitsche_Ematrix(0,0)=N[0];
            //     Kvp_nitsche_Ematrix(1,0)=N[0];
            //     Kvp_nitsche_Ematrix(2,0)=N[1];
            //     Kvp_nitsche_Ematrix(3,0)=N[1];
            //     Kvp_nitsche_Ematrix(4,0)=N[2];
            //     Kvp_nitsche_Ematrix(5,0)=N[2];
            //
            //     Kvp_nitsche_Fmatrix(0,0)=N[0]*normal(0)+N[0]*normal(1);
            //     Kvp_nitsche_Fmatrix(0,1)=N[1]*normal(0)+N[1]*normal(1);
            //     Kvp_nitsche_Fmatrix(0,2)=N[2]*normal(0)+N[2]*normal(1);

            Kvp_nitsche_Dmatrix(0, 0) = N[0] * normal(0);
            Kvp_nitsche_Dmatrix(0, 1) = N[1] * normal(0);
            Kvp_nitsche_Dmatrix(0, 2) = N[2] * normal(0);
            Kvp_nitsche_Dmatrix(1, 0) = N[0] * normal(1);
            Kvp_nitsche_Dmatrix(1, 1) = N[1] * normal(1);
            Kvp_nitsche_Dmatrix(1, 2) = N[2] * normal(1);

            Kvp_nitsche = prod(Kvu_nitsche_Cmatrix, Kvp_nitsche_Dmatrix);

            for (int i = 0; i < (TDim + 1); i++)
            {
                for (int j = 0; j < (TDim + 1); j++)
                {
                    for (int k = 0; k < TDim; k++)
                    {
                        Kvp_nitsche_large(i * (TDim + 1) + k, j * (TDim + 1) + TDim) += Kvp_nitsche(i * TDim + k, j);
                    }
                }
            }

            Kvp_nitsche_large = prod(Kvp_nitsche_large, trans(Tdublicate_neighbor));
            //     Kvp_nitsche_large = prod(Kvp_nitsche_large,trans(Tdublicate));
            Kvp_nitsche_large = prod(Tdublicate, Kvp_nitsche_large);
            //Kvp_nitsche calculation-end

            //Kuv_nitsche calculation-start
            //Term-3
            Kuv_nitsche = trans(Kvu_nitsche);

            for (int i = 0; i < (TDim + 1); i++)
            {
                for (int j = 0; j < (TDim + 1); j++)
                {
                    for (int k = 0; k < TDim; k++)
                    {
                        for (int l = 0; l < TDim; l++)
                        {
                            Kuv_nitsche_self_large(i * (TDim + 1) + k, j * (TDim + 1) + l) += Kuv_nitsche(i * (TDim) + k, j * (TDim) + l);
                            Kuv_nitsche_neighbor_large(i * (TDim + 1) + k, j * (TDim + 1) + l) -= Kuv_nitsche(i * (TDim) + k, j * (TDim) + l);
                        }
                    }
                }
            }

            Kuv_nitsche_self_large = prod(Kuv_nitsche_self_large, trans(Tdublicate));
            Kuv_nitsche_self_large = prod(Tdublicate_neighbor, Kuv_nitsche_self_large);
            //     Kuv_nitsche_self_large = prod(Tdublicate, Kuv_nitsche_self_large);

            Kuv_nitsche_neighbor_large = prod(Kuv_nitsche_neighbor_large, trans(Tdublicate_neighbor));
            Kuv_nitsche_neighbor_large = prod(Tdublicate_neighbor, Kuv_nitsche_neighbor_large);
            //     Kuv_nitsche_neighbor_large = prod(Tdublicate, Kuv_nitsche_neighbor_large);
            //Kuv_nitsche calculation-end

            //Kpv_nitsche calculation-start
            Kpv_nitsche = trans(Kvp_nitsche);

            for (int i = 0; i < (TDim + 1); i++)
            {
                for (int j = 0; j < (TDim + 1); j++)
                {
                    for (int k = 0; k < TDim; k++)
                    {
                        Kpv_nitsche_self_large(i * (TDim + 1) + TDim, j * (TDim + 1) + k) += Kpv_nitsche(i, j * TDim + k);
                        Kpv_nitsche_neighbor_large(i * (TDim + 1) + TDim, j * (TDim + 1) + k) -= Kpv_nitsche(i, j * TDim + k);
                    }
                }
            }

            Kpv_nitsche_self_large = prod(Kpv_nitsche_self_large, trans(Tdublicate));
            Kpv_nitsche_self_large = prod(Tdublicate_neighbor, Kpv_nitsche_self_large);
            //     Kpv_nitsche_self_large = prod(Tdublicate, Kpv_nitsche_self_large);

            Kpv_nitsche_neighbor_large = prod(Kpv_nitsche_neighbor_large, trans(Tdublicate_neighbor));
            Kpv_nitsche_neighbor_large = prod(Tdublicate_neighbor, Kpv_nitsche_neighbor_large);
            //     Kpv_nitsche_neighbor_large = prod(Tdublicate, Kpv_nitsche_neighbor_large);
            //Kpv_nitsche calculation-end

            //Kvu_nitsche_alpha_self and Kvu_nitsche_alpha_neighbor calculation-start
            Kvu_nitsche_alpha_self = prod(Kvu_nitsche_Cmatrix, trans(Kvu_nitsche_Cmatrix));

            for (int i = 0; i < (TDim + 1); i++)
            {
                for (int j = 0; j < (TDim + 1); j++)
                {
                    for (int k = 0; k < TDim; k++)
                    {
                        for (int l = 0; l < TDim; l++)
                        {
                            Kvu_nitsche_alpha_self_large(i * (TDim + 1) + k, j * (TDim + 1) + l) += (1.0 / (alph * mElemSize)) * Kvu_nitsche_alpha_self(i * (TDim) + k, j * (TDim) + l);
                            Kvu_nitsche_alpha_neighbor_large(i * (TDim + 1) + k, j * (TDim + 1) + l) -= (1.0 / (alph * mElemSize)) * Kvu_nitsche_alpha_self(i * (TDim) + k, j * (TDim) + l);
                        }
                    }
                }
            }

            Kvu_nitsche_alpha_self_large = prod(Kvu_nitsche_alpha_self_large, trans(Tdublicate));
            Kvu_nitsche_alpha_self_large = prod(Tdublicate, Kvu_nitsche_alpha_self_large);

            Kvu_nitsche_alpha_neighbor_large = prod(Kvu_nitsche_alpha_neighbor_large, trans(Tdublicate_neighbor));
            Kvu_nitsche_alpha_neighbor_large = prod(Tdublicate, Kvu_nitsche_alpha_neighbor_large);

            for (int i = 0; i < 18; i++)
            {
                for (int j = 0; j < 18; j++)
                {
                    lhs(i, j) += Kvu_nitsche_large(i, j);
                    //        KRATOS_WATCH(Kvu_nitsche_large);
                    lhs(i, j) += Kvp_nitsche_large(i, j);
                    //        KRATOS_WATCH(Kvp_nitsche_large);

                    lhs(i, j) += Kuv_nitsche_self_large(i, j);
                    //        KRATOS_WATCH(Kuv_nitsche_self_large);
                    lhs(i, j) += Kuv_nitsche_neighbor_large(i, j);
                    //        KRATOS_WATCH(Kuv_nitsche_neighbor_large);

                    lhs(i, j) += Kpv_nitsche_self_large(i, j);
                    //        KRATOS_WATCH(Kpv_nitsche_self_large);
                    lhs(i, j) += Kpv_nitsche_neighbor_large(i, j);
                    //        KRATOS_WATCH(Kpv_nitsche_neighbor_large);

                    lhs(i, j) += Kvu_nitsche_alpha_self_large(i, j);
                    //        KRATOS_WATCH(Kvu_nitsche_alpha_self_large);
                    lhs(i, j) += Kvu_nitsche_alpha_neighbor_large(i, j);
                    //        KRATOS_WATCH(Kvu_nitsche_alpha_neighbor_large);
                }
            }
        }
    }

} // Namespace Kratos
