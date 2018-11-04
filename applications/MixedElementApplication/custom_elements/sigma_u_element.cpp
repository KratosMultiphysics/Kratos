//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/sigma_u_element.h"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "mixedelement_application.h"
#include "utilities/geometry_utilities.h"
//#include <omp.h>

//#define SYMM_FORM

namespace Kratos
{
//************************************************************************************
//************************************************************************************

SigmaUElement::SigmaUElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

SigmaUElement::SigmaUElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

Element::Pointer SigmaUElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SigmaUElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

SigmaUElement::~SigmaUElement()
{
}

std::string SigmaUElement::Info() const
{
    std::stringstream buffer;
    buffer << "Linear Element" << std::endl;
    return buffer.str();
}


//************************************************************************************
//************************************************************************************

void SigmaUElement::Initialize()
{
    KRATOS_TRY

    unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int nintegration_points = GetGeometry().size(); //we will use a nodal integration rule
//        KRATOS_WATCH("ln114")
    if (dim == 2)
    {
        boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX;
        array_1d<double, 3 > N;

        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, mArea0);
        mDN_DX.resize(3, 2, false);
        noalias(mDN_DX) = DN_DX;
    }
    else
    {
        boost::numeric::ublas::bounded_matrix<double, 4,3 > DN_DX;
        array_1d<double, 4 > N;

        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, mArea0);
        mDN_DX.resize(4,3, false);
        noalias(mDN_DX) = DN_DX;    
    }


    //compute and save original area and shape functions

    //Constitutive Law initialisation
    if (mConstitutiveLawVector.size() != nintegration_points)
    {
        mConstitutiveLawVector.resize(nintegration_points);
    }
    InitializeMaterial();
//        KRATOS_WATCH(mConstitutiveLawVector)
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SigmaUElement::CalculateAll(MatrixType& rLeftHandSideMatrix,
                                 VectorType& rRightHandSideVector,
                                 ProcessInfo& rCurrentProcessInfo,
                                 bool CalculateStiffnessMatrixFlag,
                                 bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    const unsigned int nnodes = GetGeometry().size();

    unsigned int StrainSize;
    if (dim == 2)
        StrainSize = 3;
    else
        StrainSize = 6;

    //resizing LHS and RHS to the correct size
    const unsigned int var_block_size = dim + StrainSize;
    const unsigned int MatSize = var_block_size*nnodes;
    if(rLeftHandSideMatrix.size1() != MatSize)
        rLeftHandSideMatrix.resize(MatSize,MatSize);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize,MatSize);

    if(rRightHandSideVector.size() != MatSize)
        rRightHandSideVector.resize(MatSize);
    noalias(rRightHandSideVector) = ZeroVector(MatSize);

    //compute B (it is constant over the whole element)
    Matrix B(StrainSize,dim*nnodes);
    CalculateB(B, mDN_DX, StrainSize);
//KRATOS_WATCH(StrainSize)
    Vector discontinuous_strain(StrainSize);
    Vector u_list(dim*nnodes);
    GetDisplacements(u_list,dim);
    noalias(discontinuous_strain) = prod(B,u_list);

//KRATOS_WATCH(discontinuous_strain)
    //nodal integration weight is used
    const double weight = 1.0/static_cast<double>(dim+1);
    Matrix C(StrainSize,StrainSize);
    Matrix Cavg(StrainSize,StrainSize,0.0);
    Vector eps_h(StrainSize);
    Vector avg_eps_h(StrainSize,0.0);
//        Vector stress_vector(StrainSize);
    Matrix block11(StrainSize*nnodes,StrainSize*nnodes,0.0);
    Matrix block12(StrainSize*nnodes,dim*nnodes,0.0);
    Matrix block12_aux(StrainSize,dim*nnodes,0.0);
    Matrix block21(dim*nnodes,StrainSize*nnodes,0.0);
    Matrix block22(dim*nnodes,dim*nnodes);
    Vector N(nnodes);
    Matrix F;
    std::vector< Vector > stresses_vector(nnodes);
    Vector stress_discontinuous(StrainSize,0.0);
    Vector aux_stress_discontinuous(StrainSize);
    
    std::vector< Vector > nodal_disc_stress_vector(nnodes);
    std::vector< Vector > nodal_eps_vector(nnodes);
    
    Matrix identity(StrainSize,StrainSize);
    identity = ZeroMatrix(StrainSize,StrainSize);
    for(unsigned int k=0; k<StrainSize; k++)
      identity(k,k) = 1.0;
    


    for(unsigned int igauss=0; igauss<GetGeometry().size(); igauss++)
    {
        //get nodal strain
        GetNodalVariable(eps_h, igauss, dim);
	
	nodal_eps_vector[igauss].resize(StrainSize,false);
	noalias(nodal_eps_vector[igauss]) = eps_h;

        noalias(avg_eps_h ) += weight*eps_h;

        noalias(N)=ZeroVector(nnodes);
        N[igauss] = 1.0;

        //compute continuous stress
        stresses_vector[igauss].resize(StrainSize,false);
        mConstitutiveLawVector[igauss]->CalculateMaterialResponse(eps_h,F,stresses_vector[igauss],C,rCurrentProcessInfo,GetProperties(),GetGeometry(),N,true,1,true);

        //average constitutitve law on the center point
        noalias(Cavg) += weight * C;
	
	//auxiliary stress computed from discontinuous strain but with nodal C
	nodal_disc_stress_vector[igauss].resize(StrainSize,false);
	noalias(nodal_disc_stress_vector[igauss]) = prod(identity,discontinuous_strain);

#ifndef SYMM_FORM
        //compute block 11
        for(unsigned int k=0; k<StrainSize; k++)
            for(unsigned int l=0; l<StrainSize; l++)
                block11(igauss*StrainSize+k,igauss*StrainSize+l) += mArea0*weight*identity(k,l);
#else
        //compute block 11
        for(unsigned int k=0; k<StrainSize; k++)
            for(unsigned int l=0; l<StrainSize; l++)
                block11(igauss*StrainSize+k,igauss*StrainSize+l) += mArea0*weight*C(k,l);
		
#endif
	//compute s-u block (12)
        noalias(block12_aux) = mArea0*weight * prod(C,B);

#ifndef SYMM_FORM
        //write block12_aux into block12
        for(unsigned int k=0; k<StrainSize; k++)
            for(unsigned int l=0; l<dim*nnodes; l++)
                block12(igauss*StrainSize+k,l) = mArea0*weight*B(k,l);
#else
        //write block12_aux into block12
        for(unsigned int k=0; k<StrainSize; k++)
            for(unsigned int l=0; l<dim*nnodes; l++)
                block12(igauss*StrainSize+k,l) = block12_aux(k,l);
		
#endif	    
        //write block12_aux into block12
        for(unsigned int k=0; k<StrainSize; k++)
            for(unsigned int l=0; l<dim*nnodes; l++)
                block21(l,igauss*StrainSize+k) = block12_aux(k,l);	    
    }

    //compute discontinuous stressavg_eps_h
//        noalias(stress_discontinuous) = prod(Cavg,avg_eps_h);
    noalias(stress_discontinuous) = prod(Cavg,discontinuous_strain);
////        mpDiscontinuousConstitutiveLaw->CalculateMaterialResponse(discontinuous_strain,F,stress_discontinuous,C,rCurrentProcessInfo,GetProperties(),GetGeometry(),N,true,1,true);
////noalias(Cavg) = C;
//        KRATOS_WATCH(234)
//     noalias(block21) = trans(block12);
    
//KRATOS_WATCH(Cavg);
    //compute block 22 (note that Cavg is used here)
    Cavg *= mArea0;
    Matrix tmp = prod(Cavg,B);
    noalias(block22) = prod(trans(B), tmp);

    //hystorical variables are saved according to avg_eps_h instead of discontinuous_strain
//        mpDiscontinuousConstitutiveLaw->CalculateMaterialResponse(avg_eps_h,F,aux_stress_discontinuous,C,rCurrentProcessInfo,GetProperties(),GetGeometry(),N,true,1,true);

    //compute tau
    double tau=0.01;

//        double L = rCurrentProcessInfo[DIAMETER];
//        double A = GetGeometry().Area();
//        const double he = sqrt(2.0*A);
//        const double cepsilon = 50.0;
//        double tau = cepsilon*he/L;

    //assemble blocks (multiplication by tau is done here)
    for(unsigned int i=0; i<GetGeometry().size(); i++)
    {
        int i_sigma_block_start = i*var_block_size;
        int i_u_block_start = i_sigma_block_start + StrainSize;
        for(unsigned int j=0; j<GetGeometry().size(); j++)
        {
            int j_sigma_block_start = j*var_block_size;
            int j_u_block_start = j_sigma_block_start + StrainSize;

            //assembling simga rows
            for(unsigned int k=0; k<StrainSize; k++)
            {
                //assembling sigma columns
                for(unsigned int l=0; l<StrainSize; l++)
                    rLeftHandSideMatrix(i_sigma_block_start+k,j_sigma_block_start+l) = -(1-tau)*block11(i*StrainSize+k,j*StrainSize+l);

                //assembling u columns
                for(unsigned int l=0; l<dim; l++)
                    rLeftHandSideMatrix(i_sigma_block_start+k,j_u_block_start+l) = (1-tau)*block12(i*StrainSize+k,j*dim+l);
            }

            //assembling u rows
            for(unsigned int k=0; k<dim; k++)
            {
                //assembling sigma columns
                for(unsigned int l=0; l<StrainSize; l++)
                    rLeftHandSideMatrix(i_u_block_start+k,j_sigma_block_start+l) = (1-tau)*block21(i*dim+k,j*StrainSize+l);

                //assembling u columns
                for(unsigned int l=0; l<dim; l++)
                    rLeftHandSideMatrix(i_u_block_start+k,j_u_block_start+l) = tau*block22(i*dim+k,j*dim+l);
            }
        }
    }

    //finish the formation of the residua
    Vector tmp1(dim*nnodes);
    Vector stress_stabilized = tau*stress_discontinuous;
    for(unsigned int i=0; i<GetGeometry().size(); i++)
        noalias(stress_stabilized) += ((1.0-tau)*weight)*stresses_vector[i] ;
    stress_stabilized *= mArea0;

//        Vector tmp1(dim*nnodes);
//        Vector stress_stabilized(StrainSize,0.0);
//        for(unsigned int i=0; i<GetGeometry().size(); i++)
//            noalias(stress_stabilized) += ((1.0-tau)*weight)*stresses_vector[i] ;
//        stress_stabilized *= mArea0;

    noalias(tmp1) = -prod(trans(B),stress_stabilized);

//        KRATOS_WATCH(tmp1)
//        KRATOS_WATCH(tau*prod(trans(B),stress_discontinuous))
//        KRATOS_WATCH(tau* prod(block22,u_list))
//        Vector all_stresses(nnodes*StrainSize);
//
//        for(unsigned int i=0; i<GetGeometry().size(); i++)
//        {
//            noalias(tmp1) -= ((1.0-tau)) * prod(block21,all_stresses);
//        }

    //GetS(all_stresses,dim);
    //noalias(tmp1) -= ((1.0-tau)) * prod(block21,all_stresses);

    for(unsigned int i=0; i<GetGeometry().size(); i++)
    {
        int i_sigma_block_start = i*var_block_size;
        int i_u_block_start = i_sigma_block_start + StrainSize;
        const double coeff = mArea0*weight*(1.0-tau);
        for(unsigned int l=0; l<StrainSize; l++)
        {
#ifndef SYMM_FORM
            rRightHandSideVector[i_sigma_block_start+l] = coeff*(nodal_eps_vector[i][l] - discontinuous_strain[l] /*stress_discontinuous[l]*/ );
#else
	    rRightHandSideVector[i_sigma_block_start+l] = coeff*(stresses_vector[i][l] - stress_discontinuous[l] /*stress_discontinuous[l]*/ );    
#endif
        }

        for(unsigned int l=0; l<dim; l++)
            rRightHandSideVector[i_u_block_start+l] = /*body_force[l]*mArea0*/ tmp1[i*dim+l];
    }

//        Vector tmp1(MatSize);
//        GetValuesVector(tmp1,0);
//
//        Vector aux = prod(rLeftHandSideMatrix,tmp1);
//         for(unsigned int i=0; i<GetGeometry().size(); i++)
//        {
//            int i_sigma_block_start = i*var_block_size;
//            int i_u_block_start = i_sigma_block_start + StrainSize;
//            for(unsigned int l=0; l<dim; l++)
//                aux[i_u_block_start+l] = 0;
//        }
//        noalias(rRightHandSideVector) -= aux;

//        KRATOS_WATCH(rLeftHandSideMatrix);
//                KRATOS_WATCH(rRightHandSideVector);

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************

void SigmaUElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo & rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
}

//************************************************************************************
//************************************************************************************

void SigmaUElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo & rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);


}



////************************************************************************************
////************************************************************************************

void SigmaUElement::InitializeSolutionStep(ProcessInfo & CurrentProcessInfo)
{
//        double factor = 1.0 / static_cast<double> (geom.size());
    Vector N(GetGeometry().size());

    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
    {
        //nodal integration rule is used
        for (unsigned int j = 0; j < N.size(); j++)
        {
            if (i == j) N[j] = 1.0;
            else N[j] = 0.0;
        }
        mConstitutiveLawVector[i]->InitializeSolutionStep(GetProperties(),
                GetGeometry(), N,
                CurrentProcessInfo);
    }

//        for (unsigned int j = 0; j < N.size(); j++) N[j] = 1.0/3.0;
//        mpDiscontinuousConstitutiveLaw->InitializeSolutionStep(GetProperties(),GetGeometry(),N,CurrentProcessInfo);
}

////************************************************************************************
////************************************************************************************

void SigmaUElement::FinalizeSolutionStep(ProcessInfo & CurrentProcessInfo)
{
    Vector N(GetGeometry().size());
    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
    {
        for (unsigned int j = 0; j < N.size(); j++)
        {
            if (i == j) N[j] = 1.0;
            else N[j] = 0.0;
        }

        mConstitutiveLawVector[i]->FinalizeSolutionStep(GetProperties(),
                GetGeometry(),
                N,
                CurrentProcessInfo);
    }



//        for (unsigned int j = 0; j < N.size(); j++) N[j] = 1.0/3.0;
//        mpDiscontinuousConstitutiveLaw->FinalizeSolutionStep(GetProperties(),GetGeometry(),N,CurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************

void SigmaUElement::InitializeMaterial()
{
    KRATOS_TRY

    Vector N(GetGeometry().size());

    if (GetProperties()[CONSTITUTIVE_LAW] != NULL)
    {
        for (unsigned int i = 0; i < GetGeometry().size(); i++)
        {
            //nodal integration rule is used
            for (unsigned int j = 0; j < N.size(); j++)
            {
                if (i == j) N[j] = 1.0;
                else N[j] = 0.0;
            }

            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();

            mConstitutiveLawVector[i]->InitializeMaterial(GetProperties(), GetGeometry(), N);
        }

//            for (unsigned int j = 0; j < N.size(); j++) N[j] = 1.0/3.0;
//            mpDiscontinuousConstitutiveLaw= GetProperties()[CONSTITUTIVE_LAW]->Clone();
//            mpDiscontinuousConstitutiveLaw->InitializeMaterial(GetProperties(), GetGeometry(), N);

    }
    else
        KRATOS_THROW_ERROR(std::logic_error, "a constitutive law needs to be specified for the element with ID ", this->Id())

        KRATOS_CATCH("")
    }

//************************************************************************************
//************************************************************************************

void SigmaUElement::CalculateB(
    Matrix& B,
    Matrix& DN_DX,
    unsigned int StrainSize)
{
    KRATOS_TRY
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();


    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        unsigned int index = dimension * i;

        if (dimension == 2)
        {
            B(0, index + 0) = DN_DX(i, 0);
            B(0, index + 1) = 0.0;
            B(1, index + 0) = 0.0;
            B(1, index + 1) = DN_DX(i, 1);
            B(2, index + 0) = DN_DX(i, 1);
            B(2, index + 1) = DN_DX(i, 0);
        }
        else
        {
	    //xx
            B(0, index + 0) = DN_DX(i, 0);
            B(0, index + 1) = 0.0;
            B(0, index + 2) = 0.0;
	    
	    //yy
            B(1, index + 0) = 0.0;
            B(1, index + 1) = DN_DX(i, 1);
	    B(1, index + 2) = 0.0;
	    
	    //zz
            B(2, index + 0) = 0.0;
            B(2, index + 1) = 0.0;
	    B(2, index + 2) = DN_DX(i, 2);	    
	    
	    //eps_xy
            B(3, index + 0) = DN_DX(i, 1);
            B(3, index + 1) = DN_DX(i, 0);
	    B(3, index + 2) = 0.0;
	    
	    //eps_ xz
            B(4, index + 0) = DN_DX(i, 2);
            B(4, index + 1) = 0.0;
	    B(4, index + 2) = DN_DX(i, 0);
	    
	    //eps yz
            B(5, index + 0) = 0.0;
            B(5, index + 1) = DN_DX(i, 2);
	    B(5, index + 2) = DN_DX(i, 1);	    
        }
    }

    KRATOS_CATCH("")
}



//************************************************************************************
//************************************************************************************

void SigmaUElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo & CurrentProcessInfo)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();

    if (dim == 2)
    {
        unsigned int block_size = 5;
        unsigned int MatSize = number_of_nodes * block_size;
        if (rResult.size() != MatSize) rResult.resize(MatSize, false);

        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            unsigned int index = i * block_size;
            rResult[index]     = GetGeometry()[i].GetDof(SX).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(SY).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(SXY).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 4] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        }
    }
    else
    {
        unsigned int block_size = 9;
        unsigned int MatSize = number_of_nodes * block_size;
        if (rResult.size() != MatSize) rResult.resize(MatSize, false);

        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            unsigned int index = i * block_size;
            rResult[index]     = GetGeometry()[i].GetDof(SX).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(SY).EquationId();
	    rResult[index + 2] = GetGeometry()[i].GetDof(SZ).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof(SXY).EquationId();
	    rResult[index + 4] = GetGeometry()[i].GetDof(SXZ).EquationId();
	    rResult[index + 5] = GetGeometry()[i].GetDof(SYZ).EquationId();
            rResult[index + 6] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 7] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
	    rResult[index + 8] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
        }
    }
}

//************************************************************************************
//************************************************************************************

void SigmaUElement::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo & CurrentProcessInfo)
{
    ElementalDofList.resize(0);
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();

    if (dim == 2)
    {
        for (unsigned int i = 0; i < GetGeometry().size(); i++)
        {
            ElementalDofList.push_back(GetGeometry()[i].pGetDof(SX));
            ElementalDofList.push_back(GetGeometry()[i].pGetDof(SY));
            ElementalDofList.push_back(GetGeometry()[i].pGetDof(SXY));
            ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
        }
    }
    else
    {
        for (unsigned int i = 0; i < GetGeometry().size(); i++)
        {
            ElementalDofList.push_back(GetGeometry()[i].pGetDof(SX));
            ElementalDofList.push_back(GetGeometry()[i].pGetDof(SY));
	    ElementalDofList.push_back(GetGeometry()[i].pGetDof(SZ));
            ElementalDofList.push_back(GetGeometry()[i].pGetDof(SXY));
	    ElementalDofList.push_back(GetGeometry()[i].pGetDof(SXZ));
	    ElementalDofList.push_back(GetGeometry()[i].pGetDof(SYZ));
            ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
	    ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }
    }
}

//************************************************************************************
//************************************************************************************

void SigmaUElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo & rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SigmaUElement::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo & rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SigmaUElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& Output, const ProcessInfo & rCurrentProcessInfo)
{
    if (Output.size() != mConstitutiveLawVector.size())
        Output.resize(mConstitutiveLawVector.size());

    for (unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++)
        Output[ii] = mConstitutiveLawVector[ii]->GetValue(rVariable, Output[ii]);
}

//************************************************************************************
//************************************************************************************

void SigmaUElement::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& Output, const ProcessInfo & rCurrentProcessInfo)
{
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();

    unsigned int StrainSize;

    if (dim == 2)
        StrainSize = 3;
    else
        StrainSize = 6;

    Vector StrainVector(StrainSize);

    if (Output.size() != mConstitutiveLawVector.size())
        Output.resize(mConstitutiveLawVector.size());

    for (unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++)
        Output[ii] = mConstitutiveLawVector[ii]->GetValue(rVariable, Output[ii]);

}

//************************************************************************************
//************************************************************************************

void SigmaUElement::CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& Output, const ProcessInfo & rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SigmaUElement::SetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo & rCurrentProcessInfo)
{
    for (unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++)
    {
        mConstitutiveLawVector[PointNumber]->SetValue(rVariable,
                rValues[PointNumber], rCurrentProcessInfo);
    }

}


//************************************************************************************
//************************************************************************************

void SigmaUElement::SetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo & rCurrentProcessInfo)
{
    for (unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++)
    {
        mConstitutiveLawVector[PointNumber]->SetValue(rVariable,
                rValues[PointNumber], rCurrentProcessInfo);
    }

}

//************************************************************************************
//************************************************************************************

void SigmaUElement::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo & rCurrentProcessInfo)
{
    if (rValues.size() != mConstitutiveLawVector.size())
        rValues.resize(mConstitutiveLawVector.size(), false);

    for (unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++)
        rValues[ii] = mConstitutiveLawVector[ii]->GetValue(rVariable, rValues[ii]);
}


//************************************************************************************
//************************************************************************************

void SigmaUElement::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo & rCurrentProcessInfo)
{

    for (unsigned int PointNumber = 0;
            PointNumber < mConstitutiveLawVector.size();
            PointNumber++)
    {
        rValues[PointNumber] =
            mConstitutiveLawVector[PointNumber]->GetValue(rVariable, rValues[PointNumber]);
    }

}

//************************************************************************************
//************************************************************************************

void SigmaUElement::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues, const ProcessInfo & rCurrentProcessInfo)
{
    CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************

void SigmaUElement::GetValuesVector(Vector& values, int Step)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();

    if (dim == 2)
    {
        unsigned int block_size = 5;
        unsigned int MatSize = number_of_nodes * block_size;
        if (values.size() != MatSize) values.resize(MatSize, false);

        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            unsigned int index = i * block_size;
            values[index] = GetGeometry()[i].FastGetSolutionStepValue(SX, Step);
            values[index + 1] = GetGeometry()[i].FastGetSolutionStepValue(SY, Step);
            values[index + 2] = GetGeometry()[i].FastGetSolutionStepValue(SXY, Step);
            values[index + 3] = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_X, Step);
            values[index + 4] = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Y, Step);
        }
    }
    else
    {
        unsigned int block_size = 5;
        unsigned int MatSize = number_of_nodes * block_size;
        if (values.size() != MatSize) values.resize(MatSize, false);

        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            unsigned int index = i * block_size;
            values[index] = GetGeometry()[i].FastGetSolutionStepValue(SX, Step);
            values[index + 1] = GetGeometry()[i].FastGetSolutionStepValue(SY, Step);
	    values[index + 2] = GetGeometry()[i].FastGetSolutionStepValue(SZ, Step);
	    values[index + 3] = GetGeometry()[i].FastGetSolutionStepValue(SXY, Step);
	    values[index + 4] = GetGeometry()[i].FastGetSolutionStepValue(SXZ, Step);
            values[index + 5] = GetGeometry()[i].FastGetSolutionStepValue(SYZ, Step);
            values[index + 6] = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_X, Step);
            values[index + 7] = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Y, Step);
	    values[index + 8] = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Z, Step);
        }
    }
}


//************************************************************************************
//************************************************************************************

void SigmaUElement::GetFirstDerivativesVector(Vector& values, int Step)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();

    if (dim == 2)
    {
        unsigned int block_size = 5;
        unsigned int MatSize = number_of_nodes * block_size;
        if (values.size() != MatSize) values.resize(MatSize, false);
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            unsigned int index = i * block_size;
            values[index]     = 0.0;
            values[index + 1] = 0.0;
            values[index + 2] = 0.0;
            values[index + 3] = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X, Step);
            values[index + 4] = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y, Step);
        }
    }
    else
    {
        unsigned int block_size = 9;
        unsigned int MatSize = number_of_nodes * block_size;
        if (values.size() != MatSize) values.resize(MatSize, false);
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            unsigned int index = i * block_size;
            values[index]     = 0.0;
            values[index + 1] = 0.0;
            values[index + 2] = 0.0;
	    values[index + 3] = 0.0;
	    values[index + 4] = 0.0;
	    values[index + 5] = 0.0;
	    
	    const array_1d<double,3>& vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
            values[index + 6] = vel[0];
            values[index + 7] = vel[1];
	    values[index + 8] = vel[2];
        }
    }


}

//************************************************************************************
//************************************************************************************

void SigmaUElement::GetSecondDerivativesVector(Vector& values, int Step)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();

    if (dim == 2)
    {
        unsigned int block_size = 5;
        unsigned int MatSize = number_of_nodes * block_size;
        if (values.size() != MatSize) values.resize(MatSize, false);
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            unsigned int index = i * block_size;
            values[index]     = 0.0;
            values[index + 1] = 0.0;
            values[index + 2] = 0.0;
            values[index + 3] = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION_X, Step);
            values[index + 4] = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION_Y, Step);
        }
    }
    else
    {
        unsigned int block_size = 9;
        unsigned int MatSize = number_of_nodes * block_size;
        if (values.size() != MatSize) values.resize(MatSize, false);
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            unsigned int index = i * block_size;
            values[index]     = 0.0;
            values[index + 1] = 0.0;
            values[index + 2] = 0.0;
	    values[index + 3] = 0.0;
	    values[index + 4] = 0.0;
	    values[index + 5] = 0.0;
	    
	    const array_1d<double,3>& acc = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
            values[index + 6] = acc[0];
            values[index + 7] = acc[1];
	    values[index + 8] = acc[2];
        }
    }
}

//************************************************************************************
//************************************************************************************

void SigmaUElement::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo & rCurrentProcessInfo)
{
}


//************************************************************************************
//************************************************************************************
void SigmaUElement::GetNodalVariable(Vector& value, unsigned int i, const unsigned int dim)
{
    if(dim == 2)
    {
        value[0] = GetGeometry()[i].FastGetSolutionStepValue(SX);
        value[1] = GetGeometry()[i].FastGetSolutionStepValue(SY);
        value[2] = GetGeometry()[i].FastGetSolutionStepValue(SXY);
    }
    else
    {
        value[0] = GetGeometry()[i].FastGetSolutionStepValue(SX);
        value[1] = GetGeometry()[i].FastGetSolutionStepValue(SY);
        value[2] = GetGeometry()[i].FastGetSolutionStepValue(SZ);
        value[3] = GetGeometry()[i].FastGetSolutionStepValue(SXY);
        value[4] = GetGeometry()[i].FastGetSolutionStepValue(SXZ);
        value[5] = GetGeometry()[i].FastGetSolutionStepValue(SYZ);
    }

}

//************************************************************************************
//************************************************************************************
void SigmaUElement::GetDisplacements(Vector& value, const unsigned int dim)
{
    if(value.size() != GetGeometry().size()*dim)
        value.resize(GetGeometry().size()*dim);

    unsigned int counter=0;
    for(unsigned int i=0; i<GetGeometry().size(); i++)
        for(unsigned int j=0; j<dim; j++)
        {
            const array_1d<double,3>& disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
            value[counter++] = disp[j];
        }
}

void SigmaUElement::GetS(Vector& value, const unsigned int dim)
{
    if(dim == 2)
    {
        const unsigned int StrainSize = 3;
        if(value.size() != StrainSize*GetGeometry().size())
            value.resize(StrainSize*GetGeometry().size(),false);

        unsigned int counter=0;
        for(unsigned int i=0; i<GetGeometry().size(); i++)
        {
            value[counter++] = GetGeometry()[i].FastGetSolutionStepValue(SX);
            value[counter++] = GetGeometry()[i].FastGetSolutionStepValue(SY);
            value[counter++] = GetGeometry()[i].FastGetSolutionStepValue(SXY);
        }
    }
    else
    {
        const unsigned int StrainSize = 6;
        if(value.size() != StrainSize*GetGeometry().size())
            value.resize(StrainSize*GetGeometry().size(),false);

        unsigned int counter=0;
        for(unsigned int i=0; i<GetGeometry().size(); i++)
        {
            value[counter++] = GetGeometry()[i].FastGetSolutionStepValue(SX);
            value[counter++] = GetGeometry()[i].FastGetSolutionStepValue(SY);
	    value[counter++] = GetGeometry()[i].FastGetSolutionStepValue(SZ);
            value[counter++] = GetGeometry()[i].FastGetSolutionStepValue(SXY);
            value[counter++] = GetGeometry()[i].FastGetSolutionStepValue(SXZ);
            value[counter++] = GetGeometry()[i].FastGetSolutionStepValue(SYZ);
        }
    }
}

} // Namespace Kratos


