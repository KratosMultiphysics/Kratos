//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hencky_plastic_3d_law.hpp"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"
#include "particle_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************
HenckyElasticPlastic3DLaw::HenckyElasticPlastic3DLaw()
    : HyperElastic3DLaw()
{

}



HenckyElasticPlastic3DLaw::HenckyElasticPlastic3DLaw(MPMFlowRulePointer pMPMFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : HyperElastic3DLaw( )
{

}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyElasticPlastic3DLaw::HenckyElasticPlastic3DLaw(const HenckyElasticPlastic3DLaw&  rOther)
    : HyperElastic3DLaw(rOther)
    ,mElasticLeftCauchyGreen(rOther.mElasticLeftCauchyGreen)
    ,mpYieldCriterion(rOther.mpYieldCriterion)
    ,mpHardeningLaw(rOther.mpHardeningLaw)
{
    mpMPMFlowRule       = rOther.mpMPMFlowRule->Clone();
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyElasticPlastic3DLaw::Clone() const
{
    HenckyElasticPlastic3DLaw::Pointer p_clone(new HenckyElasticPlastic3DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyElasticPlastic3DLaw::~HenckyElasticPlastic3DLaw()
{
}
//***********************HAS : DOUBLE - VECTOR - MATRIX*******************************
//************************************************************************************

bool HenckyElasticPlastic3DLaw::Has( const Variable<double>& rThisVariable )
{
    return false;
}

bool HenckyElasticPlastic3DLaw::Has( const Variable<Vector>& rThisVariable )
{
    return false;
}

bool HenckyElasticPlastic3DLaw::Has( const Variable<Matrix>& rThisVariable )
{
    return false;
}


//***********************GET VALUE: DOUBLE - VECTOR - MATRIX**************************
//************************************************************************************

double& HenckyElasticPlastic3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if (rThisVariable==PLASTIC_STRAIN)
    {
        const MPMFlowRule::InternalVariables& InternalVariables = mpMPMFlowRule->GetInternalVariables();
        rValue=InternalVariables.EquivalentPlasticStrain;
    }


    if (rThisVariable==MIU)
    {

        rValue=mPlasticRegion;
    }
    return( rValue );
}

Vector& HenckyElasticPlastic3DLaw::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    return( rValue );
}

Matrix& HenckyElasticPlastic3DLaw::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    return( rValue );
}


//***********************SET VALUE: DOUBLE - VECTOR - MATRIX**************************
//************************************************************************************


void HenckyElasticPlastic3DLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
        const ProcessInfo& rCurrentProcessInfo )
{
    //if (rThisVariable == DETERMINANT_F)
    //{
    //mDeterminantF0 = rValue;
    //}

    //if (rThisVariable == PLASTIC_STRAIN)
    //{
    //mpMPMFlowRule->SetInternalVariables().EquivalentPlasticStrain = rValue;
    //}

    //if (rThisVariable == DELTA_PLASTIC_STRAIN)
    //{
    //mpMPMFlowRule->SetInternalVariables().DeltaPlasticStrain = rValue;
    //}

}

void HenckyElasticPlastic3DLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
        const ProcessInfo& rCurrentProcessInfo )
{

}

void HenckyElasticPlastic3DLaw::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
        const ProcessInfo& rCurrentProcessInfo )
{

}


//************************************************************************************
//************************************************************************************


void HenckyElasticPlastic3DLaw::InitializeMaterial(const Properties& rProps,
        const GeometryType& rGeom,
        const Vector& rShapeFunctionsValues)
{
    mDeterminantF0                = 1;
    mInverseDeformationGradientF0 = identity_matrix<double> (3);
    mElasticLeftCauchyGreen       = identity_matrix<double> (3);

    mPlasticRegion = 0;

    mpMPMFlowRule->InitializeMaterial(mpYieldCriterion, mpHardeningLaw, rProps);
}
//************************************************************************************
//************************************************************************************


void HenckyElasticPlastic3DLaw::InitializeSolutionStep( const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry, //this is just to give the array of nodes
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo)
{
    mpHardeningLaw->SetProperties(rMaterialProperties);
}

//************************************************************************************
//************************************************************************************


void HenckyElasticPlastic3DLaw::FinalizeSolutionStep( const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry, //this is just to give the array of nodes
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo)
{

}

//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

//*****************************MATERIAL RESPONSES*************************************
//************************************************************************************


void HenckyElasticPlastic3DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
{

    //-----------------------------//

    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);

    //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const ProcessInfo& CurrentProcessInfo = rValues.GetProcessInfo();
    //const Properties& MaterialProperties   = rValues.GetMaterialProperties();

    const Matrix&   DeformationGradientF   = rValues.GetDeformationGradientF();
    const double&   DeterminantF           = rValues.GetDeterminantF();

    const GeometryType&  DomainGeometry    = rValues.GetElementGeometry ();
    const Vector&        ShapeFunctions    = rValues.GetShapeFunctionsValues ();

    Vector& StrainVector                   = rValues.GetStrainVector();
    Vector& StressVector                   = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix             = rValues.GetConstitutiveMatrix();

    //-----------------------------//

    //0.- Initialize parameters
    MaterialResponseVariables ElasticVariables;
    ElasticVariables.Identity = identity_matrix<double> ( 3 );

    ElasticVariables.SetElementGeometry(DomainGeometry);
    ElasticVariables.SetShapeFunctionsValues(ShapeFunctions);

    MPMFlowRule::RadialReturnVariables ReturnMappingVariables;
    //ReturnMappingVariables.initialize(); //it has to be called at the start
    ReturnMappingVariables.clear();

    // Initialize variables from the process information
    ReturnMappingVariables.DeltaTime = CurrentProcessInfo[DELTA_TIME];

    if(CurrentProcessInfo[IMPLEX] == 1)
        ReturnMappingVariables.Options.Set(MPMFlowRule::IMPLEX_ACTIVE,true);
    else
        ReturnMappingVariables.Options.Set(MPMFlowRule::IMPLEX_ACTIVE,false);


    //1.- Lame constants
    // const double& YoungModulus       = MaterialProperties[YOUNG_MODULUS];
    // const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];

    //ElasticVariables.LameLanda      = (YoungModulus*PoissonCoefficient)/((1.0+PoissonCoefficient)*(1.0-2.0*PoissonCoefficient));
    //ElasticVariables.LameMu         =  YoungModulus/(2.0*(1.0+PoissonCoefficient));

    // lluis
    // ReturnMappingVariables.YoungModulus          =  MaterialProperties[YOUNG_MODULUS];
    // ReturnMappingVariables.PoissonCoefficient    =  MaterialProperties[POISSON_RATIO];


    //2.-Determinant of the Total Deformation Gradient
    ElasticVariables.DeterminantF = DeterminantF;

    //3.-Compute Incremental DeformationGradient (in 3D)
    ElasticVariables.DeformationGradientF = DeformationGradientF;
    ElasticVariables.DeformationGradientF = Transform2DTo3D(ElasticVariables.DeformationGradientF);
    //std::cout<<" ElasticVariables.DeformationGradientF total "<<ElasticVariables.DeformationGradientF<<std::endl;
    //std::cout<<" mInverseDeformationGradientF0 "<<mInverseDeformationGradientF0<<std::endl;
    ElasticVariables.DeformationGradientF = prod(ElasticVariables.DeformationGradientF,mInverseDeformationGradientF0);
    //std::cout<<" ElasticVariables.DeformationGradientF partial "<<ElasticVariables.DeformationGradientF<<std::endl;
    //4.-Left Cauchy Green tensor b: (stored in the CauchyGreenMatrix)
    //mElasticLeftCauchyGreen how is saved?


    ElasticVariables.CauchyGreenMatrix = prod(mElasticLeftCauchyGreen,trans(ElasticVariables.DeformationGradientF));
    ElasticVariables.CauchyGreenMatrix = prod(ElasticVariables.DeformationGradientF,ElasticVariables.CauchyGreenMatrix);


    //5.-Almansi Strain:
    if(Options.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ))
    {
        // e= 0.5*(1-invbT*invb)
        this->CalculateAlmansiStrain(ElasticVariables.CauchyGreenMatrix, StrainVector);

    }


    //6.-Calculate Total Kirchhoff stress

    if( Options.Is(ConstitutiveLaw::COMPUTE_STRESS ) || Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {

        Matrix StressMatrix     = ZeroMatrix(3,3);
        Vector HenckyMainStrainVector = ZeroVector(3);


        //std::cout<<" mElasticLeftCauchyGreen "<<mElasticLeftCauchyGreen<<std::endl;
        //std::cout<<" ElasticVariables.DeformationGradientF "<<ElasticVariables.DeformationGradientF<<std::endl;
        //std::cout<<" ElasticVariables.CauchyGreenMatrix "<<ElasticVariables.CauchyGreenMatrix<<std::endl;

        this->CalculateHenckyMainStrain(ElasticVariables.CauchyGreenMatrix, ReturnMappingVariables, HenckyMainStrainVector);
        // mpFlowRule->CalculateReturnMapping( ReturnMappingVariables, StressMatrix, HenckyMainStrain );
        ReturnMappingVariables.StrainMatrix = ZeroMatrix(3,3);
        ReturnMappingVariables.TrialIsoStressMatrix = ZeroMatrix(3,3);
        Matrix HenckyMainStrainMatrix = ZeroMatrix(3,3);
        for (unsigned int i = 0; i<3; ++i)
            HenckyMainStrainMatrix(i,i) = HenckyMainStrainVector[i];
        //std::cout<<" HenckyMainStrainVector "<<HenckyMainStrainVector<<std::endl;
        //Matrix NewElasticLeftCauchyGreen = ElasticVariables.CauchyGreenMatrix;
        Matrix NewElasticLeftCauchyGreen = HenckyMainStrainMatrix;

        this->CalculatePrincipalStressTrial(ElasticVariables, rValues, ReturnMappingVariables,NewElasticLeftCauchyGreen, StressMatrix );
        //std::cout<<" trial rStressMatrix "<< StressMatrix<<std::endl;
        //Attention!! when I call the return mapping function NewElasticLeftCauchyGreen represents the Hencky strain in matricial form
        //Attention!! when the return mapping is finished NewElasticLeftCauchyGreen is the NEW elastic left cauchy green tensor
        //Attention!! if and only it GetElasticLeftCachyGreen is a protected member of the flow rule that I am using.
        //Otherwise a public member of the flow rule base class has to be added as in this case.
        mpMPMFlowRule->CalculateReturnMapping( ReturnMappingVariables, ElasticVariables.DeformationGradientF, StressMatrix, NewElasticLeftCauchyGreen);
        mPlasticRegion = 0;
        if( ReturnMappingVariables.Options.Is(MPMFlowRule::PLASTIC_REGION) )
        {
            mPlasticRegion = 1;
        }
        //std::cout<<" StressMatrix after return mapping"<<StressMatrix<<std::endl;
        //if(DomainGeometry[0].Id() == 2295 && DomainGeometry[1].Id() == 2315 && DomainGeometry[2].Id() == 2313)
        //{

        ////std::cout<<" MainStrain "<<MainStrain<<std::endl;
        //std::cout<<" trial state function"<<ReturnMappingVariables.TrialStateFunction<<std::endl;

        //}
        this->CorrectDomainPressure( StressMatrix, ElasticVariables);

        //mpFlowRule->UpdateInternalVariables ( ReturnMappingVariables );


        //Evaluate the new LeftCauchyGreen tensor
        //mElasticLeftCauchyGreen = mpFlowRule->GetElasticLeftCauchyGreen(ReturnMappingVariables); //
        //std::cout<<" mElasticLeftCauchyGreen after return mapping"<<mElasticLeftCauchyGreen<<std::endl;
        //Stress vector updated
        StressVector = MathUtils<double>::StressTensorToVector(StressMatrix, StressVector.size());



        //StressVector = this->SetStressMatrixToAppropiateVectorDimension(StressVector, StressMatrix);
        //std::cout<<" StressVector after return mapping"<<StressVector<<std::endl;
        //}
        //if(DomainGeometry[0].Id() == 153 && DomainGeometry[1].Id() == 145 && DomainGeometry[2].Id() == 148)
        //{

        //std::cout<<" NewElasticLeftCauchyGreen "<<NewElasticLeftCauchyGreen<<std::endl;
        //std::cout<<" StressMatrix "<<StressMatrix<<std::endl;
        //std::cout<<" StressVector "<<StressVector<<std::endl;

        //}

        ////8.-Calculate Constitutive Matrix related to Total Kirchhoff stress
        //if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
        //{
        //initialize constitutive tensors
        ConstitutiveMatrix.clear();

        //Matrix ElastoPlasticTangentMatrix = ZeroMatrix(6,6);
        //Matrix ElastoPlasticTangentMatrix = ZeroMatrix(3,3);
        Matrix AuxConstitutiveMatrix = ZeroMatrix(6,6);
        //ElasticVariables.CauchyGreenMatrix = ElasticVariables.Identity;

        //if (ReturnMappingVariables.Options.Is(FlowRule::PLASTIC_REGION)  ) {
        double alfa = 0.0;
        //in the flow rule I distinguish if I am in plasticity or not
        //this->CalculateElastoPlasticTangentMatrix( ReturnMappingVariables, NewElasticLeftCauchyGreen, rAlpha, ElastoPlasticTangentMatrix, ElasticVariables);

        this->CalculateElastoPlasticTangentMatrix( ReturnMappingVariables, ElasticVariables.CauchyGreenMatrix, alfa, AuxConstitutiveMatrix, ElasticVariables);
        //mpFlowRule->ComputeElastoPlasticTangentMatrix( ReturnMappingVariables,  ElasticVariables.CauchyGreenMatrix, alfa, AuxConstitutiveMatrix);


        //this->CalculateElastoPlasticTangentMatrix(ReturnMappingVariables, ElasticVariables.CauchyGreenMatrix, StressMatrix, ElastoPlasticTangentMatrix, AuxConstitutiveMatrix);
        //std::cout<<"ConstitutiveMatrix "<<ConstitutiveMatrix<<std::endl;
        //std::cout<<"ConstitutiveMatrix "<<ConstitutiveMatrix<<std::endl;
        //}
        //else  {
        //mpFlowRule->CalculatePrincipalAxisAlgorithmicTangent(ReturnMappingVariables, StressMatrix, ConstitutiveMatrix);
        //}

        ConstitutiveMatrix = this->SetConstitutiveMatrixToAppropiateDimension(ConstitutiveMatrix, AuxConstitutiveMatrix);
        //std::cout<<"here 5"<<std::endl;
        //if(DomainGeometry[0].Id() == 2304 && DomainGeometry[1].Id() == 2326 && DomainGeometry[2].Id() == 2311)
        //{


        //std::cout<<" StressVector "<<StressVector<<std::endl;
        //std::cout<<" ConstitutiveMatrix "<<ConstitutiveMatrix<<std::endl;

        //}
    }



    if( Options.Is( ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE ) )
    {
        mpMPMFlowRule->UpdateInternalVariables ( ReturnMappingVariables );
        //std::cout<<" mElasticLeftCauchyGreen before in finalize material response"<<mElasticLeftCauchyGreen<<std::endl;
        mElasticLeftCauchyGreen = mpMPMFlowRule->GetElasticLeftCauchyGreen(ReturnMappingVariables);
        //std::cout<<" mElasticLeftCauchyGreen in finalize material response"<<mElasticLeftCauchyGreen<<std::endl;
        ElasticVariables.DeformationGradientF = DeformationGradientF;
        ElasticVariables.DeformationGradientF = Transform2DTo3D(ElasticVariables.DeformationGradientF);
        MathUtils<double>::InvertMatrix( ElasticVariables.DeformationGradientF, mInverseDeformationGradientF0, mDeterminantF0);
        mDeterminantF0 = DeterminantF; //special treatment of the determinant
    }

}

void HenckyElasticPlastic3DLaw::CalculatePrincipalStressTrial(const MaterialResponseVariables & rElasticVariables, Parameters& rValues, const MPMFlowRule::RadialReturnVariables & rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Matrix& rStressMatrix)
{
    Vector MainStrain      = ZeroVector(3);
    const Properties& MaterialProperties   = rValues.GetMaterialProperties();
    //const GeometryType&  DomainGeometry =  rElasticVariables.GetElementGeometry();
    for (unsigned int i = 0; i<3; ++i)
    {
        MainStrain[i] = rNewElasticLeftCauchyGreen(i,i);
    }
    //1- calculate the elastic matrix
    Matrix ElasticMatrix = ZeroMatrix(3,3);
    const double& Young       = MaterialProperties[YOUNG_MODULUS];
    const double& Nu = MaterialProperties[POISSON_RATIO];
    double diagonal   = Young/(1.0+Nu)/(1.0-2.0*Nu) * (1.0-Nu);//
    double nodiagonal = Young/(1.0+Nu)/(1.0-2.0*Nu) * ( Nu);//



    for (unsigned int i = 0; i<3; ++i)
    {
        for (unsigned int j = 0; j<3; ++j)
        {
            if (i == j)
            {
                ElasticMatrix(i,i) = diagonal;
            }
            else
            {
                ElasticMatrix(i,j) = nodiagonal;
            }
        }
    }

    Vector PrincipalStress = ZeroVector(3);

    //Evalute the Kirchhoff principal stress
    PrincipalStress = prod(ElasticMatrix, MainStrain);

    for(unsigned int i=0; i<3; i++)
    {
        rStressMatrix(i,i) = PrincipalStress(i);
    }

}
void HenckyElasticPlastic3DLaw::GetDomainPressure( double& rPressure, const MaterialResponseVariables& rElasticVariables)
{

    rPressure = 0.0;
    const GeometryType&  DomainGeometry =  rElasticVariables.GetElementGeometry();
    const Vector& ShapeFunctionsValues  =  rElasticVariables.GetShapeFunctionsValues();

    const unsigned int number_of_nodes  =  DomainGeometry.size();

    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
        rPressure += ShapeFunctionsValues[j] * DomainGeometry[j].FastGetSolutionStepValue(PRESSURE); //NOOOOO
    }

}
void HenckyElasticPlastic3DLaw::CorrectDomainPressure( Matrix& rStressMatrix, const MaterialResponseVariables & rElasticVariables)
{

}

void HenckyElasticPlastic3DLaw::CalculateElastoPlasticTangentMatrix( const MPMFlowRule::RadialReturnVariables & rReturnMappingVariables, const Matrix& rNewElasticLeftCauchyGreen, const double& rAlpha, Matrix& rElastoPlasticTangentMatrix, const MaterialResponseVariables& rElasticVariables )
{

    mpMPMFlowRule->ComputeElastoPlasticTangentMatrix( rReturnMappingVariables,  rNewElasticLeftCauchyGreen, rAlpha, rElastoPlasticTangentMatrix);


}
//void HenckyElasticPlastic3DLaw::CalculateElastoPlasticTangentMatrix( const FlowRule::RadialReturnVariables & rReturnMappingVariables, const Matrix& rTrialElasticLeftCauchyGreen, const Matrix& StressMatrix, Matrix& rElastoPlasticTangentMatrix, Matrix& rConsistentMatrix )
//{
//Matrix CostitutiveMatrix1 = ZeroMatrix(6,6);
//double alfa = 0.0;
////1. Evalute the tangent matrix depending on the type of return mapping(surface, line or apex)
//mpFlowRule->ComputeElastoPlasticTangentMatrix( rReturnMappingVariables,  rTrialElasticLeftCauchyGreen, alfa, rElastoPlasticTangentMatrix);

////2. Evaluate the eigenbasis matrix relative to the three principal directions
//Matrix EigenBasesMatrix = ZeroMatrix(3,9);
//EigenBasesMatrix = this->CalculateEigenbases(rReturnMappingVariables, EigenBasesMatrix);

//Matrix Ma = ZeroMatrix(3,3);
//Matrix Mb = ZeroMatrix(3,3);

//for(unsigned int i = 0;i<3;i++)
//{
//for(unsigned int j = 0;j<3;j++)
//{
////3. select the two eigenbases used to build the tensor base
//for(unsigned int k = 0;k<3;k++)
//{
//for(unsigned int l = 0;l<3;l++)
//{
//int index_A = l + 3 * i;
//int index_B = l + 3 * j;
//Ma(k,l) = EigenBasesMatrix(k, index_A);
//Mb(k,l) = EigenBasesMatrix(k, index_B);
//}
//}
//Matrix TensorProductMatrix = ZeroMatrix(6,6);
////4. Tensor product between the two eigenbases
//this->MyTensorProduct(Ma, Mb, TensorProductMatrix);
////5. evaluation of the -ij contribution
//CostitutiveMatrix1 += rElastoPlasticTangentMatrix(i,j) *  TensorProductMatrix;
//}
//}

//Matrix CostitutiveMatrix2 = ZeroMatrix(6,6);
////Matrix SpinTensor = ZeroMatrix(6,6);
////Matrix SpinTensor1 = ZeroMatrix(6,6);
////Matrix SpinTensor2 = ZeroMatrix(6,6);
////Matrix SpinTensor3 = ZeroMatrix(6,6);
////Matrix SpinTensor4 = ZeroMatrix(6,6);
////Matrix SpinTensor5 = ZeroMatrix(6,6);
////Matrix IdentityMatrix = identity_matrix<double> (3);
////Vector PrincipalStresses = ZeroVector(3);
////Vector PrincipalStrains = ZeroVector(3);

////for(unsigned int i=0;i<3;i++)
////{
////PrincipalStrains(i) = rReturnMappingVariables.StrainMatrix(i,i);
////PrincipalStresses(i) = rReturnMappingVariables.TrialIsoStressMatrix(i,i);
////}


////if(PrincipalStrains(0) == 0.0 || PrincipalStrains(1) == 0.0 || PrincipalStrains(2) == 0.0)
////{

////}
////else if (PrincipalStrains(0) == PrincipalStrains(1) || PrincipalStrains(0) == PrincipalStrains(2))
////{
//////std::cout<<"the same root!"<<std::endl;
//////std::cout<<" EIGENVALUES BEFORE"<<PrincipalStrains(0)<<" "<<PrincipalStrains(1)<<" "<<PrincipalStrains(2)<<std::endl;
////PrincipalStrains(1) = 1.00000001*PrincipalStrains(0);
////PrincipalStrains(2) = 1.00000002*PrincipalStrains(0);
//////std::cout<<" EIGENVALUES AFTER"<<PrincipalStrains(0)<<" "<<PrincipalStrains(1)<<" "<<PrincipalStrains(2)<<std::endl;

//////Matrix PrincipalDirectionStress = ZeroMatrix(3,3);
//////Matrix PrincipalDirectionStrain = ZeroMatrix(3,3);


//////mpFlowRule->GetPrincipalStressAndStrain(PrincipalStresses, PrincipalStrains);

//////double tol = 1e-9;
//////int iter = 100;

//////SolidMechanicsMathUtilities<double>::EigenVectors(rTrialElasticLeftCauchyGreen, PrincipalDirectionStrain, PrincipalStrains, tol, iter);
//////SolidMechanicsMathUtilities<double>::EigenVectors(rTrialElasticLeftCauchyGreen, PrincipalDirectionStress, PrincipalStresses, tol, iter);

////double FirstInvariant = 0.0;
////double ThirdInvariant = 0.0;

////for(unsigned int i=0; i<3; i++)
////{
////FirstInvariant += PrincipalStrains(i);
////ThirdInvariant *= PrincipalStrains(i);
////}
//////Matrix TensorProductMatrix = ZeroMatrix(6,6);
////for(unsigned int i=0; i<3; i++)
////{
////for(unsigned int k = 0;k<3;k++)
////{
////for(unsigned int l = 0;l<3;l++)
////{
////int index_A = l + 3 * i;

////Ma(k,l) = EigenBasesMatrix(k, index_A);
////}
////}
////double Da = 2 * PrincipalStrains(i) * PrincipalStrains(i) - FirstInvariant * PrincipalStrains(i) + ThirdInvariant /PrincipalStrains(i);
////double Coeff3 = ThirdInvariant/(PrincipalStrains(i) * PrincipalStrains(i));
////double Coeff1 = FirstInvariant + Coeff3 - 4 * PrincipalStrains(i);
////double Coeff2 = PrincipalStrains(i) * PrincipalStrains(i);
////this->MyTensorProduct(Ma, Ma, SpinTensor1);
//////SpinTensor1 = SpinTensor1*(-Coeff2);

////this->MyTensorProduct2(rTrialElasticLeftCauchyGreen, Ma, SpinTensor2);
//////SpinTensor2 = SpinTensor2* Coeff2;

////this->MyTensorProduct2(IdentityMatrix, Ma, SpinTensor3);
////SpinTensor3 = SpinTensor3*(-Coeff1);

////this->MyTensorProduct3(IdentityMatrix, SpinTensor4);
////SpinTensor4 = SpinTensor4 * Coeff1;

////this->MyTensorProduct4(rTrialElasticLeftCauchyGreen, SpinTensor5);
////SpinTensor = -Coeff1 * SpinTensor1 + Coeff2 * SpinTensor2 + Coeff3 *(SpinTensor4-SpinTensor3) + SpinTensor5;

////CostitutiveMatrix2 += PrincipalStresses(i)* 2 * SpinTensor / Da;
////}
////}
//rConsistentMatrix = CostitutiveMatrix1 + CostitutiveMatrix2;

//}
//************************************************************************************
//************************************************************************************
Vector HenckyElasticPlastic3DLaw::SetStressMatrixToAppropiateVectorDimension(Vector& rStressVector, const Matrix& rStressMatrix )
{
    rStressVector(0) = rStressMatrix(0,0);
    rStressVector(1) = rStressMatrix(1,1);
    rStressVector(2) = rStressMatrix(2,2);
    rStressVector(3) = rStressMatrix(0,1);
    rStressVector(4) = rStressMatrix(1,2);
    rStressVector(5) = rStressMatrix(0,2);
    return rStressVector;

}

Matrix HenckyElasticPlastic3DLaw::SetConstitutiveMatrixToAppropiateDimension(Matrix& rConstitutiveMatrix, const Matrix& rElastoPlasticTangentMatrix)
{

    rConstitutiveMatrix = rElastoPlasticTangentMatrix;
    return rConstitutiveMatrix;

}

/** CalculateEigenbases is a function which calculate the matrix with the eigenbases of the
 * three principal directions
*/
Matrix HenckyElasticPlastic3DLaw::CalculateEigenbases(const MPMFlowRule::RadialReturnVariables& rReturnMappingVariables, Matrix& rEigenbasesMatrix)
{
    //1- compute eigendirections
    Vector N1 = ZeroVector(3);
    Vector N2 = ZeroVector(3);
    Vector N3 = ZeroVector(3);
    //std::cout<<" rReturnMappingVariables.MainDirections"<<rReturnMappingVariables.MainDirections<<std::endl;
    for (unsigned int i = 0; i<3; ++i)
    {
        N1(i) = rReturnMappingVariables.MainDirections(i,0);
        N2(i) = rReturnMappingVariables.MainDirections(i,1);
        N3(i) = rReturnMappingVariables.MainDirections(i,2);
    }
    //std::cout<<" N1 "<<N1<<std::endl;
    //std::cout<<" N2 "<<N2<<std::endl;
    //std::cout<<" N3 "<<N3<<std::endl;

    //2- compute the eigenbases matrix
    Matrix M1 = ZeroMatrix(3,3);
    Matrix M2 = ZeroMatrix(3,3);
    Matrix M3 = ZeroMatrix(3,3);
    M1 = MathUtils<double>::TensorProduct3(N1, N1);
    M2 = MathUtils<double>::TensorProduct3(N2, N2);
    M3 = MathUtils<double>::TensorProduct3(N3, N3);

    // EigenbasesMatrix[3x9]
    for(unsigned int i = 0; i<3; i++)
    {
        for(unsigned int j = 0; j<3; j++)
        {
            int index_j = j+3;
            int index_j2 = j+6;

            rEigenbasesMatrix(i,j) = M1(i,j);
            rEigenbasesMatrix(i, index_j) = M2(i,j);
            rEigenbasesMatrix(i, index_j2) = M3(i,j);

        }
    }
    return rEigenbasesMatrix;
}

void HenckyElasticPlastic3DLaw::MyTensorProduct ( const Matrix& rMA, const Matrix& rMB,
        Matrix& rEigenbasesProductMatrix)
{

    rEigenbasesProductMatrix.clear();

    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
            rEigenbasesProductMatrix( i, j ) = TensorComponent(rEigenbasesProductMatrix( i, j ), rMA, rMB,
                                               this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
        }

    }


}

double& HenckyElasticPlastic3DLaw::TensorComponent(double & rCabcd,
        const Matrix& rMA, const Matrix& rMB,
        const unsigned int& a, const unsigned int& b,
        const unsigned int& c, const unsigned int& d)
{

    rCabcd =(rMA(a,b)*rMB(c,d));

    return rCabcd;
}
void HenckyElasticPlastic3DLaw::MyTensorProduct2 ( const Matrix& rMA, const Matrix& rMB,
        Matrix& rEigenbasesProductMatrix)
{

    rEigenbasesProductMatrix.clear();

    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
            rEigenbasesProductMatrix( i, j ) = TensorComponent2(rEigenbasesProductMatrix( i, j ), rMA, rMB,
                                               this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
        }

    }


}

double& HenckyElasticPlastic3DLaw::TensorComponent2(double & rCabcd,
        const Matrix& rMA, const Matrix& rMB,
        const unsigned int& a, const unsigned int& b,
        const unsigned int& c, const unsigned int& d)
{

    rCabcd =(rMA(a,b)*rMB(c,d) + rMB(a,b)*rMA(c,d));

    return rCabcd;
}
void HenckyElasticPlastic3DLaw::MyTensorProduct3 ( const Matrix& rMA, Matrix& rEigenbasesProductMatrix)
{

    rEigenbasesProductMatrix.clear();

    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
            rEigenbasesProductMatrix( i, j ) = TensorComponent3(rEigenbasesProductMatrix( i, j ), rMA,
                                               this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
        }

    }


}

double& HenckyElasticPlastic3DLaw::TensorComponent3(double & rCabcd,
        const Matrix& rMA,
        const unsigned int& a, const unsigned int& b,
        const unsigned int& c, const unsigned int& d)
{

    rCabcd =rMA(a,b)*rMA(c,d) - 0.5*(( rMA(a,c)*rMA(b,d))+( rMA(a,d)*rMA(b,c)));

    return rCabcd;
}
void HenckyElasticPlastic3DLaw::MyTensorProduct4 ( const Matrix& rMA, Matrix& rEigenbasesProductMatrix)
{

    rEigenbasesProductMatrix.clear();

    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
            rEigenbasesProductMatrix( i, j ) = TensorComponent4(rEigenbasesProductMatrix( i, j ), rMA,
                                               this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
        }

    }


}

double& HenckyElasticPlastic3DLaw::TensorComponent4(double & rCabcd,
        const Matrix& rMA,
        const unsigned int& a, const unsigned int& b,
        const unsigned int& c, const unsigned int& d)
{

    rCabcd = 0.5 *(rMA(a,c)*rMA(b,d) + rMA(a,d)*rMA(b,c)) - rMA(a,b)*rMA(c,d) ;

    return rCabcd;
}

//************************************************************************************
//************************************************************************************


Vector& HenckyElasticPlastic3DLaw::GetStressVectorFromMatrix(const Matrix& rStressMatrix,
        Vector& rMainStress,
        const Matrix& rEigenVectors)
{
//   EigenVectors = mpFlowRule->GetEigenVectors();
//Ll:
    Matrix auxMatrix = ZeroMatrix(3,3);
    auxMatrix = prod( rStressMatrix, trans(rEigenVectors));
    auxMatrix = prod( (rEigenVectors), auxMatrix);

    rMainStress = ZeroVector(3);
    for (unsigned int i = 0; i<3; i++)
        rMainStress[i] = auxMatrix(i,i);

    return rMainStress;
}



//************************************************************************************
//************************************************************************************


void HenckyElasticPlastic3DLaw::CalculateHenckyMainStrain(const Matrix& rCauchyGreenMatrix,
        MPMFlowRule::RadialReturnVariables& rReturnMappingVariables,
        Vector& rMainStrain)
{
    Matrix EigenVectors  = ZeroMatrix(3,3);
    Vector EigenValues   = ZeroVector(3);
    Vector EigenValues2   = ZeroVector(3);
    double tol = 1e-9;
    int iter = 100;
    //std::cout<<" rCauchyGreenMatrix "<<rCauchyGreenMatrix<<std::endl;
    SolidMechanicsMathUtilities<double>::EigenVectors(rCauchyGreenMatrix, EigenVectors, EigenValues, tol, iter);
    //EigenValues2 = SolidMechanicsMathUtilities<double>::EigenValues(rCauchyGreenMatrix, tol, iter);
    // lluis
    rReturnMappingVariables.MainDirections     = EigenVectors;
    // rReturnMappingVariables.TrialEigenValues = EigenValues;
    //std::cout<<" EigenValues "<<EigenValues<<std::endl;
    //std::cout<<" EigenVectors "<<EigenVectors<<std::endl;
    //std::cout<<" EigenValues2 "<<EigenValues2<<std::endl;
    for (unsigned int i = 0; i<3; ++i)
        rMainStrain[i] = 0.50 * std::log(EigenValues[i]);
    //std::cout<<" rMainStrain in cl"<<rMainStrain<<std::endl;
    //std::cout<<" rReturnMappingVariables.MainDirections in cl"<<rReturnMappingVariables.MainDirections<<std::endl;
}
//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void HenckyElasticPlastic3DLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
    rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
    rFeatures.mOptions.Set( FINITE_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the strain size
    rFeatures.mStrainSize = GetStrainSize();

    //Set the spacedimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();

}


//******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW*************************
//************************************************************************************

bool HenckyElasticPlastic3DLaw::CheckParameters(Parameters& rValues)
{
    return rValues.CheckAllParameters();
}



int HenckyElasticPlastic3DLaw::Check(const Properties& rMaterialProperties,
                                     const GeometryType& rElementGeometry,
                                     const ProcessInfo& rCurrentProcessInfo)
{

    if(YOUNG_MODULUS.Key() == 0 || rMaterialProperties[YOUNG_MODULUS]<= 0.00)
    {
        KRATOS_ERROR << "YOUNG_MODULUS has Key zero or invalid value " << std::endl;
    }

    const double& nu = rMaterialProperties[POISSON_RATIO];
    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );

    if(POISSON_RATIO.Key() == 0 || check==true)
    {
        KRATOS_ERROR << "POISSON_RATIO has Key zero invalid value " << std::endl;
    }


    if(DENSITY.Key() == 0 || rMaterialProperties[DENSITY]<0.00)
    {
        KRATOS_ERROR << "DENSITY has Key zero or invalid value " << std::endl;
    }


    return 0;

}

} // namespace Kratos
