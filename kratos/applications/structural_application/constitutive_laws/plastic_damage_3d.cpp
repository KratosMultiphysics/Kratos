/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
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

THE  SOFTWARE IS  PROVIDED  "AS  variablesIS", WITHOUT  WARRANTY  OF ANY  KIND,
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
*   Last Modified by:    $Author: nelson $
*   Date:                $Date: 2008-09-03
\*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/


// System includes
#include <iostream>
#include <algorithm>

// External includes
#include<cmath>

// Project includes
#include "constitutive_laws/plastic_damage_3d.h"



namespace Kratos
{
    /**
     * TO BE TESTED!!!
     */

    PlasticDamage3D::PlasticDamage3D()
            : ConstitutiveLaw()

    {
        KRATOS_ERROR( std::logic_error, "Calling the empty constructor.", "" );
    }

    PlasticDamage3D::PlasticDamage3D(
        FluencyCriteriaPointer FluencyCriteria,
        SofteningHardeningCriteriaPointer SofteningBehavior,
        PropertiesPointer Property )
            : ConstitutiveLaw()
    {
        mpFluencyCriteria_1      = FluencyCriteria;
        mpSofteningBehavior      = SofteningBehavior;
        mpProperties             = Property;

        // mFluencyCriteria = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        //KRATOS_WATCH(*mProperties)
    }

    /**
     * TO BE TESTED!!!
     */
    PlasticDamage3D::~PlasticDamage3D()
    {
    }


    bool PlasticDamage3D::Has( const Variable<double>& rThisVariable )
    {
        return false;
    }

    bool PlasticDamage3D::Has( const Variable<Vector>& rThisVariable )
    {
        return false;
    }

    bool PlasticDamage3D::Has( const Variable<Matrix>& rThisVariable )
    {
        return false;
    }

    double& PlasticDamage3D::GetValue( const Variable<double>& rThisVariable, double& rValue )
    {
        if ( rThisVariable == DAMAGE )
        {
            return mlocal_fail_factor;
        } //mplastic_damage ;}

//    else if(rThisVariable == COHESION)
//               {return mcohesion;}
//           else if(rThisVariable == DILATANCY_ANGLE)
//               {return mdilatancy_angle*180.00/PI;}
//           else if(rThisVariable == FRICTION_INTERNAL_ANGLE)
//               {return mfriction_angle*180.00/PI;}
//           else
        return rValue;
    }


    Vector& PlasticDamage3D::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
    {
        return( rValue );
    }

    Matrix& PlasticDamage3D::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
    {
        return( rValue );
    }

    void PlasticDamage3D::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                                    const ProcessInfo& rCurrentProcessInfo )
    {
    }

    void PlasticDamage3D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                                    const ProcessInfo& rCurrentProcessInfo )
    {
    }

    void PlasticDamage3D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                                    const ProcessInfo& rCurrentProcessInfo )
    {
    }

    void PlasticDamage3D::Calculate( const Variable<Matrix >& rVariable, Matrix& rResult,
                                     const ProcessInfo& rCurrentProcessInfo )
    {
        for ( unsigned int ii = 0; ii < 3; ii++ )
            rResult( 0, ii ) = mplastic_strain( ii );
    }

    void PlasticDamage3D::Calculate( const Variable<double>& rVariable,
                                     double& Output,
                                     const ProcessInfo& rCurrentProcessInfo )
    {

        Output = sqrt( mE / mDE );
        //double local_damage =  fabs(1.00-mplastic_damage );
        //Output = sqrt(mE/mDE);
        //if(Output==0){Output = 1;}
    }

//***********************************************************************************************
//***********************************************************************************************

    void PlasticDamage3D::InitializeMaterial( const Properties& props,
            const GeometryType& geom,
            const Vector& ShapeFunctionsValues )

    {

        mE                       = ( *mpProperties )[YOUNG_MODULUS];
        mNU                      = ( *mpProperties )[POISSON_RATIO];
        mDE                      = ( *mpProperties )[DENSITY];
        mlength                  = pow( fabs( geom.Volume() ), 0.333333333333333333333 );


        mFt[0]                   = ( *mpProperties )[FT];
        mFt[1]                   = ( *mpProperties )[FT];
        mFt[2]                   = ( *mpProperties )[FT];

        mplastic_damage          = 0.00;
        mdisipation              = 0.00;
        mefective_plastic_strain = 0.00;
        mlocal_fail_factor       = 0.00;

        noalias( mplastic_strain )         = ZeroVector( 6 );
        noalias( mcurrent_plastic_strain ) = ZeroVector( 6 );

        double Gc           = ( *mpProperties )[CRUSHING_ENERGY] / mlength;
        double length_limit = 2.00 * mE * Gc / (( *mpProperties )[FC] * ( *mpProperties )[FC] );

        if ( length_limit < mlength )
        {
            std::cout << "Element length greater than permitted" << std::endl;
        }

        mpFluencyCriteria_1->InitializeMaterial( *mpProperties );


    }



//***********************************************************************************************
//***********************************************************************************************

    void PlasticDamage3D::InitializeSolutionStep( const Properties& props,
            const GeometryType& geom,
            const Vector& ShapeFunctionsValues ,
            const ProcessInfo& CurrentProcessInfo )
    {

    }


//***********************************************************************************************
//***********************************************************************************************

    void PlasticDamage3D::FinalizeSolutionStep( const Properties& props,
            const GeometryType& geom,
            const Vector& ShapeFunctionsValues ,
            const ProcessInfo& CurrentProcessInfo )
    {

        //uptating
        //KRATOS_WATCH("FINALIZE_FINALIZE_FINALIZE_FINALIZE_FINALIZE_FINALIZE_FINALIZE")
        mefective_plastic_strain =  mcurrent_efective_plastic_strain;
        noalias( mplastic_strain ) =  mcurrent_plastic_strain;
        mplastic_damage          =  mcurrent_plastic_damage;
        mdisipation              =  mcurrent_disipation;
        noalias( mFt )             =  mcurrent_Ft;
        //KRATOS_WATCH(mFt )


        array_1d<double, 4> Variables;
        Variables[0] = mFt[0];
        Variables[1] = mFt[1];
        Variables[2] = mFt[2];
        Variables[3] =  mlength;

        mpFluencyCriteria_1->Finalize();
        mpFluencyCriteria_1->UpdateVariables( Variables );



        Matrix Aux_V         = ZeroMatrix( 3, 3 );
        Matrix PlasticTensor = ZeroMatrix( 3, 3 );
        Vector Aux_P         = ZeroVector( 3 );

        PlasticTensor        = SD_MathUtils<double>::StrainVectorToTensor( mplastic_strain );

        double Euc           =  2.00 * ( *mpProperties )[FRACTURE_ENERGY] / (( *mpProperties )[FT] * mlength );
        SD_MathUtils<double>::EigenVectors( PlasticTensor, Aux_V, Aux_P, 1E-9, 100 );
        double Ef = ( *max_element( Aux_P.begin(), Aux_P.end() ) );
        mlocal_fail_factor =  Ef / Euc;
        //if (mlocal_fail_factor > 1.00) {mlocal_fail_factor = 1.00; }

    }


//***********************************************************************************************
//***********************************************************************************************


    void PlasticDamage3D::CalculateElasticMatrix( boost::numeric::ublas::bounded_matrix<double, 6, 6>& C )
    {

//setting up material matrix
        double c1 = mE / (( 1.00 + mNU ) * ( 1.00 - 2.00 * mNU ) );
        double c2 = c1 * ( 1.00 - mNU );
        double c3 = c1 * mNU;
        double c4 = c1 * 0.5 * ( 1.00 - 2.00 * mNU );
//filling material matrix
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


//***********************************************************************************************
//***********************************************************************************************

    void PlasticDamage3D::CalculateElasticStress( const array_1d<double, 6>& Strain, array_1d<double, 6>& Stress )
    {

        /*boost::numeric::ublas::bounded_matrix<double,4,4> C;
        CalculateElasticMatrix(C);
        noalias(Stress) = prod(C, Strain);
        KRATOS_WATCH(Stress)
        */
///* Owen
        double G          = 0.5 * mE / ( 1.00 + mNU );
        double K          = mE / ( 3.00 * ( 1.00 - 2.00 * mNU ) );
        double vol_strain = Strain[0] + Strain[1] + Strain[2];
        double pt         = K * vol_strain;

        Stress[0] = 2.00 * G * ( Strain[0] - vol_strain / 3.00 ) + pt;
        Stress[1] = 2.00 * G * ( Strain[1] - vol_strain / 3.00 ) + pt;
        Stress[2] = 2.00 * G * ( Strain[2] - vol_strain / 3.00 ) + pt;
        Stress[3] = G * ( Strain[3] );
        Stress[4] = G * ( Strain[4] );
        Stress[5] = G * ( Strain[5] );

    }


//***********************************************************************************************
//***********************************************************************************************
    void PlasticDamage3D::CalculateMaterialResponse( const Vector& StrainVector,
            const Matrix& DeformationGradient,
            Vector& StressVector,
            Matrix& AlgorithmicTangent,
            const ProcessInfo& CurrentProcessInfo,
            const Properties& props,
            const GeometryType& geom,
            const Vector& ShapeFunctionsValues,
            bool CalculateStresses,
            int CalculateTangent,
            bool SaveInternalVariables )
    {
        if ( CalculateStresses )
            CalculateStress( StrainVector, StressVector );

        if ( CalculateTangent != 0 )
            CalculateStressAndTangentMatrix( StressVector, StrainVector, AlgorithmicTangent );
    }


    void PlasticDamage3D::CalculateConstitutiveMatrix( const Vector& StrainVector, Matrix& ConstitutiveMatrix )
    {

        ConstitutiveMatrix.resize( 6, 6 );
        Matrix& C = ConstitutiveMatrix;
//setting up material matrix
        double c1 = mE / (( 1.00 + mNU ) * ( 1.00 - 2.00 * mNU ) );
        double c2 = c1 * ( 1.00 - mNU );
        double c3 = c1 * mNU;
        double c4 = c1 * 0.5 * ( 1.00 - 2.00 * mNU );
//filling material matrix
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


//***********************************************************************************************
//***********************************************************************************************

    void PlasticDamage3D::CalculateStress( const Vector& StrainVector,
                                           Vector& StressVector )
    {

        double ElasticDomain_1               = 0.00;
        double toler                         = 1E-6;
        array_1d<double, 6> ElasticStrain     = ZeroVector( 6 );
        array_1d<double, 6> ElasticStress     = ZeroVector( 6 );
        StressVector                         = ZeroVector( 6 );


        ///* The elastic strain trial
        noalias( ElasticStrain ) = StrainVector - mcurrent_plastic_strain;

        ///*calculating elastic stress trial
        CalculateElasticStress( ElasticStrain, ElasticStress );


        ///* Comprobado criterio de fluencia
        mpFluencyCriteria_1->CalculateEquivalentUniaxialStress( ElasticStress, ElasticDomain_1 );

        ///* dominio elastico

        if ( ElasticDomain_1 <= toler )
        {
            noalias( mcurrent_Ft )             = mFt;
            mcurrent_efective_plastic_strain = mefective_plastic_strain;
            noalias( mcurrent_plastic_strain ) = mplastic_strain;
            mcurrent_plastic_damage          = mplastic_damage;
            mComputeTangentMatrix            = false;
        }

        else
        {
            mComputeTangentMatrix = false;
            ///* Spectral Descomposition
            Compute_Principal_Stress( ElasticStress );
            array_1d<double, 3 > PrincipalStress;
            array_1d<double, 3 > Sigma;
            array_1d<unsigned int, 3 > Order;
            Vector delta_gamma;


            ///*
            noalias( PrincipalStress ) = mPrincipalStress;

            ///* Guardando el orden de los cambios
            IdentifyMaximunAndMinumumPrincipalStres_CalculateOrder( PrincipalStress , Order );

            ///* Regresando ala superficie de fluencia
            mpFluencyCriteria_1->ReturnMapping( ElasticStress, ElasticStrain, delta_gamma, Sigma );

            ///* Calculando las tensiones elsaticas
            AssembleUpdateStressAndStrainTensor( Sigma, StrainVector, Order, ElasticStrain, ElasticStress );
            //KRATOS_WATCH(Sigma)

            ///* General evoluation of accumulated hardening
            mcurrent_efective_plastic_strain = mefective_plastic_strain + norm_1( delta_gamma );

            ///* aculated Von mises plastic
            //double aux_var = sqrt(delta_gamma_a*delta_gamma_a + delta_gamma_b * delta_gamma_b + delta_gamma_c * delta_gamma_c);
            //mcurrent_efective_plastic_strain = mefective_plastic_strain + (2.00 / sqrt(3.00) ) * aux_var;

            ///* updating current variables
            Vector Result( 3 );
            mpFluencyCriteria_1->GetValue( Result );
            noalias( mcurrent_Ft ) = Result;


        }


        StressVector[0] = ElasticStress[0];

        StressVector[1] = ElasticStress[1];
        StressVector[2] = ElasticStress[2];
        StressVector[3] = ElasticStress[3];
        StressVector[4] = ElasticStress[4];
        StressVector[5] = ElasticStress[5];

    }


    void PlasticDamage3D::CalculateStressAndTangentMatrix( Vector& StressVector,
            const Vector& StrainVector,
            Matrix& algorithmicTangent )
    {

        CalculateConstitutiveMatrix( StrainVector, algorithmicTangent );

    }




    void PlasticDamage3D::CalculateCauchyStresses(
        Vector& rCauchy_StressVector,
        const Matrix& rF,
        const Vector& rPK2_StressVector,
        const Vector& rGreenLagrangeStrainVector )
    {

    }

    void PlasticDamage3D::UpdateMaterial( const Vector& StrainVector,
                                          const Properties& props,
                                          const GeometryType& geom,
                                          const Vector& ShapeFunctionsValues,
                                          const ProcessInfo& CurrentProcessInfo )
    {
        noalias( mcurrent_plastic_strain ) = mplastic_strain;
        mcurrent_plastic_damage          = mplastic_damage;
        mcurrent_disipation              = mdisipation;
        mcurrent_efective_plastic_strain = mefective_plastic_strain;
        noalias( mcurrent_Ft )             = mFt;

        array_1d<double, 4> Variables;
        Variables[0] = mFt[0];
        Variables[1] = mFt[1];
        Variables[2] = mFt[2];
        Variables[3] =  mlength;

        mpFluencyCriteria_1->UpdateVariables( Variables );


    }




    void PlasticDamage3D::Compute_Principal_Stress( const Vector& StressVector )
    {

        int    iter = 50;
        double zero = 1.0E-9;
        Matrix StressTensor    = ZeroMatrix( 3, 3 );
        Matrix EigenVectors    = ZeroMatrix( 3, 3 );
        mPrincipalStress       = ZeroVector( 3 );

        mpFluencyCriteria_1->State_Tensor( StressVector, StressTensor );
        SD_MathUtils<double>::EigenVectors( StressTensor, EigenVectors, mPrincipalStress, zero, iter );

        mEigenVectors[0] = ZeroVector( 3 );
        mEigenVectors[1] = ZeroVector( 3 );
        mEigenVectors[2] = ZeroVector( 3 );

        mEigenVectors[0][0] = EigenVectors( 0, 0 );
        mEigenVectors[0][1] = EigenVectors( 0, 1 );
        mEigenVectors[0][2] = EigenVectors( 0, 2 );

        mEigenVectors[1][0] = EigenVectors( 1, 0 );
        mEigenVectors[1][1] = EigenVectors( 1, 1 );
        mEigenVectors[1][2] = EigenVectors( 1, 2 );

        mEigenVectors[2][0] = EigenVectors( 2, 0 );
        mEigenVectors[2][1] = EigenVectors( 2, 1 );
        mEigenVectors[2][2] = EigenVectors( 2, 2 );

    }

    void PlasticDamage3D::IdentifyMaximunAndMinumumPrincipalStres_CalculateOrder( array_1d<double, 3 >& PrincipalStress , array_1d<unsigned int, 3>& order )
    {

        order[0] = 0;
        order[2] = 0;
        PrincipalStress[0] = mPrincipalStress[order[0]];
        PrincipalStress[2] = mPrincipalStress[order[2]];

        for ( unsigned int i = 1; i < 3; i++ )
        {
            if ( mPrincipalStress[i] >= PrincipalStress[0] )
            {
                order[0] = i;
                PrincipalStress[0] = mPrincipalStress[order[0]];
            }

            if ( mPrincipalStress[i] <= PrincipalStress[2] )
            {
                order[2] = i;
                PrincipalStress[2] = mPrincipalStress[order[2]];
            }
        }

        if ( order[0] != 0 && order[2] != 0 )
        {
            order[1] = 0;
        }

        if ( order[0] != 1 && order[2] != 1 )
        {
            order[1] = 1;
        }

        if ( order[0] != 2 && order[2] != 2 )
        {
            order[1] = 2;
        }

        PrincipalStress[1] = mPrincipalStress[order[1]];

        if ( PrincipalStress[0] == PrincipalStress[1] )
        {
            order[0] = 0;
            order[1] = 1;
            order[2] = 2;
        }

        ///* Evitando redondeos matematicos absurdos
        if ( fabs( PrincipalStress[0] ) < 1E-3 )
        {
            PrincipalStress[0] = 0.00;
        }

        if ( fabs( PrincipalStress[1] ) < 1E-3 )
        {
            PrincipalStress[1] = 0.00;
        }

        if ( fabs( PrincipalStress[2] ) < 1E-3 )
        {
            PrincipalStress[2] = 0.00;
        }

    }



    void PlasticDamage3D::AssembleUpdateStressAndStrainTensor(
        const array_1d<double, 3>& Sigma,
        const array_1d<double, 6>& Total_Strain,
        const array_1d<unsigned int, 3>& order,
        array_1d<double, 6>& ElasticStrain,
        array_1d<double, 6>& ElasticStress )
    {

        Matrix StressTensor      = ZeroMatrix( 3, 3 );
        Matrix DesviatoricTensor = ZeroMatrix( 3, 3 );
        Matrix VolumnetricTensor = ZeroMatrix( 3, 3 );
        Matrix StrainTensor      = ZeroMatrix( 3, 3 );
        const identity_matrix<double> I( 3 );

//KRATOS_WATCH(order)

///* Updating the  elastic stress tensor

        for ( unsigned int i = 0; i < 3; i++ )
        {
            noalias( StressTensor ) = StressTensor + Sigma[i] * Matrix( outer_prod( mEigenVectors[order[i]],  mEigenVectors[order[i]] ) );
        }


        ElasticStress[0] = StressTensor( 0, 0 );

        ElasticStress[1] = StressTensor( 1, 1 );
        ElasticStress[2] = StressTensor( 2, 2 );
        ElasticStress[3] = StressTensor( 0, 1 );
        ElasticStress[4] = StressTensor( 1, 2 );
        ElasticStress[5] = StressTensor( 0, 2 );


///* Updating the elastic tensor
        double p  = ( ElasticStress[0] + ElasticStress[1] + ElasticStress[2] ) / 3.00;
        double G  = 0.5 * mE / ( 1.00 + mNU );
        double K  = mE / ( 3.00 * ( 1.00 - 2.00 * mNU ) );
//noalias(DesviatoricTensor) = StressTensor - p * I;
//noalias(StrainTensor)      = (1.00/(2.00 * G ) ) * DesviatoricTensor + p/(3.00 * K) * I;

        ElasticStrain[0] = ( ElasticStress[0] - p ) / ( 2.00 * G ) + p / ( 3.00 * K );
        ElasticStrain[1] = ( ElasticStress[1] - p ) / ( 2.00 * G ) + p / ( 3.00 * K );
        ElasticStrain[2] = ( ElasticStress[2] - p ) / ( 2.00 * G ) + p / ( 3.00 * K );
        ElasticStrain[3] =  ElasticStress[3] / G;
        ElasticStrain[4] =  ElasticStress[4] / G;
        ElasticStrain[5] =  ElasticStress[5] / G;

        noalias( mcurrent_plastic_strain ) = Total_Strain - ElasticStrain;

//Updated_Internal_Variables(ElasticStress, mcurrent_plastic_strain);////
//Tensile_Fracture_Model(ElasticStress,  mcurrent_plastic_strain);


    }

    void  PlasticDamage3D::Tensile_Fracture_Model( const Vector& Stress, const Vector& plastic_strain )
    {
// Anisotropic Rotating crack  Model
        const int    iter      = 50;
        const double zero      = 1.0E-9;
        double& ft = ( *mpProperties )[FT];
        double& Gf = ( *mpProperties )[FRACTURE_ENERGY];

        double critical_strain = 2.00 * Gf / ( mlength * ft ); // WARNING = Unidades en Kg - cm

        array_1d<double, 3> H =  ZeroVector( 3 );

//array_1d<double,3> Rotating;
        Vector principal_strain_current = ZeroVector( 3 );
        Vector principal_strain_old     = ZeroVector( 3 );
        Matrix EigenVectors             = ZeroMatrix( 3, 3 );
        Matrix StrainTensor             = ZeroMatrix( 3, 3 );

/// plastic strain = current plastic strain
        StrainTensor = SD_MathUtils<double>::StrainVectorToTensor( plastic_strain );
        SD_MathUtils<double>::EigenVectors( StrainTensor, EigenVectors, principal_strain_current, zero, iter );
        sort( principal_strain_current.begin(), principal_strain_current.end() );
        reverse( principal_strain_current.begin(), principal_strain_current.end() );


///mplastic_strain = the old strain before do the updtate
        StrainTensor = SD_MathUtils<double>::StrainVectorToTensor( mplastic_strain );
        SD_MathUtils<double>::EigenVectors( StrainTensor, EigenVectors, principal_strain_old, zero, iter );
        sort( principal_strain_old.begin(), principal_strain_old.end() );
        reverse( principal_strain_old.begin(), principal_strain_old.end() );



        array_1d<double, 3> Delta_Strain                 = ZeroVector( 3 );
        noalias( Delta_Strain ) = principal_strain_current - principal_strain_old;

        for ( unsigned int i = 0; i < 3; i++ )
        {
            if ( Delta_Strain[i] <= 0.00 )
            {
                Delta_Strain[i] = 0.00;
            }

            if ( principal_strain_current[i] <= critical_strain )
            {
                H[i]   = mlength * ft * ft / ( 2.00 * Gf );
            }
        }



        mcurrent_Ft[0] = mFt[0] -  H[0] * Delta_Strain[0];

        mcurrent_Ft[1] = mFt[1] -  H[1] * Delta_Strain[1];
        mcurrent_Ft[2] = mFt[2] -  H[2] * Delta_Strain[2];

        if ( mcurrent_Ft[0] <= 0.00 )
        {
            mcurrent_Ft[0] = 0.00;
        }

        if ( mcurrent_Ft[1] <= 0.00 )
        {
            mcurrent_Ft[1] = 0.00;
        }

        if ( mcurrent_Ft[2] <= 0.00 )
        {
            mcurrent_Ft[2] = 0.00;
        }


        mpFluencyCriteria_1->UpdateVariables( mcurrent_Ft );



    }




}

