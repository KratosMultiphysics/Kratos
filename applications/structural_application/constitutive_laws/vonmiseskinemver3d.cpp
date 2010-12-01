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
*   Last Modified by:    $Author: janosch $
*   Date:                $Date: 2008-10-23 12:22:22 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes

#include "includes/define.h"
#include "constitutive_laws/von_mises_kinemVer_3d.h"
#include "includes/constitutive_law.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "nonlinear_materials_application.h"
#include "includes/properties.h"

namespace Kratos
{
    /**
     * TO BE TESTED!!!
     */
    VonMisesKinemVer3D::VonMisesKinemVer3D()
            : ConstitutiveLaw< Node<3> >()
    {
    }

    /**
     * TO BE TESTED!!!
     */
    VonMisesKinemVer3D::~VonMisesKinemVer3D()
    {
    }

    bool VonMisesKinemVer3D::Has( const Variable<double>& rThisVariable )
    {
        return false;
    }

    bool VonMisesKinemVer3D::Has( const Variable<Vector>& rThisVariable )
    {
        if ( rThisVariable == INSITU_STRESS )
            return true;

        return false;
    }

    bool VonMisesKinemVer3D::Has( const Variable<Matrix>& rThisVariable )
    {
        return false;
    }

    double VonMisesKinemVer3D::GetValue( const Variable<double>& rThisVariable )
    {
        KRATOS_ERROR( std::logic_error, "Vector Variable case not considered" , "" );
    }

    Vector VonMisesKinemVer3D::GetValue( const Variable<Vector>& rThisVariable )
    {
        if ( rThisVariable == INSITU_STRESS )
            return mInSituStress;

        if ( rThisVariable == INTERNAL_VARIABLES )
        {
            Vector dummy = ZeroVector( 8 );
            dummy[0] = f_history;
            dummy[1] = mPlasticStrainVector( 0 );
            dummy[2] = mPlasticStrainVector( 1 );
            dummy[3] = mPlasticStrainVector( 2 );
            dummy[4] = mCurrentStress( 0 );
            dummy[5] = mCurrentStress( 1 );
            dummy[6] = mCurrentStress( 2 );
            dummy[7] = mCurrentAlpha;


            return( dummy );
        }

    }

    Matrix VonMisesKinemVer3D::GetValue( const Variable<Matrix>& rThisVariable )
    {
        KRATOS_ERROR( std::logic_error, "Vector Variable case not considered", "" );
    }

    void VonMisesKinemVer3D::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                                       const ProcessInfo& rCurrentProcessInfo )
    {
    }

    void VonMisesKinemVer3D::SetValue( const Variable<array_1d<double, 3> >& rThisVariable,
                                       const array_1d<double, 3>& rValue,
                                       const ProcessInfo& rCurrentProcessInfo )
    {
    }

    void VonMisesKinemVer3D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                                       const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rThisVariable == INSITU_STRESS )
        {
            mInSituStress = rValue;
        }
    }

    void VonMisesKinemVer3D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                                       const ProcessInfo& rCurrentProcessInfo )
    {
    }

    /**
     * TO BE TESTED!!!
     */
    void VonMisesKinemVer3D::InitializeMaterial( const Properties& props,
            const GeometryType& geom,
            const Vector& ShapeFunctionsValues )
    {
        f_history = 0.0;
        mCurrentStress = ZeroVector( 6 );
        mInSituStress = ZeroVector( 6 );
        mPlasticStrainVector = ZeroVector( 6 );
        mCurrentPlasticStrainVector = ZeroVector( 6 );
        mAlpha = 0.0;
        mCurrentAlpha = 0.0;
        mRho = ZeroVector( 6 );
        mCurrentRho = ZeroVector( 6 );

        mCtangent = ZeroMatrix( 6, 6 );
        //Kompression modulus and Shear modulus
        mK = props[MATERIAL_PARAMETERS][0] / ( 3.0 * ( 1.0 - 2.0 * props[MATERIAL_PARAMETERS][1] ) );
        mMU = props[MATERIAL_PARAMETERS][0] / ( 2.0 * ( 1.0 + props[MATERIAL_PARAMETERS][1] ) );

        // initial yielding parameter Y0 and hardening modulus H
        mSigmaY0 = props[MATERIAL_PARAMETERS][7];
        mH = props[MATERIAL_PARAMETERS][8];
        mBeta = props[MATERIAL_PARAMETERS][20];

        //deviatoric unity matrix
        mIdev = ZeroMatrix( 6, 6 );

        for ( int i = 0; i < 6; i++ )
            mIdev( i, i ) = 1.0;

        for ( int i = 0; i < 3; i++ )
            for ( int j = 0; j < 3; j++ )
                mIdev( i, j ) -= 1.0 / 3.0;


    }

    void VonMisesKinemVer3D::InitializeSolutionStep( const Properties& props,
            const GeometryType& geom, //this is just to give the array of nodes
            const Vector& ShapeFunctionsValues ,
            const ProcessInfo& CurrentProcessInfo )
    {
    }

    void VonMisesKinemVer3D::FinalizeSolutionStep( const Properties& props,
            const GeometryType& geom, //this is just to give the array of nodes
            const Vector& ShapeFunctionsValues ,
            const ProcessInfo& CurrentProcessInfo )
    {
//         std::cout <<"******************** in FinalizeSolutionStep ********************"<< std::endl;

        if ( CurrentProcessInfo[CALCULATE_INSITU_STRESS] )
        {
            mInSituStress -= mCurrentStress;
            //SetValue( INSITU_STRESS, mInSituStress, CurrentProcessInfo );
        }

        //update plastic strains
        mPlasticStrainVector = mCurrentPlasticStrainVector;

        //update alpha
        mAlpha = mCurrentAlpha;

        //update rho
        mRho = mCurrentRho;
    }

    void VonMisesKinemVer3D::UpdateMaterial( const Vector& StrainVector,
            const Properties& props,
            const GeometryType& geom,
            const Vector& ShapeFunctionsValues,
            const ProcessInfo& CurrentProcessInfo )
    {

        if ( geom.WorkingSpaceDimension() != 3 )
        {
            KRATOS_ERROR( std::logic_error, "This constitutive law is defined for 3D only!" , "" );
        }

        //set up elastic matrix
        mCelastic = ZeroMatrix( 6, 6 );

        //set up tangent matrix
        mCtangent = ZeroMatrix( 6, 6 );

        //update current plastic strains
        mCurrentPlasticStrainVector = mPlasticStrainVector;

        //update current alpha
        mCurrentAlpha = mAlpha;

        //update current rho
        mCurrentRho = mRho;

        //calculate elastic matrix
        CalculateElasticMatrix( mCelastic, mK, mMU );

        //calculate elastic stress
        mCurrentStress = prod( mCelastic, ( StrainVector - mPlasticStrainVector ) );

        //calculate deviatoric stress
        Vector SigmaDev = mCurrentStress;

        for ( int i = 0; i < 3; i++ )
        {
            SigmaDev[i] -= ( mCurrentStress[0] + mCurrentStress[1] + mCurrentStress[2] ) / 3.0;
        }

        //calculate eta
        Vector Eta = SigmaDev - mRho;

        //calculate norm of deviatoric stress
        double Etanorm = 0.0;

        for ( int i = 0; i < 3; i++ )
        {
            Etanorm += Eta[i] * Eta[i];
            Etanorm += 2.0 * Eta[i+3] * Eta[i+3];
        }

        Etanorm = sqrt( Etanorm );

        //calculate yield function
        double f = ( Etanorm ) - sqrt( 2.0 / 3.0 ) * ( mSigmaY0 + mBeta * mH * mAlpha );


        if ( f > 0.0 )
        {
            //calculate norm vector / projection vector
            Vector N = Eta / Etanorm;

            //calculate plastic multiplier
            double dgamma = ( Etanorm - sqrt( 2.0 / 3.0 ) * ( mSigmaY0 + mBeta * mH * mAlpha ) ) / ( 2.0 * mMU );

            //calculate increment of the hardening parameter alpha
            double dAlpha = sqrt( 2.0 / 3.0 ) * dgamma;

            //calculate increment of the kinematic hardening vector rho
            Vector dRho = sqrt( 2.0 / 3.0 ) * ( 1.0 - mBeta ) * mH * dgamma * N;

            //update plastic strains
            mCurrentPlasticStrainVector = mPlasticStrainVector + dgamma * N;

            //update alpha
            mCurrentAlpha = mAlpha + dAlpha;

            //update eta
            mCurrentRho = mRho + dRho;

            //update stress
            mCurrentStress = prod( mCelastic, ( StrainVector - mCurrentPlasticStrainVector ) );

            //update tangent matrix
            mCtangent = mCelastic - 2.0 * mMU * (( 1.0 / ( 1.0 + (( mBeta * mH ) / ( 3.0 * mMU ) ) ) ) - ( 2.0 * mMU * dgamma / Etanorm ) ) * outer_prod( N, N );
            mCtangent = mCtangent - ( 4.0 * mMU * mMU * dgamma / Etanorm ) * mIdev;
        }
        else
        {
            mCtangent = mCelastic;
        }

    }


    void VonMisesKinemVer3D::CalculateElasticMatrix( Matrix& C, const double K, const double MU )
    {
        C = ZeroMatrix( 6, 6 );
        //setting up material matrix

        double c1 = K - ( 2.0 * MU / 3.0 );
        double c2 = 2.0 * MU;
        double c3 = MU;
        //filling material matrix
        C( 0, 0 ) = c1 + c2;
        C( 0, 1 ) = c1;
        C( 0, 2 ) = c1;
        C( 0, 3 ) = 0.0;
        C( 0, 4 ) = 0.0;
        C( 0, 5 ) = 0.0;
        C( 1, 0 ) = c1;
        C( 1, 1 ) = c1 + c2;
        C( 1, 2 ) = c1;
        C( 1, 3 ) = 0.0;
        C( 1, 4 ) = 0.0;
        C( 1, 5 ) = 0.0;
        C( 2, 0 ) = c1;
        C( 2, 1 ) = c1;
        C( 2, 2 ) = c1 + c2;
        C( 2, 3 ) = 0.0;
        C( 2, 4 ) = 0.0;
        C( 2, 5 ) = 0.0;
        C( 3, 0 ) = 0.0;
        C( 3, 1 ) = 0.0;
        C( 3, 2 ) = 0.0;
        C( 3, 3 ) = c3;
        C( 3, 4 ) = 0.0;
        C( 3, 5 ) = 0.0;
        C( 4, 0 ) = 0.0;
        C( 4, 1 ) = 0.0;
        C( 4, 2 ) = 0.0;
        C( 4, 3 ) = 0.0;
        C( 4, 4 ) = c3;
        C( 4, 5 ) = 0.0;
        C( 5, 0 ) = 0.0;
        C( 5, 1 ) = 0.0;
        C( 5, 2 ) = 0.0;
        C( 5, 3 ) = 0.0;
        C( 5, 4 ) = 0.0;
        C( 5, 5 ) = c3;

    }

    /**
     * TO BE TESTED!!!
     */
    void VonMisesKinemVer3D::CalculateStress( const Vector& StrainVector, Vector& StressVector )
    {

        if ( StressVector.size() != 6 )
        {
            StressVector.resize( 6 );
        }

        noalias( StressVector ) = mCurrentStress;

    }


    void VonMisesKinemVer3D::CalculateConstitutiveMatrix( const Vector& StrainVector, Matrix& rResult )
    {
        noalias( rResult ) = mCtangent;

    }


    void VonMisesKinemVer3D::CalculateCauchyStresses( Vector& rCauchy_StressVector,
            const Matrix& rF,
            Vector& rPK2_StressVector, const Vector& GreenLagrangeStrainVector )
    {

        Matrix S = MathUtils<double>::StressVectorToTensor( rPK2_StressVector );
        double J = MathUtils<double>::Det3( rF );

        boost::numeric::ublas::bounded_matrix<double, 3, 3> temp;
        boost::numeric::ublas::bounded_matrix<double, 3, 3> aux;

        noalias( temp ) = prod( rF, S );
        noalias( aux ) = prod( mstemp, trans( rF ) );
        aux *= J;

        if ( rCauchy_StressVector.size() != 6 )
            rCauchy_StressVector.resize( 6 );

        rCauchy_StressVector[0] = aux( 0, 0 );

        rCauchy_StressVector[1] = aux( 1, 1 );

        rCauchy_StressVector[2] = aux( 2, 2 );

        rCauchy_StressVector[3] = aux( 1, 2 );

        rCauchy_StressVector[4] = aux( 1, 3 );

        rCauchy_StressVector[5] = aux( 2, 3 );

    }
} // Namespace Kratos
