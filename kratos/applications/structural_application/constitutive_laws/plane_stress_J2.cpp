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
 *   Last Modified by:    $Author: kazem $
 *   Date:                $Date: 2009-01-16 10:50:24 $
 *   Revision:            $Revision: 1.14 $
 *
 * ***********************************************************/


// System includes 
#include <iostream>

// External includes 
#include<cmath>

// Project includes 

#include "includes/define.h"
#include "constitutive_laws/plane_stress_J2.h"

#include "includes/constitutive_law.h"

#include "utilities/math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "structural_application.h"
#include "includes/properties.h"

namespace Kratos
{

    /**
     *	TO BE TESTED!!!
     */
    PlaneStressJ2::PlaneStressJ2()
    : ConstitutiveLaw< Node < 3 > >()
    {
        noalias(mOldPlasticStrain) = ZeroVector(3);
        noalias(mbeta_old) = ZeroVector(3);
        malpha_old = 0.0;

        malpha_current = 0.0;
        noalias(mCurrentPlasticStrain) = ZeroVector(3);
        noalias(mbeta_n1) = ZeroVector(3);


    }

    /**
     *	TO BE TESTED!!!
     */
    PlaneStressJ2::~PlaneStressJ2() { }

    bool PlaneStressJ2::Has(const Variable<double>& rThisVariable)
    {
        return false;
    }

    bool PlaneStressJ2::Has(const Variable<Vector>& rThisVariable)
    {
        return false;
    }

    bool PlaneStressJ2::Has(const Variable<Matrix>& rThisVariable)
    {
        return false;
    }

    double PlaneStressJ2::GetValue(const Variable<double>& rThisVariable)
    {
        return 0;
        //KRATOS_ERROR(std::logic_error, "double Variable case not considered", "");
    }

    Vector PlaneStressJ2::GetValue(const Variable<Vector>& rThisVariable)
    {
        KRATOS_ERROR(std::logic_error, "Vector Variable case not considered", "");
    }

    Matrix PlaneStressJ2::GetValue(const Variable<Matrix>& rThisVariable)
    {
        Matrix temp(1, 3);
        noalias(row(temp, 0)) = mOldPlasticStrain;
        return temp;
    }

    void PlaneStressJ2::SetValue(const Variable<double>& rThisVariable, const double& rValue,
            const ProcessInfo& rCurrentProcessInfo) { }

    void PlaneStressJ2::SetValue(const Variable<Vector>& rThisVariable, const Vector& rValue,
            const ProcessInfo& rCurrentProcessInfo) { }

    void PlaneStressJ2::SetValue(const Variable<Matrix>& rThisVariable, const Matrix& rValue,
            const ProcessInfo& rCurrentProcessInfo) { }

    void PlaneStressJ2::Calculate(const Variable<Matrix >& rVariable, Matrix& rResult,
            const ProcessInfo& rCurrentProcessInfo) { }

    /**
     *	TO BE TESTED!!!
     */
    void PlaneStressJ2::InitializeMaterial(const Properties& props,
            const GeometryType& geom,
            const Vector& ShapeFunctionsValues) { }

    /**
     *	TO BE TESTED!!!
     */
    void PlaneStressJ2::CalculateElasticMatrix(boost::numeric::ublas::bounded_matrix<double, 3, 3 > & C)
    {
        double c1 = mE / (1.00 - mNU * mNU);
        double c2 = c1 * mNU;
        double c3 = 0.5 * mE / (1.0 + mNU);

        C(0, 0) = c1;
        C(0, 1) = c2;
        C(0, 2) = 0.0;
        C(1, 0) = c2;
        C(1, 1) = c1;
        C(1, 2) = 0.0;
        C(2, 0) = 0.0;
        C(2, 1) = 0.0;
        C(2, 2) = c3;
    }

    void PlaneStressJ2::CalculateInverseElasticMatrix(boost::numeric::ublas::bounded_matrix<double, 3, 3 > & Cinv)
    {
        double c1 = 1.0 / mE;
        double c2 = -mNU * c1;
        double c3 = (1 + mNU) / (0.5 * mE);

        Cinv(0, 0) = c1;
        Cinv(0, 1) = c2;
        Cinv(0, 2) = 0.0;
        Cinv(1, 0) = c2;
        Cinv(1, 1) = c1;
        Cinv(1, 2) = 0.0;
        Cinv(2, 0) = 0.0;
        Cinv(2, 1) = 0.0;
        Cinv(2, 2) = c3;
    }

    void PlaneStressJ2::CalculateP(boost::numeric::ublas::bounded_matrix<double, 3, 3 > & P)
    {
        P(0, 0) = 0.6666666666666667;
        P(0, 1) = -0.3333333333333333;
        P(0, 2) = 0.0;
        P(1, 0) = -0.3333333333333333;
        P(1, 1) = 0.6666666666666667;
        P(1, 2) = 0.0;
        P(2, 0) = 0.0;
        P(2, 1) = 0.0;
        P(2, 2) = 2.0;
    }

    void PlaneStressJ2::CalculateMaterialResponse(
            const Vector& StrainVector,
            Vector& StressVector,
            Matrix& algorithmicTangent,
            bool calculate_stress_flag,
            bool calculate_tangent_flag,
            bool save_internal_variables
            )
    {
        KRATOS_TRY

                double theta = 0.0;

        //resize output quantities
        if (StressVector.size() != 3 && calculate_stress_flag == true) StressVector.resize(3, false);
        if (algorithmicTangent.size1() != 3 && calculate_tangent_flag == true) algorithmicTangent.resize(3, 3, false);


        array_1d<double, 3 > elastic_strain;
        noalias(elastic_strain) = StrainVector;
        noalias(elastic_strain) -= mOldPlasticStrain;
        //KRATOS_WATCH(StrainVector);
        //KRATOS_WATCH(mOldPlasticStrain);

        boost::numeric::ublas::bounded_matrix<double, 3, 3 > C, Cinv, P;
        CalculateElasticMatrix(C);
        CalculateInverseElasticMatrix(Cinv);
        CalculateP(P);

//                KRATOS_WATCH(C);

        array_1d<double, 3 > s_trial = prod(C, elastic_strain);
        array_1d<double, 3 > xi_trial = s_trial;
        noalias(s_trial) -= mbeta_old;

        //                KRATOS_WATCH(compute_f_trial(xi_trial, malpha_old));

//        double fbar2 = fbar_2(0.0, xi_trial);
//        double fbar_value = sqrt(fbar_2(0.0, xi_trial));
//        double r2_value = R_2(0.0, fbar_value, malpha_old);
//        double aaa = fbar_value - sqrt(r2_value);
        //                KRATOS_WATCH(sqrt(r2_value));
        //                KRATOS_WATCH(fbar_value);
        //                KRATOS_WATCH(sqrt(r2_value));
        //                KRATOS_WATCH(aaa);

        double H1 = (1.0 - mtheta) * mH;
        ////        KRATOS_WATCH(xi_trial);
        ////        KRATOS_WATCH(mbeta_old)
        if (compute_f_trial(xi_trial, malpha_old) < 0) //elastic case
        {
            if (calculate_stress_flag == true) noalias(StressVector) = s_trial;
            if (calculate_tangent_flag == true) noalias(algorithmicTangent) = C;
            //note that in this case internal variables are not modified!

        } else
        {
            //algorithm copied identically from the Simo Hughes, BOX3.3, pag 130
            double dgamma = ComputeDGamma(xi_trial, malpha_old);
            double ccc = 0.6666666666666667 * dgamma*H1;

            //            KRATOS_WATCH(dgamma);
            //            KRATOS_WATCH(xi_trial);
            //            KRATOS_WATCH(malpha_old);

            //calculate XImat
            //note that here i reuse the C as auxiliary variable as i don't need it anymore
            boost::numeric::ublas::bounded_matrix<double, 3, 3 > XImat;
            noalias(C) = Cinv;
            noalias(C) += (dgamma / (1.0 + ccc)) * P;
            double detC;
            InvertMatrix(C, XImat, detC);
            //            KRATOS_WATCH(XImat);


            array_1d<double, 3 > aux, xi;
            noalias(aux) = prod(Cinv, xi_trial);
            noalias(xi) = prod(XImat, aux);
            xi /= (1.0 + ccc);

            //            KRATOS_WATCH(compute_f_trial(xi, malpha_old));

            noalias(mbeta_n1) = mbeta_old;
            noalias(mbeta_n1) += ccc*xi;

            array_1d<double, 3 > stress;
            noalias(stress) = xi;
            noalias(stress) += mbeta_n1;
            if (calculate_stress_flag == true) noalias(StressVector) = s_trial;

            malpha_current = malpha_old + sqrt(0.6666666666666667) * dgamma * sqrt(fbar_2(dgamma, xi));


            //KRATOS_WATCH(StressVector);
            noalias(aux) = prod(P, xi);

            noalias(mCurrentPlasticStrain) = mOldPlasticStrain;
            noalias(mCurrentPlasticStrain) += dgamma*aux;



            if (calculate_tangent_flag == true)
            {

                //compute tangent
                array_1d<double, 3 > XPXi = prod(XImat, aux);

                double K1_n1 = theta*mH; //msigma_y + theta * mH*alpha;
                //                double K1_n1 = msigma_y + theta * mH*alpha;
                double theta1 = 1.0 + 0.6666666666666667 * H1*dgamma;
                double theta2 = 1.0 - 0.6666666666666667 * K1_n1*dgamma;
                double beta_val = inner_prod(xi, aux);
                beta_val *= 0.6666666666666667 * (theta1 / theta2) * (K1_n1 * theta1 + H1 * theta2);

                double denom = inner_prod(aux, XPXi);
                denom += beta_val;
                denom = sqrt(denom);
                XPXi /= denom;

                noalias(algorithmicTangent) = XImat;
                noalias(algorithmicTangent) -= outer_prod(XPXi, XPXi);
//                KRATOS_WATCH(XImat);
//                KRATOS_WATCH(XPXi);
//
//                                KRATOS_WATCH(algorithmicTangent);
            }
            //KRATOS_WATCH(algorithmicTangent);



        }

        KRATOS_CATCH("")
    }

    //***********************************************************************************************
    //***********************************************************************************************

    void PlaneStressJ2::UpdateMaterial(const Vector& StrainVector,
            const Properties& props,
            const GeometryType& geom,
            const Vector& ShapeFunctionsValues,
            const ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        mE = props.GetValue(YOUNG_MODULUS);
        mNU = props.GetValue(POISSON_RATIO);
        msigma_y = props.GetValue(YIELD_STRESS);
        mH = 0.0;
        mtheta = 0.0;

        KRATOS_CATCH("");
    }

    //***********************************************************************************************
    //***********************************************************************************************

    void PlaneStressJ2::InvertMatrix(const boost::numeric::ublas::bounded_matrix<double, 3, 3 > & InputMatrix,
            boost::numeric::ublas::bounded_matrix<double, 3, 3 > & InvertedMatrix,
            double& InputMatrixDet) //VERIFIED!!!
    {
        //filling the inverted matrix with the algebraic complements
        //first column
        InvertedMatrix(0, 0) = InputMatrix(1, 1) * InputMatrix(2, 2) - InputMatrix(1, 2) * InputMatrix(2, 1);
        InvertedMatrix(1, 0) = -InputMatrix(1, 0) * InputMatrix(2, 2) + InputMatrix(1, 2) * InputMatrix(2, 0);
        InvertedMatrix(2, 0) = InputMatrix(1, 0) * InputMatrix(2, 1) - InputMatrix(1, 1) * InputMatrix(2, 0);

        //second column
        InvertedMatrix(0, 1) = -InputMatrix(0, 1) * InputMatrix(2, 2) + InputMatrix(0, 2) * InputMatrix(2, 1);
        InvertedMatrix(1, 1) = InputMatrix(0, 0) * InputMatrix(2, 2) - InputMatrix(0, 2) * InputMatrix(2, 0);
        InvertedMatrix(2, 1) = -InputMatrix(0, 0) * InputMatrix(2, 1) + InputMatrix(0, 1) * InputMatrix(2, 0);

        //third column
        InvertedMatrix(0, 2) = InputMatrix(0, 1) * InputMatrix(1, 2) - InputMatrix(0, 2) * InputMatrix(1, 1);
        InvertedMatrix(1, 2) = -InputMatrix(0, 0) * InputMatrix(1, 2) + InputMatrix(0, 2) * InputMatrix(1, 0);
        InvertedMatrix(2, 2) = InputMatrix(0, 0) * InputMatrix(1, 1) - InputMatrix(0, 1) * InputMatrix(1, 0);

        //calculation of determinant (of the input matrix)
        InputMatrixDet = InputMatrix(0, 0) * InvertedMatrix(0, 0) + InputMatrix(0, 1) * InvertedMatrix(1, 0) + InputMatrix(0, 2) * InvertedMatrix(2, 0);

        //finalizing the calculation of the inverted matrix
        InvertedMatrix /= InputMatrixDet;
    }

    //**********************************************************************
    //**********************************************************************

    void PlaneStressJ2::Calculate(const Variable<double>& rVariable,
            double& Output,
            const ProcessInfo& rCurrentProcessInfo) { }

    double PlaneStressJ2::fbar_2(const double dgamma, const array_1d<double, 3 > & xi)
    {
        double mu = 0.5 * mE / (1.0 - mNU);
        double H1 = (1.0 - mtheta) * mH;
        double coeff = mE / (3.0 * (1.0 - mNU));
        double value = pow(xi[0] + xi[1], 2)*0.16666666666667 / pow((1.0 + dgamma * (coeff + 0.666666666666667 * H1)), 2);
        value += (0.5 * pow(xi[0] - xi[1], 2) + 2.0 * xi[2] * xi[2]) / pow((1.0 + dgamma * (2.0 * mu + 0.666666666666667 * H1)), 2);

        return value;
    }

    double PlaneStressJ2::R_2(const double dgamma, const double fbar_value, const double alpha)
    {
        double KK = K(alpha + 0.81649658092773 * dgamma * fbar_value);
        return 0.3333333333333 * KK * KK;
    }

    double PlaneStressJ2::f2(const double dgamma, const array_1d<double, 3 > & xi, const double alpha)
    {
        double fbar2 = fbar_2(dgamma, xi);
        double fbar_value = sqrt(fbar_2(dgamma, xi));
        double value = -R_2(dgamma, fbar_value, alpha);
        value += 0.5 * fbar2;
        return value;
    }

    double PlaneStressJ2::ComputeDGamma(const array_1d<double, 3 > & xi, const double alpha)
    {
        double x0 = 0.0;
        double g0 = f2(x0, xi, alpha);

        double x1 = 1e-14;
        double g1 = f2(x1, xi, alpha);

        double tol = 1e-9 * g0;
        unsigned it = 0;
        double increment = 1.0;
        const unsigned int itmax = 1000;
        while (fabs(g1) > tol && fabs(increment)>1e-14 && it++ < itmax)
        {
            increment = -g1 * (x1 - x0) / (g1 - g0);
            x0 = x1;
            g0 = g1;
            x1 += increment;

            g1 = f2(x1, xi, alpha);

//                        std::cout << increment << " " << fabs(g1) << " " << tol  <<std::endl;
        }

        if (it == itmax)
            std::cout << "reached itmax iterations" << std::endl;

        return x1;
    }

    double PlaneStressJ2::K(const double alpha)
    {
        return (msigma_y + mtheta * mH * alpha);
    }

    double PlaneStressJ2::compute_f_trial(const array_1d<double, 3 > & xi_trial, const double alpha)
    {
        array_1d<double, 3 > eta;
        eta[0] = 0.70710678118655 * (xi_trial[0] + xi_trial[1]);
        eta[1] = 0.70710678118655 * (-xi_trial[0] + xi_trial[1]);
        eta[2] = xi_trial[2];

        double f = eta[0]*0.3333333333333333333 * (2.0 * eta[0] - eta[1]);
        f += eta[1]*0.3333333333333333333 * (-eta[0] + 2.0 * eta[1]);
        f += 2.0 * eta[2] * eta[2];
        f = sqrt(f);

        f -= 0.81649658092773 * K(alpha);

        return f;

    }

    void PlaneStressJ2::FinalizeSolutionStep(const Properties& props,
            const GeometryType& geom,
            const Vector& ShapeFunctionsValues,
            const ProcessInfo& CurrentProcessInfo)
    {
        malpha_old = malpha_current;
        noalias(mOldPlasticStrain) = mCurrentPlasticStrain;
        noalias(mbeta_old) = mbeta_n1;

    }

    void PlaneStressJ2::CalculateStress( const Vector& StrainVector,
                                          Vector& StressVector)
    {
        bool calculate_stress_flag = true;
        bool calculate_tangent_flag = false;
        bool save_internal_variables = false;
        Matrix algorithmicTangent(0,0);
        CalculateMaterialResponse(StrainVector,StressVector,algorithmicTangent,calculate_stress_flag,calculate_tangent_flag,save_internal_variables);
    }

        void PlaneStressJ2::CalculateConstitutiveMatrix( const Vector& StrainVector,
                    Matrix& ElasticityTensor )
    {
        bool calculate_stress_flag = false;
        bool calculate_tangent_flag = true;
        bool save_internal_variables = false;
        Vector StressVector(0);
        CalculateMaterialResponse(StrainVector,StressVector,ElasticityTensor,calculate_stress_flag,calculate_tangent_flag,save_internal_variables);
    }

        void PlaneStressJ2::CalculateStressAndTangentMatrix(Vector& StressVector,
                    const Vector& StrainVector,
                    Matrix& algorithmicTangent)
    {
        bool calculate_stress_flag = true;
        bool calculate_tangent_flag = true;
        bool save_internal_variables = false;
        CalculateMaterialResponse(StrainVector,StressVector,algorithmicTangent,calculate_stress_flag,calculate_tangent_flag,save_internal_variables);
    }


} // Namespace Kratos
