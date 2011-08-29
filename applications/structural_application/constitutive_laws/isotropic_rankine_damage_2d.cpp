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
#include "constitutive_laws/isotropic_rankine_damage_2d.h"

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
    IsotropicRankineDamage2D::IsotropicRankineDamage2D()
    : ConstitutiveLaw()
    {
    }

    /**
     *	TO BE TESTED!!!
     */
    IsotropicRankineDamage2D::~IsotropicRankineDamage2D()
    {
    }

    bool IsotropicRankineDamage2D::Has(const Variable<double>& rThisVariable)
    {
        return false;
    }

    bool IsotropicRankineDamage2D::Has(const Variable<Vector>& rThisVariable)
    {
        return false;
    }

    bool IsotropicRankineDamage2D::Has(const Variable<Matrix>& rThisVariable)
    {
        return false;
    }

    void IsotropicRankineDamage2D::SetValue(const Variable<double>& rThisVariable, const double& rValue,
            const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void IsotropicRankineDamage2D::SetValue(const Variable<Vector>& rThisVariable, const Vector& rValue,
            const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void IsotropicRankineDamage2D::SetValue(const Variable<Matrix>& rThisVariable, const Matrix& rValue,
            const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void IsotropicRankineDamage2D::Calculate(const Variable<Matrix >& rVariable, Matrix& rResult,
            const ProcessInfo& rCurrentProcessInfo)
    {
    }

    double& IsotropicRankineDamage2D::GetValue(const Variable<double>& rThisVariable, double& rValue)
    {
        if(rThisVariable==DAMAGE)
        {
//            std::cout << this << " md " << md << std::endl;
                rValue = md;
        }

        return rValue;

    }

    /**
     *	TO BE TESTED!!!
     */
    void IsotropicRankineDamage2D::InitializeMaterial(const Properties& props,
            const GeometryType& geom,
            const Vector& ShapeFunctionsValues)
    {
        mrold = props[YIELD_STRESS];
        mr = mrold;
    }

    /**
     *	TO BE TESTED!!!
     */
    void IsotropicRankineDamage2D::CalculateElasticMatrix(Matrix& C, const double E, const double NU)
    {
        double c1 = E / (1.00 - NU * NU);
        double c2 = c1 * NU;
        double c3 = 0.5 * E / (1 + NU);

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


    //**********************************************************************

    void IsotropicRankineDamage2D::CalculateCauchyStresses(
            Vector& rCauchy_StressVector,
            const Matrix& rF,
            const Vector& rPK2_StressVector,
            const Vector& rGreenLagrangeStrainVector)
    {
        KRATOS_ERROR(std::logic_error,"method not yet implemented","")
    }

    //**********************************************************************
    //**********************************************************************

    void IsotropicRankineDamage2D::Calculate(const Variable<double>& rVariable,
            double& Output,
            const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void IsotropicRankineDamage2D::CalculateMaterialResponse(const Vector& StrainVector,
            const Matrix& DeformationGradient,
            Vector& StressVector,
            Matrix& AlgorithmicTangent,
            const ProcessInfo& CurrentProcessInfo,
            const Properties& props,
            const GeometryType& geom,
            const Vector& ShapeFunctionsValues,
            bool CalculateStresses,
            int CalculateTangent,
            bool SaveInternalVariables)
    {

        //get material parameters
        const double E  = props[YOUNG_MODULUS];
        const double NU = props[POISSON_RATIO];
        const double Gf = props[FRACTURE_ENERGY]; //***************************
        const double sigma0 = props[YIELD_STRESS]; //***************************
        const double r0 = sigma0; //***************************
        
        //compute elastic stress
        Matrix Cel(3,3);
        Vector stress(3);
        CalculateElasticMatrix(Cel, E, NU);
        noalias(stress) = prod(Cel,StrainVector);

        //compute stress max eigenvalue
        const double sigma_m = 0.5*(stress[0]+  stress[1]);
        const double R = sqrt(stress[2]*stress[2] + 0.25*(stress[0]-stress[1])*(stress[0]-stress[1]) );
        double smax = sigma_m + R;

        //compute tau
        double tau = smax;
        if(tau < 0) tau=0.0;

        //compute actualized damage indicator
        double r = mrold;
        if(tau > r) r=tau;
//        KRATOS_WATCH(tau);
//        KRATOS_WATCH(r);
//        KRATOS_WATCH(r0);
//        KRATOS_WATCH(mr);
//        KRATOS_WATCH(mrold);



        //compute element lenght
        double A = geom.Area();
        const double he = sqrt(2.0*A);

        const double lch=2.0*he;
        const double ls =2.0*E*Gf/(sigma0*sigma0);
        
        double Hs = lch/(ls -lch);
        if(Hs < 0.0) Hs=1.0e10;

//         KRATOS_WATCH(Hs);

        double d=0.0;
        if(r>=r0)
            d = 1.0 - r0/r*exp(-2.0*Hs*(r-r0)/r0);
        if(d > 0.999)
            d = 0.999;

//        KRATOS_WATCH(d);
        

        //write outputs as needed
        if(CalculateStresses == true)
            noalias(StressVector) = (1.0-d)*stress;
        if(CalculateTangent == 1)
            noalias(AlgorithmicTangent) = (1.0-d)*Cel;

        //keep trace of internal variables if allowed, otherwise reset them
        if(SaveInternalVariables == false)
            mr = mrold;
        else
        {        
            mr = r;

            //save variable for output
            md=d;
        }

//        KRATOS_WATCH(md);
//std::cout << this  << " inside calculate md " << md << std::endl;

    }

    /// Turn back information as a string.
    std::string IsotropicRankineDamage2D::Info() const
    {
        std::stringstream buffer;
        buffer << "Isotropic Rankine Damage 2D" << std::endl;
        return buffer.str();
    }







} // Namespace Kratos
