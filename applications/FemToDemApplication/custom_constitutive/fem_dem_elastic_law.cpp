//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez

// System includes

// External includes

// Project includes
#include "custom_constitutive/fem_dem_elastic_law.hpp"

#include "fem_to_dem_application_variables.h"

namespace Kratos
{
    //******************************CONSTRUCTOR*******************************************
//************************************************************************************

FemDemElasticLaw::FemDemElasticLaw()
    : LinearElastic3DLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

FemDemElasticLaw::FemDemElasticLaw(const FemDemElasticLaw& rOther)
    : LinearElastic3DLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer FemDemElasticLaw::Clone() const
{
    return Kratos::make_shared<FemDemElasticLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

FemDemElasticLaw::~FemDemElasticLaw()
{
}

double& FemDemElasticLaw::CalculateValue(Parameters& rParameterValues,
                                         const Variable<double>& rThisVariable,
                                         double& rValue) 
{
    if (rThisVariable == DAMAGE_ELEMENT) {
        return rValue = mDamage;
    } else if(rThisVariable == STRESS_THRESHOLD) {
        return rValue = mThreshold;
    }
}

double& FemDemElasticLaw::GetValue(const Variable<double>& rThisVariable,
                                   double& rValue) 
{
    return LinearElastic3DLaw::GetValue(rThisVariable, rValue);
}

Vector& FemDemElasticLaw::GetValue(const Variable<Vector>& rThisVariable,
                                   Vector& rValue) 
{
	return rValue;
}

Matrix& FemDemElasticLaw::GetValue(const Variable<Matrix>& rThisVariable,
                                   Matrix& rValue)
{
	return rValue;
}

void FemDemElasticLaw::SetValue(const Variable<double>& rThisVariable,
                                const double& rValue,
                                const ProcessInfo& rCurrentProcessInfo) 
{
    if (rThisVariable == DETERMINANT_F) {
        mDeterminantF0 = rValue;
    } else if(rThisVariable == DAMAGE_ELEMENT) {
        mDamage = rValue;
    } else if (rThisVariable == STRESS_THRESHOLD) {
        mThreshold = rValue;
    }
}

void FemDemElasticLaw::SetValue(const Variable<Vector>& rThisVariable,
                                const Vector& rValue,
                                const ProcessInfo& rCurrentProcessInfo) 
{
    LinearElastic3DLaw::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

void FemDemElasticLaw::SetValue(const Variable<Matrix>& rThisVariable,
                                const Matrix& rValue,
                                const ProcessInfo& rCurrentProcessInfo) 
{
    LinearElastic3DLaw::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

} // Namespace Kratos