//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/newtonian_2d_law.h"

#include "fluid_dynamics_application_variables.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

Newtonian2DLaw::Newtonian2DLaw()
    : ConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

Newtonian2DLaw::Newtonian2DLaw(const Newtonian2DLaw& rOther)
    : ConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer Newtonian2DLaw::Clone() const
{
    return Kratos::make_shared<Newtonian2DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

Newtonian2DLaw::~Newtonian2DLaw()
{
}


ConstitutiveLaw::SizeType Newtonian2DLaw::WorkingSpaceDimension() {
    return 2;
}

ConstitutiveLaw::SizeType Newtonian2DLaw::GetStrainSize() {
    return 3;
}


void  Newtonian2DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{

    //-----------------------------//

    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    //CheckParameters(rValues);

    //b.- Get Values to compute the constitutive law:
    Flags &Options = rValues.GetOptions();

    const Properties& MaterialProperties  = rValues.GetMaterialProperties();

    Vector& S = rValues.GetStrainVector(); //using the short name S to reduce the lenght of the expressions
    Vector& StressVector = rValues.GetStressVector();

    //-----------------------------//

    //1.- Lame constants
    const double mu = MaterialProperties[DYNAMIC_VISCOSITY];
    const double trS = S[0]+S[1];
    const double eps_vol = trS/3.0;

    //computation of stress
    StressVector[0] = 2.0*mu*(S[0] - eps_vol);
    StressVector[1] = 2.0*mu*(S[1] - eps_vol);
    StressVector[2] = mu*S[2];

    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
        Matrix& C = rValues.GetConstitutiveMatrix();

        noalias(C) = ZeroMatrix(3,3);

        C(0,0) = 4.0/3.0*mu;
        C(0,1) = -2.0/3.0*mu;
        C(1,0) = -2.0/3.0*mu;
        C(1,1) = 4.0/3.0*mu;
        C(2,2) = mu;

    }

    this->mViscosity = mu;

}


//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void Newtonian2DLaw::GetLawFeatures(Features& rFeatures)
{
    	//Set the type of law
	rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
	rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

	//Set the strain size
	rFeatures.mStrainSize = 3;

	//Set the spacedimension
	rFeatures.mSpaceDimension = 2;

}

//******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW*************************
//************************************************************************************

// bool Newtonian2DLaw::CheckParameters(Parameters& rValues)
// {
//     return rValues.CheckAllParameters();
// }    bool& GetValue(const Variable<bool>& rThisVariable, bool& rValue) override;

int Newtonian2DLaw::Check(const Properties& rMaterialProperties,
                          const GeometryType& rElementGeometry,
                          const ProcessInfo& rCurrentProcessInfo)
{

    if(DYNAMIC_VISCOSITY.Key() == 0 || rMaterialProperties[DYNAMIC_VISCOSITY]<= 0.00)
        KRATOS_THROW_ERROR( std::invalid_argument,"DYNAMIC_VISCOSITY has Key zero or invalid value ", "" )

    return 0;

}

bool& Newtonian2DLaw::GetValue(const Variable<bool>& rThisVariable,bool& rValue) {}
int& Newtonian2DLaw::GetValue(const Variable<int>& rThisVariable,int& rValue) {}

double& Newtonian2DLaw::GetValue(const Variable<double>& rThisVariable, double& rValue) {
    return this->mViscosity;
}

Vector& Newtonian2DLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue) {}
Matrix& Newtonian2DLaw::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue) {}
array_1d<double, 3 > & Newtonian2DLaw::GetValue(const Variable<array_1d<double, 3 > >& rVariable,array_1d<double, 3 > & rValue) {}
array_1d<double, 6 > & Newtonian2DLaw::GetValue(const Variable<array_1d<double, 6 > >& rVariable, array_1d<double, 6 > & rValue) {}


} // Namespace Kratos
