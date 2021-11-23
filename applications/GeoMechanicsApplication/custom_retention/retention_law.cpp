// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

/* System includes */
#include <iostream>

/* External includes */

/* Project includes */
#include "custom_retention/saturated_law.h"

namespace Kratos
{

//-------------------------------------------------------------------------------------------------
RetentionLaw::RetentionLaw()
{
    // nothing
}

//-------------------------------------------------------------------------------------------------
RetentionLaw::~RetentionLaw()
{
    // nothing
}

//-------------------------------------------------------------------------------------------------
RetentionLaw::Pointer RetentionLaw::Clone() const
{
    return Kratos::make_shared<RetentionLaw>(*this);
}

//-------------------------------------------------------------------------------------------------
bool RetentionLaw::Has(const Variable<bool> &rThisVariable)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling Has<bool> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
bool RetentionLaw::Has(const Variable<int> &rThisVariable)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling Has<int> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
bool RetentionLaw::Has(const Variable<double> &rThisVariable)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling Has<double> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
bool RetentionLaw::Has(const Variable<Vector> &rThisVariable)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling Has<Vector> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
bool RetentionLaw::Has(const Variable<Matrix> &rThisVariable)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling Has<Matrix> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
bool RetentionLaw::Has(const Variable<array_1d<double, 3>> &rThisVariable)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling Has<array_1d<double, 3>> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
bool RetentionLaw::Has(const Variable<array_1d<double, 6>> &rThisVariable)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling Has<array_1d<double, 6>> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
bool& RetentionLaw::GetValue(const Variable<bool> &rThisVariable, bool &rValue)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling GetValue<bool> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
int& RetentionLaw::GetValue(const Variable<int> &rThisVariable, int &rValue)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling GetValue<int> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
double& RetentionLaw::GetValue(const Variable<double> &rThisVariable, double &rValue)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling GetValue<double> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
Vector& RetentionLaw::GetValue(const Variable<Vector> &rThisVariable, Vector &rValue)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling GetValue<Vector> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
Matrix& RetentionLaw::GetValue(const Variable<Matrix> &rThisVariable, Matrix &rValue)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling GetValue<Matrix> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
array_1d<double, 3>& RetentionLaw::
    GetValue(const Variable<array_1d<double, 3>> &rThisVariable,
             array_1d<double, 3> &rValue)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling GetValue<array_1d<double, 3>> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
array_1d<double, 6>& RetentionLaw::
    GetValue(const Variable<array_1d<double, 6>> &rThisVariable,
             array_1d<double, 6> &rValue)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling GetValue<array_1d<double, 6>> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
void RetentionLaw::SetValue(const Variable<bool> &rVariable,
                            const bool &Value,
                            const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling SetValue<bool> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
void RetentionLaw::SetValue(const Variable<int> &rVariable,
                            const int &Value,
                            const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR <<"calling SetValue<int> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
void RetentionLaw::SetValue(const Variable<double> &rVariable,
                            const double &rValue,
                            const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling SetValue<double> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
void RetentionLaw::SetValue(const Variable<Vector> &rVariable,
                            const Vector &rValue,
                            const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling SetValue<Vector> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
void RetentionLaw::SetValue(const Variable<Matrix> &rVariable,
                            const Matrix &rValue,
                            const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling SetValue<Matrix> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
void RetentionLaw::SetValue(const Variable<array_1d<double, 3>> &rVariable,
                            const array_1d<double, 3> &rValue,
                            const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling SetValue<array_1d<double, 3>> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
void RetentionLaw::SetValue(const Variable<array_1d<double, 6>> &rVariable,
                            const array_1d<double, 6> &rValue,
                            const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling SetValue<array_1d<double, 6>> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
bool& RetentionLaw::CalculateValue(Parameters &rParameters,
                                   const Variable<bool> &rThisVariable,
                                   bool &rValue)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling CalculateValue<bool> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
int& RetentionLaw::CalculateValue(Parameters &rParameters,
                                  const Variable<int> &rThisVariable,
                                  int &rValue)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling CalculateValue<int> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
double& RetentionLaw::CalculateValue(Parameters &rParameters,
                                     const Variable<double> &rThisVariable,
                                     double &rValue)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling CalculateValue<double> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
Vector& RetentionLaw::CalculateValue(Parameters &rParameters,
                                     const Variable<Vector> &rThisVariable,
                                     Vector &rValue)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling CalculateValue<Vector> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
Matrix& RetentionLaw::CalculateValue(Parameters &rParameters,
                                    const Variable<Matrix> &rThisVariable,
                                    Matrix &rValue)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling CalculateValue<Matrix> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
array_1d<double, 3>& RetentionLaw::
    CalculateValue(Parameters &rParameters,
                   const Variable<array_1d<double, 3>> &rVariable,
                   array_1d<double, 3> &rValue)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling CalculateValue<double, 3> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
array_1d<double, 6>& RetentionLaw::
    CalculateValue(Parameters &rParameters,
                   const Variable<array_1d<double, 6>> &rVariable,
                   array_1d<double, 6> &rValue)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling CalculateValue<double, 6> function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
double RetentionLaw::CalculateSaturation(Parameters &rParameters)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling CalculateSaturation function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
double RetentionLaw::CalculateEffectiveSaturation(Parameters &rParameters)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling CalculateEffectiveSaturation function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
double RetentionLaw::CalculateDerivativeOfSaturation(Parameters &rParameters)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling CalculateDerivativeOfSaturation function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
double RetentionLaw::CalculateRelativePermeability(Parameters &rParameters)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling CalculateRelativePermeability function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
double RetentionLaw::CalculateBishopCoefficient(Parameters &rParameters)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling CalculateBishopCoefficient function from base class ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------
void RetentionLaw::
    InitializeMaterial(const Properties& rMaterialProperties,
                       const GeometryType& rElementGeometry,
                       const Vector& rShapeFunctionsValues)
{
    // nothing
}

//-------------------------------------------------------------------------------------------------
void RetentionLaw::Initialize(Parameters &rParameters)
{
    // nothing
}

//-------------------------------------------------------------------------------------------------
void RetentionLaw::InitializeSolutionStep(Parameters &rParameters)
{
    // nothing
}

//-------------------------------------------------------------------------------------------------
void RetentionLaw::FinalizeSolutionStep(Parameters &rParameters)
{
    // nothing
}

//-------------------------------------------------------------------------------------------------
void RetentionLaw::Finalize(Parameters &rParameters)
{
    // nothing
}

//-------------------------------------------------------------------------------------------------
void RetentionLaw::ResetMaterial(const Properties &rMaterialProperties,
                                 const GeometryType &rElementGeometry,
                                 const Vector &rShapeFunctionsValues)
{
    // nothing
}

//-------------------------------------------------------------------------------------------------
int RetentionLaw::Check(const Properties &rMaterialProperties,
                        const ProcessInfo &rCurrentProcessInfo)
{
    // nothing
    return 1;
}

//-------------------------------------------------------------------------------------------------
void RetentionLaw::save(Serializer& rSerializer) const
{
    // there is no member variables to be saved
}

//-------------------------------------------------------------------------------------------------
void RetentionLaw::load(Serializer& rSerializer)
{
    // there is no member variables to be loaded
}

} /* namespace Kratos.*/
