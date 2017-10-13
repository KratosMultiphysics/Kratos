/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
 *   Last Modified by:    $Author:   Massimo Petracca$
 *   Date:                $Date:     30-10-2013$
 *   Revision:            $Revision: 1.0$
 *
 * ***********************************************************/

#include "includes/dummy_constitutive_law.h"


namespace Kratos
{

    DummyConstitutiveLaw::DummyConstitutiveLaw() 
        : ConstitutiveLaw()
    {
    }
 
    ConstitutiveLaw::Pointer DummyConstitutiveLaw::Clone() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }

    DummyConstitutiveLaw::SizeType DummyConstitutiveLaw::WorkingSpaceDimension()
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }

    DummyConstitutiveLaw::SizeType DummyConstitutiveLaw::GetStrainSize()
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }

    bool DummyConstitutiveLaw::Has(const Variable<double>& rThisVariable)
    {
        return false;
    }
    
    bool DummyConstitutiveLaw::Has(const Variable<Vector>& rThisVariable)
    {
        return false;
    }

    bool DummyConstitutiveLaw::Has(const Variable<Matrix>& rThisVariable)
    {
        return false;
    }

    bool DummyConstitutiveLaw::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
    {
        return false;
    }
    
    bool DummyConstitutiveLaw::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
    {
        return false;
    }

    double& DummyConstitutiveLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
    {
        return rValue;
    }

    Vector& DummyConstitutiveLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
    {
        return rValue;
    }

    Matrix& DummyConstitutiveLaw::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
    {
        return rValue;
    }

    array_1d<double, 3 > & DummyConstitutiveLaw::GetValue(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & rValue)
    {
        return rValue;
    }

    array_1d<double, 6 > & DummyConstitutiveLaw::GetValue(const Variable<array_1d<double, 6 > >& rVariable, array_1d<double, 6 > & rValue)
    {
        return rValue;
    }

     void DummyConstitutiveLaw::SetValue(const Variable<double>& rVariable,
                                         const double& rValue,
                                         const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }

    void DummyConstitutiveLaw::SetValue(const Variable<Vector >& rVariable,
                                        const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }

    void DummyConstitutiveLaw::SetValue(const Variable<Matrix >& rVariable,
                                        const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }

    void DummyConstitutiveLaw::SetValue(const Variable<array_1d<double, 3 > >& rVariable,
                                         const array_1d<double, 3 > & rValue,
                                         const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }

    void DummyConstitutiveLaw::SetValue(const Variable<array_1d<double, 6 > >& rVariable,
                                        const array_1d<double, 6 > & rValue,
                                        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }

    bool DummyConstitutiveLaw::ValidateInput(const Properties& rMaterialProperties)
    {
        return false;
    }

    DummyConstitutiveLaw::StrainMeasure DummyConstitutiveLaw::GetStrainMeasure()
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }
    
    DummyConstitutiveLaw::StressMeasure DummyConstitutiveLaw::GetStressMeasure()
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }

    bool DummyConstitutiveLaw::IsIncremental()
    {
        return false;
    }

    void DummyConstitutiveLaw::InitializeMaterial(const Properties& rMaterialProperties,
                                                  const GeometryType& rElementGeometry,
                                                  const Vector& rShapeFunctionsValues)
    {
    }

    void DummyConstitutiveLaw::InitializeSolutionStep(const Properties& rMaterialProperties,
                                                      const GeometryType& rElementGeometry,
                                                      const Vector& rShapeFunctionsValues,
                                                      const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void DummyConstitutiveLaw::FinalizeSolutionStep(const Properties& rMaterialProperties,
                                                    const GeometryType& rElementGeometry,
                                                    const Vector& rShapeFunctionsValues,
                                                    const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void DummyConstitutiveLaw::InitializeNonLinearIteration(const Properties& rMaterialProperties,
                                                            const GeometryType& rElementGeometry,
                                                            const Vector& rShapeFunctionsValues,
                                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }

    void DummyConstitutiveLaw::FinalizeNonLinearIteration(const Properties& rMaterialProperties,
                                                          const GeometryType& rElementGeometry,
                                                          const Vector& rShapeFunctionsValues,
                                                          const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }

    void DummyConstitutiveLaw::CalculateMaterialResponsePK1 (Parameters& rValues)
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }

    void DummyConstitutiveLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }

    void DummyConstitutiveLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }

    void DummyConstitutiveLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }

    void DummyConstitutiveLaw::FinalizeMaterialResponsePK1 (Parameters& rValues)
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }

    void DummyConstitutiveLaw::FinalizeMaterialResponsePK2 (Parameters& rValues)
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }

    void DummyConstitutiveLaw::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }

    void DummyConstitutiveLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }

    void DummyConstitutiveLaw::ResetMaterial(const Properties& rMaterialProperties,
                                             const GeometryType& rElementGeometry,
                                             const Vector& rShapeFunctionsValues)
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }

    void DummyConstitutiveLaw::GetLawFeatures(Features& rFeatures)
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }

    int DummyConstitutiveLaw::Check(const Properties& rMaterialProperties,
                                    const GeometryType& rElementGeometry,
                                    const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        return 0;
        KRATOS_CATCH("");
    }

    //*** OUTDATED METHODS: ***//

    void DummyConstitutiveLaw::CalculateMaterialResponse(const Vector& StrainVector,
                                                         const Matrix& DeformationGradient,
                                                         Vector& StressVector,
                                                         Matrix& AlgorithmicTangent,
                                                         const ProcessInfo& rCurrentProcessInfo,
                                                         const Properties& rMaterialProperties,
                                                         const GeometryType& rElementGeometry,
                                                         const Vector& rShapeFunctionsValues,
                                                         bool CalculateStresses,
                                                         int CalculateTangent,
                                                         bool SaveInternalVariables)
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }

    void DummyConstitutiveLaw::CalculateVolumetricResponse(const double VolumetricStrain,
                                                           const Matrix& DeformationGradient,
                                                           double& VolumetricStress,
                                                           double& AlgorithmicBulk,
                                                           const ProcessInfo& rCurrentProcessInfo,
                                                           const Properties& rMaterialProperties,
                                                           const GeometryType& rElementGeometry,
                                                           const Vector& rShapeFunctionsValues,
                                                           bool CalculateStresses,
                                                           int CalculateTangent,
                                                           bool SaveInternalVariables)
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }

    void DummyConstitutiveLaw::CalculateDeviatoricResponse(const Vector& StrainVector,
                                                           const Matrix& DeformationGradient,
                                                           Vector& StressVector,
                                                           Matrix& AlgorithmicTangent,
                                                           const ProcessInfo& rCurrentProcessInfo,
                                                           const Properties& rMaterialProperties,
                                                           const GeometryType& rElementGeometry,
                                                           const Vector& rShapeFunctionsValues,
                                                           bool CalculateStresses,
                                                           int CalculateTangent,
                                                           bool SaveInternalVariables)
    {
        KRATOS_THROW_ERROR(std::logic_error, "DUMMY C.LAW TEMPLATE!", "");
    }

    void DummyConstitutiveLaw::CalculateCauchyStresses(Vector& Cauchy_StressVector,
                                                       const Matrix& F,
                                                       const Vector& PK2_StressVector,
                                                       const Vector& GreenLagrangeStrainVector)
    {
    }




} /* namespace Kratos.*/
