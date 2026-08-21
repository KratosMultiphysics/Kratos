#pragma once

#include "includes/model_part.h"
#include "includes/define.h"
#include "includes/constitutive_law.h"

namespace Kratos
{

class KRATOS_API(IGA_APPLICATION) ComputeInterfaceTractionShell3pUtility
{
public:

    using SizeType = std::size_t;
    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(ComputeInterfaceTractionShell3pUtility);

    ComputeInterfaceTractionShell3pUtility() = default;

    virtual ~ComputeInterfaceTractionShell3pUtility() = default;

    static void ComputeAndSetInterfaceTraction(
        ModelPart& rModelPart,
        const Variable<array_1d<double,3>>& rTractionVariable);

private:

    enum class ConfigurationType
    {
        Reference,
        Current
    };

    struct KinematicVariables
    {
        KinematicVariables(const SizeType WorkingSpaceDimension);

        array_1d<double,3> a1;
        array_1d<double,3> a2;
        array_1d<double,3> a3;
        array_1d<double,3> a3_tilde;

        array_1d<double,3> t;
        array_1d<double,3> n;

        array_1d<double,2> n_contravariant;

        array_1d<double,3> a_ab_covariant;

        double dA;
    };

    struct ConstitutiveVariables
    {
        ConstitutiveVariables(const SizeType VoigtSize);

        Vector StrainVector;
        Vector StressVector;
        Matrix ConstitutiveMatrix;
    };

    static void CalculateKinematics(
        const Condition& rCondition,
        IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables,
        const Matrix& rShapeFunctionGradientValues,
        const ConfigurationType Configuration);

    static void CalculateTransformation(
        const KinematicVariables& rKinematicVariables,
        Matrix& rT,
        Matrix& rT_hat);

    static void CalculateConstitutiveVariables(
        const Condition& rCondition,
        const ProcessInfo& rCurrentProcessInfo,
        const IndexType IntegrationPointIndex,
        const Matrix& rT,
        const array_1d<double,3>& rReferenceMetric,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rConstitutiveVariables);

    static void CalculateTraction(
        array_1d<double,3>& rTraction,
        const Matrix& rT_hat,
        const array_1d<double,2>& rNContravariantVector,
        const KinematicVariables& rActualKinematic,
        const ConstitutiveVariables& rConstitutiveVariables);
};

} // namespace Kratos