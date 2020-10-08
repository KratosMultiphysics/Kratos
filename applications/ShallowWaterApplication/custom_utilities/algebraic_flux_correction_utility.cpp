//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "algebraic_flux_correction_utility.h"
#include "shallow_water_application_variables.h"

namespace Kratos
{
    AlgebraicFluxCorrectionUtility::AlgebraicFluxCorrectionUtility(
        ModelPart& rModelPart,
        Parameters ThisParameters) :
        mrModelPart(rModelPart)
    {
        ThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
        mRebuildLevel = ThisParameters["rebuild_level"].GetInt();

        size_t number_of_elements = mrModelPart.NumberOfElements();
        mHighOrderValues.reserve(number_of_elements);
        mLowOrderValues.reserve(number_of_elements);
        mAlgebraicFluxCorrections.reserve(number_of_elements);
        mElementalDofs.reserve(number_of_elements);

        GetElementalDofList();
    }

    void AlgebraicFluxCorrectionUtility::SetProcessInfoHighOrderFlags()
    {
        ProcessInfo& r_process_info = mrModelPart.GetProcessInfo();
        r_process_info.GetValue(LUMPED_MASS_FACTOR) = 1.0;
    }

    void AlgebraicFluxCorrectionUtility::SetProcessInfoLowOrderFlags()
    {
        ProcessInfo& r_process_info = mrModelPart.GetProcessInfo();
        r_process_info.GetValue(LUMPED_MASS_FACTOR) = 0.0;
    }

    void AlgebraicFluxCorrectionUtility::GetElementalHighOrderValues()
    {
        #pragma omp parallel for
        for (size_t i = 0; i < mrModelPart.NumberOfElements(); ++i)
        {
            (mrModelPart.ElementsBegin() + i)->GetValuesVector(*(mHighOrderValues.begin() + i), 0);
        }
    }

    void AlgebraicFluxCorrectionUtility::GetElementalLowOrderValues()
    {
        #pragma omp parallel for
        for (size_t i = 0; i < mrModelPart.NumberOfElements(); ++i)
        {
            (mrModelPart.ElementsBegin() + i)->GetValuesVector(*(mLowOrderValues.begin() + i), 0);
        }
    }

    void AlgebraicFluxCorrectionUtility::GetElementalPreviousValues()
    {
        #pragma omp parallel for
        for (size_t i = 0; i < mrModelPart.NumberOfElements(); ++i)
        {
            (mrModelPart.ElementsBegin() + i)->GetValuesVector(*(mPreviousValues.begin() + i), 0);
        }
    }

    void AlgebraicFluxCorrectionUtility::ComputeElementalAlgebraicFluxCorrections()
    {
        #pragma omp parallel for
        for (size_t i = 0; i < mrModelPart.NumberOfElements(); ++i)
        {
            auto it_corrections = mAlgebraicFluxCorrections.begin() + i;
            auto it_ho_values = mHighOrderValues.begin() + i;
            auto it_lo_values = mLowOrderValues.begin() + i;
            *it_corrections = *it_ho_values - *it_lo_values;
        }
    }

    void AlgebraicFluxCorrectionUtility::AssembleCorrections()
    {
        #pragma omp parallel for
        for (size_t i = 0; i < mrModelPart.NumberOfElements(); ++i)
        {
            Vector& corrections = *(mAlgebraicFluxCorrections.begin() + i);
            DofsVectorType dofs = *(mElementalDofs.begin() + i);
            for (size_t d = 0; d < dofs.size(); ++d)
            {
                #pragma omp critical
                (*dofs[d])(0) += corrections[d];
            }
        }
    }

    const Parameters AlgebraicFluxCorrectionUtility::GetDefaultParameters() const
    {
        const Parameters default_parameters = Parameters(R"(
        {
            "name"          : "algebraic_flux_correction_utility",
            "rebuild_level" : 0
        })");
        return default_parameters;
    }

    void AlgebraicFluxCorrectionUtility::GetElementalDofList()
    {
        const ProcessInfo& r_process_info = mrModelPart.GetProcessInfo();
        #pragma omp parallel for
        for (size_t i = 0; i < mrModelPart.NumberOfElements(); ++i)
        {
            (mrModelPart.ElementsBegin() + i)->GetDofList(*(mElementalDofs.begin() + i), r_process_info);
        }
    }
}
