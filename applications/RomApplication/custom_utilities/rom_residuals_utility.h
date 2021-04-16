//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//  Kratos default license: kratos/license.txt
//
//  Main authors:    RAUL BRAVO
//

#if !defined( ROM_RESIDUALS_UTILITY_H_INCLUDED )
#define  ROM_RESIDUALS_UTILITY_H_INCLUDED

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "spaces/ublas_space.h"

/* Application includes */
#include "rom_application_variables.h"

namespace Kratos
{
    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme<SparseSpaceType, LocalSpaceType> BaseSchemeType;

    // This utility returns the converged residuals projected onto the ROM basis Phi.
    class RomResidualsUtility
    {
        public:

        KRATOS_CLASS_POINTER_DEFINITION(RomResidualsUtility);

        RomResidualsUtility(
        ModelPart& rModelPart,
        Parameters ThisParameters,
        BaseSchemeType::Pointer pScheme
        ): mpModelPart(rModelPart), mpScheme(pScheme){
        // Validate default parameters
        Parameters default_parameters = Parameters(R"(
        {
            "nodal_unknowns" : [],
            "number_of_rom_dofs" : 10
        })" );

        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        mNodalVariablesNames = ThisParameters["nodal_unknowns"].GetStringArray();

        mNodalDofs = mNodalVariablesNames.size();
        mRomDofs = ThisParameters["number_of_rom_dofs"].GetInt();

        // Setting up mapping: VARIABLE_KEY --> CORRECT_ROW_IN_BASIS
        for(int k=0; k<mNodalDofs; k++){
            if(KratosComponents<Variable<double>>::Has(mNodalVariablesNames[k]))
            {
                const auto& var = KratosComponents<Variable<double>>::Get(mNodalVariablesNames[k]);
                MapPhi[var.Key()] = k;
            }
            else
                KRATOS_ERROR << "variable \""<< mNodalVariablesNames[k] << "\" not valid" << std::endl;

        }
    }

        ~RomResidualsUtility()= default;


        void GetPhiElemental(
            Matrix &PhiElemental,
            const Element::DofsVectorType &dofs,
            const Element::GeometryType &geom)
        {
            const auto *pcurrent_rom_nodal_basis = &(geom[0].GetValue(ROM_BASIS));
            int counter = 0;
            for(unsigned int k = 0; k < dofs.size(); ++k){
                auto variable_key = dofs[k]->GetVariable().Key();
                if(k==0)
                    pcurrent_rom_nodal_basis = &(geom[counter].GetValue(ROM_BASIS));
                else if(dofs[k]->Id() != dofs[k-1]->Id()){
                    counter++;
                    pcurrent_rom_nodal_basis = &(geom[counter].GetValue(ROM_BASIS));
                }
                if (dofs[k]->IsFixed())
                    noalias(row(PhiElemental, k)) = ZeroVector(PhiElemental.size2());
                else
                    noalias(row(PhiElemental, k)) = row(*pcurrent_rom_nodal_basis, MapPhi[variable_key]);
            }
        }


        Matrix Calculate()
        {
            // Getting the number of elements and conditions from the model
            const int nelements = static_cast<int>(mpModelPart.Elements().size());
            const int nconditions = static_cast<int>(mpModelPart.Conditions().size());

            const auto& CurrentProcessInfo = mpModelPart.GetProcessInfo();
            const auto el_begin = mpModelPart.ElementsBegin();
            const auto cond_begin = mpModelPart.ConditionsBegin();

            //contributions to the system
            Matrix LHS_Contribution = ZeroMatrix(0, 0);
            Vector RHS_Contribution = ZeroVector(0);

            //vector containing the localization in the system of the different terms
            Element::EquationIdVectorType EquationId;
            Matrix MatrixResiduals( (nelements + nconditions), mRomDofs); // Matrix of reduced residuals.
            Matrix PhiElemental;
            #pragma omp parallel firstprivate(nelements, nconditions, LHS_Contribution, RHS_Contribution, EquationId, PhiElemental, el_begin, cond_begin)
            {
                #pragma omp for nowait
                for (int k = 0; k < nelements; k++){
                    auto it_el = el_begin + k;
                    //detect if the element is active or not. If the user did not make any choice the element is active by default
                    bool element_is_active = true;
                    if ((it_el)->IsDefined(ACTIVE))
                        element_is_active = (it_el)->Is(ACTIVE);
                    if (element_is_active){
                        //calculate elemental contribution
                        mpScheme->CalculateSystemContributions(*it_el, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                        Element::DofsVectorType dofs;
                        it_el->GetDofList(dofs, CurrentProcessInfo);
                        //assemble the elemental contribution - here is where the ROM acts
                        //compute the elemental reduction matrix PhiElemental
                        const auto& geom = it_el->GetGeometry();
                        if(PhiElemental.size1() != dofs.size() || PhiElemental.size2() != mRomDofs)
                            PhiElemental.resize(dofs.size(), mRomDofs,false);
                        GetPhiElemental(PhiElemental, dofs, geom);
                        noalias(row(MatrixResiduals, k)) = prod(trans(PhiElemental), RHS_Contribution); // The size of the residual will vary only when using more ROM modes, one row per condition
                    }

                }

                #pragma omp for nowait
                for (int k = 0; k < nconditions;  k++){
                    ModelPart::ConditionsContainerType::iterator it = cond_begin + k;
                    //detect if the condition is active or not. If the user did not make any choice the condition is active by default
                    bool condition_is_active = true;
                    if ((it)->IsDefined(ACTIVE))
                        condition_is_active = (it)->Is(ACTIVE);
                    if (condition_is_active){
                        Condition::DofsVectorType dofs;
                        it->GetDofList(dofs, CurrentProcessInfo);
                        //calculate elemental contribution
                        mpScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                        //assemble the elemental contribution - here is where the ROM acts
                        //compute the elemental reduction matrix PhiElemental
                        const auto& geom = it->GetGeometry();
                        if(PhiElemental.size1() != dofs.size() || PhiElemental.size2() != mRomDofs)
                            PhiElemental.resize(dofs.size(), mRomDofs,false);
                        GetPhiElemental(PhiElemental, dofs, geom);
                        noalias(row(MatrixResiduals, k+nelements)) = prod(trans(PhiElemental), RHS_Contribution); // The size of the residual will vary only when using more ROM modes, one row per condition
                    }
                }
            }
        return MatrixResiduals;
        }

        protected:
            std::vector< std::string > mNodalVariablesNames;
            int mNodalDofs;
            unsigned int mRomDofs;
            BaseSchemeType::Pointer mpScheme;
            ModelPart& mpModelPart;
            std::unordered_map<Kratos::VariableData::KeyType,int> MapPhi;
        };



} // namespace Kratos



#endif // ROM_RESIDUALS_UTILITY_H_INCLUDED  defined