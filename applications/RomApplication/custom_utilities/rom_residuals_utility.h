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

#if !defined( GET_ROM_RESIDUALS_H_INCLUDED )
#define  GET_ROM_RESIDUALS_H_INCLUDED

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/strategies/solving_strategy.h"

/* Application includes */
#include "rom_application_variables.h"

namespace Kratos
{    
    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme<SparseSpaceType, LocalSpaceType> BaseSchemeType;

    // This utility returns the converged residuals projected onto the ROM basis Phi.
    class GetRomResiduals
    {
        public:
        
        KRATOS_CLASS_POINTER_DEFINITION(GetRomResiduals);
        
        GetRomResiduals(
        ModelPart& rModelPart,
        Parameters ThisParameters,
        BaseSchemeType::Pointer pScheme             // Need to define a scheme in a simple way     
        ): mpModelPart(rModelPart), mpScheme(pScheme)
    {
        // Validate default parameters
        Parameters default_parameters = Parameters(R"(
        {
            "nodal_unknowns" : [],
            "number_of_rom_dofs" : 10
        })" );

        ThisParameters.ValidateAndAssignDefaults(default_parameters);
        
        mNodalVariablesNames = ThisParameters["nodal_unknowns"].GetStringArray();        
        //Need to read the type of the variable and optain its size, incorrectly done here
        mNodalDofs = mNodalVariablesNames.size();
        mRomDofs = ThisParameters["number_of_rom_dofs"].GetInt();
        //this->mpModelPart = rModelPart;
        //this->mpScheme = pScheme;
    }

        ~GetRomResiduals()= default;


        Matrix Calculate() 
        {
            // Getting the elements from the model
            const int nelements = static_cast<int>(mpModelPart.Elements().size());

            // Getting the array of the conditions
            const int nconditions = static_cast<int>(mpModelPart.Conditions().size());
            auto& CurrentProcessInfo = mpModelPart.GetProcessInfo();
            auto el_begin = mpModelPart.ElementsBegin();
            auto cond_begin = mpModelPart.ConditionsBegin();

            //contributions to the system
            Matrix LHS_Contribution = ZeroMatrix(0, 0);
            Vector RHS_Contribution = ZeroVector(0);

            //vector containing the localization in the system of the different
            //terms
            Element::EquationIdVectorType EquationId;
            Matrix MatrixResiduals( (nelements + nconditions), mRomDofs);   // Matrix of reduced residuals.
            
            for (int k = 0; k < nelements; k++)
            {
                auto it_el = el_begin + k;
                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool element_is_active = true;
                if ((it_el)->IsDefined(ACTIVE))
                    element_is_active = (it_el)->Is(ACTIVE);              
                
                if (element_is_active)
                {   //calculate elemental contribution
                    mpScheme->CalculateSystemContributions(*(it_el.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);                    
                    Element::DofsVectorType dofs;
                    it_el->GetDofList(dofs, CurrentProcessInfo);
                    //assemble the elemental contribution - here is where the ROM acts
                    //compute the elemental reduction matrix PhiElemental
                    const auto& geom = it_el->GetGeometry();
                    Matrix PhiElemental(geom.size()*mNodalDofs, mRomDofs);
                    Vector ResidualReduced(mRomDofs); // The size of the residual will vary only when using more ROM modes, one row per element

                    for(unsigned int i=0; i<geom.size(); ++i)
                    {
                        const Matrix& rom_nodal_basis = geom[i].GetValue(ROM_BASIS);
                        for(unsigned int k=0; k<rom_nodal_basis.size1(); ++k)
                        {
                            if (dofs[i*mNodalDofs + k]->IsFixed())
                                row(PhiElemental, i*mNodalDofs + k) = ZeroVector(PhiElemental.size2());
                            else
                                row(PhiElemental, i*mNodalDofs+k) = row(rom_nodal_basis,k);
                        }
                    }
                    ResidualReduced = prod(trans(PhiElemental), RHS_Contribution);
                    row(MatrixResiduals, k) = ResidualReduced;
                    // clean local elemental me overridemory
                    mpScheme->CleanMemory(*(it_el.base()));                
                }
                
            }

            // #pragma omp for  schedule(guided , 512)
            for (int k = 0; k < nconditions;  k++)
            {
                ModelPart::ConditionsContainerType::iterator it = cond_begin + k;
                //detect if the condition is active or not. If the user did not make any choice the condition
                //is active by default
                bool condition_is_active = true;
                if ((it)->IsDefined(ACTIVE))
                    condition_is_active = (it)->Is(ACTIVE);

                if (condition_is_active)
                {
                    Condition::DofsVectorType dofs;
                    it->GetDofList(dofs, CurrentProcessInfo);
                    //calculate elemental contribution
                    mpScheme->Condition_CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                    //assemble the elemental contribution - here is where the ROM acts
                    //compute the elemental reduction matrix PhiElemental
                    const auto& r_geom = it->GetGeometry();
                    Matrix PhiElemental(r_geom.size()*mNodalDofs, mRomDofs);
                    Vector ResidualReduced(mRomDofs); // The size of the residual will vary only when using more ROM modes, one row per condition

                    for(unsigned int i=0; i<r_geom.size(); ++i)
                    {
                        const Matrix& rom_nodal_basis = r_geom[i].GetValue(ROM_BASIS);
                        for(unsigned int k=0; k<rom_nodal_basis.size1(); ++k)
                        {
                            if (dofs[i*mNodalDofs + k]->IsFixed())
                                row(PhiElemental, i*mNodalDofs + k) = ZeroVector(PhiElemental.size2());
                            else
                                row(PhiElemental, i*mNodalDofs+k) = row(rom_nodal_basis,k);
                        }
                    }                                       
                    ResidualReduced = prod(trans(PhiElemental), RHS_Contribution);
                    row(MatrixResiduals, k+nelements) = ResidualReduced;

                    // clean local elemental memory
                    mpScheme->CleanMemory(*(it.base()));
                }

            }
        return MatrixResiduals;        
        }


        protected:
        
            std::vector< std::string > mNodalVariablesNames;
            int mNodalDofs;
            int mRomDofs;
            BaseSchemeType::Pointer mpScheme;
            ModelPart& mpModelPart;

        };



} // namespace Kratos



#endif // GET_ROM_RESIDUALS_H_INCLUDED  defined