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

#if !defined( RESIDUALS_UTILITY_H_INCLUDED )
#define  RESIDUALS_UTILITY_H_INCLUDED

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/strategies/solving_strategy.h"

/* Include pybind to convert residual matrices to numpy  */
#include <pybind11/pybind11.h>

namespace Kratos
{
    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme<SparseSpaceType, LocalSpaceType> BaseSchemeType;

    // This utility returns the converged residuals projected onto the ROM basis Phi.
    class GetResiduals
    {
        public:

        KRATOS_CLASS_POINTER_DEFINITION(GetResiduals);

        GetResiduals(
        ModelPart& rModelPart,
        BaseSchemeType::Pointer pScheme             // Need to define a scheme in a simple way
        ): mpModelPart(rModelPart), mpScheme(pScheme)
        {};

        ~GetResiduals()= default;


        //  ########################## WORKING FOR CONDIDIONAL RESIDUAL ###########################
        Matrix ConditionalResidual()
        {
            // Getting the array of the conditions
            const int nconditions = static_cast<int>(mpModelPart.Conditions().size());
            auto& CurrentProcessInfo = mpModelPart.GetProcessInfo();
            auto cond_begin = mpModelPart.ConditionsBegin();

            //contributions to the system
            Matrix LHS_Contribution = ZeroMatrix(0, 0);
            Vector RHS_Contribution = ZeroVector(0);

            int LargestNodes = 0;
            //finding the maximum number of nodes in conditions
            for (int k = 0; k < nconditions; k++)
            {
                auto it = cond_begin + k;
                int NumNodes = it->GetGeometry().size();
                if (NumNodes > LargestNodes)
                {
                    LargestNodes = NumNodes;
                }
            }
            KRATOS_WATCH(LargestNodes)

            Element::EquationIdVectorType EquationId;
            Matrix MatrixResiduals = ZeroMatrix( nconditions , (LargestNodes * 1));   // Matrix of reduced residuals.

            for (int k = 0; k < nconditions; k++)
            {
                auto it = cond_begin + k;
                int NumDofs = it->GetGeometry().size() * 1;
                mpScheme->Condition_CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                for (unsigned int i = 0; i < NumDofs; i++)
                    MatrixResiduals(k,i) = RHS_Contribution(i);
                // clean local elemental me overridemory
                mpScheme->CleanMemory(*(it.base()));
            }
        KRATOS_WATCH(MatrixResiduals)
        return MatrixResiduals;        
        }

        //  ########################## WORKING FOR ELEMENTAL RESIDUAL ###########################
        Matrix ElementalResidual()
        {
            // Getting the elements from the model
            const int nelements = static_cast<int>(mpModelPart.Elements().size());
            auto& CurrentProcessInfo = mpModelPart.GetProcessInfo();
            auto el_begin = mpModelPart.ElementsBegin();

            //contributions to the system
            Matrix LHS_Contribution = ZeroMatrix(0, 0);
            Vector RHS_Contribution = ZeroVector(0);

            int LargestNodes = 0;
            //finding the maximum number of nodes in elements
            for (int k = 0; k < nelements; k++)
            {
                auto it_el = el_begin + k;
                int NumNodes = it_el->GetGeometry().size();
                if (NumNodes > LargestNodes)
                {
                    LargestNodes = NumNodes;
                }
            }
            KRATOS_WATCH(LargestNodes)

            Element::EquationIdVectorType EquationId;
            Matrix MatrixResiduals = ZeroMatrix( nelements , (LargestNodes * 1 ));   // Matrix of reduced residuals.

            for (int k = 0; k < nelements; k++)
            {
                auto it_el = el_begin + k;
                int NumDofs = it_el->GetGeometry().size() * 1 ;
                mpScheme->CalculateSystemContributions(*(it_el.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                for (unsigned int i = 0; i < NumDofs; i++)
                    MatrixResiduals(k,i) = RHS_Contribution(i);
                // clean local elemental me overridemory
                mpScheme->CleanMemory(*(it_el.base()));
            }
        KRATOS_WATCH(MatrixResiduals)
        return MatrixResiduals;
        
        }


        void convert_to_numpy(const Matrix & input, pybind11::object obj)
        {
            PyObject* pobj = obj.ptr();
            Py_buffer pybuf;
            PyObject_GetBuffer(pobj, &pybuf, PyBUF_SIMPLE);
            void *buf = pybuf.buf;
            double *p = (double*)buf;
            Py_XDECREF(pobj);

			unsigned int n_rows = input.size1();
			unsigned int n_cols = input.size2();
            for (unsigned int i = 0; i < n_rows; i++)
            {
                for (unsigned int j = 0; j < n_cols; j++)
                {
                    p[i*n_cols+j] = input(i,j);
                }
            }
        }

        protected:

            std::vector< std::string > mNodalVariablesNames;
            int mNodalDofs;
            int mRomDofs;
            BaseSchemeType::Pointer mpScheme;
            ModelPart& mpModelPart;

        };



} // namespace Kratos



#endif // RESIDUALS_UTILITY_H_INCLUDED  defined