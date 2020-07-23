// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//

// System includes


// External includes


// Project includes
#include "includes/checks.h"

#include "adjoint_semi_analytic_point_load_condition.h"
#include "structural_mechanics_application_variables.h"
#include "custom_conditions/point_load_condition.h"

namespace Kratos
{
    template <class TPrimalCondition>
    void AdjointSemiAnalyticPointLoadCondition<TPrimalCondition>::CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const SizeType number_of_nodes = this->GetGeometry().size();
        const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();
        const SizeType mat_size = number_of_nodes * dimension;

        if( rDesignVariable == POINT_LOAD )
        {
            if ((rOutput.size1() != mat_size) || (rOutput.size2() != mat_size))
                rOutput.resize(mat_size, mat_size, false);

            noalias(rOutput) = ZeroMatrix(mat_size,mat_size);
            for(IndexType i = 0; i < mat_size; ++i)
                rOutput(i,i) = 1.0;
        }
        else if( rDesignVariable == SHAPE_SENSITIVITY )
        {
            rOutput = ZeroMatrix(mat_size, mat_size);
        }
        else
        {
            rOutput = ZeroMatrix(0, mat_size);
        }

        KRATOS_CATCH( "" )
    }

    // TODO find out what to do with KRATOS_API
    template class AdjointSemiAnalyticPointLoadCondition<PointLoadCondition>;

} // Namespace Kratos


