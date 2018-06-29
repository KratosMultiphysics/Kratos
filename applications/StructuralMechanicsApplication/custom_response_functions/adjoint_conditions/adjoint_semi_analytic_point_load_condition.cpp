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
#include "adjoint_semi_analytic_point_load_condition.h"
#include "structural_mechanics_application_variables.h"
#include "includes/checks.h"

namespace Kratos
{
    //******************************* CONSTRUCTOR ****************************************
    //************************************************************************************

    AdjointSemiAnalyticPointLoadCondition::AdjointSemiAnalyticPointLoadCondition(Condition::Pointer pPrimalCondition )
        : AdjointSemiAnalyticBaseCondition(pPrimalCondition )
    {
    }

    //********************************* CREATE *******************************************
    //************************************************************************************

    Condition::Pointer AdjointSemiAnalyticPointLoadCondition::Create(Condition::Pointer pPrimalCondition ) const
    {
        return Kratos::make_shared<AdjointSemiAnalyticPointLoadCondition>( pPrimalCondition );
    }

    //******************************* DESTRUCTOR *****************************************
    //************************************************************************************

    AdjointSemiAnalyticPointLoadCondition::~AdjointSemiAnalyticPointLoadCondition()
    {
    }

    //************************************************************************************
    //************************************************************************************

    void AdjointSemiAnalyticPointLoadCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        if (rResult.size() != dim * number_of_nodes)
        {
            rResult.resize(dim*number_of_nodes,false);
        }

        const unsigned int pos = this->GetGeometry()[0].GetDofPosition(ADJOINT_DISPLACEMENT_X);

        if(dim == 2)
        {
            for (unsigned int i = 0; i < number_of_nodes; ++i)
            {
                const unsigned int index = i * 2;
                rResult[index    ] = GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_X,pos    ).EquationId();
                rResult[index + 1] = GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Y,pos + 1).EquationId();
            }
        }
        else
        {
            for (unsigned int i = 0; i < number_of_nodes; ++i)
            {
                const unsigned int index = i * 3;
                rResult[index    ] = GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_X,pos    ).EquationId();
                rResult[index + 1] = GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Y,pos + 1).EquationId();
                rResult[index + 2] = GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Z,pos + 2).EquationId();
            }
        }
        KRATOS_CATCH("")
    }

    //***********************************************************************
    //***********************************************************************

    void AdjointSemiAnalyticPointLoadCondition::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim =  GetGeometry().WorkingSpaceDimension();
        ElementalDofList.resize(0);
        ElementalDofList.reserve(dim * number_of_nodes);

        if(dim == 2)
        {
            for (unsigned int i = 0; i < number_of_nodes; ++i)
            {
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_X));
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Y));
            }
        }
        else
        {
            for (unsigned int i = 0; i < number_of_nodes; ++i)
            {
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_X));
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Y));
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Z));
            }
        }

        KRATOS_CATCH("")
    }

    //***********************************************************************
    //***********************************************************************

    void AdjointSemiAnalyticPointLoadCondition::GetValuesVector(Vector& rValues, int Step)
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const unsigned int mat_size = number_of_nodes * dim;

        if (rValues.size() != mat_size)
        {
            rValues.resize(mat_size, false);
        }

        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            const array_1d<double, 3 > & Displacement = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_DISPLACEMENT, Step);
            unsigned int index = i * dim;
            for(unsigned int k = 0; k < dim; ++k)
            {
                rValues[index + k] = Displacement[k];
            }
        }
    }

    //************************************************************************************
    //************************************************************************************

    void AdjointSemiAnalyticPointLoadCondition::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        KRATOS_ERROR << "There is no scalar design variable available!" << std::endl;

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    void AdjointSemiAnalyticPointLoadCondition::CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int mat_size = number_of_nodes * dimension;

        rOutput = ZeroMatrix(mat_size,mat_size);

        if( rDesignVariable == POINT_LOAD )
        {
            const auto backup = mpPrimalCondition->GetValue(POINT_LOAD);
            mpPrimalCondition->SetValue(POINT_LOAD, this->GetValue(POINT_LOAD));

            Vector RHS;
            ProcessInfo copy_process_info = rCurrentProcessInfo;

            mpPrimalCondition->CalculateRightHandSide(RHS, copy_process_info);

            // Check if intensities of the point load are still available after the replacement process
            const double numerical_limit = std::numeric_limits<double>::epsilon();
            const double load_norm = MathUtils<double>::Norm(RHS);
            KRATOS_ERROR_IF(load_norm < numerical_limit) << "The norm of the point load is zero. Something went wrong!" << std::endl;

            int k = 0;
            for(unsigned int i = 0; i < RHS.size(); ++i)
            {
                rOutput(k,i)=(RHS[i]>0) ? 1.0 : -1.0;
                k++;
            }

            mpPrimalCondition->SetValue(POINT_LOAD, backup);
        }
        else
            rOutput.clear();

        KRATOS_CATCH( "" )
    }

    //***********************************************************************
    //***********************************************************************

    int AdjointSemiAnalyticPointLoadCondition::Check( const ProcessInfo& rCurrentProcessInfo )
    {
        // verify that the variables are correctly initialized
        KRATOS_CHECK_VARIABLE_KEY(ADJOINT_DISPLACEMENT);
        KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);

        // Check dofs
        GeometryType& r_geom = GetGeometry();
        for (unsigned int i = 0; i < r_geom.size(); i++)
        {
            auto& r_node = r_geom[i];
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_DISPLACEMENT, r_node);

            KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_X, r_node);
            KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Y, r_node);
            KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Z, r_node);
        }

        return 0;
    }

} // Namespace Kratos


