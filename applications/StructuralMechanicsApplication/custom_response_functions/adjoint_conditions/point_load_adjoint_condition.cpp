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
#include "includes/define.h"
#include "point_load_adjoint_condition.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"
#include "includes/checks.h"

namespace Kratos
{
    //******************************* CONSTRUCTOR ****************************************
    //************************************************************************************

    PointLoadAdjointCondition::PointLoadAdjointCondition( IndexType NewId, GeometryType::Pointer pGeometry )
        : PointLoadCondition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    PointLoadAdjointCondition::PointLoadAdjointCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : PointLoadCondition( NewId, pGeometry, pProperties )
    {
    }

    //********************************* CREATE *******************************************
    //************************************************************************************

    Condition::Pointer PointLoadAdjointCondition::Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const
    {
        return Kratos::make_shared<PointLoadAdjointCondition>(NewId, pGeom, pProperties);
    }

    //************************************************************************************
    //************************************************************************************

    Condition::Pointer PointLoadAdjointCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
    {
        return Kratos::make_shared<PointLoadAdjointCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
    }

    //******************************* DESTRUCTOR *****************************************
    //************************************************************************************

    PointLoadAdjointCondition::~PointLoadAdjointCondition()
    {
    }

      //************************************************************************************
    //************************************************************************************

    void PointLoadAdjointCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
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
    void PointLoadAdjointCondition::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
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

    void PointLoadAdjointCondition::GetValuesVector(Vector& rValues, int Step)
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

    void PointLoadAdjointCondition::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                       ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int mat_size = number_of_nodes * dimension;

        if ( rLeftHandSideMatrix.size1() != mat_size )
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size );

        KRATOS_CATCH( "" )
    }

    //***********************************************************************
    //***********************************************************************

     //************************************************************************************
    //************************************************************************************

    void PointLoadAdjointCondition::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        KRATOS_ERROR << "There is no scalar design variable availible!" << std::endl; 

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    void PointLoadAdjointCondition::CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
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
            Vector RHS;
            Matrix dummy_LHS;
            ProcessInfo copy_process_info = rCurrentProcessInfo;

            PointLoadCondition::CalculateAll(dummy_LHS, RHS, copy_process_info, false, true);

            int k = 0;
            for(unsigned int i = 0; i < RHS.size(); ++i)
            {
                if( RHS[i] > 1e-12 || RHS[i] < -1e-12 )
                    rOutput(k, i) = 1.0; // or is it sign dependent?
                k++;
            }
        }
        else if( rDesignVariable == SHAPE_SENSITIVITY )
        {
            rOutput.clear();
        }
        else
            rOutput.clear();

        KRATOS_CATCH( "" )
    }

    //***********************************************************************
    //***********************************************************************

    int PointLoadAdjointCondition::Check( const ProcessInfo& rCurrentProcessInfo )
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


