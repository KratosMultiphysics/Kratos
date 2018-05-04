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
#include "surface_load_adjoint_condition_3d.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"
#include "includes/checks.h"

namespace Kratos
{
    //******************************* CONSTRUCTOR ****************************************
    //************************************************************************************

    SurfaceLoadAdjointCondition3D::SurfaceLoadAdjointCondition3D()
    {
    }

    //***********************************************************************************
    //***********************************************************************************

    SurfaceLoadAdjointCondition3D::SurfaceLoadAdjointCondition3D(
        IndexType NewId,
        GeometryType::Pointer pGeometry
        )
        : SurfaceLoadCondition3D(NewId, pGeometry)
    {
    }

    //***********************************************************************************
    //***********************************************************************************

    SurfaceLoadAdjointCondition3D::SurfaceLoadAdjointCondition3D(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
        )
        : SurfaceLoadCondition3D(NewId, pGeometry, pProperties)
    {
    }

    //********************************* CREATE *******************************************
    //************************************************************************************

    Condition::Pointer SurfaceLoadAdjointCondition3D::Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const
    {
        return Kratos::make_shared<SurfaceLoadAdjointCondition3D>(NewId, pGeom, pProperties);
    }

    //***********************************************************************************
    //***********************************************************************************

    Condition::Pointer SurfaceLoadAdjointCondition3D::Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const
    {
        return Kratos::make_shared<SurfaceLoadAdjointCondition3D>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

       //************************************************************************************
    //************************************************************************************

    void SurfaceLoadAdjointCondition3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
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
    void SurfaceLoadAdjointCondition3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
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

    void SurfaceLoadAdjointCondition3D::GetValuesVector(Vector& rValues, int Step)
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const unsigned int MatSize = number_of_nodes * dim;

        if (rValues.size() != MatSize)
        {
            rValues.resize(MatSize, false);
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

    void SurfaceLoadAdjointCondition3D::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                       ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int MatSize = number_of_nodes * dimension;

        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize );
        // return zero vector because surface load has no influence on adjoint problem
        KRATOS_CATCH( "" )
    }

    //***********************************************************************
    //***********************************************************************

     //************************************************************************************
    //************************************************************************************

    void SurfaceLoadAdjointCondition3D::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        KRATOS_ERROR << "Chosen design variable not availible for Surface Load Condition!" << std::endl;

        KRATOS_CATCH( "" )

    }

    //************************************************************************************
    //************************************************************************************

    void SurfaceLoadAdjointCondition3D::CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        /*if( this->Has(rDesignVariable) )
        {

        }
        else */if( rDesignVariable == SHAPE_SENSITIVITY )
        {

            // Define working variables
            Vector RHS_undist;
            Vector RHS_dist;
            Matrix dummy_LHS;
            ProcessInfo copy_process_info = rCurrentProcessInfo;

            // Get disturbance measure
            double delta = this->GetValue(DISTURBANCE_MEASURE);

            const int number_of_nodes = GetGeometry().PointsNumber();
            const int dimension = this->GetGeometry().WorkingSpaceDimension();

            // compute RHS before disturbing
            this->CalculateAll(dummy_LHS, RHS_undist, copy_process_info, false, true);

            rOutput.resize(dimension * number_of_nodes, RHS_undist.size());

            //TODO: look that this works also for parallel computing
            for(int j = 0; j < number_of_nodes; j++)
            {
                //begin: derive w.r.t. x-coordinate---------------------------------------------------

                // disturb the design variable
                this->GetGeometry()[j].X0() += delta;

                // compute RHS after disturbance
                this->CalculateAll(dummy_LHS, RHS_dist, copy_process_info, false, true);

                //compute derivative of RHS w.r.t. design variable with finite differences
                noalias(RHS_dist) -= RHS_undist;
                RHS_dist /= delta;
                for(unsigned int i = 0; i < RHS_dist.size(); i++)
                    rOutput( (0 + j*dimension), i) = RHS_dist[i];

                // Reset pertubed vector
                RHS_dist = Vector(0);

                // undisturb the design variable
                this->GetGeometry()[j].X0() -= delta;

                //end: derive w.r.t. x-coordinate-----------------------------------------------------

                //begin: derive w.r.t. y-coordinate---------------------------------------------------

                // disturb the design variable
                this->GetGeometry()[j].Y0() += delta;

                // compute RHS after disturbance
                this->CalculateAll(dummy_LHS, RHS_dist, copy_process_info, false, true);

                //compute derivative of RHS w.r.t. design variable with finite differences
                noalias(RHS_dist) -= RHS_undist;
                RHS_dist /= delta;
                for(unsigned int i = 0; i < RHS_dist.size(); i++)
                    rOutput((1 + j*dimension),i) = RHS_dist[i];

                // Reset pertubed vector
                RHS_dist = Vector(0);

                // undisturb the design variable
                this->GetGeometry()[j].Y0() -= delta;

                //end: derive w.r.t. y-coordinate-----------------------------------------------------

                //begin: derive w.r.t. z-coordinate---------------------------------------------------

                // disturb the design variable
                this->GetGeometry()[j].Z0() += delta;

                // compute RHS after disturbance
                this->CalculateAll(dummy_LHS, RHS_dist, copy_process_info, false, true);

                //compute derivative of RHS w.r.t. design variable with finite differences
                noalias(RHS_dist) -= RHS_undist;
                RHS_dist /= delta;
                for(unsigned int i = 0; i < RHS_dist.size(); i++)
                    rOutput((2 + j*dimension),i) = RHS_dist[i];

                // Reset pertubed vector
                RHS_dist = Vector(0);

                // undisturb the design variable
                this->GetGeometry()[j].Z0() -= delta;

                // Compute RHS again in order to ensure that changed member variables get back their origin values
                this->CalculateAll(dummy_LHS, RHS_dist, copy_process_info, false, true);

                //end: derive w.r.t. z-coordinate-----------------------------------------------------

            }// end loop over element nodes

        }
        else
            KRATOS_ERROR << "Chosen design variable not availible for Surface Load Condition!" << std::endl;

        KRATOS_CATCH( "" )
    }

    //***********************************************************************
    //***********************************************************************

    int SurfaceLoadAdjointCondition3D::Check( const ProcessInfo& rCurrentProcessInfo )
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


    //******************************* DESTRUCTOR *****************************************
    //************************************************************************************

    SurfaceLoadAdjointCondition3D::~SurfaceLoadAdjointCondition3D()
    {
    }

    //***********************************************************************************
    //***********************************************************************************

} // Namespace Kratos.
