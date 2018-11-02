//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

// Project includes

// Application includes
#include "custom_conditions/surface_load_from_DEM_condition_3d.h"


namespace Kratos
{
    //******************************* CONSTRUCTOR ****************************************
    //************************************************************************************

    SurfaceLoadFromDEMCondition3D::SurfaceLoadFromDEMCondition3D()
    {
    }

    //***********************************************************************************
    //***********************************************************************************

    SurfaceLoadFromDEMCondition3D::SurfaceLoadFromDEMCondition3D(
        IndexType NewId,
        GeometryType::Pointer pGeometry
        )
        : SurfaceLoadCondition3D(NewId, pGeometry)
    {
    }

    //***********************************************************************************
    //***********************************************************************************

    SurfaceLoadFromDEMCondition3D::SurfaceLoadFromDEMCondition3D(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
        )
        : SurfaceLoadCondition3D(NewId, pGeometry, pProperties)
    {
    }

    //********************************* CREATE *******************************************
    //************************************************************************************

    Condition::Pointer SurfaceLoadFromDEMCondition3D::Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const
    {
        return Kratos::make_shared<SurfaceLoadFromDEMCondition3D>(NewId, pGeom, pProperties);
    }

    //***********************************************************************************
    //***********************************************************************************

    Condition::Pointer SurfaceLoadFromDEMCondition3D::Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const
    {
        return Kratos::make_shared<SurfaceLoadFromDEMCondition3D>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    //******************************* DESTRUCTOR *****************************************
    //************************************************************************************

    SurfaceLoadFromDEMCondition3D::~SurfaceLoadFromDEMCondition3D()
    {
    }

    //***********************************************************************************
    //***********************************************************************************

    void SurfaceLoadFromDEMCondition3D::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
        )
    {
        KRATOS_TRY;

        const GeometryType& Geometry = this->GetGeometry();

        const unsigned int number_of_nodes = Geometry.size();
        const unsigned int mat_size = number_of_nodes * 3;

        //Resizing as needed the LHS
        if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
        {
            if (rLeftHandSideMatrix.size1() != mat_size)
            {
                rLeftHandSideMatrix.resize(mat_size, mat_size, false);
            }

            noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); //resetting LHS
        }

        // Resizing as needed the RHS
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            if (rRightHandSideVector.size() != mat_size)
            {
                rRightHandSideVector.resize(mat_size, false);
            }

            rRightHandSideVector = ZeroVector(mat_size); //resetting RHS
        }

        // Reading integration points and local gradients
        IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(GetGeometry());
        const GeometryType::IntegrationPointsArrayType& integration_points = Geometry.IntegrationPoints(integration_method);
        const Matrix& Ncontainer = Geometry.ShapeFunctionsValues(integration_method);

        // Calculating actual jacobian
        GeometryType::JacobiansType J;
        J = Geometry.Jacobian(J,integration_method);

        // Vector with a loading applied to the elemnt
        array_1d<double, 3 > surface_load = ZeroVector(3);
        if( this->Has( DEM_SURFACE_LOAD ) )
        {
            noalias(surface_load) = this->GetValue( DEM_SURFACE_LOAD );
        }

        for (unsigned int point_number = 0; point_number < integration_points.size(); point_number++)
        {
            const double det_j = MathUtils<double>::GeneralizedDet(J[point_number]);
            const double integration_weight = GetIntegrationWeight(integration_points, point_number, det_j);
            auto& N = row(Ncontainer, point_number);

            //generic load on gauss point
            array_1d<double, 3> gauss_load = surface_load;
            for (unsigned int ii = 0; ii < number_of_nodes; ++ii)
            {
                if( Geometry[ii].SolutionStepsDataHas( DEM_SURFACE_LOAD ) )
                {
                    noalias(gauss_load) += N[ii]*Geometry[ii].FastGetSolutionStepValue( DEM_SURFACE_LOAD );
                }
            }

            for (unsigned int ii = 0; ii < number_of_nodes; ++ii)
            {
                const unsigned int base = ii * 3;
                for(unsigned int k = 0; k < 3; ++k)
                {
                    rRightHandSideVector[base+k] += integration_weight * N[ii] * gauss_load[k];
                }
            }
        }

        KRATOS_CATCH("")
    }

} // Namespace Kratos.
