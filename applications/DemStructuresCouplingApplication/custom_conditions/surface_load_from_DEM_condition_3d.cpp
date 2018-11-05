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

    GeometryData::IntegrationMethod SurfaceLoadFromDEMCondition3D::GetIntegrationMethod()
    {
        return GeometryData::GI_GAUSS_2;
        //return this->GetGeometry().GetDefaultIntegrationMethod();
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
        GeometryData::IntegrationMethod integration_method = GetIntegrationMethod();
        const GeometryType::IntegrationPointsArrayType& integration_points = Geometry.IntegrationPoints(integration_method);
        const Matrix& Ncontainer = Geometry.ShapeFunctionsValues(integration_method);

        // Calculating actual jacobian
        GeometryType::JacobiansType J;
        J = Geometry.Jacobian(J,integration_method);

        // Vector with a loading applied to the elemnt
        array_1d<double, 3 > surface_load;

        for (unsigned int point_number = 0; point_number < integration_points.size(); point_number++)
        {
            const double det_j = MathUtils<double>::GeneralizedDet(J[point_number]);
            const double integration_weight = GetIntegrationWeight(integration_points, point_number, det_j);

            //generic load on gauss point
            this->InterpolateSurfaceLoad(surface_load, Ncontainer, number_of_nodes, point_number);

            for (unsigned int i = 0; i < number_of_nodes; ++i)
            {
                const unsigned int base = i * 3;
                for(unsigned int k = 0; k < 3; ++k)
                {
                    rRightHandSideVector[base+k] += integration_weight * Ncontainer(point_number,i) * surface_load[k];
                }
            }
        }

        KRATOS_CATCH("")
    }

    //***********************************************************************************
    //***********************************************************************************

    void SurfaceLoadFromDEMCondition3D::InterpolateSurfaceLoad(array_1d<double,3>& r_surface_load,
                                                                const Matrix& n_container,
                                                                const unsigned int& number_of_nodes,
                                                                const unsigned int& g_point)
    {
        const GeometryType& Geometry = this->GetGeometry();

        //generic load on gauss point
        noalias(r_surface_load) = ZeroVector(3);

        for (unsigned int i = 0; i < number_of_nodes; ++i) {
            
            if (Geometry[i].SolutionStepsDataHas(DEM_SURFACE_LOAD)) {
                
                noalias(r_surface_load) += n_container(g_point,i) * Geometry[i].FastGetSolutionStepValue(DEM_SURFACE_LOAD);
            }
        }
    }
} // Namespace Kratos
