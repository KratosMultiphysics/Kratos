//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes
#include <iostream>
#include <fstream>

// External includes

// Project includes
#include "master_slave_process.h"

namespace Kratos
{

MasterSlaveProcess::MasterSlaveProcess(
    Model& rModel,
    Parameters ThisParameters)
    : Process()
    , mrModel(rModel)
    , mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());

    std::string constraints_model_part_name = mThisParameters["constraints_model_part_name"].GetString();
    mpConstraintsModelPart = &rModel.GetModelPart(constraints_model_part_name);

    std::string analysis_model_part_name = mThisParameters["analysis_model_part_name"].GetString();
    mpAnalysisModelPart = &rModel.GetModelPart(analysis_model_part_name);

    std::vector<std::string> constraints_list = mThisParameters["variables"].GetStringArray();
    mConstraintsList = constraints_list;

    bool master_slave_flag = mThisParameters["master_slave_flag"].GetBool();
    mMasterSlaveFlag = master_slave_flag;
}

void MasterSlaveProcess::ExecuteBeforeSolutionLoop()
{
    //0) Trick: create quadrature point geometries (TO DO)
    bool is_the_first_geometry_master = mMasterSlaveFlag; 

    auto& r_brep_curve_on_surface_master = mpConstraintsModelPart->GeometriesBegin()->GetGeometryPart(!is_the_first_geometry_master);
    auto& r_brep_curve_on_surface_slave = mpConstraintsModelPart->GeometriesBegin()->GetGeometryPart(is_the_first_geometry_master); 
    
    GeometriesArrayType slave_quadrature_points;
    IntegrationInfo integration_info = r_brep_curve_on_surface_slave.GetDefaultIntegrationInfo();
    r_brep_curve_on_surface_slave.CreateQuadraturePointGeometries(slave_quadrature_points, 2, integration_info);
      
    //1) Get node_id slave
    std::vector<IndexType> node_id;
    for (auto& geom : slave_quadrature_points) {
        auto& r_N = geom.ShapeFunctionsValues();

        for (IndexType i = 0; i<r_N.size2();++i)
        {
            if(r_N(0,i) > 1e-6)
            {
                mpConstraintsModelPart->AddNode(geom.pGetPoint(i));
            }
        }
    }
    
    for (auto& node : mpConstraintsModelPart->Nodes()) {
        node_id.push_back(node.Id());
    }

    //2) Construct and apply the constraints for each individual slave node
    const auto& p_nurbs_surface_master = mpConstraintsModelPart->GeometriesBegin()->GetGeometryPart(!is_the_first_geometry_master).pGetGeometryPart(std::numeric_limits<IndexType>::max());
    const auto& p_nurbs_surface_slave = mpConstraintsModelPart->GeometriesBegin()->GetGeometryPart(is_the_first_geometry_master).pGetGeometryPart(std::numeric_limits<IndexType>::max());

    std::unordered_map<int, std::unordered_map<int, double>> master_slave_constraints;

    for (auto& id : node_id) {
        CoordinatesArrayType local_coords_master = ZeroVector(3);
        CoordinatesArrayType local_coords_slave = ZeroVector(3); 
        CoordinatesArrayType global_coords_master = ZeroVector(3);
        CoordinatesArrayType global_coords_slave = ZeroVector(3); 

        //a) ProjectionPointGlobalToLocalSpace of the slave (TO DO: Check!)
        int success_slave = p_nurbs_surface_slave->ProjectionPointGlobalToLocalSpace(
            mpConstraintsModelPart->GetNode(id), local_coords_slave);

        p_nurbs_surface_slave->GlobalCoordinates(global_coords_slave, local_coords_slave);

        CurveTessellation<PointerVector<Node>> curve_tesselation;
        curve_tesselation.Tessellate(*(p_nurbs_surface_master.get()), 1e-2, p_nurbs_surface_master->PolynomialDegree(0));

        curve_tesselation.GetClosestPoint(
            global_coords_slave,
            global_coords_master,
            local_coords_master);

        p_nurbs_surface_master->ProjectionPointGlobalToLocalSpace(
            global_coords_slave,
            local_coords_master);

        p_nurbs_surface_master->GlobalCoordinates(global_coords_master, local_coords_master);

        //b) Evaluate the shape function and obtain the corresponding control points via master's quadrature points

        GeometriesArrayType master_quadrature_points;
        IntegrationPointsArrayType master_integration_points(1);
        master_integration_points[0] = IntegrationPoint<3>(local_coords_master, 1.0); //dummy weight

        IntegrationInfo integration_info_master = p_nurbs_surface_master->GetDefaultIntegrationInfo();
        p_nurbs_surface_master->CreateQuadraturePointGeometries(
                        master_quadrature_points, 2, master_integration_points, integration_info_master);

        //c) Data Structure (TO DO: make it more efficient)
        std::unordered_map<int, double> master_constraints;

        for (auto& geom : master_quadrature_points)
        {
            auto& r_N = geom.ShapeFunctionsValues();

            for (IndexType i = 0; i<r_N.size2();++i)
            {
                if(r_N(0,i) > 1e-6)
                {
                    master_constraints.emplace(geom.pGetPoint(i)->Id(), r_N(0,i));
                }
            }
        }
        master_slave_constraints.emplace(id, master_constraints);
    }

    //// Print the master slave constraints for debugging
    // for (const auto& [slave_id, constraints] : master_slave_constraints) {
    //     std::cout << "slave_id: " << slave_id << ":\n";
    //     for (const auto& [master_id, value] : constraints) {
    //         std::cout << "master_id: " << master_id << " => " << value << '\n';
    //     }
    // }

    //3) Version 1 : Linear master slave constraints
    Matrix relation_matrix;
    std::vector<IndexType> slave_ids;
    std::vector<IndexType> master_ids;
    GetRelationMatrix(master_slave_constraints, slave_ids, master_ids, relation_matrix);

    IndexType constraint_index = 55; //TO DO

    for (IndexType i = 0; i < mConstraintsList.size(); i++) 
    {
        LinearMasterSlaveConstraints(slave_ids, master_ids, relation_matrix, KratosComponents<Variable<double>>::Get(mConstraintsList[i]), constraint_index);
        constraint_index++;
    }

    //4) Version 2 : Linear multi freedom constraints (TO DO)
}

void MasterSlaveProcess::GetRelationMatrix(
    const std::unordered_map<int, std::unordered_map<int, double>>& MasterSlaveConstraints,
    std::vector<IndexType>& SlaveIds,
    std::vector<IndexType>& MasterIds,
    Matrix& RelationMatrix
)
{
    // a) collect and sort unique slave and master ids
    std::set<int> slave_set, master_set;
    for (const auto& [slave_id, constraints] : MasterSlaveConstraints) {
        slave_set.insert(slave_id); 
        for (const auto& [master_id, value] : constraints) {
            master_set.insert(master_id);
        }
    }
    std::vector<IndexType> slaves(slave_set.begin(), slave_set.end());
    std::vector<IndexType> masters(master_set.begin(), master_set.end());

    SlaveIds = slaves;
    MasterIds = masters;

    // b) index maps (master id -> column, slave id -> row)
    std::unordered_map<IndexType,IndexType> row, col;
    for (IndexType r = 0; r < slaves.size(); ++r)
    {
        row[slaves[r]] = r;
    } 
    for (IndexType c = 0; c < masters.size(); ++c) 
    {
        col[masters[c]] = c;   
    }

    // c) build the relation matrix
    RelationMatrix = ZeroMatrix(slaves.size(), masters.size());
    for (const auto& [slave_id, constraints] : MasterSlaveConstraints) {
        auto r = row[slave_id];
        for (const auto& [master_id, value] : constraints) {
            auto c = col[master_id];
            RelationMatrix(r,c) = value;
        }
    }
}

void MasterSlaveProcess::LinearMasterSlaveConstraints(
    const std::vector<IndexType>& SlaveIds,
    const std::vector<IndexType>& MasterIds,
    const Matrix& RelationMatrix,
    const Variable<double>& ConstraintVariable,
    const IndexType ConstraintIndex
)
{
    std::vector<Dof<double>*> dofs_slave;
    for(IndexType slave_id = 0; slave_id < SlaveIds.size(); slave_id++)
    {
        dofs_slave.push_back(mpAnalysisModelPart->GetNode(SlaveIds[slave_id]).pGetDof(ConstraintVariable));
    }

    std::vector<Dof<double>*> dofs_master;
    for(IndexType master_id = 0; master_id < MasterIds.size(); master_id++)
    {
        dofs_master.push_back(mpAnalysisModelPart->GetNode(MasterIds[master_id]).pGetDof(ConstraintVariable));
    }

    Vector constraint_vector = ZeroVector(RelationMatrix.size2());

    mpAnalysisModelPart->AddMasterSlaveConstraint(LinearMasterSlaveConstraint::Pointer(
        new LinearMasterSlaveConstraint(
            ConstraintIndex,
            dofs_master,
            dofs_slave,
            RelationMatrix,
            constraint_vector
        )
    ));
}


const Parameters MasterSlaveProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "constraints_model_part_name"          : "",
        "analysis_model_part_name"             : "",
        "variables"                            : [],
        "master_slave_flag"                    : true
    })" );
    return default_parameters;
}

} // namespace Kratos
