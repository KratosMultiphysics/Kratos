//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "director_utilities.h"
#include "factories/linear_solver_factory.h"


namespace Kratos
{
    DirectorUtilities::DirectorUtilities(
        ModelPart& rModelPart,
        Parameters JsonParameters)
        : mrModelPart(rModelPart)
        , mParameters(JsonParameters)
    {
    }

    void DirectorUtilities::ComputeDirectors()
    {
        Vector brep_ids = mParameters["brep_ids"].GetVector();
        for (IndexType i = 0; i < brep_ids.size(); ++i) {
            auto p_geom = mrModelPart.GetRootModelPart().pGetGeometry((IndexType)brep_ids[i]);
            KRATOS_WATCH(p_geom->Id())
        }

        auto elementList = mrModelPart.GetRootModelPart().pElements();

        size_t elenum = mrModelPart.NumberOfElements();
        size_t number_of_control_points = mrModelPart.NumberOfNodes();

        SparseMatrixType NTN(number_of_control_points, number_of_control_points, number_of_control_points*4); //inital guess how much non-zero are there
        DenseVectorType directorAtIntgrationPoints(number_of_control_points, 3);

        SparseSpaceType::SetToZero(NTN);
        SparseSpaceType::SetToZero(directorAtIntgrationPoints);
        Matrix Nele;
        Matrix NTNele;
        Matrix RhsEle;

        for (SizeType i = 0; i < elenum; ++i)
        {
            const auto& ele = mrModelPart.pGetElement(i);
            const auto& eleGeometry = ele->GetGeometry();

            const SizeType numNodes = eleGeometry.size();
            RhsEle = zero_matrix<double>(numNodes,3);
            NTNele = zero_matrix<double>(numNodes);

            Nele = eleGeometry.ShapeFunctionsValues();
            const SizeType r_number_of_integration_points = ele->GetGeometry().IntegrationPointsNumber();
            for (SizeType iP = 0; iP < r_number_of_integration_points; ++iP)
            {
                const Vector& Nip = row(Nele, iP);
                const array_1d<double, 3> A3 = ele->GetGeometry().UnitNormal(iP);

                NTNele += outer_prod(Nip, Nip);

                for(SizeType inodes=0; inodes< numNodes;++inodes)
                   row( RhsEle,inodes) += Nip[inodes]*trans(A3);

            }
            for (SizeType inodes = 0; inodes < numNodes; ++inodes)
            {
                row(directorAtIntgrationPoints, eleGeometry[inodes].GetId()) += row(RhsEle, inodes);

                //how to  effiecently use ublas_space AssembleLHS
                for (SizeType jnodes = 0; inodes < numNodes; ++jnodes)
                    NTN(eleGeometry[inodes].GetId(), eleGeometry[jnodes].GetId()) += NTNele(inodes, jnodes);

            }
        }
        
        Parameters solver_parameters(mParameters["linear_solver_settings"]);
        if (!solver_parameters.Has("solver_type")) solver_parameters.AddString("solver_type", "skyline_lu_factorization");

    DenseVectorType nodalDirectors(number_of_control_points, 3);
        auto solver = LinearSolverFactory<SparseSpaceType, LocalSpaceType>().Create(solver_parameters);

        solver->Solve(NTN, nodalDirectors, directorAtIntgrationPoints);

        for (SizeType i = 0; i < number_of_control_points; ++i)
            //mrModelPart.GetRootModelPart().pGetGeometry()[i].set //How to send solution back to nodes?
            mrModelPart.Nodes()[i].SetValue(DIRECTOR, row(nodalDirectors, i);


        KRATOS_WATCH(solution)

    }



}  // namespace Kratos.