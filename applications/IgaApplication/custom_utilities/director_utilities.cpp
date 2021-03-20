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
#include "iga_application_variables.h"
#include "custom_elements/shell_5p_element.h"

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
            auto& r_geometry = mrModelPart.GetRootModelPart().GetGeometry((IndexType)brep_ids[i]);

            size_t number_of_control_points = r_geometry.size();

            std::vector<IndexType> eqID(number_of_control_points);
            for (SizeType inodes = 0; inodes < number_of_control_points; ++inodes)
                eqID[inodes] = r_geometry[inodes].GetId();

            SparseMatrixType NTN(number_of_control_points, number_of_control_points, number_of_control_points*4); //inital guess how much non-zero are there
            // Can't we solve each patch independently? -> much more efficient
            Matrix directorAtIntgrationPoints{ ZeroMatrix(number_of_control_points, 3) };

            SparseSpaceType::SetToZero(NTN);
            Matrix NTNele;
            Matrix RhsEle;

            RhsEle = zero_matrix<double>(number_of_control_points,3);
            NTNele = zero_matrix<double>(number_of_control_points);

            PointerVector<Geometry<Node<3>>> quad_points;
            r_geometry.CreateQuadraturePointGeometries(quad_points, 3);
            for (SizeType iP = 0; iP < quad_points.size(); ++iP)
            {
                const Matrix& r_N = quad_points[iP].ShapeFunctionsValues();
                const array_1d<double, 3> A3 = quad_points[iP].UnitNormal(iP);

                NTNele += outer_prod(row(r_N,0), row(r_N, 0));

                for(SizeType inodes=0; inodes< number_of_control_points;++inodes)
                   row( RhsEle,inodes) += r_N(0, inodes)*trans(A3);
            }

            for (SizeType inodes = 0; inodes < number_of_control_points; ++inodes)
                row(directorAtIntgrationPoints, eqID[inodes]) += row(RhsEle, inodes);

            SparseSpaceType::AssembleLHS(NTN, NTNele, eqID);

            std::cout << "asdasdads" << std::endl;
        
            Parameters solver_parameters(mParameters["linear_solver_settings"]);
            if (!solver_parameters.Has("solver_type")) solver_parameters.AddString("solver_type", "skyline_lu_factorization");

            Matrix nodalDirectors(number_of_control_points, 3);
            auto solver = LinearSolverFactory<SparseSpaceType, LocalSpaceType>().Create(solver_parameters);
            std::cout << "asdasdads" << std::endl;
           // DenseVectorType nodalDirectorsVec(number_of_control_points);
           // DenseVectorType directorAtIntgrationPointsVec(number_of_control_points);
            solver->Solve(NTN, nodalDirectors, directorAtIntgrationPoints);
           // for (SizeType i = 0; i<3; ++i)
           // { 
            //   solver->Solve(NTN, nodalDirectorsVec, directorAtIntgrationPointsVec);
           // }
            std::cout << "asdasdads" << brep_ids.size()<< std::endl;

            std::cout << nodalDirectors << number_of_control_points<< std::endl;

            for (SizeType i = 0; i < number_of_control_points; ++i)
            {
                std::cout << row(nodalDirectors, i) << std::endl;
                r_geometry[i].SetValue(DIRECTOR, row(nodalDirectors, i));
                r_geometry[i].SetValue(DIRECTORLENGTH, norm_2(row(nodalDirectors, i)));
                KRATOS_WATCH(Shell5pElement::TangentSpaceFromStereographicProjection(row(nodalDirectors, i)))
                r_geometry[i].SetValue(DIRECTORTANGENTSPACE, Shell5pElement::TangentSpaceFromStereographicProjection(row(nodalDirectors, i)));
            }
            KRATOS_WATCH(nodalDirectors)
        }
    }
}  // namespace Kratos.
