//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//

// System includes

// External includes

// Project includes
#include "custom_utilities/mpm_particle_generator_utility.h"
#include "custom_utilities/particle_mechanics_math_utilities.h"


namespace Kratos
{
namespace MPMParticleGeneratorUtility
{

/***********************************************************************************/
/***********************************************************************************/
    

    Matrix MP16ShapeFunctions()
    {
        const double Na1 = 0.33333333333333;
        const double Nb1 = 0.45929258829272;
        const double Nb2 = 0.08141482341455;
        const double Nc1 = 0.17056930775176;
        const double Nc2 = 0.65886138449648;

        const double Nd1 = 0.05054722831703;
        const double Nd2 = 0.89890554336594;

        const double Ne1 = 0.26311282963464;
        const double Ne2 = 0.72849239295540;
        const double Ne3 = 0.00839477740996;

        BoundedMatrix<double,16,3> MP_shape_functions;

        MP_shape_functions(0,0) = Na1;
        MP_shape_functions(0,1) = Na1;
        MP_shape_functions(0,2) = Na1;

        MP_shape_functions(1,0) = Nb1;
        MP_shape_functions(1,1) = Nb1;
        MP_shape_functions(1,2) = Nb2;

        MP_shape_functions(2,0) = Nb1;
        MP_shape_functions(2,1) = Nb2;
        MP_shape_functions(2,2) = Nb1;

        MP_shape_functions(3,0) = Nb2;
        MP_shape_functions(3,1) = Nb1;
        MP_shape_functions(3,2) = Nb1;

        MP_shape_functions(4,0) = Nc1;
        MP_shape_functions(4,1) = Nc1;
        MP_shape_functions(4,2) = Nc2;

        MP_shape_functions(5,0) = Nc1;
        MP_shape_functions(5,1) = Nc2;
        MP_shape_functions(5,2) = Nc1;

        MP_shape_functions(6,0) = Nc2;
        MP_shape_functions(6,1) = Nc1;
        MP_shape_functions(6,2) = Nc1;

        MP_shape_functions(7,0) = Nd1;
        MP_shape_functions(7,1) = Nd1;
        MP_shape_functions(7,2) = Nd2;

        MP_shape_functions(8,0) = Nd1;
        MP_shape_functions(8,1) = Nd2;
        MP_shape_functions(8,2) = Nd1;

        MP_shape_functions(9,0) = Nd2;
        MP_shape_functions(9,1) = Nd1;
        MP_shape_functions(9,2) = Nd1;

        MP_shape_functions(10,0) = Ne1;
        MP_shape_functions(10,1) = Ne2;
        MP_shape_functions(10,2) = Ne3;

        MP_shape_functions(11,0) = Ne2;
        MP_shape_functions(11,1) = Ne3;
        MP_shape_functions(11,2) = Ne1;

        MP_shape_functions(12,0) = Ne3;
        MP_shape_functions(12,1) = Ne1;
        MP_shape_functions(12,2) = Ne2;

        MP_shape_functions(13,0) = Ne2;
        MP_shape_functions(13,1) = Ne1;
        MP_shape_functions(13,2) = Ne3;

        MP_shape_functions(14,0) = Ne1;
        MP_shape_functions(14,1) = Ne3;
        MP_shape_functions(14,2) = Ne2;

        MP_shape_functions(15,0) = Ne3;
        MP_shape_functions(15,1) = Ne2;
        MP_shape_functions(15,2) = Ne1;

        //MP_shape_functions = [(Na1, Na1, Na1),(Nb1, Nb1, Nb2),(Nb1, Nb2, Nb1),(Nb2, Nb1, Nb1),
        //                    (Nc1, Nc1, Nc2),(Nc1, Nc2, Nc1),(Nc2, Nc1, Nc1),(Nd1, Nd1, Nd2),
        //                    (Nd1, Nd2, Nd1),(Nd2, Nd1, Nd1),(Ne1, Ne2, Ne3),(Ne2, Ne3, Ne1),
        //                    (Ne3, Ne1, Ne2),(Ne2, Ne1, Ne3),(Ne1, Ne3, Ne2),(Ne3, Ne2, Ne1)];

        return MP_shape_functions;
    }

/***********************************************************************************/
/***********************************************************************************/

    Matrix MP33ShapeFunctions()
    {
        const double Na2 = 0.02356522045239;
        const double Na1 = 0.488217389773805;

        const double Nb2 = 0.120551215411079;
        const double Nb1 = 0.43972439229446;

        const double Nc2 = 0.457579229975768;
        const double Nc1 = 0.271210385012116;

        const double Nd2 = 0.744847708916828;
        const double Nd1 = 0.127576145541586;

        const double Ne2 = 0.957365299093579;
        const double Ne1 = 0.021317350453210;

        const double Nf1 = 0.115343494534698;
        const double Nf2 = 0.275713269685514;
        const double Nf3 = 0.608943235779788;

        const double Ng1 = 0.022838332222257;
        const double Ng2 = 0.281325580989940;
        const double Ng3 = 0.695836086787803;

        const double Nh1 = 0.025734050548330;
        const double Nh2 = 0.116251915907597;
        const double Nh3 = 0.858014033544073;

        BoundedMatrix<double,33,3> MP_shape_functions;

        MP_shape_functions(0,0) = Na1;
        MP_shape_functions(0,1) = Na1;
        MP_shape_functions(0,2) = Na2;

        MP_shape_functions(1,0) = Na1;
        MP_shape_functions(1,1) = Na2;
        MP_shape_functions(1,2) = Na1;

        MP_shape_functions(2,0) = Na2;
        MP_shape_functions(2,1) = Na1;
        MP_shape_functions(2,2) = Na1;


        MP_shape_functions(3,0) = Nb1;
        MP_shape_functions(3,1) = Nb1;
        MP_shape_functions(3,2) = Nb2;

        MP_shape_functions(4,0) = Nb1;
        MP_shape_functions(4,1) = Nb2;
        MP_shape_functions(4,2) = Nb1;

        MP_shape_functions(5,0) = Nb2;
        MP_shape_functions(5,1) = Nb1;
        MP_shape_functions(5,2) = Nb1;

        MP_shape_functions(6,0) = Nc1;
        MP_shape_functions(6,1) = Nc1;
        MP_shape_functions(6,2) = Nc2;

        MP_shape_functions(7,0) = Nc1;
        MP_shape_functions(7,1) = Nc2;
        MP_shape_functions(7,2) = Nc1;

        MP_shape_functions(8,0) = Nc2;
        MP_shape_functions(8,1) = Nc1;
        MP_shape_functions(8,2) = Nc1;

        MP_shape_functions(9,0) = Nd1;
        MP_shape_functions(9,1) = Nd1;
        MP_shape_functions(9,2) = Nd2;

        MP_shape_functions(10,0) = Nd1;
        MP_shape_functions(10,1) = Nd2;
        MP_shape_functions(10,2) = Nd1;

        MP_shape_functions(11,0) = Nd2;
        MP_shape_functions(11,1) = Nd1;
        MP_shape_functions(11,2) = Nd1;

        MP_shape_functions(12,0) = Ne1;
        MP_shape_functions(12,1) = Ne1;
        MP_shape_functions(12,2) = Ne2;

        MP_shape_functions(13,0) = Ne1;
        MP_shape_functions(13,1) = Ne2;
        MP_shape_functions(13,2) = Ne1;

        MP_shape_functions(14,0) = Ne2;
        MP_shape_functions(14,1) = Ne1;
        MP_shape_functions(14,2) = Ne1;

        MP_shape_functions(15,0) = Nf1;
        MP_shape_functions(15,1) = Nf2;
        MP_shape_functions(15,2) = Nf3;

        MP_shape_functions(16,0) = Nf2;
        MP_shape_functions(16,1) = Nf3;
        MP_shape_functions(16,2) = Nf1;

        MP_shape_functions(17,0) = Nf3;
        MP_shape_functions(17,1) = Nf1;
        MP_shape_functions(17,2) = Nf2;

        MP_shape_functions(18,0) = Nf2;
        MP_shape_functions(18,1) = Nf1;
        MP_shape_functions(18,2) = Nf3;

        MP_shape_functions(19,0) = Nf1;
        MP_shape_functions(19,1) = Nf3;
        MP_shape_functions(19,2) = Nf2;

        MP_shape_functions(20,0) = Nf3;
        MP_shape_functions(20,1) = Nf2;
        MP_shape_functions(20,2) = Nf1;

        MP_shape_functions(21,0) = Ng1;
        MP_shape_functions(21,1) = Ng2;
        MP_shape_functions(21,2) = Ng3;

        MP_shape_functions(22,0) = Ng2;
        MP_shape_functions(22,1) = Ng3;
        MP_shape_functions(22,2) = Ng1;

        MP_shape_functions(23,0) = Ng3;
        MP_shape_functions(23,1) = Ng1;
        MP_shape_functions(23,2) = Ng2;

        MP_shape_functions(24,0) = Ng2;
        MP_shape_functions(24,1) = Ng1;
        MP_shape_functions(24,2) = Ng3;

        MP_shape_functions(25,0) = Ng1;
        MP_shape_functions(25,1) = Ng3;
        MP_shape_functions(25,2) = Ng2;

        MP_shape_functions(26,0) = Ng3;
        MP_shape_functions(26,1) = Ng2;
        MP_shape_functions(26,2) = Ng1;

        MP_shape_functions(27,0) = Nh1;
        MP_shape_functions(27,1) = Nh2;
        MP_shape_functions(27,2) = Nh3;

        MP_shape_functions(28,0) = Nh2;
        MP_shape_functions(28,1) = Nh3;
        MP_shape_functions(28,2) = Nh1;

        MP_shape_functions(29,0) = Nh3;
        MP_shape_functions(29,1) = Nh1;
        MP_shape_functions(29,2) = Nh2;

        MP_shape_functions(30,0) = Nh2;
        MP_shape_functions(30,1) = Nh1;
        MP_shape_functions(30,2) = Nh3;

        MP_shape_functions(31,0) = Nh1;
        MP_shape_functions(31,1) = Nh3;
        MP_shape_functions(31,2) = Nh2;

        MP_shape_functions(32,0) = Nh3;
        MP_shape_functions(32,1) = Nh2;
        MP_shape_functions(32,2) = Nh1;

        return MP_shape_functions;
    }

    void GetIntegrationPointVolumes(const GeometryType& rGeom, const IntegrationMethod IntegrationMethod, Vector& rIntVolumes)
    {
        auto int_points = rGeom.IntegrationPoints(IntegrationMethod);
        if (rIntVolumes.size() != int_points.size()) rIntVolumes.resize(int_points.size(),false);
        DenseVector<Matrix> jac_vec(int_points.size());
        rGeom.Jacobian(jac_vec, IntegrationMethod);
        for (size_t i = 0; i < int_points.size(); ++i) {
            rIntVolumes[i] = MathUtils<double>::DetMat(jac_vec[i]) * int_points[i].Weight();
        }
    }

    void GetIntegrationPointArea(const GeometryType& rGeom, const IntegrationMethod IntegrationMethod, Vector& rIntVolumes)
    {
        const double area = rGeom.Area();
        auto int_points = rGeom.IntegrationPoints(IntegrationMethod);
        if (rIntVolumes.size() != int_points.size()) rIntVolumes.resize(int_points.size(),false);
        for (size_t i = 0; i < int_points.size(); ++i) {
            rIntVolumes[i] = area * 0.5 * int_points[i].Weight();
        }
    }

    void DetermineIntegrationMethodAndShapeFunctionValues(const GeometryType& rGeom, const SizeType ParticlesPerElement,
        IntegrationMethod& rIntegrationMethod, Matrix& rN, bool& IsEqualVolumes)
    {
        const GeometryData::KratosGeometryType geo_type = rGeom.GetGeometryType();
        const SizeType domain_size = rGeom.WorkingSpaceDimension();

        if (geo_type == GeometryData::Kratos_Tetrahedra3D4 || geo_type == GeometryData::Kratos_Triangle2D3)
        {
            switch (ParticlesPerElement)
            {
            case 1:
                rIntegrationMethod = GeometryData::GI_GAUSS_1;
                break;
            case 3:
                rIntegrationMethod = GeometryData::GI_GAUSS_2;
                break;
            case 6:
                rIntegrationMethod = GeometryData::GI_GAUSS_4;
                break;
            case 12:
                rIntegrationMethod = GeometryData::GI_GAUSS_5;
                break;
            case 16:
                if (domain_size == 2) {
                    IsEqualVolumes = true;
                    KRATOS_INFO("MPMParticleGeneratorUtility") << "WARNING: "
                        << "16 particles per triangle element is only valid for undistorted triangles." << std::endl;
                    rN = MP16ShapeFunctions();
                    break;
                }
            case 33:
                if (domain_size == 2) {
                    IsEqualVolumes = true;
                    KRATOS_INFO("MPMParticleGeneratorUtility") << "WARNING: "
                        << "33 particles per triangle element is only valid for undistorted triangles." << std::endl;
                    rN = MP33ShapeFunctions();
                    break;
                }
            default:
                rIntegrationMethod = GeometryData::GI_GAUSS_2; // default to 3 particles per tri

                std::string warning_msg = "The input number of PARTICLES_PER_ELEMENT: " + std::to_string(ParticlesPerElement);
                warning_msg += " is not available for Triangular" + std::to_string(domain_size) + "D.\n";
                warning_msg += "Available options are: 1, 3, 6, 12, 16 (only 2D), and 33 (only 2D).\n";
                warning_msg += "The default number of particle: 3 is currently assumed.";
                KRATOS_INFO("MPMParticleGeneratorUtility") << "WARNING: " << warning_msg << std::endl;
                break;
            }
        }
        else if (geo_type == GeometryData::Kratos_Hexahedra3D8 || geo_type == GeometryData::Kratos_Quadrilateral2D4)
        {
            switch (ParticlesPerElement)
            {
            case 1:
                rIntegrationMethod = GeometryData::GI_GAUSS_1;
                break;
            case 4:
                rIntegrationMethod = GeometryData::GI_GAUSS_2;
                break;
            case 9:
                rIntegrationMethod = GeometryData::GI_GAUSS_3;
                break;
            case 16:
                rIntegrationMethod = GeometryData::GI_GAUSS_4;
                break;
            default:
                rIntegrationMethod = GeometryData::GI_GAUSS_2; // default to 4 particles per quad

                std::string warning_msg = "The input number of PARTICLES_PER_ELEMENT: " + std::to_string(ParticlesPerElement);
                warning_msg += " is not available for Quadrilateral" + std::to_string(domain_size) + "D.\n";
                warning_msg += "Available options are: 1, 4, 9, 16.\n";
                warning_msg += "The default number of particle: 4 is currently assumed.";
                KRATOS_INFO("MPMParticleGeneratorUtility") << "WARNING: " << warning_msg << std::endl;
                break;
            }
        }

        // Get shape function values
        if (!IsEqualVolumes) rN = rGeom.ShapeFunctionsValues(rIntegrationMethod);
    }

    void DetermineConditionIntegrationMethodAndShapeFunctionValues(const GeometryType& rGeom, const SizeType ParticlesPerCondition,
        IntegrationMethod& rIntegrationMethod, Matrix& rN, bool& IsEqualVolumes)
    {
        const GeometryData::KratosGeometryType geo_type = rGeom.GetGeometryType();
        const SizeType domain_size = rGeom.WorkingSpaceDimension();

        if (geo_type == GeometryData::Kratos_Point2D  || geo_type == GeometryData::Kratos_Point3D)
        {
            switch (ParticlesPerCondition)
            {
                case 0: // Default case
                    IsEqualVolumes = true;
                    rN = ZeroMatrix(1,1);
                    break;
                case 1: // Only nodal
                    IsEqualVolumes = true;
                    rN = ZeroMatrix(1,1);
                    break;
                default:
                    std::string warning_msg = "The input number of PARTICLES_PER_CONDITION: " + std::to_string(ParticlesPerCondition);
                    warning_msg += " is not available for Point" + std::to_string(domain_size) + "D.\n";
                    warning_msg += "Available option is: 1 (default).\n";
                    warning_msg += "The default number of particle: 1 is currently assumed.";
                    KRATOS_INFO("MPMParticleGeneratorUtility") << "WARNING: " << warning_msg << std::endl;
                    break;
            }

            

        }
        else if (geo_type == GeometryData::Kratos_Line2D2  || geo_type == GeometryData::Kratos_Line3D2)
        {
            switch (ParticlesPerCondition)
            {
            case 1:
                rIntegrationMethod = GeometryData::GI_GAUSS_1;
                break;
            case 2:
                rIntegrationMethod = GeometryData::GI_GAUSS_2;
                break;
            case 3:
                rIntegrationMethod = GeometryData::GI_GAUSS_3;
                break;
            case 4:
                rIntegrationMethod = GeometryData::GI_GAUSS_4;
                break;
            case 5:
                rIntegrationMethod = GeometryData::GI_GAUSS_5;
                break;
            default:
                std::string warning_msg = "The input number of PARTICLES_PER_CONDITION: " + std::to_string(ParticlesPerCondition);
                warning_msg += " is not available for Line" + std::to_string(domain_size) + "D.\n";
                warning_msg += "Available options are: 1 (default), 2, 3, 4, 5.\n";
                warning_msg += "The default number of particle: 1 is currently assumed.";
                KRATOS_INFO("MPMParticleGeneratorUtility") << "WARNING: " << warning_msg << std::endl;
                break;
            }
            

        }
        else if (geo_type == GeometryData::Kratos_Triangle3D3)
        {
            switch (ParticlesPerCondition)
            {
            case 1:
                rIntegrationMethod = GeometryData::GI_GAUSS_1;
                break;
            case 3:
                rIntegrationMethod = GeometryData::GI_GAUSS_2;
                break;
            case 6:
                rIntegrationMethod = GeometryData::GI_GAUSS_4;
                break;
            case 12:
                rIntegrationMethod = GeometryData::GI_GAUSS_5;
                break;
            case 16:
                IsEqualVolumes = true;
                KRATOS_INFO("MPMParticleGeneratorUtility") << "WARNING: "
                    << "16 particles per triangle element is only valid for undistorted triangles." << std::endl;
                rN = MP16ShapeFunctions();
                break;
            case 33:
                IsEqualVolumes = true;
                KRATOS_INFO("MPMParticleGeneratorUtility") << "WARNING: "
                    << "33 particles per triangle element is only valid for undistorted triangles." << std::endl;
                rN = MP33ShapeFunctions();
                break;
            default:
                std::string warning_msg = "The input number of PARTICLES_PER_CONDITION: " + std::to_string(ParticlesPerCondition);
                warning_msg += " is not available for Triangular" + std::to_string(domain_size) + "D.\n";
                warning_msg += "Available options are: 1 (default), 3, 6, 12, 16 and 33.\n";
                warning_msg += "The default number of particle: 1 is currently assumed.";
                KRATOS_INFO("MPMParticleGeneratorUtility") << "WARNING: " << warning_msg << std::endl;
                break;
            }
            
        }
        else if (geo_type == GeometryData::Kratos_Quadrilateral3D4)
        {
            switch (ParticlesPerCondition)
            {
            case 1:
                rIntegrationMethod = GeometryData::GI_GAUSS_1;
                break;
            case 4:
                rIntegrationMethod = GeometryData::GI_GAUSS_2;
                break;
            case 9:
                rIntegrationMethod = GeometryData::GI_GAUSS_3;
                break;
            case 16:
                rIntegrationMethod = GeometryData::GI_GAUSS_4;
                break;
            default:
                std::string warning_msg = "The input number of PARTICLES_PER_CONDITION: " + std::to_string(ParticlesPerCondition);
                warning_msg += " is not available for Triangular" + std::to_string(domain_size) + "D.\n";
                warning_msg += "Available options are: 1 (default), 4, 9 and 16.\n";
                warning_msg += "The default number of particle: 1 is currently assumed.";
                KRATOS_INFO("MPMParticleGeneratorUtility") << "WARNING: " << warning_msg << std::endl;
                break;
            }
           
        }

        // Get shape function values
        if (!IsEqualVolumes) rN = rGeom.ShapeFunctionsValues(rIntegrationMethod);
    }

} // end namespace MPMParticleGeneratorUtility
} // end namespace Kratos



