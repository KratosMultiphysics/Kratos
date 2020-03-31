//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//                   Vicente Mataix Ferrandiz
//
//

// System includes
#include <set>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/geometry.h"
#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

    /// Auxiliar check functions (from geometry_tester.h)
    /// - All this functions should probably me moved somewhere else.

    /** Gets the corresponding string of the integration method provided.
     * Gets the corresponding string of the integration method provided.
     * @param  ThisGeometry       Geometry that is used for nothing
     * @param  ThisMethod Input Integration method
     * @return            String with the name of the input integration method
     */
    std::string GetIntegrationName(
        GeometryType& ThisGeometry,
        GeometryType::IntegrationMethod ThisMethod
        )
    {
      switch(ThisMethod) {
        case GeometryData::GI_GAUSS_1 :
          return "GI_GAUSS_1";
        case GeometryData::GI_GAUSS_2 :
          return "GI_GAUSS_2";
        case GeometryData::GI_GAUSS_3 :
          return "GI_GAUSS_3";
        case GeometryData::GI_GAUSS_4 :
          return "GI_GAUSS_4";
        case GeometryData::GI_GAUSS_5 :
          return "GI_GAUSS_5";
        case GeometryData::GI_EXTENDED_GAUSS_1 :
          return "GI_EXTENDED_GAUSS_1";
        case GeometryData::GI_EXTENDED_GAUSS_2 :
          return "GI_EXTENDED_GAUSS_2";
        case GeometryData::GI_EXTENDED_GAUSS_3 :
          return "GI_EXTENDED_GAUSS_3";
        case GeometryData::GI_EXTENDED_GAUSS_4 :
          return "GI_EXTENDED_GAUSS_4";
        case GeometryData::GI_EXTENDED_GAUSS_5 :
          return "GI_EXTENDED_GAUSS_5";
        case GeometryData::NumberOfIntegrationMethods :
          return "NumberOfIntegrationMethods";
      };

      return "UnknownIntegrationMethod";
    }

    /** Gets the corresponding string of the geometry name.
     * Gets the corresponding string of the geometry name.
     * @param  ThisGeometry Input Geometry
     * @return      String corresponding to the name of the input geometry
     */
    std::string GetGeometryName(GeometryType& ThisGeometry)
    {
      GeometryData::KratosGeometryType geom_type = ThisGeometry.GetGeometryType();

      switch(geom_type) {
        case GeometryData::Kratos_generic_type :
          return "Kratos_generic_type";
        case GeometryData::Kratos_Hexahedra3D20 :
          return "Kratos_Hexahedra3D20";
        case GeometryData::Kratos_Hexahedra3D27 :
          return "Kratos_Hexahedra3D27";
        case GeometryData::Kratos_Hexahedra3D8 :
          return "Kratos_Hexahedra3D8";
        case GeometryData::Kratos_Prism3D15 :
          return "Kratos_Prism3D15";
        case GeometryData::Kratos_Prism3D6 :
          return "Kratos_Prism3D6";
        case GeometryData::Kratos_Quadrilateral2D4 :
          return "Kratos_Quadrilateral2D4";
        case GeometryData::Kratos_Quadrilateral2D8 :
          return "Kratos_Quadrilateral2D8";
        case GeometryData::Kratos_Quadrilateral2D9 :
          return "Kratos_Quadrilateral2D9";
        case GeometryData::Kratos_Quadrilateral3D4 :
          return "Kratos_Quadrilateral3D4";
        case GeometryData::Kratos_Quadrilateral3D8 :
          return "Kratos_Quadrilateral3D8";
        case GeometryData::Kratos_Quadrilateral3D9 :
          return "Kratos_Quadrilateral3D9";
        case GeometryData::Kratos_Tetrahedra3D10 :
          return "Kratos_Tetrahedra3D10";
        case GeometryData::Kratos_Tetrahedra3D4 :
          return "Kratos_Tetrahedra3D4";
        case GeometryData::Kratos_Triangle2D3 :
          return "Kratos_Triangle3D3";
        case GeometryData::Kratos_Triangle2D6 :
          return "Kratos_Triangle2D6";
        case GeometryData::Kratos_Triangle3D3 :
          return "Kratos_Triangle3D3";
        case GeometryData::Kratos_Triangle3D6 :
          return "Kratos_Triangle3D6";
        case GeometryData::Kratos_Line2D2 :
          return "Kratos_Line2D2";
        case GeometryData::Kratos_Line2D3 :
          return "Kratos_Line2D3";
        case GeometryData::Kratos_Line3D2 :
          return "Kratos_Line3D2";
        case GeometryData::Kratos_Line3D3 :
          return "Kratos_Line3D3";
        case GeometryData::Kratos_Point2D :
          return "Kratos_Point2D";
        case GeometryData::Kratos_Point3D :
          return "Kratos_Point3D";
        case GeometryData::Kratos_Sphere3D1 :
          return "Kratos_Sphere3D1";
      };

      return "UnknownGeometry";
    }

    /** Computes the linear strain matrix.
     * Computes the linear strain matrix which is useful to verify that
     * a constant strain can be correctly reproduced
     * @param B               [description]
     * @param DN_DX           [description]
     * @param NumberOfNodes   Number of nodes of the geometry
     * @param Dimension       Dimension (i.e. 1, 2 or 3)
     */
    void CalculateB(
        Matrix& B,
        Matrix& DN_DX,
        const SizeType NumberOfNodes,
        const SizeType Dimension
        )
    {
        if(Dimension == 2) {
                B.resize(3, 2*NumberOfNodes, false);
        } else {
                B.resize(6, 3*NumberOfNodes, false);
        }

        for(SizeType i = 0; i < NumberOfNodes; i++) {
            SizeType index = Dimension * i;

            if(Dimension == 2) {
                B(0, index + 0) = DN_DX(i, 0);
                B(0, index + 1) = 0.0;

                B(1, index + 0) = 0.0;
                B(1, index + 1) = DN_DX(i, 1);

                B(2, index + 0) = DN_DX(i, 1);
                B(2, index + 1) = DN_DX(i, 0);
            } else {
                B(0, index + 0) = DN_DX(i, 0);
                B(0, index + 1) = 0.0;
                B(0, index + 2) = 0.0;

                B(1, index + 0) = 0.0;
                B(1, index + 1) = DN_DX(i, 1);
                B(1, index + 2) = 0.0;

                B(2, index + 0) = 0.0;
                B(2, index + 1) = 0.0;
                B(2, index + 2) = DN_DX(i, 2);

                B(3, index + 0) = DN_DX(i, 1);
                B(3, index + 1) = DN_DX(i, 0);
                B(3, index + 2) = 0.0;

                B(4, index + 0) = 0.0;
                B(4, index + 1) = DN_DX(i, 2);
                B(4, index + 2) = DN_DX(i, 1);

                B(5, index + 0) = DN_DX(i, 2);
                B(5, index + 1) = 0.0;
                B(5, index + 2) = DN_DX(i, 0);
            }
        }
    }

    /** Verifies the area of the ThisGeometryetry using the integration method.
     * Verifies the area of the ThisGeometryetry using the integration method.
     * @param  ThisGeometry           Geometry to be tested
     * @param  ThisMethod     Integration method used
     * @param  reference_area Expected area
     * @param  error_msg      Buffer to write the error message
     * @return                Area claculated using the selected integration method.
     */
    double CalculateAreaByIntegration(
        GeometryType& ThisGeometry,
        GeometryType::IntegrationMethod ThisMethod
        )
    {
        double area = 0.0;

        if(ThisGeometry.WorkingSpaceDimension() != ThisGeometry.LocalSpaceDimension()) {
            KRATOS_ERROR << "VerifyStrainExactness can not be used if LocalSpaceDimension and WorkingSpaceDimension do not coincide --> ThisGeometryetry is " << GetGeometryName(ThisGeometry) << std::endl;
        }

        const Element::GeometryType::IntegrationPointsArrayType& integration_points = ThisGeometry.IntegrationPoints( ThisMethod );

        if(integration_points.size() == 0) {
            KRATOS_ERROR << "Geometry Type = " << GetGeometryName(ThisGeometry) << " - IntegrationMethod = " << GetIntegrationName(ThisGeometry,ThisMethod) << " -- the integration method is not supported " << std::endl;
        }

        // Resizing jacobian inverses containers
        Matrix InvJ0(ThisGeometry.WorkingSpaceDimension(), ThisGeometry.WorkingSpaceDimension());

        Element::GeometryType::JacobiansType J0;
        J0 = ThisGeometry.Jacobian( J0, ThisMethod );

        Vector determinants;
        ThisGeometry.DeterminantOfJacobian(determinants, ThisMethod);

        for(SizeType PointNumber = 0; PointNumber < integration_points.size(); PointNumber++) {
            const double IntegrationWeight = integration_points[PointNumber].Weight();

            // Calculating and storing inverse of the jacobian and the parameters needed
            double DetJ0 = MathUtils<double>::Det( J0[PointNumber] );

            KRATOS_ERROR_IF( std::abs(determinants[PointNumber] - DetJ0)/std::abs(DetJ0) > 1e-13) << "Geometry Type = " << GetGeometryName(ThisGeometry) << " - IntegrationMethod = " << GetIntegrationName(ThisGeometry,ThisMethod) << " --> " << " determinant as computed from DeterminantOfJacobian does not match the value computed by taking the determinant of J "  << std::endl;

            // Calculating the total area
            area += DetJ0 * IntegrationWeight;
        }

        return area;
    }

    /** Verifies that a displacement field produces the expected strain distribution.
     * Verifies that a displacement field which varies linearly in space, produces the expected strain distribution.
     * This shall be considered a test for shape function derivatives
     * @param ThisGeometry       Geometry to be tested
     * @param ThisMethod Integration method used
     * @param error_msg  Buffer to write the error message
     */
    void VerifyStrainExactness(
        GeometryType& ThisGeometry,
        GeometryType::IntegrationMethod ThisMethod
        )
    {
        const Element::GeometryType::IntegrationPointsArrayType& integration_points = ThisGeometry.IntegrationPoints( ThisMethod );
        const SizeType NumberOfNodes = ThisGeometry.PointsNumber();
        const SizeType dim = ThisGeometry.WorkingSpaceDimension();

        if(dim != ThisGeometry.LocalSpaceDimension()) {
            KRATOS_THROW_ERROR(std::logic_error,"VerifyStrainExactness can not be used if LocalSpaceDimension and WorkingSpaceDimension do not coincide ",GetGeometryName(ThisGeometry) );
        }

        if(integration_points.size() == 0) {
            KRATOS_ERROR << "Geometry Type = " << GetGeometryName(ThisGeometry) << " - IntegrationMethod = " << GetIntegrationName(ThisGeometry,ThisMethod) << " -- the integration method is not supported " << std::endl;
        } else {
            // Charlie: Shouldn't this be initialized for dim = 1?
            SizeType strain_size = 1;

            if(dim == 2) {
                strain_size = 3;
            } else {
                strain_size = 6;
            }

            // Definition of the expected strain
            Matrix MatrixA(dim,dim);
            Vector VectorB(dim);

            for(SizeType i=0; i<dim; i++) {
                VectorB[i]=i*i+0.567; //arbitrary values
                for(SizeType j=0; j<dim; j++) {
                    MatrixA(i,j)=i*j + 0.12345; //initialization fo the values of this matrix is arbitrary
                }
            }

            Vector expected_strain(strain_size);

            if(dim == 2) {
                expected_strain[0] = MatrixA(0,0);
                expected_strain[1] = MatrixA(1,1);
                expected_strain[2] = MatrixA(0,1)+MatrixA(1,0);
            } else {
                expected_strain[0] = MatrixA(0,0);
                expected_strain[1] = MatrixA(1,1);
                expected_strain[2] = MatrixA(2,2);
                expected_strain[3] = MatrixA(0,1)+MatrixA(1,0);
                expected_strain[4] = MatrixA(1,2)+MatrixA(2,1);
                expected_strain[5] = MatrixA(0,2)+MatrixA(2,0);
            }

            // Resizing jacobian inverses containers
            Matrix InvJ0(dim,dim);
            double DetJ0;
            Matrix B;
            Matrix DN_DX;
            Vector displacements(dim*NumberOfNodes);

            const Element::GeometryType::ShapeFunctionsGradientsType& DN_De = ThisGeometry.ShapeFunctionsLocalGradients( ThisMethod );
            const Matrix& Ncontainer = ThisGeometry.ShapeFunctionsValues( ThisMethod );

            Element::GeometryType::JacobiansType J0;
            ThisGeometry.Jacobian( J0, ThisMethod );

            // Untested functions to be tested
            Element::GeometryType::ShapeFunctionsGradientsType DN_DX_ThisGeometry;
            ThisGeometry.ShapeFunctionsIntegrationPointsGradients( DN_DX_ThisGeometry, ThisMethod );

            Element::GeometryType::JacobiansType Jinv;
            ThisGeometry.InverseOfJacobian(Jinv, ThisMethod);

            bool succesful = true;
            for(SizeType PointNumber = 0; PointNumber < integration_points.size(); PointNumber++) {
                // Check that shape functions sum to 1
                double sum = 0.0;
                for(SizeType k = 0; k<NumberOfNodes; k++) {
                    sum += Ncontainer(PointNumber,k);
                }

                if(std::abs(sum-1.0)>1e-14) {
                    KRATOS_ERROR << "Geometry Type = " << GetGeometryName(ThisGeometry) << " - IntegrationMethod = " << GetIntegrationName(ThisGeometry,ThisMethod) << " --> " << " error: shape functions do not sum to 1 on gauss point" << std::endl;
                }

                // Calculating and storing inverse of the jacobian and the parameters needed
                MathUtils<double>::InvertMatrix( J0[PointNumber], InvJ0, DetJ0 );
                DN_DX  = prod( DN_De[PointNumber], InvJ0 );

                // Check that the shape function gradients as obtained from the ThisGeometryety match what is obtained here starting from the local_gradients
                if(norm_frobenius(DN_DX_ThisGeometry[PointNumber] - DN_DX)/norm_frobenius(DN_DX) > 1e-13) {
                    KRATOS_ERROR << "Geometry Type = " << GetGeometryName(ThisGeometry) << " - IntegrationMethod = " << GetIntegrationName(ThisGeometry,ThisMethod) << " -->  " << std::endl;
                    KRATOS_ERROR << "     error: shape function gradients are wrongly calculated in function ShapeFunctionsIntegrationPointsGradients: DN_DX_ThisGeometry " << DN_DX_ThisGeometry[PointNumber] << " vs " << DN_DX << std::endl;
                    KRATOS_ERROR << " norm_frobenius(DN_DX_ThisGeometry[PointNumber] - DN_DX)/norm_frobenius(DN_DX) = " << norm_frobenius(DN_DX_ThisGeometry[PointNumber] - DN_DX)/norm_frobenius(DN_DX) <<std::endl;
                }

                if(norm_frobenius(Jinv[PointNumber] - InvJ0)/norm_frobenius(InvJ0) > 1e-13) {
                    KRATOS_ERROR << "Geometry Type = " << GetGeometryName(ThisGeometry) << " - IntegrationMethod = " << GetIntegrationName(ThisGeometry,ThisMethod) << " --> " << std::endl;
                    KRATOS_ERROR << "     error: shape function gradients are wrongly calculated in function ShapeFunctionsIntegrationPointsGradients: DN_DX_ThisGeometry " << DN_DX_ThisGeometry[PointNumber] << " vs " << DN_DX << std::endl;
                    KRATOS_ERROR << " norm_frobenius(Jinv[PointNumber] - InvJ0)/norm_frobenius(InvJ0) = " << norm_frobenius(Jinv[PointNumber] - InvJ0)/norm_frobenius(InvJ0) <<std::endl;
                }

                CalculateB(B, DN_DX, NumberOfNodes, dim);

                // Calculate a displacement_field which varies linearly in the space
                for(SizeType i=0; i<NumberOfNodes; i++) {
                    const array_1d<double,3>& coords = ThisGeometry[i].Coordinates();
                    Vector disp(dim);

                    for(SizeType k=0; k<dim; k++) {
                        disp[k] = VectorB[k];
                        for(SizeType l=0; l<dim; l++) {
                            disp[k] += MatrixA(k,l)*coords[l] ;
                        }
                    }

                    // Vector disp = prod(MatrixA,) + VectorB;
                    for(SizeType k=0; k<dim; k++) {
                        displacements[i*dim+k] = disp[k];
                    }
                }

                Vector strain = prod(B,displacements);
                Vector strain_err = strain-expected_strain;

                if( norm_2(strain_err)/norm_2(expected_strain) < 1e-14) {
                    //do nothing
                } else {
                    succesful = false;
                    KRATOS_ERROR << "Geometry Type = " << GetGeometryName(ThisGeometry) << " - IntegrationMethod = " << GetIntegrationName(ThisGeometry,ThisMethod) << " --> " << " error: expected strain found was not correctly recovered on gauss point. recovered strain = " << strain << " expected value "  << expected_strain << std::endl;
                }
            }

            KRATOS_CHECK(succesful);
        }
    }

    /// Test self assigned geometry Id
    KRATOS_TEST_CASE_IN_SUITE(GeometryIdSelfAssigned, KratosCoreGeometriesFastSuite) {
        auto this_geometry = Geometry<Point>();

        KRATOS_CHECK_IS_FALSE(this_geometry.IsIdGeneratedFromString());
        KRATOS_CHECK(this_geometry.IsIdSelfAssigned());

        this_geometry.SetId(2);
        KRATOS_CHECK_IS_FALSE(this_geometry.IsIdGeneratedFromString());
        KRATOS_CHECK_IS_FALSE(this_geometry.IsIdSelfAssigned());

        this_geometry.SetId("ThisGeometry");
        KRATOS_CHECK(this_geometry.IsIdGeneratedFromString());
        KRATOS_CHECK_IS_FALSE(this_geometry.IsIdSelfAssigned());
    }

    /// Test geometry Id with name
    KRATOS_TEST_CASE_IN_SUITE(GeometryName, KratosCoreGeometriesFastSuite) {
        auto this_geometry = Geometry<Point>("Geometry1");

        KRATOS_CHECK(this_geometry.IsIdGeneratedFromString());
        KRATOS_CHECK_IS_FALSE(this_geometry.IsIdSelfAssigned());
        KRATOS_CHECK_EQUAL(this_geometry.Id(), Geometry<Point>::GenerateId("Geometry1"));
    }

    /// Test geometry Id
    KRATOS_TEST_CASE_IN_SUITE(GeometryId, KratosCoreGeometriesFastSuite) {
        auto this_geometry = Geometry<Point>(1);

        KRATOS_CHECK_IS_FALSE(this_geometry.IsIdGeneratedFromString());
        KRATOS_CHECK_IS_FALSE(this_geometry.IsIdSelfAssigned());
        KRATOS_CHECK_EQUAL(this_geometry.Id(), 1);

        // Check for higher Id.
        auto this_geometry_2 = Geometry<Point>(717);
        KRATOS_CHECK_IS_FALSE(this_geometry_2.IsIdGeneratedFromString());
        KRATOS_CHECK_IS_FALSE(this_geometry_2.IsIdSelfAssigned());
        KRATOS_CHECK_EQUAL(this_geometry_2.Id(), 717);
    }
} // namespace Testing.
} // namespace Kratos.
