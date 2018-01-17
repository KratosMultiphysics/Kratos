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
//
//

// System includes
#include <set>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/geometry.h"
#include "tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

    /// Auxiliar check functions (from geometry_tester.h)
    /// - All this functions should probably me moved somewhere else.

    /** Gets the corresponding string of the integration method provided.
     * Gets the corresponding string of the integration method provided.
     * @param  geom       Geometry that is used for nothing
     * @param  ThisMethod Input Integration method
     * @return            String with the name of the input integration method
     */
    std::string GetIntegrationName(Geometry<Node<3>>& geom, Geometry<Node<3>>::IntegrationMethod ThisMethod) {
      switch(ThisMethod) {
        case GeometryData::GI_GAUSS_1 :
          return std::string("GI_GAUSS_1");
        case GeometryData::GI_GAUSS_2 :
          return std::string("GI_GAUSS_2");
        case GeometryData::GI_GAUSS_3 :
          return std::string("GI_GAUSS_3");
        case GeometryData::GI_GAUSS_4 :
          return std::string("GI_GAUSS_4");
        case GeometryData::GI_GAUSS_5 :
          return std::string("GI_GAUSS_5");
        case GeometryData::GI_EXTENDED_GAUSS_1 :
          return std::string("GI_EXTENDED_GAUSS_1");
        case GeometryData::GI_EXTENDED_GAUSS_2 :
          return std::string("GI_EXTENDED_GAUSS_2");
        case GeometryData::GI_EXTENDED_GAUSS_3 :
          return std::string("GI_EXTENDED_GAUSS_3");
        case GeometryData::GI_EXTENDED_GAUSS_4 :
          return std::string("GI_EXTENDED_GAUSS_4");
        case GeometryData::GI_EXTENDED_GAUSS_5 :
          return std::string("GI_EXTENDED_GAUSS_5");
        case GeometryData::NumberOfIntegrationMethods :
          return std::string("NumberOfIntegrationMethods");
      };

      return std::string("UnknownIntegrationMethod");
    }

    /** Gets the corresponding string of the geometry name.
     * Gets the corresponding string of the geometry name.
     * @param  geom Input Geometry
     * @return      String corresponding to the name of the input geometry
     */
    std::string GetGeometryName(Geometry<Node<3>>& geom) {
      GeometryData::KratosGeometryType geom_type = geom.GetGeometryType();

      switch(geom_type) {
        case GeometryData::Kratos_generic_type :
          return std::string("Kratos_generic_type");
        case GeometryData::Kratos_Hexahedra3D20 :
          return std::string("Kratos_Hexahedra3D20");
        case GeometryData::Kratos_Hexahedra3D27 :
          return std::string("Kratos_Hexahedra3D27");
        case GeometryData::Kratos_Hexahedra3D8 :
          return std::string("Kratos_Hexahedra3D8");
        case GeometryData::Kratos_Prism3D15 :
          return std::string("Kratos_Prism3D15");
        case GeometryData::Kratos_Prism3D6 :
          return std::string("Kratos_Prism3D6");
        case GeometryData::Kratos_Quadrilateral2D4 :
          return std::string("Kratos_Quadrilateral2D4");
        case GeometryData::Kratos_Quadrilateral2D8 :
          return std::string("Kratos_Quadrilateral2D8");
        case GeometryData::Kratos_Quadrilateral2D9 :
          return std::string("Kratos_Quadrilateral2D9");
        case GeometryData::Kratos_Quadrilateral3D4 :
          return std::string("Kratos_Quadrilateral3D4");
        case GeometryData::Kratos_Quadrilateral3D8 :
          return std::string("Kratos_Quadrilateral3D8");
        case GeometryData::Kratos_Quadrilateral3D9 :
          return std::string("Kratos_Quadrilateral3D9");
        case GeometryData::Kratos_Tetrahedra3D10 :
          return std::string("Kratos_Tetrahedra3D10");
        case GeometryData::Kratos_Tetrahedra3D4 :
          return std::string("Kratos_Tetrahedra3D4");
        case GeometryData::Kratos_Triangle2D3 :
          return std::string("Kratos_Triangle3D3");
        case GeometryData::Kratos_Triangle2D6 :
          return std::string("Kratos_Triangle2D6");
        case GeometryData::Kratos_Triangle3D3 :
          return std::string("Kratos_Triangle3D3");
        case GeometryData::Kratos_Triangle3D6 :
          return std::string("Kratos_Triangle3D6");
        case GeometryData::Kratos_Line2D2 :
          return std::string("Kratos_Line2D2");
        case GeometryData::Kratos_Line2D3 :
          return std::string("Kratos_Line2D3");
        case GeometryData::Kratos_Line3D2 :
          return std::string("Kratos_Line3D2");
        case GeometryData::Kratos_Line3D3 :
          return std::string("Kratos_Line3D3");
        case GeometryData::Kratos_Point2D :
          return std::string("Kratos_Point2D");
        case GeometryData::Kratos_Point3D :
          return std::string("Kratos_Point3D");
        case GeometryData::Kratos_Sphere3D1 :
          return std::string("Kratos_Sphere3D1");
      };

      return std::string("UnknownGeometry");
    }

    /** Computes the linear strain matrix.
     * Computes the linear strain matrix which is useful to verify that
     * a constant strain can be correctly reproduced
     * @param B               [description]
     * @param DN_DX           [description]
     * @param number_of_nodes Nuber of nodes of the geometry
     * @param dimension       Dimension (i.e. 1, 2 or 3)
     */
    void CalculateB(Matrix& B, Matrix& DN_DX, const unsigned int number_of_nodes, const unsigned int dimension) {
      if(dimension == 2) {
        B.resize(3, 2*number_of_nodes, false);
      } else {
        B.resize(6, 3*number_of_nodes, false);
      }

      for(unsigned int i = 0; i < number_of_nodes; i++) {
        unsigned int index = dimension * i;

        if(dimension == 2) {
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

    /** Verifies the area of the geometry using the integration method.
     * Verifies the area of the geometry using the integration method.
     * @param  geom           Geometry to be tested
     * @param  ThisMethod     Integration method used
     * @param  reference_area Expected area
     * @param  error_msg      Buffer to write the error message
     * @return                Area claculated using the selected integration method.
     */
    double CalculateAreaByIntegration(Geometry<Node<3>>& geom, Geometry<Node<3> >::IntegrationMethod ThisMethod) {
      double area = 0.0;

      if(geom.WorkingSpaceDimension() != geom.LocalSpaceDimension()) {
        KRATOS_ERROR << "VerifyStrainExactness can not be used if LocalSpaceDimension and WorkingSpaceDimension do not coincide --> geometry is " << GetGeometryName(geom) << std::endl;
      }

      const Element::GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints( ThisMethod );

      if(integration_points.size() == 0) {
        KRATOS_ERROR << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " -- the integration method is not supported " << std::endl;
      }

      // Resizing jacobian inverses containers
      Matrix InvJ0(geom.WorkingSpaceDimension(), geom.WorkingSpaceDimension());

      Element::GeometryType::JacobiansType J0;
      J0 = geom.Jacobian( J0, ThisMethod );

      Vector determinants;
      geom.DeterminantOfJacobian(determinants, ThisMethod);

      for(unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++) {
        const double IntegrationWeight = integration_points[PointNumber].Weight();

        // Calculating and storing inverse of the jacobian and the parameters needed
        double DetJ0 = MathUtils<double>::Det( J0[PointNumber] );

        if( std::abs(determinants[PointNumber] - DetJ0)/std::abs(DetJ0) > 1e-13) {
          KRATOS_ERROR << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " --> " << " determinant as computed from DeterminantOfJacobian does not match the value computed by taking the determinant of J "  << std::endl;
        }

        // Calculating the total area
        area += DetJ0 * IntegrationWeight;
      }

      return area;
    }

    /** Verifies that a displacement field produces the expected strain distribution.
     * Verifies that a displacement field which varies linearly in space, produces the expected strain distribution.
     * This shall be considered a test for shape function derivatives
     * @param geom       Geometry to be tested
     * @param ThisMethod Integration method used
     * @param error_msg  Buffer to write the error message
     */
    void VerifyStrainExactness(Geometry<Node<3>>& geom,  Geometry<Node<3> >::IntegrationMethod ThisMethod) {

      const Element::GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints( ThisMethod );
      const unsigned int number_of_nodes = geom.PointsNumber();
      const unsigned int dim = geom.WorkingSpaceDimension();

      if(dim != geom.LocalSpaceDimension()) {
        KRATOS_THROW_ERROR(std::logic_error,"VerifyStrainExactness can not be used if LocalSpaceDimension and WorkingSpaceDimension do not coincide ",GetGeometryName(geom) );
      }

      if(integration_points.size() == 0) {
        KRATOS_ERROR << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " -- the integration method is not supported " << std::endl;
      } else {
        // Charlie: Shouldn't this be initialized for dim = 1?
        unsigned int strain_size = 1;

        if(dim == 2) {
          strain_size = 3;
        } else {
          strain_size = 6;
        }

        // Definition of the expected strain
        Matrix MatrixA(dim,dim);
        Vector VectorB(dim);

        for(unsigned int i=0; i<dim; i++) {
          VectorB[i]=i*i+0.567; //arbitrary values
          for(unsigned int j=0; j<dim; j++) {
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
        Vector displacements(dim*number_of_nodes);

        const Element::GeometryType::ShapeFunctionsGradientsType& DN_De = geom.ShapeFunctionsLocalGradients( ThisMethod );
        const Matrix& Ncontainer = geom.ShapeFunctionsValues( ThisMethod );

        Element::GeometryType::JacobiansType J0;
        geom.Jacobian( J0, ThisMethod );

        // Untested functions to be tested
        Element::GeometryType::ShapeFunctionsGradientsType DN_DX_geom;
        geom.ShapeFunctionsIntegrationPointsGradients( DN_DX_geom, ThisMethod );

        Element::GeometryType::JacobiansType Jinv;
        geom.InverseOfJacobian(Jinv, ThisMethod);

        bool succesful = true;
        for(unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++) {
          // Check that shape functions sum to 1
          double sum = 0.0;
          for(unsigned int k = 0; k<number_of_nodes; k++) {
            sum += Ncontainer(PointNumber,k);
          }

          if(std::abs(sum-1.0)>1e-14) {
            KRATOS_ERROR << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " --> " << " error: shape functions do not sum to 1 on gauss point" << std::endl;
          }

          // Calculating and storing inverse of the jacobian and the parameters needed
          MathUtils<double>::InvertMatrix( J0[PointNumber], InvJ0, DetJ0 );
          DN_DX  = prod( DN_De[PointNumber], InvJ0 );

          // Check that the shape function gradients as obtained from the geomety match what is obtained here starting from the local_gradients
          if(norm_frobenius(DN_DX_geom[PointNumber] - DN_DX)/norm_frobenius(DN_DX) > 1e-13) {
            KRATOS_ERROR << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " -->  " << std::endl;
            KRATOS_ERROR << "     error: shape function gradients are wrongly calculated in function ShapeFunctionsIntegrationPointsGradients: DN_DX_geom " << DN_DX_geom[PointNumber] << " vs " << DN_DX << std::endl;
            KRATOS_ERROR << " norm_frobenius(DN_DX_geom[PointNumber] - DN_DX)/norm_frobenius(DN_DX) = " << norm_frobenius(DN_DX_geom[PointNumber] - DN_DX)/norm_frobenius(DN_DX) <<std::endl;
          }

          if(norm_frobenius(Jinv[PointNumber] - InvJ0)/norm_frobenius(InvJ0) > 1e-13) {
            KRATOS_ERROR << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " --> " << std::endl;
            KRATOS_ERROR << "     error: shape function gradients are wrongly calculated in function ShapeFunctionsIntegrationPointsGradients: DN_DX_geom " << DN_DX_geom[PointNumber] << " vs " << DN_DX << std::endl;
            KRATOS_ERROR << " norm_frobenius(Jinv[PointNumber] - InvJ0)/norm_frobenius(InvJ0) = " << norm_frobenius(Jinv[PointNumber] - InvJ0)/norm_frobenius(InvJ0) <<std::endl;
          }

          CalculateB(B, DN_DX, number_of_nodes, dim);

          // Calculate a displacement_field which varies linearly in the space
          for(unsigned int i=0; i<number_of_nodes; i++) {
            const array_1d<double,3>& coords = geom[i].Coordinates();
            Vector disp(dim);

            for(unsigned int k=0; k<dim; k++) {
              disp[k] = VectorB[k];
              for(unsigned int l=0; l<dim; l++) {
                disp[k] += MatrixA(k,l)*coords[l] ;
              }
            }

            // Vector disp = prod(MatrixA,) + VectorB;
            for(unsigned int k=0; k<dim; k++) {
              displacements[i*dim+k] = disp[k];
            }
          }

          Vector strain = prod(B,displacements);
          Vector strain_err = strain-expected_strain;

          if( norm_2(strain_err)/norm_2(expected_strain) < 1e-14) {
            //do nothing
          } else {
            succesful = false;
            KRATOS_ERROR << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " --> " << " error: expected strain found was not correctly recovered on gauss point. recovered strain = " << strain << " expected value "  << expected_strain << std::endl;
          }
        }

        if(succesful == true) {
          // Charlie: Do nothing
          // std::cout << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " --> " << " Strain Calculation Test: OK "  << std::endl;
        }
      }
    }

} // namespace Testing.
} // namespace Kratos.
