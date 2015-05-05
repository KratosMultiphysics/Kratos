//
//   Project Name:        Kratos
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes
#include "includes/define.h"
#include "custom_elements/nurbs_poisson_2d.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{



    //************************************************************************************
    //************************************************************************************
    NurbsPoisson2D::NurbsPoisson2D(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************
    NurbsPoisson2D::NurbsPoisson2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {

    }

    Element::Pointer NurbsPoisson2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new NurbsPoisson2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    NurbsPoisson2D::~NurbsPoisson2D()
    {
    }

    //************************************************************************************
    //************************************************************************************
    void NurbsPoisson2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        GeometryType& rGeom = this->GetGeometry();
        const SizeType NumNodes = rGeom.PointsNumber();
        const SizeType Dim = rGeom.WorkingSpaceDimension();
        KRATOS_WATCH(NumNodes)
        KRATOS_WATCH(Dim)

        // Check sizes and initialize
        if( rLeftHandSideMatrix.size1() != NumNodes )
            rLeftHandSideMatrix.resize(NumNodes,NumNodes);

        rLeftHandSideMatrix = ZeroMatrix(NumNodes,NumNodes);

        if( rRightHandSideVector.size() != NumNodes )
            rRightHandSideVector.resize(NumNodes);

        rRightHandSideVector = ZeroVector(NumNodes);

        // Shape functions and integration points
        ShapeFunctionDerivativesArrayType DN_DX;
        Matrix NContainer;
        VectorType GaussWeights;
        this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);
        const unsigned int NumGauss = GaussWeights.size();
        KRATOS_WATCH(DN_DX)
        KRATOS_WATCH(GaussWeights)

        // Loop on integration points
        for (unsigned int g = 0; g < NumGauss; g++)
        {
            const double GaussWeight = GaussWeights[g];
//            const ShapeFunctionsType& N = row(NContainer,g);
            const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];

            // Evaluate required variables at the integration point
            double Conductivity = this->GetProperties().GetValue(CONDUCTIVITY);

            // Add convection, stabilization and RHS contributions to the local system equation
            for (SizeType i = 0; i < NumNodes; ++i)
            {
                // LHS contribution
                for (SizeType j = 0; j < NumNodes; ++j)
                {
                    double Lij = 0.0;
                    for (SizeType d = 0; d < Dim; ++d)
                        Lij += rDN_DX(i,d) * rDN_DX(j,d);

                    rLeftHandSideMatrix(i,j) += GaussWeight * Lij * Conductivity;
                }

                // RHS contribution
                // Source term to be added

            }
        }
        KRATOS_WATCH(this->GetProperties().GetValue(CONDUCTIVITY))
        KRATOS_WATCH(rLeftHandSideMatrix);

        // Residual formulation
        Vector Values = ZeroVector(NumNodes);
        for (unsigned int i = 0; i < NumNodes; i++)
            Values[i] = rGeom[i].FastGetSolutionStepValue(TEMPERATURE);
        KRATOS_WATCH(Values);


        rRightHandSideVector -= prod(rLeftHandSideMatrix,Values);

        KRATOS_WATCH(rRightHandSideVector);

        EquationIdVectorType EqId;
        EquationIdVector(EqId,rCurrentProcessInfo);
        std::cout << "EqIds: ";
        for (unsigned int i = 0; i < NumNodes; i++)
               std::cout << EqId[i] << " (" << GetGeometry()[i].Id() << ")" << " " ;
        std::cout << std::endl;
        //KRATOS_WATCH(EqId);
        double Tgauss = 0.0;
        double DTDXgauss = 0.0;
        double DTDYgauss = 0.0;
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            Tgauss += Values[i]*NContainer(0,i);
            DTDXgauss += Values[i]*DN_DX[0](i,0);
            DTDYgauss += Values[i]*DN_DX[0](i,1);
        }

        KRATOS_WATCH(Tgauss);
        KRATOS_WATCH(DTDXgauss);
        KRATOS_WATCH(DTDYgauss);


        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************
    void NurbsPoisson2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
    }

    //************************************************************************************
    //************************************************************************************
    // this subroutine calculates the nodal contributions for the explicit steps of the
    // fractional step procedure
    void NurbsPoisson2D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************
    void NurbsPoisson2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
    {
        unsigned int number_of_nodes = GetGeometry().PointsNumber();
        if(rResult.size() != number_of_nodes)
            rResult.resize(number_of_nodes,false);

        for (unsigned int i=0;i<number_of_nodes;i++)
                rResult[i] = GetGeometry()[i].GetDof(TEMPERATURE).EquationId();
    }

    //************************************************************************************
    //************************************************************************************
      void NurbsPoisson2D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_WATCH("aaaaaaa")
        unsigned int number_of_nodes = GetGeometry().PointsNumber();
        if(ElementalDofList.size() != number_of_nodes)
            ElementalDofList.resize(number_of_nodes);

        for (unsigned int i=0;i<number_of_nodes;i++)
        {
            ElementalDofList[i] = GetGeometry()[i].pGetDof(TEMPERATURE);
            KRATOS_WATCH(ElementalDofList[i]);
        }

    }


      void NurbsPoisson2D::CalculateGeometryData(ShapeFunctionDerivativesArrayType &rDN_DX,
                                                 Matrix &ShapeFunctionsValues,
                                                 Vector &rGaussWeights)
      {
          const GeometryType& rGeom = this->GetGeometry();
          Vector DetJ;

          rGeom.ShapeFunctionsIntegrationPointsGradients(rDN_DX,
                                                         DetJ,
                                                         GeometryData::GI_GAUSS_3,
                                                         ShapeFunctionsValues);
//          rGeom.ShapeFunctionsIntegrationPointsGradients(rDN_DX,
//                                                         DetJ,
//                                                         GeometryData::GI_GAUSS_3);



//          NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_3);
          const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_3);

//          KRATOS_WATCH(IntegrationPoints.size());
          rGaussWeights.resize(rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_3),false);

          double AreaGuess = 0.0;
          double TemperatureInsideGeom(0);
          double sumSF(0),sumDervY(0);
          for (unsigned int g = 0; g < rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_3); g++)
          {
              rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
              AreaGuess += rGaussWeights[g];
//              KRATOS_WATCH(IntegrationPoints[g].Weight());
//              KRATOS_WATCH(AreaGuess);

                TemperatureInsideGeom = 0;
                sumSF=0;
                sumDervY=0;
              for (unsigned int a=0;a<rGeom.PointsNumber();a++)
              {
//                KRATOS_WATCH(rGeom(a)->GetSolutionStepValue(TEMPERATURE,0));
//                  KRATOS_WATCH(rGeom(a)->Id());
//                  KRATOS_WATCH(rDN_DX[g](a,0));
//                  KRATOS_WATCH(rDN_DX[g](a,1));
                  sumSF += rDN_DX[g](a,1)*rGeom(a)->X();
                  sumDervY += rDN_DX[g](a,1)*rGeom(a)->Y();
               // KRATOS_WATCH(ShapeFunctionsValues(g,a));
                  TemperatureInsideGeom += rGeom(a)->GetSolutionStepValue(TEMPERATURE,0) * ShapeFunctionsValues(g,a);
                }

//              KRATOS_WATCH(TemperatureInsideGeom);
//                KRATOS_WATCH(sumSF);
//                KRATOS_WATCH(sumDervY);
//                KRATOS_WATCH(AreaGuess);
//                      for (int i=0;i<8;i++)
//              {
//                      std::cout<<i+1<<". ShapeFunctionValue =" <<NContainer(g,i)<<std::endl;
//              }
          }


          /*
          const GeometryType& rGeom = this->GetGeometry();
          const SizeType NumNodes = rGeom.PointsNumber();
          const unsigned int NumGauss = rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2);

          // Initialize arrays to proper size
          rDN_DX.resize(NumGauss);
          rDetJ.resize(NumGauss);

          const GeometryType::ShapeFunctionsGradientsType& DN_De = rGeom.ShapeFunctionsLocalGradients( GeometryData::GI_GAUSS_2 );

          // Temporary container for inverse of J
          Matrix InvJ;

          GeometryType::JacobiansType J;
          rGeom.Jacobian( J, GeometryData::GI_GAUSS_2 );

          for (unsigned int g = 0; g < NumGauss; g++)
          {
              // calculate inverse of the jacobian and its determinant
              MathUtils<double>::InvertMatrix( J[g], InvJ, rDetJ[g] );

              // calculate the shape function derivatives in global coordinates
              rDN_DX[g].resize(NumNodes,TDim);
              noalias( rDN_DX[g] ) = prod( DN_De[g], InvJ );
          }*/
      }


} // Namespace Kratos
