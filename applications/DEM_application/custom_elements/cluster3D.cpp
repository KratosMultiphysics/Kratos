//
//   Project Name:        Kratos
//   Last Modified by:    $Author: M.Santasusana $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "cluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application.h"
#include "includes/kratos_flags.h"



namespace Kratos
{
     // using namespace GeometryFunctions;

      Cluster3D::Cluster3D()
        : Element()
        
        {

        }
    
        
      Cluster3D::Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
        
        {
      
        }

      Cluster3D::Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
      : Element(NewId, pGeometry, pProperties)
      {
       
      }

      Cluster3D::Cluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
      : Element(NewId, ThisNodes)
      {
          
      }

      Element::Pointer Cluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
      {
           return Element::Pointer(new Cluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));

      }      

      /// Destructor.
      Cluster3D::~Cluster3D(){}

      void Cluster3D::Initialize()
      {
          KRATOS_TRY 
          
            KRATOS_WATCH("HOLA SOY UN CLUSTER")

          KRATOS_CATCH( "" )
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void Cluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector,
                ProcessInfo& rCurrentProcessInfo){}
  
      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void Cluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo){}

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void Cluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
      {
 
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void Cluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo){}

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void Cluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
      {
          KRATOS_TRY

          ElementalDofList.resize(0);

//           for (unsigned int i = 0; i < GetGeometry().size(); i++){
//               ElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_X));
//               ElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Y));
//               if (GetGeometry().WorkingSpaceDimension() == 3){
//                   ElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Z));
//               }
// 
//               ElementalDofList.push_back(GetGeometry()[i].pGetDof(ANGULAR_VELOCITY_X));
//               ElementalDofList.push_back(GetGeometry()[i].pGetDof(ANGULAR_VELOCITY_Y));
//               if (GetGeometry().WorkingSpaceDimension() == 3){
//                   ElementalDofList.push_back(GetGeometry()[i].pGetDof(ANGULAR_VELOCITY_Z));
//               }
// 
//           }

          KRATOS_CATCH("")
      }


      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void Cluster3D::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
      {
          KRATOS_TRY


          KRATOS_CATCH("")

      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void Cluster3D::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
      {
        

      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void Cluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
      {
          
        KRATOS_TRY

        

          KRATOS_CATCH("")

      }// Calculate

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void Cluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo)
      {

        
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void Cluster3D::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
      void Cluster3D::Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}

                                              

}  // namespace Kratos.

