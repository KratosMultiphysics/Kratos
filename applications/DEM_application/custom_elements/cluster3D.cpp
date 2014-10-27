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
#include <stdlib.h>

// External includes

// Project includes
#include "includes/define.h"
#include "cluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"


namespace Kratos {
        
    // using namespace GeometryFunctions;

    Cluster3D::Cluster3D() : Element() {}
            
      
    Cluster3D::Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}
      
      
    Cluster3D::Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

      
    Cluster3D::Cluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Element(NewId, ThisNodes) {}

      
    Element::Pointer Cluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Element::Pointer(new Cluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    Cluster3D::~Cluster3D() {}

      
    void Cluster3D::Initialize() {}
        
      
    void Cluster3D::CreateParticles(ParticleCreatorDestructor::Pointer p_creator_destructor, ModelPart& dem_model_part) {
          
        KRATOS_TRY 
          
        KRATOS_WATCH("CREATING THE PARTICLES OF THE CLUSTER..........................")
                
        //unsigned int& max_Id = p_creator_destructor->mMaxNodeId; 
        
	//bool mBallsModelPartHasSphericity(false);        
        //bool mBallsModelPartHasRotation(false);
       
        std::string ElementNameString = "SphericParticle3D";
                
        const Element& r_reference_element = KratosComponents<Element>::Get(ElementNameString); //crea una spheric particle y la guarda como ref a un Element
        
        Node<3>& lele = GetGeometry()[0]; //NODE
    
        //PropertiesProxy* p_fast_properties = NULL;
       
        array_1d<double, 3> pepito;
        
        double r, s, t, R;
        
        for (int i = 1; i <= 250; i++) {
            
            //KRATOS_WATCH("RANDOM NUMBERS ARE: *************************************")
                    
            r = ((double) rand() / (RAND_MAX)) + 1;
            s = ((double) rand() / (RAND_MAX)) + 1;
            t = r - s;
            //KRATOS_WATCH(t)
            pepito[0] = lele.Coordinates()[0] + t;
            
            r = ((double) rand() / (RAND_MAX)) + 1;
            s = ((double) rand() / (RAND_MAX)) + 1;
            t = r - s;
            //KRATOS_WATCH(t)
            pepito[1] = lele.Coordinates()[1] + t;
            
            r = ((double) rand() / (RAND_MAX)) + 1;
            s = ((double) rand() / (RAND_MAX)) + 1;
            t = r - s;
            //KRATOS_WATCH(t)
            pepito[2] = lele.Coordinates()[2] + t;
            
            R = 0.5 * (double) rand() / (RAND_MAX);
            //KRATOS_WATCH(R)
                
            p_creator_destructor->ElementCreatorForClusters(dem_model_part,
            40000 + 3*i, R, pepito, 100.0, this->pGetProperties(), r_reference_element);
        }
                
        //KRATOS_WATCH(*this->pGetProperties())
        //KRATOS_WATCH(pepito)
        //KRATOS_WATCH(max_Id)
        //max_Id++;
        KRATOS_CATCH("")
        
    }
      

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Cluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Cluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Cluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Cluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Cluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) {
        
        KRATOS_TRY

        ElementalDofList.resize(0);

//      for (unsigned int i = 0; i < GetGeometry().size(); i++){
//          ElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_X));
//          ElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Y));
//          if (GetGeometry().WorkingSpaceDimension() == 3) {
//              ElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Z));
//          }
// 
//          ElementalDofList.push_back(GetGeometry()[i].pGetDof(ANGULAR_VELOCITY_X));
//          ElementalDofList.push_back(GetGeometry()[i].pGetDof(ANGULAR_VELOCITY_Y));
//          if (GetGeometry().WorkingSpaceDimension() == 3) {
//              ElementalDofList.push_back(GetGeometry()[i].pGetDof(ANGULAR_VELOCITY_Z));
//          }
// 
//      }

        KRATOS_CATCH("")
        
    }


    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Cluster3D::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Cluster3D::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void Cluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Cluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Cluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
    void Cluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}
    

}  // namespace Kratos.

