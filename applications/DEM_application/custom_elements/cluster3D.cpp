//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Salva $
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
#include "includes/variables.h"

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

      
    void Cluster3D::Initialize() {
    
        KRATOS_TRY

        array_1d<double, 3> velocity                  = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
        array_1d<double, 3> totalforces               = GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES);
        double sqrt_of_mass                           = GetGeometry()[0].FastGetSolutionStepValue(SQRT_OF_MASS);
        mSqrtOfRealMass                               = sqrt_of_mass;
        array_1d<double, 3> angularvelocity           = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
        array_1d<double, 3> particlemoment            = GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT);
        array_1d<double, 3> principalmomentsofinertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        array_1d<double, 3> eulerangles               = GetGeometry()[0].FastGetSolutionStepValue(EULER_ANGLES);
                
        KRATOS_CATCH("")
    
    }
        
      
    void Cluster3D::CreateParticles(ParticleCreatorDestructor::Pointer p_creator_destructor, ModelPart& dem_model_part) {
          
        KRATOS_TRY 
          
        std::string ElementNameString = "SphericParticle3D";
            
        //We now create a spheric particle and keep it as a reference to an Element
        const Element& r_reference_element = KratosComponents<Element>::Get(ElementNameString);
        
        Node<3>& lele = GetGeometry()[0]; //NODE
        
        array_1d<double, 3> pepito;
        
        double r, s, t, R;
        
        for (int i = 1; i <= 250; i++) {
                    
            r = ((double) rand() / (RAND_MAX)) + 1;
            s = ((double) rand() / (RAND_MAX)) + 1;
            t = r - s;
            pepito[0]  = lele.Coordinates()[0] + t;
            
            r = ((double) rand() / (RAND_MAX)) + 1;
            s = ((double) rand() / (RAND_MAX)) + 1;
            t = r - s;
            pepito[1]  = lele.Coordinates()[1] + t;
            
            r = ((double) rand() / (RAND_MAX)) + 1;
            s = ((double) rand() / (RAND_MAX)) + 1;
            t = r - s;
            pepito[2]  = lele.Coordinates()[2] + t;
            
            R = 0.5 * (double) rand() / (RAND_MAX);
             
            p_creator_destructor->ElementCreatorForClusters(dem_model_part,
            i, R, pepito, 100.0, this->pGetProperties(), r_reference_element);
        }
                
        KRATOS_CATCH("")
        
    }
      
    void Cluster3D::UpdatePositionOfSpheres(double RotationMatrix[3][3]) {}
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
    
    double Cluster3D::GetSqrtOfRealMass()                                                       { return mSqrtOfRealMass; }
    //mListOfCoordinates         = GetClusterCoordinates();
    //mListOfRadii               = GetClusterRadii();
    

}  // namespace Kratos.

