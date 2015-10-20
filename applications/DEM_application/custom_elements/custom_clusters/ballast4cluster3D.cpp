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

// External includes

// Project includes
#include "includes/define.h"
#include "ballast4cluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    Ballast4Cluster3D::Ballast4Cluster3D() : Cluster3D() {}
            
      
    Ballast4Cluster3D::Ballast4Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
      
      
    Ballast4Cluster3D::Ballast4Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}

      
    Ballast4Cluster3D::Ballast4Cluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}

      
    Element::Pointer Ballast4Cluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Cluster3D::Pointer(new Ballast4Cluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    Ballast4Cluster3D::~Ballast4Cluster3D() {}
      
    
    void Ballast4Cluster3D::CustomInitialize() {
        
        int number_of_spheres = 20;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        // 1.183 (in meters) was the medium diameter of the rock in GiD (rock05.gid file)
        // so we should multiply every size that follow by the inverse of that number, 0.8453,
        // to obtain a 'unity' rock.
        // We then have to multiply again everything by 'cl' to obtain the desired dimensions
        // to adjust to the characteristic length given
                
        cl *= 0.8453;
                
//         double a, b, c;
//         
//         a = 0.572749; // this is the original semi-axis (1.1455/2)
//         b = a; c = b;
        
        mListOfCoordinates[ 0][0]= 0.231256; mListOfCoordinates[ 0][1]=-0.074296; mListOfCoordinates[ 0][2]=-0.068691;
        mListOfCoordinates[ 1][0]=-0.419411; mListOfCoordinates[ 1][1]= 0.135353; mListOfCoordinates[ 1][2]=-0.408928;
        mListOfCoordinates[ 2][0]= 0.146031; mListOfCoordinates[ 2][1]= 0.275717; mListOfCoordinates[ 2][2]=-0.179708;
        mListOfCoordinates[ 3][0]= 0.139874; mListOfCoordinates[ 3][1]= 0.087744; mListOfCoordinates[ 3][2]=-0.112625;
        mListOfCoordinates[ 4][0]=-0.417929; mListOfCoordinates[ 4][1]=-0.098586; mListOfCoordinates[ 4][2]=-0.343806;
        mListOfCoordinates[ 5][0]= 0.334970; mListOfCoordinates[ 5][1]= 0.162886; mListOfCoordinates[ 5][2]=-0.058594;
        mListOfCoordinates[ 6][0]= 0.103533; mListOfCoordinates[ 6][1]= 0.236640; mListOfCoordinates[ 6][2]=-0.035744;
        mListOfCoordinates[ 7][0]= 0.056260; mListOfCoordinates[ 7][1]=-0.228340; mListOfCoordinates[ 7][2]=-0.005488;
        mListOfCoordinates[ 8][0]=-0.371725; mListOfCoordinates[ 8][1]= 0.223389; mListOfCoordinates[ 8][2]=-0.453942;
        mListOfCoordinates[ 9][0]= 0.111432; mListOfCoordinates[ 9][1]=-0.073221; mListOfCoordinates[ 9][2]=-0.127622;
        mListOfCoordinates[10][0]=-0.133260; mListOfCoordinates[10][1]=-0.538192; mListOfCoordinates[10][2]= 0.010085;
        mListOfCoordinates[11][0]= 0.093703; mListOfCoordinates[11][1]=-0.099668; mListOfCoordinates[11][2]=-0.343238;
        mListOfCoordinates[12][0]=-0.088034; mListOfCoordinates[12][1]=-0.327130; mListOfCoordinates[12][2]=-0.036903;
        mListOfCoordinates[13][0]=-0.057795; mListOfCoordinates[13][1]=-0.228695; mListOfCoordinates[13][2]=-0.207391;
        mListOfCoordinates[14][0]= 0.154340; mListOfCoordinates[14][1]= 0.234559; mListOfCoordinates[14][2]=-0.322492;
        mListOfCoordinates[15][0]= 0.166286; mListOfCoordinates[15][1]=-0.084923; mListOfCoordinates[15][2]=-0.101633;
        mListOfCoordinates[16][0]= 0.239770; mListOfCoordinates[16][1]=-0.014129; mListOfCoordinates[16][2]=-0.074446;
        mListOfCoordinates[17][0]=-0.108138; mListOfCoordinates[17][1]=-0.023480; mListOfCoordinates[17][2]=-0.177239;
        mListOfCoordinates[18][0]=-0.051764; mListOfCoordinates[18][1]= 0.010696; mListOfCoordinates[18][2]=-0.152177;
        mListOfCoordinates[19][0]= 0.088039; mListOfCoordinates[19][1]= 0.133217; mListOfCoordinates[19][2]=-0.102246;


        for (int i = 0; i < number_of_spheres; i++) { 
            mListOfCoordinates[i][0] *= cl;  mListOfCoordinates[i][1] *= cl;  mListOfCoordinates[i][2] *= cl;
        }
        
        mListOfRadii[ 0] = 0.444471;
        mListOfRadii[ 1] = 0.370286;
        mListOfRadii[ 2] = 0.363593;
        mListOfRadii[ 3] = 0.445970;
        mListOfRadii[ 4] = 0.243500;
        mListOfRadii[ 5] = 0.359937;
        mListOfRadii[ 6] = 0.345782;
        mListOfRadii[ 7] = 0.348147;
        mListOfRadii[ 8] = 0.280183;
        mListOfRadii[ 9] = 0.426531;
        mListOfRadii[10] = 0.095537;
        mListOfRadii[11] = 0.346096;
        mListOfRadii[12] = 0.369412;
        mListOfRadii[13] = 0.378931;
        mListOfRadii[14] = 0.326994;
        mListOfRadii[15] = 0.418579;
        mListOfRadii[16] = 0.380625;
        mListOfRadii[17] = 0.428503;
        mListOfRadii[18] = 0.408325;
        mListOfRadii[19] = 0.409263;
        
        for (int i = 0; i < number_of_spheres; i++) { mListOfRadii[i]= mListOfRadii[i] * cl; }
       
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * 0.910205 * cl * cl * cl;
                
        array_1d<double, 3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = cluster_mass * cl * cl * 0.128633293;
        base_principal_moments_of_inertia[1] = cluster_mass * cl * cl * 0.155345640;
        base_principal_moments_of_inertia[2] = cluster_mass * cl * cl * 0.189666335; 
  
    }     
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void Ballast4Cluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
    void Ballast4Cluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}
    double Ballast4Cluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

