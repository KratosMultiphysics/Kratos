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
#include "cluster3D.h"
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
        
        mListOfCoordinates[ 0][0] = 0.109377;mListOfCoordinates[ 0][1] =-0.009710;mListOfCoordinates[ 0][2] =-0.281412;
        mListOfCoordinates[ 1][0] =-0.142170;mListOfCoordinates[ 1][1] =-0.054794;mListOfCoordinates[ 1][2] = 0.352726;
        mListOfCoordinates[ 2][0] = 0.184175;mListOfCoordinates[ 2][1] = 0.239158;mListOfCoordinates[ 2][2] =-0.031211;
        mListOfCoordinates[ 3][0] = 0.097582;mListOfCoordinates[ 3][1] = 0.104613;mListOfCoordinates[ 3][2] =-0.131214;
        mListOfCoordinates[ 4][0] =-0.230406;mListOfCoordinates[ 4][1] =-0.230793;mListOfCoordinates[ 4][2] = 0.226359;
        mListOfCoordinates[ 5][0] = 0.232129;mListOfCoordinates[ 5][1] = 0.217145;mListOfCoordinates[ 5][2] =-0.245093;
        mListOfCoordinates[ 6][0] = 0.046670;mListOfCoordinates[ 6][1] = 0.273874;mListOfCoordinates[ 6][2] =-0.054866;
        mListOfCoordinates[ 7][0] =-0.102738;mListOfCoordinates[ 7][1] =-0.131159;mListOfCoordinates[ 7][2] =-0.246952;
        mListOfCoordinates[ 8][0] =-0.056793;mListOfCoordinates[ 8][1] = 0.005332;mListOfCoordinates[ 8][2] = 0.371649;
        mListOfCoordinates[ 9][0] = 0.053185;mListOfCoordinates[ 9][1] =-0.049706;mListOfCoordinates[ 9][2] =-0.183108;
        mListOfCoordinates[10][0] =-0.325390;mListOfCoordinates[10][1] =-0.417720;mListOfCoordinates[10][2] =-0.263742;
        mListOfCoordinates[11][0] = 0.177439;mListOfCoordinates[11][1] =-0.185131;mListOfCoordinates[11][2] =-0.130760;
        mListOfCoordinates[12][0] =-0.216448;mListOfCoordinates[12][1] =-0.248704;mListOfCoordinates[12][2] =-0.185033;
        mListOfCoordinates[13][0] =-0.059379;mListOfCoordinates[13][1] =-0.244949;mListOfCoordinates[13][2] =-0.118607;
        mListOfCoordinates[14][0] = 0.277810;mListOfCoordinates[14][1] = 0.130129;mListOfCoordinates[14][2] =-0.021583;
        mListOfCoordinates[15][0] = 0.077362;mListOfCoordinates[15][1] =-0.041935;mListOfCoordinates[15][2] =-0.233118;
        mListOfCoordinates[16][0] = 0.131849;mListOfCoordinates[16][1] = 0.041977;mListOfCoordinates[16][2] =-0.257780;
        mListOfCoordinates[17][0] =-0.079260;mListOfCoordinates[17][1] =-0.050416;mListOfCoordinates[17][2] = 0.005213;
        mListOfCoordinates[18][0] =-0.044223;mListOfCoordinates[18][1] =-0.001916;mListOfCoordinates[18][2] =-0.024154;
        mListOfCoordinates[19][0] = 0.058239;mListOfCoordinates[19][1] = 0.145923;mListOfCoordinates[19][2] =-0.076332;

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
        
        base_principal_moments_of_inertia[0] = cluster_mass * cl * cl * 0.167992017;
        base_principal_moments_of_inertia[1] = cluster_mass * cl * cl * 0.164636030;
        base_principal_moments_of_inertia[2] = cluster_mass * cl * cl * 0.141015686; 
  
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

