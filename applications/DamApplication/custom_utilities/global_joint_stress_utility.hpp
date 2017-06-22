//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Lorenzo Gracia
//
//

#if !defined(KRATOS_GLOBAL_JOINT_STRESS_UTILITIES )
#define  KRATOS_GLOBAL_JOINT_STRESS_UTILITIES

#include <cmath>

// Project includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "utilities/openmp_utils.h"
#include "utilities/math_utils.h"
#include "includes/element.h"

// Application includes
#include "dam_application_variables.h"

// How to use it?
// The aim of this utility is to obtain the tangential and normal forces at some plane of interest in the joint element.
// It can be called from MainKratos.py following next commmands.
//
// Example of use from python:
// import global_joint_stress_utility
// global_joint_stress_utility.GlobalJoinStresstUtility(main_model_part.GetSubModelPart("Parts_Parts_Auto3")).ComputingGlobalStress()
//
// The input plane must be introduce in the global_joint_stress_utility.py

namespace Kratos
{
    
class GlobalJointStressUtility
{

public:
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    GlobalJointStressUtility(ModelPart& model_part,
                            Parameters rParameters
                            ) : mr_model_part(model_part)
    {
        KRATOS_TRY

        // Getting values of the interested plane
        mpmid0_0 = rParameters["pmid0_0"].GetDouble();
        mpmid0_1 = rParameters["pmid0_1"].GetDouble();
        mpmid0_2 = rParameters["pmid0_2"].GetDouble();
        mpmid1_0 = rParameters["pmid1_0"].GetDouble();
        mpmid1_1 = rParameters["pmid1_1"].GetDouble();
        mpmid1_2 = rParameters["pmid1_2"].GetDouble();
        mpmid2_0 = rParameters["pmid2_0"].GetDouble();
        mpmid2_1 = rParameters["pmid2_1"].GetDouble();
        mpmid2_2 = rParameters["pmid2_2"].GetDouble();

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------
    
    /// Destructor
    ~GlobalJointStressUtility() {}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Transforming local stress vector in global coordinates 
    void ComputingGlobalStress()
    {
        const int nelements = static_cast<int>(mr_model_part.Elements().size());
        ModelPart::ElementsContainerType::iterator el_begin = mr_model_part.ElementsBegin();
        GeometryData::IntegrationMethod MyIntegrationMethod;
        const ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
        boost::numeric::ublas::bounded_matrix<double,3,3> RotationPlane;

        //Definition of the plane
        array_1d<double, 3> pmid0;
        array_1d<double, 3> pmid1;
        array_1d<double, 3> pmid2;

        // Coordinates of interest plane
        pmid0[0]= mpmid0_0;
        pmid0[1]= mpmid0_1;
        pmid0[2]= mpmid0_2;

        pmid1[0]= mpmid1_0;
        pmid1[1]= mpmid1_1;
        pmid1[2]= mpmid1_2;

        pmid2[0]= mpmid2_0;
        pmid2[1]= mpmid2_1;
        pmid2[2]= mpmid2_2;

        this->CalculateRotationMatrix(RotationPlane,pmid0,pmid1,pmid2);
        array_1d<double,3> VectorForceinPlane;
        array_1d<double,3> GlobalElementVectorForce = ZeroVector(3);

        for(int k = 0; k<nelements; k++)
        {
            ModelPart::ElementsContainerType::iterator it = el_begin + k;
            Element::GeometryType& Geom = it->GetGeometry();
            boost::numeric::ublas::bounded_matrix<double,3,3> RotationMatrix;

            //Define mid-plane points for prism_interface_3d_6
            noalias(pmid0) = 0.5 * (Geom.GetPoint( 0 ) + Geom.GetPoint( 3 ));
            noalias(pmid1) = 0.5 * (Geom.GetPoint( 1 ) + Geom.GetPoint( 4 ));
            noalias(pmid2) = 0.5 * (Geom.GetPoint( 2 ) + Geom.GetPoint( 5 ));

            this->CalculateRotationMatrix(RotationMatrix,pmid0,pmid1,pmid2);
            MyIntegrationMethod = it->GetIntegrationMethod();
            const Element::GeometryType::IntegrationPointsArrayType& IntegrationPoints = Geom.IntegrationPoints(MyIntegrationMethod);
            unsigned int NumGPoints = IntegrationPoints.size();
            std::vector<array_1d<double,3>> LocalStressVector;
            array_1d<double,3> LocalElementStress = ZeroVector(3);
            array_1d<double,3> LocalElementVectorForce;
            it->GetValueOnIntegrationPoints(LOCAL_STRESS_VECTOR,LocalStressVector,CurrentProcessInfo);
            
            for(unsigned int GPoint=0; GPoint<NumGPoints; GPoint++)
            {
                noalias(LocalElementStress) += LocalStressVector[GPoint];
            }

            // Computing area at mid plane
            double InvNumGP = 1.0/static_cast<double>(NumGPoints);
            LocalElementStress[0] *= InvNumGP;
            LocalElementStress[1] *= InvNumGP;
            LocalElementStress[2] *= InvNumGP;
            double Area;
            this->AreaMidPlane(Area,pmid0,pmid1,pmid2);
            noalias(LocalElementVectorForce) = LocalElementStress*Area;
            noalias(GlobalElementVectorForce) += prod(trans(RotationMatrix),LocalElementVectorForce);
        }

        noalias(VectorForceinPlane) = prod(RotationPlane,GlobalElementVectorForce);
        double TangentialForce = sqrt(VectorForceinPlane[0]*VectorForceinPlane[0] + VectorForceinPlane[1]*VectorForceinPlane[1]);
        
        std::cout<< " Tangential Force (N) "<<TangentialForce<<std::endl;
        std::cout<< " Normal Force (N) "<<fabs(VectorForceinPlane[2])<<std::endl;
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

        /// Member Variables
        ModelPart& mr_model_part;
        double mpmid0_0;
        double mpmid0_1;
        double mpmid0_2;
        double mpmid1_0;
        double mpmid1_1;
        double mpmid1_2;
        double mpmid2_0;
        double mpmid2_1;
        double mpmid2_2;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    //Currently is just working for prism element

   void CalculateRotationMatrix(boost::numeric::ublas::bounded_matrix<double,3,3>& rRotationMatrix, array_1d<double, 3>& pmid0, array_1d<double, 3>& pmid1, array_1d<double, 3>& pmid2 )
    {
        KRATOS_TRY
        
        //Unitary vector in local x direction
        array_1d<double, 3> Vx;
        noalias(Vx) = pmid1 - pmid0;
        double inv_norm_x = 1.0/norm_2(Vx);
        Vx[0] *= inv_norm_x;
        Vx[1] *= inv_norm_x;
        Vx[2] *= inv_norm_x;
        
        //Unitary vector in local z direction
        array_1d<double, 3> Vy;
        noalias(Vy) = pmid2 - pmid0;
        array_1d<double, 3> Vz;
        MathUtils<double>::CrossProduct(Vz, Vx, Vy);
        double inv_norm_z = 1.0/norm_2(Vz);
        Vz[0] *= inv_norm_z;
        Vz[1] *= inv_norm_z;
        Vz[2] *= inv_norm_z;
            
        //Unitary vector in local y direction
        MathUtils<double>::CrossProduct( Vy, Vz, Vx);
    
        //Rotation Matrix
        rRotationMatrix(0,0) = Vx[0];
        rRotationMatrix(0,1) = Vx[1];
        rRotationMatrix(0,2) = Vx[2];
    
        rRotationMatrix(1,0) = Vy[0];
        rRotationMatrix(1,1) = Vy[1];
        rRotationMatrix(1,2) = Vy[2];
    
        rRotationMatrix(2,0) = Vz[0];
        rRotationMatrix(2,1) = Vz[1];
        rRotationMatrix(2,2) = Vz[2];
    
        KRATOS_CATCH( "" )
        
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    void AreaMidPlane(double& rArea, array_1d<double, 3>& pmid0, array_1d<double, 3>& pmid1, array_1d<double, 3>& pmid2 )
    {
        KRATOS_TRY
        
        //Vector declarations
        array_1d<double, 3> Vx;
        array_1d<double, 3> Vy;
        array_1d<double, 3> Vz;

        // Computing distances and area
        noalias(Vx) = pmid1 - pmid0;        
        noalias(Vy) = pmid2 - pmid0;
        MathUtils<double>::CrossProduct(Vz, Vx, Vy);
        rArea = sqrt(Vz[0]*Vz[0]+ Vz[1]*Vz[1]+Vz[2]*Vz[2])/2.0;

        KRATOS_CATCH( "" )
        
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

};//Class

} /* namespace Kratos.*/

#endif /* KRATOS_GLOBAL_JOINT_STRESS_UTILITIES defined */