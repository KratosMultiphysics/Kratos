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

#if !defined(KRATOS_GLOBAL_JOINT_STRESS_UTILITIES)
#define KRATOS_GLOBAL_JOINT_STRESS_UTILITIES

// System includes
#include <cmath>

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
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
// global_joint_stress_utility.GlobalJoinStressUtility(main_model_part.GetSubModelPart("Parts_Parts_Auto3")).ComputingGlobalStress()
//
// The input plane must be introduce in the global_joint_stress_utility.py

namespace Kratos
{

class GlobalJointStressUtility
{

  public:
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    GlobalJointStressUtility(ModelPart &model_part,
                             Parameters rParameters) : mr_model_part(model_part)
    {
        KRATOS_TRY

        // Getting values of the interested plane
        m_plane_coordinates[0] = rParameters["pmid0_0"].GetDouble();
        m_plane_coordinates[1] = rParameters["pmid0_1"].GetDouble();
        m_plane_coordinates[2] = rParameters["pmid0_2"].GetDouble();
        m_plane_coordinates[3] = rParameters["pmid1_0"].GetDouble();
        m_plane_coordinates[4] = rParameters["pmid1_1"].GetDouble();
        m_plane_coordinates[5] = rParameters["pmid1_2"].GetDouble();
        m_plane_coordinates[6] = rParameters["pmid2_0"].GetDouble();
        m_plane_coordinates[7] = rParameters["pmid2_1"].GetDouble();
        m_plane_coordinates[8] = rParameters["pmid2_2"].GetDouble();

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
        const ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
        BoundedMatrix<double, 3, 3> RotationPlane;

        //Definition of the plane
        array_1d<double, 3> point_mid0;
        array_1d<double, 3> point_mid1;
        array_1d<double, 3> point_mid2;

        // Coordinates of interest plane
        point_mid0[0] = m_plane_coordinates[0];
        point_mid0[1] = m_plane_coordinates[1];
        point_mid0[2] = m_plane_coordinates[2];

        point_mid1[0] = m_plane_coordinates[3];
        point_mid1[1] = m_plane_coordinates[4];
        point_mid1[2] = m_plane_coordinates[5];

        point_mid2[0] = m_plane_coordinates[6];
        point_mid2[1] = m_plane_coordinates[7];
        point_mid2[2] = m_plane_coordinates[8];

        this->CalculateRotationMatrix(RotationPlane, point_mid0, point_mid1, point_mid2);
        array_1d<double, 3> VectorForceinPlane;
        array_1d<double, 3> GlobalElementVectorForce = ZeroVector(3);

        for (int k = 0; k < nelements; k++)
        {
            ModelPart::ElementsContainerType::iterator it = el_begin + k;
            Element::GeometryType &rGeom = it->GetGeometry();
            BoundedMatrix<double, 3, 3> RotationMatrix;

            //Define mid-plane points for prism_interface_3d_6
            noalias(point_mid0) = 0.5 * (rGeom.GetPoint(0) + rGeom.GetPoint(3));
            noalias(point_mid1) = 0.5 * (rGeom.GetPoint(1) + rGeom.GetPoint(4));
            noalias(point_mid2) = 0.5 * (rGeom.GetPoint(2) + rGeom.GetPoint(5));

            this->CalculateRotationMatrix(RotationMatrix, point_mid0, point_mid1, point_mid2);
            MyIntegrationMethod = it->GetIntegrationMethod();
            const Element::GeometryType::IntegrationPointsArrayType &IntegrationPoints = rGeom.IntegrationPoints(MyIntegrationMethod);
            unsigned int NumGPoints = IntegrationPoints.size();
            std::vector<array_1d<double, 3>> LocalStressVector;
            array_1d<double, 3> LocalElementStress = ZeroVector(3);
            array_1d<double, 3> LocalElementVectorForce;
            it->GetValueOnIntegrationPoints(LOCAL_STRESS_VECTOR, LocalStressVector, CurrentProcessInfo);

            for (unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
            {
                noalias(LocalElementStress) += LocalStressVector[GPoint];
            }

            // Computing area at mid plane
            double InvNumGP = 1.0 / static_cast<double>(NumGPoints);
            LocalElementStress[0] *= InvNumGP;
            LocalElementStress[1] *= InvNumGP;
            LocalElementStress[2] *= InvNumGP;
            double Area;
            this->AreaMidPlane(Area, point_mid0, point_mid1, point_mid2);
            noalias(LocalElementVectorForce) = LocalElementStress * Area;
            noalias(GlobalElementVectorForce) += prod(trans(RotationMatrix), LocalElementVectorForce);
        }

        noalias(VectorForceinPlane) = prod(RotationPlane, GlobalElementVectorForce);
        const double TangentialForce = sqrt(VectorForceinPlane[0] * VectorForceinPlane[0] + VectorForceinPlane[1] * VectorForceinPlane[1]);

        std::cout << " Tangential Force (N) " << TangentialForce << std::endl;
        std::cout << " Normal Force (N) " << fabs(VectorForceinPlane[2]) << std::endl;
    }

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  protected:
    /// Member Variables
    ModelPart &mr_model_part;
    array_1d<double, 9> m_plane_coordinates;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    //Currently is just working for prism element

    void CalculateRotationMatrix(BoundedMatrix<double, 3, 3> &rRotationMatrix, array_1d<double, 3> &point_mid0, array_1d<double, 3> &point_mid1, array_1d<double, 3> &point_mid2)
    {
        KRATOS_TRY

        //Unitary vector in local x direction
        array_1d<double, 3> Vx;
        noalias(Vx) = point_mid1 - point_mid0;
        double inv_norm_x = 1.0 / norm_2(Vx);
        Vx[0] *= inv_norm_x;
        Vx[1] *= inv_norm_x;
        Vx[2] *= inv_norm_x;

        //Unitary vector in local z direction
        array_1d<double, 3> Vy;
        noalias(Vy) = point_mid2 - point_mid0;
        array_1d<double, 3> Vz;
        MathUtils<double>::CrossProduct(Vz, Vx, Vy);
        double inv_norm_z = 1.0 / norm_2(Vz);
        Vz[0] *= inv_norm_z;
        Vz[1] *= inv_norm_z;
        Vz[2] *= inv_norm_z;

        //Unitary vector in local y direction
        MathUtils<double>::CrossProduct(Vy, Vz, Vx);

        //Rotation Matrix
        rRotationMatrix(0, 0) = Vx[0];
        rRotationMatrix(0, 1) = Vx[1];
        rRotationMatrix(0, 2) = Vx[2];

        rRotationMatrix(1, 0) = Vy[0];
        rRotationMatrix(1, 1) = Vy[1];
        rRotationMatrix(1, 2) = Vy[2];

        rRotationMatrix(2, 0) = Vz[0];
        rRotationMatrix(2, 1) = Vz[1];
        rRotationMatrix(2, 2) = Vz[2];

        KRATOS_CATCH("")
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AreaMidPlane(double &rArea, array_1d<double, 3> &point_mid0, array_1d<double, 3> &point_mid1, array_1d<double, 3> &point_mid2)
    {
        KRATOS_TRY

        //Vector declarations
        array_1d<double, 3> Vx;
        array_1d<double, 3> Vy;
        array_1d<double, 3> Vz;

        // Computing distances and area
        noalias(Vx) = point_mid1 - point_mid0;
        noalias(Vy) = point_mid2 - point_mid0;
        MathUtils<double>::CrossProduct(Vz, Vx, Vy);
        rArea = sqrt(Vz[0] * Vz[0] + Vz[1] * Vz[1] + Vz[2] * Vz[2]) / 2.0;

        KRATOS_CATCH("")
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

}; //Class

} /* namespace Kratos.*/

#endif /* KRATOS_GLOBAL_JOINT_STRESS_UTILITIES defined */