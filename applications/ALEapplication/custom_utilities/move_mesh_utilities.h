//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Andreas Winterstein (a.winterstein@tum.de)
//


#if !defined(KRATOS_MESHMOVING_UTILITIES_H_INCLUDED )
#define  KRATOS_MESHMOVING_UTILITIES_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
#include "ale_application.h"
#include "ale_application.h"

namespace Kratos
{

class MoveMeshUtilities
{
public:

    MoveMeshUtilities() {};
    ~MoveMeshUtilities() {};

    void BDF_MoveMesh(unsigned int time_order, ModelPart& rModelPart)
    {
        KRATOS_TRY
        //calculating time integration coefficients
        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        double dt = CurrentProcessInfo[DELTA_TIME];

        Vector BDFcoeffs(time_order+1);

        if(time_order == 2)
        {
            if(rModelPart.GetBufferSize() < 3)
                KRATOS_THROW_ERROR(std::logic_error,"insufficient buffer size for BDF2","")

                BDFcoeffs[0] =	1.5 / dt;	//coefficient for step n+1
            BDFcoeffs[1] =	-2.0 / dt;//coefficient for step n
            BDFcoeffs[2] =	0.5 / dt;//coefficient for step n-1
        }
        else
        {
            BDFcoeffs[0] =	1.0 / dt;	//coefficient for step n+1
            BDFcoeffs[1] =	-1.0 / dt;//coefficient for step n
        }

        //update nodal coordinates
        array_1d<double,3> mesh_vel;
        for(ModelPart::NodesContainerType::iterator i = rModelPart.NodesBegin();
                i!=rModelPart.NodesEnd(); i++)
        {
            const array_1d<double,3>& disp = i->FastGetSolutionStepValue(DISPLACEMENT);
            i->X() = i->X0() + disp[0];
            i->Y() = i->Y0() + disp[1];
            i->Z() = i->Z0() + disp[2];

            //calculating the mesh velocity
            noalias(mesh_vel) = BDFcoeffs[0] * disp;
            for(unsigned int step=1; step<time_order+1; step++)
                noalias(mesh_vel) += BDFcoeffs[step]*i->FastGetSolutionStepValue(DISPLACEMENT,step);

            //saving the mesh velocity
            noalias(i->FastGetSolutionStepValue(MESH_VELOCITY)) = mesh_vel;
        }


        KRATOS_CATCH("");
    }






private:

};



}  // namespace Kratos.

#endif // KRATOS_MESHMOVING_UTILITIES_H_INCLUDED  defined
