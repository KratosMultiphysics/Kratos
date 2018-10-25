//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#ifndef KRATOS_COMPUTE_LIFT_LEVEL_SET_PROCESS_H
#define KRATOS_COMPUTE_LIFT_LEVEL_SET_PROCESS_H


#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "processes/process.h"
#include "geometries/geometry.h"
#include "utilities/geometry_utilities.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"
#include "utilities/divide_triangle_2d_3.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"

#include <string>
#include <iostream>
#include <sstream>

#include <boost/functional/hash.hpp> //TODO: remove this dependence when Kratos has en internal one
#include <unordered_map> //TODO: remove this dependence when Kratos has en internal one
#include <utility>

namespace Kratos
{

class ComputeLiftLevelSetProcess: public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(ComputeLiftLevelSetProcess);

    typedef ModelPart::ElementType ElementType;
    typedef ModelPart::ConditionType ConditionType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for ComputeLiftLevelSetProcess Process
//     ComputeLiftLevelSetProcess(ModelPart& rModelPart,
//                      KratosParameters& parameters
//                     ):
//         Process(),
//         mrModelPart(rModelPart),
//         mrOptions(Flags()),
//         mrParameters(parameters)
//     {
//     }
    /// Constructor for ComputeLiftLevelSetProcess Process
    ComputeLiftLevelSetProcess(ModelPart& rModelPart,
                Vector& rResultForce                             
                    ):
        Process(),
        mrModelPart(rModelPart),
        mrResultForce(rResultForce)
    {
    }

    /// Destructor.
    ~ComputeLiftLevelSetProcess() override {}


    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{


    /// Check elements to make sure that their jacobian is positive and conditions to ensure that their face normals point outwards
    void Execute() override
    {
        KRATOS_TRY;
        std::cout<<"Compute Lift Level Set Process Start"<< std::endl;
        double Cl,Cd,Rz;

        for(auto it=mrModelPart.ElementsBegin(); it!=mrModelPart.ElementsEnd(); ++it)
        {
            if (it->Is(BOUNDARY)){     
                auto geom = it->GetGeometry();
                const unsigned int NumNodes=geom.size();     
                array_1d<double,3> elemental_distance;
                Geometry< Node<3> >::Pointer pgeom = it-> pGetGeometry();
    
                for(unsigned int i_node = 0; i_node<NumNodes; i_node++)
                    elemental_distance[i_node] = geom[i_node].GetSolutionStepValue(LEVEL_SET_DISTANCE);               
                
                const Vector& r_elemental_distances=elemental_distance;
                
                Triangle2D3ModifiedShapeFunctions triangle_shape_functions(pgeom, r_elemental_distances);

                
                array_1d<double,3>  gp_wall = ZeroVector(3);
                Matrix positive_side_interface_sh_func;
                ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_interface_gradients;
                Vector positive_side_interface_weights;
                triangle_shape_functions.ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
                    positive_side_interface_sh_func,
                    positive_side_sh_func_interface_gradients,
                    positive_side_interface_weights,
                    GeometryData::GI_GAUSS_1);


                std::vector<Vector> cut_unit_normal;
                std::vector<Vector> gp_cut;
                triangle_shape_functions.ComputePositiveSideInterfaceAreaNormals(cut_unit_normal,GeometryData::GI_GAUSS_1);

                const unsigned int n_gauss = positive_side_interface_weights.size();
            
                std::vector<double> cp;
                
       
                for (unsigned int i_gauss=0;i_gauss<n_gauss;i_gauss++){
                    for (unsigned int i_node=0;i_node<NumNodes;i_node++){
                        array_1d<double,3> coord=geom[i_node].Coordinates();
                        for (unsigned int i_dim=0;i_dim<3;i_dim++){                    
                            gp_wall[i_dim] += coord[i_dim]*positive_side_interface_sh_func(i_gauss,i_node);//*positive_side_interface_weights(i_gauss);
                        }
                    }
                }
                // geom.PrintData(std::cout);
                // std::cout<<"gp_wall0: "<<gp_wall[0]<<std::endl;
                // std::cout<<"gp_wall1: "<<gp_wall[1]<<std::endl;
                // std::cout<<"gp_wall2: "<<gp_wall[2]<<std::endl;
                // std::cout<<std::endl;

                it->GetValueOnIntegrationPoints(PRESSURE,cp,mrModelPart.GetProcessInfo());

                double cpressure=cp[0];
                               
                Cl += cpressure*cut_unit_normal[0][1];
                Cd += cpressure*cut_unit_normal[0][0];
                Rz += cpressure*cut_unit_normal[0][2];

          
                it->SetValue(PRESSURE,cp[0]);
                it->SetValue(NORMAL,cut_unit_normal[0]);
                it->SetValue(BODY_FORCE,gp_wall); //provisional name for cp application point
                
                
            }
        }
        mrResultForce[0]=Cd;
        mrResultForce[1]=Cl;
        mrResultForce[2]=Rz;
                 
        KRATOS_CATCH("");
    }



    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ComputeLiftLevelSetProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ComputeLiftLevelSetProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        this->PrintInfo(rOStream);
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}


private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    Vector& mrResultForce;
    Flags mrOptions;


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ComputeLiftLevelSetProcess& operator=(ComputeLiftLevelSetProcess const& rOther);

    /// Copy constructor.
    ComputeLiftLevelSetProcess(ComputeLiftLevelSetProcess const& rOther);


    ///@}

}; // Class Process

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ComputeLiftLevelSetProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ComputeLiftLevelSetProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}



} // namespace Kratos


#endif // KRATOS_ComputeLift_LEVEL_SET_PROCESS_H
