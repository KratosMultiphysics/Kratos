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

#ifndef KRATOS_KUTTA_CONDITION_PROCESS_H
#define KRATOS_KUTTA_CONDITION_PROCESS_H


#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "processes/process.h"
#include "geometries/geometry.h"
#include "compressible_potential_flow_application_variables.h"
#include "utilities/geometry_utilities.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"

#include <string>
#include <iostream>
#include <sstream>

#include <boost/functional/hash.hpp> //TODO: remove this dependence when Kratos has en internal one
#include <unordered_map> //TODO: remove this dependence when Kratos has en internal one
#include <utility>

namespace Kratos
{

class KuttaConditionProcess: public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(KuttaConditionProcess);

    typedef ModelPart::ElementType ElementType;
    typedef ModelPart::ConditionType ConditionType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for KuttaConditionProcess Process
//     KuttaConditionProcess(ModelPart& rModelPart,
//                      KratosParameters& parameters
//                     ):
//         Process(),
//         mrModelPart(rModelPart),
//         mrOptions(Flags()),
//         mrParameters(parameters)
//     {
//     }
    /// Constructor for KuttaConditionProcess Process
    KuttaConditionProcess(ModelPart& rModelPart
                    ):
        Process(),
        mrModelPart(rModelPart)
    {
    }

    /// Destructor.
    ~KuttaConditionProcess() override {}


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

        if(mrModelPart.Elements().size() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "the number of elements in the domain is zero. kutta condition can not be applied","")

        using normal_map_type = std::unordered_map<std::size_t, array_1d<double,3> >;
        
        for(auto it=mrModelPart.NodesBegin(); it!=mrModelPart.NodesEnd(); ++it)
            if(it->Is(STRUCTURE)) it->Set(VISITED,false);

        PointerVector<Element> kutta_elements;
        for(auto it=mrModelPart.ElementsBegin(); it!=mrModelPart.ElementsEnd(); ++it)
        {
            auto geom = it->GetGeometry();
            
            for(unsigned int i=0; i<geom.size(); ++i)
            {
                if(geom[i].Is(STRUCTURE))
                {
                    kutta_elements.push_back(*(it.base()));
                    break;
                }
            }
        }
        
        //compute one normal per each kutta node
        normal_map_type normal_map;
        for(auto it=kutta_elements.begin(); it!=kutta_elements.end(); ++it)
        {
            if(it->Is(MARKER))
            {
                //compute normal
                const Vector& elemental_distances = it->GetValue(ELEMENTAL_DISTANCES);
                
                double vol;
                array_1d<double,4> N;
                bounded_matrix<double,4,3> DN_DX;
                GeometryUtils::CalculateGeometryData(it->GetGeometry(), DN_DX, N, vol);
                
                //get a relevant norm
                array_1d<double,3> n = prod(trans(DN_DX),elemental_distances);
                n/=(norm_2(n)+1e-30);

                //
                auto geom = it->GetGeometry();
                for(unsigned int i=0; i<geom.size(); ++i)
                {
                    if(geom[i].Is(STRUCTURE))
                    {
                        geom[i].Set(VISITED,true);
                        normal_map[geom[i].Id()] = n;;
                    }
                }
            }
        }
        
        //now compute one elemntal distance per each element in the kutta list
        for(auto it=kutta_elements.begin(); it!=kutta_elements.end(); ++it)
        {
            it->Set(MARKER,false);
            Vector& elemental_distances = it->GetValue(ELEMENTAL_DISTANCES);
            
            array_1d<double,3> n;
            auto geom = it->GetGeometry();
            int kutta_id = -1;
            for(unsigned int i=0; i<geom.size(); ++i)
            {
                if(geom[i].Is(STRUCTURE) && geom[i].Is(VISITED))
                {
                    kutta_id = i;
                    n = normal_map[geom[i].Id()];
                    break;
                }
            }
            
            if(kutta_id == -1) //no node was found with the normal correctly calculated
            {
                it->Set(ACTIVE,false);
            }
            else
            {
                   
                
                array_1d<double,3> x0 = geom[kutta_id].Coordinates();
                double dkutta = -1e-3;
                for(int i=0; i<static_cast<int>(geom.size()); ++i)
                {
                    if(kutta_id != i)
                    {
                        
                        array_1d<double,3> v = geom[i].Coordinates() - x0;
                        elemental_distances[i] = inner_prod(n,v);
                    }
                    else
                        elemental_distances[i] = dkutta;
                }
                
                
                unsigned int npos = 0;
                unsigned int nneg = 0;
                for(unsigned int i=0; i<geom.size(); ++i)
                    if(elemental_distances[i] >= 0)
                        npos++;
                    else
                        nneg++;
                    
                if(nneg > 0 && npos >0)
                    it->Set(MARKER,true);
            }
            
        }       
        
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
        return "KuttaConditionProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "KuttaConditionProcess";
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
    Flags mrOptions;


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KuttaConditionProcess& operator=(KuttaConditionProcess const& rOther);

    /// Copy constructor.
    KuttaConditionProcess(KuttaConditionProcess const& rOther);


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
                                  KuttaConditionProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const KuttaConditionProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}



} // namespace Kratos


#endif // KRATOS_KUTTA_CONDITION_PROCESS_H
