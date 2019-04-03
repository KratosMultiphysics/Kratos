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

#ifndef KRATOS_VERIFY_WATERTIGHTNESS_PROCESS_H
#define KRATOS_VERIFY_WATERTIGHTNESS_PROCESS_H


#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "processes/process.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"
#include "includes/key_hash.h"

#include <string>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <utility>

namespace Kratos
{

/** This function verifies that the skin has no holes nor overlapped geometries
 * this is accomplished by storing all of the edges in the model in a hash map and
 * verifying that no edge appears more than twice (which would imply an overlap)
 * or less than once (which would imply a gap in the skin)
 *
 * in the case such condition is violated an error is thrown
 */
class CheckSkinProcess: public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(CheckSkinProcess);

    typedef ModelPart::ElementType ElementType;
    typedef ModelPart::ConditionType ConditionType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for CheckSkinProcess Process
//     CheckSkinProcess(ModelPart& rModelPart,
//                      KratosParameters& parameters
//                     ):
//         Process(),
//         mrModelPart(rModelPart),
//         mrOptions(Flags()),
//         mrParameters(parameters)
//     {
//     }
    /// Constructor for CheckSkinProcess Process
    CheckSkinProcess(ModelPart& rModelPart,
                     Flags options
                    ):
        Process(),
        mrModelPart(rModelPart),
        mrOptions(options)
    {
    }

    /// Destructor.
    ~CheckSkinProcess() override {}


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

        if(mrModelPart.Conditions().size() == 0 && mrModelPart.Elements().size() != 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "the number of conditions is zero and the number of elements is not, hence the skin can not envelope the domain","")

        typedef std::unordered_map<DenseVector<unsigned int>, unsigned int, KeyHasherRange<DenseVector<unsigned int>>, KeyComparorRange<DenseVector<unsigned int>> > hashmap;
        hashmap edge_map;

        DenseVector<unsigned int> ids(2);

        //add 1 to the counter for every edge find in the model part
        for (ModelPart::ConditionIterator itCond = mrModelPart.ConditionsBegin(); itCond != mrModelPart.ConditionsEnd(); itCond++)
        {
            Element::GeometryType::GeometriesArrayType edges = itCond->GetGeometry().Edges();

            for(unsigned int edge=0; edge<edges.size(); edge++)
            {
                for(unsigned int i=0; i<edges[edge].size(); i++)
                {
                    ids[i] = edges[edge][i].Id();
                }

                //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
                std::sort(ids.begin(), ids.end());

                edge_map[ids] += 1;
            }
        }


        //now loop over the entire edge map.
        //all values shall have a value of 2
        //if that is not the case throw an error
        std::stringstream buffer;

        for(auto it=edge_map.begin(); it!=edge_map.end(); it++)
        {
//             std::cout << it->first << " " << it->second << std::endl;
            if(it->second > 2)
            {
                buffer << std::endl << " the edge between nodes " << it->first[0] << " and " << it->first[1] << std::endl;
                buffer << " belongs to an overlapping condition " << std::endl;
                KRATOS_THROW_ERROR(std::invalid_argument,"ERROR OVERLAPPING CONDITIONS IN SKIN FOUND : ", buffer.str());
            }
            else if(it->second < 2)
            {
                buffer << std::endl << " the edge between nodes " << it->first[0] << " and " << it->first[1] << std::endl;
                buffer << " only appears once, hence it belongs to a non watertight boundary " << std::endl;
                KRATOS_THROW_ERROR(std::invalid_argument,"ERROR NON CLOSED SKIN ", buffer.str());
            }
        }

        std::cout << "checked " << edge_map.size() << " edges in the skin. No gap or overlap found " << std::endl;


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
        return "CheckSkinProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CheckSkinProcess";
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
    CheckSkinProcess& operator=(CheckSkinProcess const& rOther);

    /// Copy constructor.
    CheckSkinProcess(CheckSkinProcess const& rOther);


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
                                  CheckSkinProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const CheckSkinProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}



} // namespace Kratos


#endif // KRATOS_VERIFY_WATERTIGHTNESS_PROCESS_H
