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

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "includes/key_hash.h"

namespace Kratos
{

/** 
 * @class CheckSkinProcess
 * @ingroup KratosCore
 * @brief This function verifies that the skin has no holes nor overlapped geometries this is accomplished by storing all of the edges in the model in a hash map and verifying that no edge appears more than twice (which would imply an overlap)
 * or less than once (which would imply a gap in the skin)
 * @details In the case such condition is violated an error is thrown
 * @author Riccardo Rossi
 */
class CheckSkinProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(CheckSkinProcess);

    /// The index type definition
    typedef std::size_t IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for CheckSkinProcess Process
    CheckSkinProcess(
        ModelPart& rModelPart,
        Flags Options
        ): Process(Options),
           mrModelPart(rModelPart)
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


    /**
     * @brief Check elements to make sure that their jacobian is positive and conditions to ensure that their face normals point outwards
     */
    void Execute() override
    {
        KRATOS_TRY;
            
        KRATOS_ERROR_IF(mrModelPart.Conditions().size() == 0 && mrModelPart.Elements().size() != 0) << "The number of conditions is zero and the number of elements is not, hence the skin can not envelope the domain" << std::endl;

        typedef std::unordered_map<DenseVector<IndexType>, IndexType, KeyHasherRange<DenseVector<IndexType>>, KeyComparorRange<DenseVector<IndexType>> > hashmap;
        hashmap edge_map;

        DenseVector<IndexType> ids(2);

        // Add 1 to the counter for every edge find in the model part
        for (auto& r_cond : mrModelPart.Conditions()) {
            const auto edges = r_cond.GetGeometry().GenerateEdges();

            for(IndexType edge=0; edge<edges.size(); edge++) {
                for(IndexType i=0; i<edges[edge].size(); i++) {
                    ids[i] = edges[edge][i].Id();
                }

                //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
                std::sort(ids.begin(), ids.end());

                edge_map[ids] += 1;
            }
        }

        // Now loop over the entire edge map.
        // All values shall have a value of 2
        // If that is not the case throw an error
        for(auto it=edge_map.begin(); it!=edge_map.end(); ++it) {
            if(it->second > 2) {
                KRATOS_ERROR << "ERROR OVERLAPPING CONDITIONS IN SKIN FOUND : " << std::endl << "The edge between nodes " << it->first[0] << " and " << it->first[1] << std::endl << " belongs to an overlapping condition " << std::endl;
            } else if(it->second < 2) {
                KRATOS_ERROR << "ERROR NON CLOSED SKIN " << std::endl << "The edge between nodes " << it->first[0] << " and " << it->first[1] << std::endl << " only appears once, hence it belongs to a non watertight boundary " << std::endl;
            }
        }

        KRATOS_INFO("CheckSkinProcess") << "Checked " << edge_map.size() << " edges in the skin. No gap or overlap found " << std::endl;

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
