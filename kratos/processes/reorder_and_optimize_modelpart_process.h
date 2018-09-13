//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

#if !defined(KRATOS_REORDER_AND_OPTIMIZE_MODELPART_PROCESS_H_INCLUDED )
#define  KRATOS_REORDER_AND_OPTIMIZE_MODELPART_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"


namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(KRATOS_CORE) ReorderAndOptimizeModelPartProcess : public Process
{
public:
    using GeometryType = Geometry<Point >;
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ReorderAndOptimizeModelPartProcess
    KRATOS_CLASS_POINTER_DEFINITION(ReorderAndOptimizeModelPartProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor is deleted.
    ReorderAndOptimizeModelPartProcess() = delete;

    /// Constructor to be used. Takes the geometry to be meshed and ModelPart to be filled
	ReorderAndOptimizeModelPartProcess(ModelPart& rModelPart, Parameters settings);

    /// The object is not copyable.
    ReorderAndOptimizeModelPartProcess(ReorderAndOptimizeModelPartProcess const& rOther) = delete;

    /// Destructor.
	~ReorderAndOptimizeModelPartProcess() override {};

    ///@}
    ///@name Operators
    ///@{

    /// It is not assignable.
    ReorderAndOptimizeModelPartProcess& operator=(ReorderAndOptimizeModelPartProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void Execute() override;


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
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    ModelPart& mrModelPart;


    ///@}
    ///@name Private Operations
    ///@{
    void ActualizeSubModelPart(ModelPart& subpart);
    void OptimizeOrdering();
    void ReorderElements();

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{

    template <bool reverse = false>
    struct CuthillMcKee
    {
        template <class Matrix>
        static void get(const Matrix &A, std::vector<int>& invperm)
        {
            const int n = A.size1();
            std::vector<int> perm(n);

            /* The data structure used to sort and traverse the level sets:
             *
             * The current level set is currentLevelSet;
             * In this level set, there are nodes with degrees from 0 (not really
             * useful) to maxDegreeInCurrentLevelSet.
             * firstWithDegree[i] points to a node with degree i, or to -1 if it
             * does not exist. nextSameDegree[firstWithDegree[i]] points to the
             * second node with that degree, etc.
             * While the level set is being traversed, the structure for the next
             * level set is generated; nMDICLS will be the next
             * maxDegreeInCurrentLevelSet and nFirstWithDegree will be
             * firstWithDegree.
             */
            int initialNode = 0; // node to start search
            int maxDegree   = 0;

            std::vector<int> degree(n);
            std::vector<int> levelSet(n, 0);
            std::vector<int> nextSameDegree(n, -1);

            for(int i = 0; i < n; ++i)
            {
                degree[i] = A.index1_data()[i+1] - A.index1_data()[i]; //backend::row_nonzeros(A, i);
                maxDegree = std::max(maxDegree, degree[i]);
            }

            std::vector<int> firstWithDegree(maxDegree + 1, -1);
            std::vector<int> nFirstWithDegree(maxDegree + 1);

            // Initialize the first level set, made up by initialNode alone
            perm[0] = initialNode;
            int currentLevelSet = 1;
            levelSet[initialNode] = currentLevelSet;
            int maxDegreeInCurrentLevelSet = degree[initialNode];
            firstWithDegree[maxDegreeInCurrentLevelSet] = initialNode;

            // Main loop
            for (int next = 1; next < n; )
            {
                int nMDICLS = 0;
                std::fill(nFirstWithDegree.begin(), nFirstWithDegree.end(), -1);
                bool empty = true; // used to detect different connected components

                int firstVal  = reverse ? maxDegreeInCurrentLevelSet : 0;
                int finalVal  = reverse ? -1 : maxDegreeInCurrentLevelSet + 1;
                int increment = reverse ? -1 : 1;

                for(int soughtDegree = firstVal; soughtDegree != finalVal; soughtDegree += increment)
                {
                    int node = firstWithDegree[soughtDegree];
                    while (node > 0)
                    {
                        // Visit neighbors
                        for(auto a = A.index1_data()[node] /*backend::row_begin(A, node)*/; a<A.index1_data()[node+1]; ++a)
                        {
                            int c = A.index2_data()[a]; //a.col();
                            if (levelSet[c] == 0)
                            {
                                levelSet[c] = currentLevelSet + 1;
                                perm[next] = c;
                                ++next;
                                empty = false; // this level set is not empty
                                nextSameDegree[c] = nFirstWithDegree[degree[c]];
                                nFirstWithDegree[degree[c]] = c;
                                nMDICLS = std::max(nMDICLS, degree[c]);
                            }
                        }
                        node = nextSameDegree[node];
                    }
                }

                ++currentLevelSet;
                maxDegreeInCurrentLevelSet = nMDICLS;
                for(int i = 0; i <= nMDICLS; ++i)
                    firstWithDegree[i] = nFirstWithDegree[i];

                if (empty)
                {
                    // The graph contains another connected component that we
                    // cannot reach.  Search for a node that has not yet been
                    // included in a level set, and start exploring from it.
                    for(int i = 0; i < n; ++i)
                    {
                        if (levelSet[i] == 0)
                        {
                            perm[next] = i;
                            ++next;
                            levelSet[i] = currentLevelSet;
                            maxDegreeInCurrentLevelSet = degree[i];
                            firstWithDegree[maxDegreeInCurrentLevelSet] = i;
                            break;
                        }
                    }
                }
            }

            //computing the inverse permutation
            #pragma omp parallel for
            for(int i = 0; i < n; ++i) invperm[perm[i]] = i;

        }
    };

    ///@}
    ///@name Un accessible methods
    ///@{



    ///@}

}; // Class ReorderAndOptimizeModelPartProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ReorderAndOptimizeModelPartProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ReorderAndOptimizeModelPartProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_REORDER_AND_OPTIMIZE_MODELPART_PROCESS_H_INCLUDED  defined
