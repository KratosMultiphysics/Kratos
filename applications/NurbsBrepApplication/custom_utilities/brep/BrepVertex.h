//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//


#if !defined(KRATOS_BREP_VERTEX_H_INCLUDED )
#define  KRATOS_BREP_VERTEX_H_INCLUDED


// Project includes
#include "nurbs_brep_application.h"
#include "nurbs_brep_application_variables.h"


namespace Kratos
{
///@name Kratos Classes
///@{
/// Edge in global space.
/** Detail class definition.
*/
class BrepVertex : public IndexedObject, public Flags
{
public:
	/* Used to descibe the topology of vertices. */
	struct Topology
	{
		int face_id;
		int trim_index;

		Topology(const int& rFaceId, const int& rTrimIndex)
		{
			face_id = rFaceId;
			trim_index = rTrimIndex;
		}
	};

    ///@name Life Cycle 
    ///@{ 
    Topology GetVertexInformation(const int& rTopologyIndex);

    /// Constructor.
    BrepVertex(unsigned int VertexId,
		std::vector<Topology>& BrepVertexTopologyVector,
		int control_point_id,
		Vector& Point);

    /// Destructor.
    virtual ~BrepVertex();

    ///@} 
protected:

private:
	///@name Member Variables
	///@{ 
	// topology 
	std::vector<Topology> m_brep_vertex_topology_vector;

	//3d point
	int    m_control_point_id;
	Vector m_point;
	///@}    
}; // Class BrepVertex 

}  // namespace Kratos.

#endif // KRATOS_BREP_VERTEX_H_INCLUDED  defined