#if !defined(KRATOS_BREP_MODEL_H_INCLUDED )
#define  KRATOS_BREP_MODEL_H_INCLUDED


// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <math.h>
#include <vector>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "BrepFace.h"
#include "BrepEdge.h"
#include "BrepVertex.h"

// ==============================================================================

namespace Kratos
{
  ///@name Kratos Classes
  ///@{
  /// Holder of all Brep entities.
  /** This class contains the full Brep description needed to describe all CAD-geometries.
  */
	class BrepModel : public IndexedObject, public Flags
	{
	public:
		///@name Type Definitions
		///@{
		KRATOS_CLASS_POINTER_DEFINITION(BrepModel);

		///@}
		///@name Life Cycle 
		///@{ 
		std::vector<BrepFace>&   GetFaceVector();
		std::vector<BrepEdge>&   GetEdgeVector();
		std::vector<BrepVertex>& GetVertexVector();

		//TODO: you need to give reading access to your internals through the Calculate function
		/// Constructor.
		BrepModel(unsigned int& brep_id,
			std::vector<BrepFace>& faces,
			std::vector<BrepEdge>& edges,
			std::vector<BrepVertex>& vertices
		);

		/// Destructor.
		virtual ~BrepModel();
		///@} 
	protected:

	private:
		///@name Member Variables
		///@{ 
		std::vector<BrepFace>   m_brep_faces;
		std::vector<BrepEdge>   m_brep_edges;
		std::vector<BrepVertex> m_brep_vertices;
		///@}    
	}; // Class BrepModel 

}  // namespace Kratos.

#endif // KRATOS_BREP_MODEL_H_INCLUDED  defined