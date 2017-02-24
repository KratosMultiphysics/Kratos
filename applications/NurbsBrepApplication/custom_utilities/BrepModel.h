#ifndef BREP_MODEL_H
#define BREP_MODEL_H

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
#include "Face.h"
#include "Edge.h"

// ==============================================================================

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

	//typedef std::vector<Vertex> VerticesVector;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.

 */

class BrepModel
{
public:
	///@name Type Definitions
	///@{
	typedef std::vector<Face> FacesVector;
	typedef std::vector<Edge> EdgesVector;



	///@}

	/// Pointer definition of BrepModel
	//    KRATOS_CLASS_POINTER_DEFINITION[BrepModel];

	/// Default constructor.
	BrepModel(FacesVector faces, EdgesVector edges)
	: m_faces(faces),
	  m_edges(edges)
	{
	}

	/// Destructor.
	virtual ~BrepModel()
	{
	}

	// ==============================================================================
	/// Turn back information as a string.
	virtual std::string Info() const
	{
		return "BrepModel";
	}

	// ==============================================================================
	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "BrepModel";
	}

	// ==============================================================================
	/// Print object's data.
	virtual void PrintData(std::ostream &rOStream) const
	{
	}


private:
	// ==============================================================================
	// Initialized by class constructor
	// ==============================================================================
	FacesVector m_faces;
	EdgesVector m_edges;
	//VerticesVector m_vertices;

	// ==============================================================================
	// General working arrays
	// ==============================================================================
	/// Assignment operator.
	//      BrepModel& operator=[BrepModel const& rOther];

	/// Copy constructor.
	//      BrepModel[BrepModel const& rOther];

}; // Class BrepModel

} // namespace Kratos.

#endif // BREP_MODEL_H
