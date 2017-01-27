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

// System includes
#include <limits>


// External includes


// Project includes
#include "includes/exception.h"
#include "processes/measure_mesh_quality_process.h"

namespace Kratos
{

	MeasureMeshQualityProcess::MeasureMeshQualityProcess(ModelPart& rModelPart, std::size_t Dimension)
		: mrModelPart(rModelPart)
		, mDimension(Dimension)
		, mNumberOfInvertedElements(0)
		, mNumberOfSlivers(0)
		, mMinArea(std::numeric_limits<double>::max())
		, mMaxArea(std::numeric_limits<double>::lowest())
		, mMinH(std::numeric_limits<double>::max())
		, mMinAngle(std::numeric_limits<double>::max())
	{

	}

	MeasureMeshQualityProcess::~MeasureMeshQualityProcess()
	{

	}

	void MeasureMeshQualityProcess::Execute()
	{
		KRATOS_TRY

		ResetMeasures();

		if (mDimension == 2)
			for (ModelPart::ElementIterator i_element = mrModelPart.ElementsBegin(); i_element != mrModelPart.ElementsEnd(); i_element++)
			{
				double area = i_element->GetGeometry().Area();
				if (area < 0.00)
					mNumberOfInvertedElements++;
				if (area < mMinArea)
					mMinArea = area;
				if (area > mMaxArea)
					mMaxArea = area;
			}
		else
			KRATOS_ERROR << "Un supported dimension " << mDimension;

		KRATOS_CATCH("");
	}

	/// Turn back information as a string.
	std::string MeasureMeshQualityProcess::Info() const
	{
		return "MeasureMeshQualityProcess";
	}

	/// Print information about this object.
	void MeasureMeshQualityProcess::PrintInfo(std::ostream& rOStream) const
	{
		rOStream << Info();
	}

	/// Print object's data.
	void MeasureMeshQualityProcess::PrintData(std::ostream& rOStream) const
	{
		rOStream << "Number of Inverted Elements : " << mNumberOfInvertedElements << std::endl;
		rOStream << "Min Area                    : " << mMinArea << std::endl;
		rOStream << "Max Area                    : " << mMaxArea << std::endl;
	}

	void MeasureMeshQualityProcess::ResetMeasures()
	{
		mNumberOfInvertedElements = 0;
		mNumberOfSlivers = 0;
		mMinArea = std::numeric_limits<double>::max();
		mMaxArea = std::numeric_limits<double>::lowest();
		mMinH = std::numeric_limits<double>::max();
		mMinAngle = std::numeric_limits<double>::max();
	}

}  // namespace Kratos.
