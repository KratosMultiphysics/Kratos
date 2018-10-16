//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Davand
//

// System includes


// External includes


// Project includes
#include "processes/find_intersected_geometrical_objects_process.h"


namespace Kratos
{

	FindIntersectedGeometricalObjectsProcess::FindIntersectedGeometricalObjectsProcess(ModelPart& rPart1, ModelPart& rPart2)
		: mrModelPart1(rPart1), mrModelPart2(rPart2)
	{
	}

	void FindIntersectedGeometricalObjectsProcess::Initialize()
	{
		GenerateOctree();
	}

	void FindIntersectedGeometricalObjectsProcess::FindIntersectedSkinObjects(std::vector<PointerVector<GeometricalObject>>& rResults)
	{
		const std::size_t number_of_elements = mrModelPart1.NumberOfElements();
		auto& r_elements = mrModelPart1.ElementsArray();
		std::vector<OctreeType::cell_type*> leaves;

		rResults.resize(number_of_elements);
		for (std::size_t i = 0; i < number_of_elements; i++) {
			auto p_element_1 = r_elements[i];
			leaves.clear();
			mOctree.GetIntersectedLeaves(p_element_1, leaves);
			FindIntersectedSkinObjects(*p_element_1, leaves, rResults[i]);
		}
	}

	void FindIntersectedGeometricalObjectsProcess::FindIntersections()
	{
		this->FindIntersectedSkinObjects(mIntersectedObjects);
	}

	std::vector<PointerVector<GeometricalObject>>& FindIntersectedGeometricalObjectsProcess::GetIntersections()
	{
		return mIntersectedObjects;
	}

	ModelPart& FindIntersectedGeometricalObjectsProcess::GetModelPart1()
	{
		return mrModelPart1;
	}

	OctreeBinary<OctreeBinaryCell<Internals::DistanceSpatialContainersConfigure>>* FindIntersectedGeometricalObjectsProcess::GetOctreePointer()
	{
		return &mOctree;
	}

	void FindIntersectedGeometricalObjectsProcess::Clear()
	{
		mIntersectedObjects.clear();
	}

	void FindIntersectedGeometricalObjectsProcess::Execute()
	{
		GenerateOctree();

		std::vector<OctreeType::cell_type*> leaves;
		const int number_of_elements = mrModelPart1.NumberOfElements();

		#pragma omp parallel for private(leaves)
		for (int i = 0; i < number_of_elements; i++)
		{
			auto p_element_1 = mrModelPart1.ElementsBegin() + i;
			leaves.clear();
			mOctree.GetIntersectedLeaves(*(p_element_1.base()), leaves);
			MarkIfIntersected(**(p_element_1.base()), leaves);
		}
	}

	/// Turn back information as a string.
	std::string FindIntersectedGeometricalObjectsProcess::Info() const {
		return "FindIntersectedGeometricalObjectsProcess";
	}

	/// Print information about this object.
	void FindIntersectedGeometricalObjectsProcess::PrintInfo(std::ostream& rOStream) const {
		rOStream << Info();
	}

	/// Print object's data.
	void FindIntersectedGeometricalObjectsProcess::PrintData(std::ostream& rOStream) const {

	}

	void FindIntersectedGeometricalObjectsProcess::GenerateOctree() {
		this->SetOctreeBoundingBox();

		// Adding mrModelPart2 to the octree
		for (auto i_node = mrModelPart2.NodesBegin(); i_node != mrModelPart2.NodesEnd(); i_node++) {
#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it 
			mOctree.Insert(i_node->Coordinates().data());

#else
			mOctree.Insert(i_node->Coordinates().data().data());
#endif // ifdef KRATOS_USE_AMATRIX
		}

		for (auto i_element = mrModelPart2.ElementsBegin(); i_element != mrModelPart2.ElementsEnd(); i_element++) {
			mOctree.Insert(*(i_element).base());
		}
	}

	void  FindIntersectedGeometricalObjectsProcess::SetOctreeBoundingBox() {
		Point low(mrModelPart1.NodesBegin()->Coordinates());
		Point high(mrModelPart1.NodesBegin()->Coordinates());

		// loop over all nodes in first modelpart
		for (auto i_node = mrModelPart1.NodesBegin(); i_node != mrModelPart1.NodesEnd(); i_node++) {
			const array_1d<double,3> &r_coordinates = i_node->Coordinates();
			for (int i = 0; i < 3; i++) {
				low[i] = r_coordinates[i] < low[i] ? r_coordinates[i] : low[i];
				high[i] = r_coordinates[i] > high[i] ? r_coordinates[i] : high[i];
			}
		}

		// loop over all skin nodes
		for (auto i_node = mrModelPart2.NodesBegin(); i_node != mrModelPart2.NodesEnd(); i_node++) {
			const array_1d<double,3>& r_coordinates = i_node->Coordinates();
			for (int i = 0; i < 3; i++) {
				low[i] = r_coordinates[i] < low[i] ? r_coordinates[i] : low[i];
				high[i] = r_coordinates[i] > high[i] ? r_coordinates[i] : high[i];
			}
		}

    for(int i = 0 ; i < 3; i++){
			low[i] -= std::numeric_limits<double>::epsilon();
			high[i] += std::numeric_limits<double>::epsilon();
		}


		// TODO: Octree needs refactoring to work with BoundingBox. Pooyan.
#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it 
	mOctree.SetBoundingBox(low.data(), high.data());
#else
	mOctree.SetBoundingBox(low.data().data(), high.data().data());
#endif // ifdef KRATOS_USE_AMATRIX
	
	}

	void  FindIntersectedGeometricalObjectsProcess::MarkIfIntersected(Element& rElement1, std::vector<OctreeType::cell_type*>& leaves) {
		for (auto p_leaf : leaves) {
			for (auto p_element_2 : *(p_leaf->pGetObjects())) {
				if (HasIntersection(rElement1.GetGeometry(),p_element_2->GetGeometry())) {
					rElement1.Set(SELECTED);
					return;
				}
			}
		}
	}

	bool FindIntersectedGeometricalObjectsProcess::HasIntersection(Element::GeometryType& rFirstGeometry, Element::GeometryType& rSecondGeometry) {
		auto faces = rFirstGeometry.Faces();
		for (auto& face : faces) {
			if (face.HasIntersection(rSecondGeometry))
				return true;
		}
		// Let check second geometry is inside the first one.
		// Considering that there are no intersection, if one point is inside all of it is inside.
		array_1d<double, 3> local_point;
		if (rFirstGeometry.IsInside(rSecondGeometry.GetPoint(0), local_point))
			return true;
		return false;
	}

	void FindIntersectedGeometricalObjectsProcess::FindIntersectedSkinObjects(Element& rElement1, std::vector<OctreeType::cell_type*>& leaves, PointerVector<GeometricalObject>& rResults) {
		for (auto p_leaf : leaves) {
			for (auto p_element_2 : *(p_leaf->pGetObjects())) {
				if (HasIntersection(rElement1.GetGeometry(), p_element_2->GetGeometry())) {
					rElement1.Set(SELECTED);
					if(std::find(rResults.ptr_begin(), rResults.ptr_end(), p_element_2) == rResults.ptr_end())
						rResults.push_back(p_element_2);
				}
			}
		}

	}

}  // namespace Kratos.
