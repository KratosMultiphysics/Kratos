//
//   Project Name:        KratosPfemFluidApplication $
//   Created by:          $Author:           AFranci $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:         August 2017 $
//   Revision:            $Revision:             0.0 $
//
//

#if !defined(KRATOS_SET_ACTIVE_FLAG_MESHER_PROCESS_H_INCLUDED)
#define KRATOS_SET_ACTIVE_FLAG_MESHER_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes

#include "spatial_containers/spatial_containers.h"

#include "custom_processes/set_active_flag_process.hpp"
#include "custom_utilities/mesher_utilities.hpp"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "utilities/math_utils.h"

///VARIABLES used:
//Data:
//StepData:
//Flags:    (checked)
//          (set)
//          (modified)
//          (reset)

namespace Kratos
{

	///@name Kratos Globals
	///@{

	///@}
	///@name Type Definitions
	///@{
	typedef ModelPart::NodesContainerType NodesContainerType;
	typedef ModelPart::ElementsContainerType ElementsContainerType;
	typedef ModelPart::MeshType::GeometryType::PointsArrayType PointsArrayType;

	typedef GlobalPointersVector<Node<3>> NodeWeakPtrVectorType;
	typedef GlobalPointersVector<Element> ElementWeakPtrVectorType;

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
	class SetActiveFlagMesherProcess
		: public SetActiveFlagProcess
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of SetActiveFlagMesherProcess
		KRATOS_CLASS_POINTER_DEFINITION(SetActiveFlagMesherProcess);

		///@}
		///@name Life Cycle
		///@{

		/// Default constructor.
		SetActiveFlagMesherProcess(ModelPart &rModelPart,
								   bool unactivePeakElements,
								   bool unactiveSliverElements,
								   int EchoLevel)
			: SetActiveFlagProcess(rModelPart, unactivePeakElements, unactiveSliverElements, EchoLevel)
		{
		}

		/// Destructor.
		virtual ~SetActiveFlagMesherProcess()
		{
		}

		void operator()()
		{
			Execute();
		}

		///@}
		///@name Operations
		///@{

		void Execute() override{

			KRATOS_TRY
#pragma omp parallel
			{
				double tolerance = 0.0000000001;
		const ProcessInfo &rCurrentProcessInfo = mrModelPart.GetProcessInfo();
		const double timeInterval = rCurrentProcessInfo[DELTA_TIME];
		const unsigned int dimension = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
		unsigned int sliversDetectedFromVolume = 0;
		unsigned int sliversDetectedFromShape = 0;
		// const unsigned int dimension = (itElem)->GetGeometry().WorkingSpaceDimension();
		ModelPart::ElementIterator ElemBegin;
		ModelPart::ElementIterator ElemEnd;
		OpenMPUtils::PartitionedIterators(mrModelPart.Elements(), ElemBegin, ElemEnd);
		double ModelPartVolume = 0;
		if (mUnactiveSliverElements == true)
		{
			MesherUtilities MesherUtils;
			ModelPartVolume = MesherUtils.ComputeModelPartVolume(mrModelPart);
		}
		for (ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
		{
			bool sliverEliminationCriteria = false;
			bool peakElementsEliminationCriteria = false;
			bool wallElementsEliminationCriteria = false;
			unsigned int numNodes = itElem->GetGeometry().size();

			// ELIMINATION CHECK FOR SLIVERS
			if (mUnactiveSliverElements == true && numNodes == (dimension + 1))
			{
				double ElementalVolume = 0;
				if (dimension == 2)
				{
					ElementalVolume = (itElem)->GetGeometry().Area();
				}
				else if (dimension == 3)
				{
					ElementalVolume = 0;
					if (itElem->GetGeometry().Dimension() == 3)
						ElementalVolume = (itElem)->GetGeometry().Volume();
				}
				else
				{
					ElementalVolume = 0;
				}
				double CriticalVolume = 0.005 * ModelPartVolume / double(mrModelPart.Elements().size());
				// if(ElementalVolume<CriticalVolume && ElementalVolume>0){
				if (fabs(ElementalVolume) < CriticalVolume)
				{
					sliverEliminationCriteria = true;
					sliversDetectedFromVolume++;
					// std::cout << "RESET ACTIVE FOR THIS SLIVER! \t";
					// std::cout << "its volume is " << ElementalVolume << " vs CriticalVolume " << CriticalVolume <<"number of elements= "<<mrModelPart.Elements().size()<<std::endl;
					// for (unsigned int i = 0; i < numNodes; i++)
					// {
					// 	itElem->GetGeometry()[i].Set(ACTIVE,true);
					// }
				}

				if (sliverEliminationCriteria == false && dimension == 3)
				{

					array_1d<double, 3> nodeA = itElem->GetGeometry()[0].Coordinates();
					array_1d<double, 3> nodeB = itElem->GetGeometry()[1].Coordinates();
					array_1d<double, 3> nodeC = itElem->GetGeometry()[2].Coordinates();
					array_1d<double, 3> nodeD = itElem->GetGeometry()[3].Coordinates();

					double a1 = 0; //slope x for plane on the first triangular face of the tetrahedra (nodes A,B,C)
					double b1 = 0; //slope y for plane on the first triangular face of the tetrahedra (nodes A,B,C)
					double c1 = 0; //slope z for plane on the first triangular face of the tetrahedra (nodes A,B,C)
					a1 = (nodeB[1] - nodeA[1]) * (nodeC[2] - nodeA[2]) - (nodeC[1] - nodeA[1]) * (nodeB[2] - nodeA[2]);
					b1 = (nodeB[2] - nodeA[2]) * (nodeC[0] - nodeA[0]) - (nodeC[2] - nodeA[2]) * (nodeB[0] - nodeA[0]);
					c1 = (nodeB[0] - nodeA[0]) * (nodeC[1] - nodeA[1]) - (nodeC[0] - nodeA[0]) * (nodeB[1] - nodeA[1]);
					double a2 = 0; //slope x for plane on the second triangular face of the tetrahedra (nodes A,B,D)
					double b2 = 0; //slope y for plane on the second triangular face of the tetrahedra (nodes A,B,D)
					double c2 = 0; //slope z for plane on the second triangular face of the tetrahedra (nodes A,B,D)
					a2 = (nodeB[1] - nodeA[1]) * (nodeD[2] - nodeA[2]) - (nodeD[1] - nodeA[1]) * (nodeB[2] - nodeA[2]);
					b2 = (nodeB[2] - nodeA[2]) * (nodeD[0] - nodeA[0]) - (nodeD[2] - nodeA[2]) * (nodeB[0] - nodeA[0]);
					c2 = (nodeB[0] - nodeA[0]) * (nodeD[1] - nodeA[1]) - (nodeD[0] - nodeA[0]) * (nodeB[1] - nodeA[1]);
					double a3 = 0; //slope x for plane on the third triangular face of the tetrahedra (nodes B,C,D)
					double b3 = 0; //slope y for plane on the third triangular face of the tetrahedra (nodes B,C,D)
					double c3 = 0; //slope z for plane on the third triangular face of the tetrahedra (nodes B,C,D)
					a3 = (nodeB[1] - nodeC[1]) * (nodeD[2] - nodeC[2]) - (nodeD[1] - nodeC[1]) * (nodeB[2] - nodeC[2]);
					b3 = (nodeB[2] - nodeC[2]) * (nodeD[0] - nodeC[0]) - (nodeD[2] - nodeC[2]) * (nodeB[0] - nodeC[0]);
					c3 = (nodeB[0] - nodeC[0]) * (nodeD[1] - nodeC[1]) - (nodeD[0] - nodeC[0]) * (nodeB[1] - nodeC[1]);
					double a4 = 0; //slope x for plane on the fourth triangular face of the tetrahedra (nodes A,C,D)
					double b4 = 0; //slope y for plane on the fourth triangular face of the tetrahedra (nodes A,C,D)
					double c4 = 0; //slope z for plane on the fourth triangular face of the tetrahedra (nodes A,C,D)
					a4 = (nodeA[1] - nodeC[1]) * (nodeD[2] - nodeC[2]) - (nodeD[1] - nodeC[1]) * (nodeA[2] - nodeC[2]);
					b4 = (nodeA[2] - nodeC[2]) * (nodeD[0] - nodeC[0]) - (nodeD[2] - nodeC[2]) * (nodeA[0] - nodeC[0]);
					c4 = (nodeA[0] - nodeC[0]) * (nodeD[1] - nodeC[1]) - (nodeD[0] - nodeC[0]) * (nodeA[1] - nodeC[1]);

					double cosAngle12 = (a1 * a2 + b1 * b2 + c1 * c2) / (sqrt(pow(a1, 2) + pow(b1, 2) + pow(c1, 2)) * sqrt(pow(a2, 2) + pow(b2, 2) + pow(c2, 2)));
					double cosAngle13 = (a1 * a3 + b1 * b3 + c1 * c3) / (sqrt(pow(a1, 2) + pow(b1, 2) + pow(c1, 2)) * sqrt(pow(a3, 2) + pow(b3, 2) + pow(c3, 2)));
					double cosAngle14 = (a1 * a4 + b1 * b4 + c1 * c4) / (sqrt(pow(a1, 2) + pow(b1, 2) + pow(c1, 2)) * sqrt(pow(a4, 2) + pow(b4, 2) + pow(c4, 2)));
					double cosAngle23 = (a3 * a2 + b3 * b2 + c3 * c2) / (sqrt(pow(a3, 2) + pow(b3, 2) + pow(c3, 2)) * sqrt(pow(a2, 2) + pow(b2, 2) + pow(c2, 2)));
					double cosAngle24 = (a4 * a2 + b4 * b2 + c4 * c2) / (sqrt(pow(a4, 2) + pow(b4, 2) + pow(c4, 2)) * sqrt(pow(a2, 2) + pow(b2, 2) + pow(c2, 2)));
					double cosAngle34 = (a4 * a3 + b4 * b3 + c4 * c3) / (sqrt(pow(a4, 2) + pow(b4, 2) + pow(c4, 2)) * sqrt(pow(a3, 2) + pow(b3, 2) + pow(c3, 2)));

					// if two faces are coplanar, I will erase the element (which is probably a sliver)
					double limit = 0.99999;
					if (fabs(cosAngle12) > limit || fabs(cosAngle13) > limit || fabs(cosAngle14) > limit || fabs(cosAngle23) > limit || fabs(cosAngle24) > limit || fabs(cosAngle34) > limit)
					{
						unsigned int fsNodes = 0;
						for (unsigned int i = 0; i < numNodes; i++)
						{
							if (itElem->GetGeometry()[i].Is(FREE_SURFACE))
							{
								fsNodes++;
							}
							NodeWeakPtrVectorType &rN = itElem->GetGeometry()[i].GetValue(NEIGHBOUR_NODES);
							unsigned int neighborNodes = rN.size();
							if (neighborNodes == numNodes)
							{
								fsNodes = 4;
								// std::cout << "ATTENTION WITH THIS, ONE NODE COULD BE ALONE IN THE LINEAR SYSTEM" << std::endl;
							}
						}
						if (fsNodes < 3)
						{
							sliverEliminationCriteria = true;
							sliversDetectedFromShape++;
						}
					}
				}
			}

			// ELIMINATION CHECK FOR PEAK ELEMENTS (those annoying elements created by pfem remeshing and placed bewteen the free-surface and the walls)
			if (mUnactivePeakElements == true && sliverEliminationCriteria == false)
			{
				double scalarProduct = 1.0;
				bool doNotErase = false;
				unsigned int elementRigidNodes = 0;
				for (unsigned int i = 0; i < numNodes; i++)
				{
					if (itElem->GetGeometry()[i].Is(RIGID) && itElem->GetGeometry()[i].IsNot(SOLID))
					{
						elementRigidNodes++;
					}
					if (itElem->GetGeometry()[i].IsNot(RIGID) && itElem->GetGeometry()[i].IsNot(FREE_SURFACE))
					{
						peakElementsEliminationCriteria = false;
						doNotErase = true;
						// break;
					}
					else if (itElem->GetGeometry()[i].Is(RIGID) && itElem->GetGeometry()[i].IsNot(SOLID) && itElem->GetGeometry()[i].Is(FREE_SURFACE) && doNotErase == false)
					{
						peakElementsEliminationCriteria = true;
						const array_1d<double, 3> &wallVelocity = itElem->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
						double normWallVelocity = norm_2(wallVelocity);
						if (normWallVelocity == 0)
						{ // up to now this is for fixed walls only
							for (unsigned int j = 0; j < numNodes; j++)
							{

								if (itElem->GetGeometry()[j].IsNot(RIGID) && itElem->GetGeometry()[j].Is(FREE_SURFACE))
								{
									Point freeSurfaceToRigidNodeVector = Point{itElem->GetGeometry()[i].Coordinates() - itElem->GetGeometry()[j].Coordinates()};
									const array_1d<double, 3> &freeSurfaceVelocity = itElem->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY);

									double freeSurfaceToRigidNodeDistance = sqrt(freeSurfaceToRigidNodeVector[0] * freeSurfaceToRigidNodeVector[0] +
																				 freeSurfaceToRigidNodeVector[1] * freeSurfaceToRigidNodeVector[1] +
																				 freeSurfaceToRigidNodeVector[2] * freeSurfaceToRigidNodeVector[2]);
									double displacementFreeSurface = timeInterval * (sqrt(freeSurfaceVelocity[0] * freeSurfaceVelocity[0] +
																						  freeSurfaceVelocity[1] * freeSurfaceVelocity[1] +
																						  freeSurfaceVelocity[2] * freeSurfaceVelocity[2]));
									if (dimension == 2)
									{
										scalarProduct = freeSurfaceToRigidNodeVector[0] * freeSurfaceVelocity[0] + freeSurfaceToRigidNodeVector[1] * freeSurfaceVelocity[1];
									}
									else if (dimension == 3)
									{
										scalarProduct = freeSurfaceToRigidNodeVector[0] * freeSurfaceVelocity[0] + freeSurfaceToRigidNodeVector[1] * freeSurfaceVelocity[1] + freeSurfaceToRigidNodeVector[2] * freeSurfaceVelocity[2];
									}
									if (scalarProduct > tolerance && displacementFreeSurface > (0.01 * freeSurfaceToRigidNodeDistance))
									{
										// if(scalarProduct>tolerance){
										peakElementsEliminationCriteria = false;
										doNotErase = true;
										break;
									}
									else
									{
										// I will not unactive the element if the free-surface node is sorrounded by rigd nodes only
										NodeWeakPtrVectorType &rN = itElem->GetGeometry()[j].GetValue(NEIGHBOUR_NODES);
										unsigned int rigidNodes = 0;
										unsigned int freeSurfaceNodes = 0;
										for (unsigned int i = 0; i < rN.size(); i++)
										{
											if (rN[i].Is(RIGID) && rN[i].IsNot(SOLID))
												rigidNodes += 1;
											if (rN[i].Is(FREE_SURFACE) && rN[i].IsNot(RIGID))
												freeSurfaceNodes += 1;
										}
										if (dimension == 2)
										{
											if (rigidNodes == rN.size())
											{
												peakElementsEliminationCriteria = false;
												doNotErase = true;
												break;
											}
										}
										else if (dimension == 3)
										{
											if (rigidNodes == rN.size() || freeSurfaceNodes == 1 || (scalarProduct > tolerance && freeSurfaceNodes < 4))
											{
												peakElementsEliminationCriteria = false;
												doNotErase = true;
												break;
											}
										}
									}
								}
							}
						}
					}
				}
				if (elementRigidNodes == numNodes)
				{
					wallElementsEliminationCriteria = true;
					Geometry<Node<3>> wallElementNodes = itElem->GetGeometry();
					this->SetPressureToIsolatedWallNodes(wallElementNodes);
				}
			}
			// ELIMINATION CHECK FOR ELEMENTS FORMED BY WALL PARTICLES ONLY (this is included for computational efficiency purpose also in the previous peak element check)
			else if (mUnactivePeakElements == false)
			{
				unsigned int elementRigidNodes = 0;
				for (unsigned int i = 0; i < numNodes; i++)
				{
					if (itElem->GetGeometry()[i].Is(RIGID) && itElem->GetGeometry()[i].IsNot(SOLID))
					{
						elementRigidNodes++;
					}
				}

				if (elementRigidNodes == numNodes)
				{
					wallElementsEliminationCriteria = true;
					Geometry<Node<3>> wallElementNodes = itElem->GetGeometry();
					this->SetPressureToIsolatedWallNodes(wallElementNodes);
				}
			}

			if (sliverEliminationCriteria == true || peakElementsEliminationCriteria == true || wallElementsEliminationCriteria == true)
			{
				(itElem)->Set(ACTIVE, false);
			}
			else
			{
				(itElem)->Set(ACTIVE, true);
			}
		}
		// if (sliversDetectedFromShape > 0)
		// {
		// 	std::cout << "I have set ACTIVE=false to " << sliversDetectedFromShape << " slivers due to shape." << std::endl;
		// }
		// if (sliversDetectedFromVolume > 0)
		// {
		// 	std::cout << "I have set ACTIVE=false to " << sliversDetectedFromVolume << " slivers due to volume." << std::endl;
		// }
	}

	KRATOS_CATCH(" ")
}; // namespace Kratos

///@}
///@name Operators
///@{

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
	return "SetActiveFlagMesherProcess";
}

/// Print information about this object.
void PrintInfo(std::ostream &rOStream) const override
{
	rOStream << "SetActiveFlagMesherProcess";
}

protected:
///@name Protected static Member Variables
///@{

///@}
///@name Protected member Variables
///@{

///@}
///@name Protected  Access
///@{

///@}
///@name Protected Inquiry
///@{

///@}
///@name Protected LifeCycle
///@{

///@}

private:
///@name Static Member Variables
///@{

///@}
///@name Member Variables
///@{

///@}
///@name Private Operators
///@{

///@}
///@name Private Operations
///@{

///@}
///@name Private  Access
///@{

///@}
///@name Private Inquiry
///@{

///@}
///@name Un accessible methods
///@{

/// Assignment operator.
SetActiveFlagMesherProcess &operator=(SetActiveFlagMesherProcess const &rOther);

/// Copy constructor.
//SetActiveFlagMesherProcess(SetActiveFlagMesherProcess const& rOther);

///@}
}
; // Class SetActiveFlagMesherProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
								SetActiveFlagMesherProcess &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
								const SetActiveFlagMesherProcess &rThis)
{
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_SET_ACTIVE_FLAG_MESHER_PROCESS_H_INCLUDED  defined
