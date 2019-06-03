//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#if !defined(KRATOS_COMPUTE_NORMALIZED_FREE_ENERGY_ON_NODES_PROCESS)
#define KRATOS_COMPUTE_NORMALIZED_FREE_ENERGY_ON_NODES_PROCESS

#include "processes/process.h"
#include "includes/model_part.h"
#include "fem_to_dem_application_variables.h"

namespace Kratos
{
class ComputeNormalizedFreeEnergyOnNodesProcess : public Process
{
  protected:
    struct NodeNormalizedFreeEnergy
    {
        int NElems;
        double NormalizedFreeEnergy;

		NodeNormalizedFreeEnergy()
        {
            NElems = 0;
            NormalizedFreeEnergy = 0.0;
        }
    };

  public:
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    /// Pointer definition of ApplyMultipointConstraintsProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeNormalizedFreeEnergyOnNodesProcess);

    typedef ModelPart::ElementsContainerType ElementsArrayType;

    // Constructor
	ComputeNormalizedFreeEnergyOnNodesProcess(ModelPart &r_model_part, Parameters ThisParameters);

    // Destructor
    ~ComputeNormalizedFreeEnergyOnNodesProcess() override = default;

    void Execute() override;

    void NormalizedFreeEnergyExtrapolation(NodeNormalizedFreeEnergy *pNodeNormalizedFreeEnergyVector);

    double CalculateNormalizedFreeEnergy(
        const Vector& rStrainVector, 
        const Vector& rStressVector, 
        const double Damage, 
        const Properties& rMatProps,
        Geometry<Node<3>>& rGeometry);

	double CalculateCharacteristicLength2D(const Geometry<Node<3>>& rGeometry);
	double CalculateCharacteristicLength3D(Geometry<Node<3>>& rGeometry);
    double ComputeTensionFactor2D(const Vector& rStressVector);
    double ComputeTensionFactor3D(const Vector& rStressVector);
    void CalculatePrincipalStresses2D(const Vector& rStressVector, Vector& rPrincipalStressVector);
    void CalculatePrincipalStresses3D(const Vector& rStressVector, Vector& rPrincipalStressVector);
    void ObtainMaximumNodeId(int &rmax_id);
    double CalculateI1Invariant(const Vector& rStressVector);
    double CalculateI2Invariant(const Vector& rStressVector);
    double CalculateI3Invariant(const Vector& rStressVector);

  protected:
    // Member Variables
    ModelPart& mrModelPart;
    unsigned int mNNodes;
    bool mComputeNormalizedFreeEnergy = false;
    bool mCorrectWithDisplacements = false;
    double mCorrectionFactor = 1.0;

}; // Class ComputeNormalizedFreeEnergyOnNodesProcess

} // namespace Kratos

#endif /* KRATOS_COMPUTE_NORMALIZED_FREE_ENERGY_ON_NODES_PROCESS defined */