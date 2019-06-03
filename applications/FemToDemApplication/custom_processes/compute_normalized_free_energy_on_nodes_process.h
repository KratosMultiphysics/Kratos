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
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class ComputeNormalizedFreeEnergyOnNodesProcess
 * @ingroup FemToDemApplication
 * @brief This class calculates the Free Energy indicator at all the nodes of the mesh. This indicator will be used to compute the metric and the remesh
 * @author Alejandro Cornejo
 */
class ComputeNormalizedFreeEnergyOnNodesProcess : public Process
{
    ///@name Type Definitions
    ///@{
    typedef ModelPart::ElementsContainerType ElementsArrayType;

    ///@}
    ///@name  Enum's
    ///@{

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

    /**
     * @brief This is the default constructor (double)
     * @param rModelPart The model part to be used
     * @param ThisParameters The input parameters
     */
	ComputeNormalizedFreeEnergyOnNodesProcess(ModelPart& rModelPart, Parameters ThisParameters);

    // Destructor
    ~ComputeNormalizedFreeEnergyOnNodesProcess() override = default;

    void Execute() override;

    /**
     * @brief This method Calculates and extrapolates the
     * free energy indicator over the nodes
     * @param pNodeNormalizedFreeEnergyVector the indicator container
     */
    void NormalizedFreeEnergyExtrapolation(NodeNormalizedFreeEnergy *pNodeNormalizedFreeEnergyVector);

    /**
     * @brief This method computes the free energy 
     * @param rStrainVector The strain vector
     * @param rStressVector The stress vector
     * @param Damage The damage variable
     * @param rMatProps The material properties
     * @param rGeometry The geometry of the element
     */
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