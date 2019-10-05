//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:         BSD License
//    Kratos default license:
//  kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#if !defined(KRATOS_GENERATE_DEM_PROCESS)
#define KRATOS_GENERATE_DEM_PROCESS


#include "includes/model_part.h"
#include "processes/process.h"
#include "fem_to_dem_application_variables.h"
#include "custom_utilities/create_and_destroy.h"
#include "processes/find_nodal_neighbours_process.h"

namespace Kratos {

typedef std::size_t SizeType;
typedef Node<3> NodeType;
typedef ModelPart::ElementsContainerType::iterator ElementIteratorType;

class GenerateDemProcess : public Process 
{
 public:

  /// Pointer definition of ApplyMultipointConstraintsProcess
  KRATOS_CLASS_POINTER_DEFINITION(GenerateDemProcess);

  typedef std::size_t IndexType;

  // Constructor
  GenerateDemProcess(ModelPart& rModelPart, ModelPart& rDemModelPart);

  // Destructor
  ~GenerateDemProcess() override = default;

  void operator()() { Execute(); }

    /**
     * @brief Creates the DEM particles
     */
  void Execute() override;

    /**
     * @brief Computes the distance between nodes
     */
  double CalculateDistanceBetweenNodes(const NodeType& Node1, const NodeType& Node2);

    /**
     * @brief This creates one particle
     */
  void CreateDEMParticle(const int Id, const array_1d<double, 3> Coordinates, 
      const Properties::Pointer pProperties, const double Radius, NodeType& rNode); 

    /**
     * @brief This returns the min value of a vector
     */
  double GetMinimumValue(const Vector& rValues);

    /**
     * @brief This returns the maximum Id of the particles
     */
  int GetMaximumDEMId();

protected:

  // Member Variables
  ModelPart& mrModelPart;
  ModelPart& mrDEMModelPart;
  ParticleCreatorDestructor mParticleCreator = ParticleCreatorDestructor();

};  // Class

}  // namespace Kratos
#endif /* KRATOS_GENERATE_DEM_PROCESS defined */