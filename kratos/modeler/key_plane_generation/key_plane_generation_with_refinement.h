//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Rubio
//

#pragma once

// System includes

// External includes

// Project includes
#include "key_plane_generation.h"

namespace Kratos {

struct RefinementZone
{
    array_1d<double,3> MinPos;
    array_1d<double,3> MaxPos;
    array_1d<double,3> Dx;

    RefinementZone(const array_1d<double,3>& rMinPos, const array_1d<double,3>& rMaxPos, const array_1d<double, 3>& rDx): MinPos(rMinPos), MaxPos(rMaxPos), Dx(rDx)
    {}
};

class KeyPlaneGenerationWithRefinement: public VoxelMesherKeyPlaneGeneration {
    typedef std::pair<array_1d<double,3>, array_1d<double,3>> BoundingBoxType;
public:
    KeyPlaneGenerationWithRefinement(VoxelMeshGeneratorModeler& rModeler, Parameters GenerationParameters):
        VoxelMesherKeyPlaneGeneration(rModeler, GenerationParameters)
    {}

    ~KeyPlaneGenerationWithRefinement() override = default;

    void ValidateParameters() override;

    void Generate() override;

protected:

    /**
     * @brief Retrieves the default parameters for the current utility.
     * @details This method provides a set of default parameters that can be used to initialize or reset the object's state.
     * @return Parameters A Parameters object containing the default values.
     */
    Parameters GetDefaultParameters() const override;

    /**
     * @brief Retrieves the default parameters for the refinement zone.
     * @details This method returns a set of parameters specifically tailored for controlling the behavior and properties of the refinement zone.
     * @return Parameters A Parameters object with default settings for the refinement zone.
     */
    Parameters GetDefaultRefinementZoneParameters();

    /**
     * @brief Calculates the bounding box for a given model part.
     * details This method computes the bounding box that encloses the specified model part, which is useful for spatial queries and optimization tasks.
     * @param rMyModelPart Reference to the ModelPart object for which the bounding box is to be calculated.
     * @return BoundingBoxType The bounding box encapsulating the given model part.
     */
    BoundingBoxType GetBoundingBox(const ModelPart& rMyModelPart);

    /**
     * @brief Checks if the provided parameters object is empty.
     * @details This method determines whether a given Parameters object lacks configuration or is considered empty.
     * @param param The Parameters object to check.
     * @return bool True if the Parameters object is empty; otherwise, false.
     */
    bool IsEmpty(Parameters param);

    /**
     * @brief Identifies refinement areas and calculates voxel sizes for each area.
     * @details This method processes a list of refinement areas and their corresponding voxel sizes based on the global voxel size configuration. It helps in refining the mesh or grid structure in specific regions.
     * @param rRefinementAreas Reference to a vector of BoundingBoxType that will be filled with the identified refinement areas.
     * @param rVoxelSizes Reference to a vector of voxel sizes, one for each refinement area.
     * @param rGlobalVoxelSize The global voxel size used as a reference for calculation.
     */
    void FindRefinements(std::vector<BoundingBoxType>& rRefinementAreas, std::vector<array_1d<double,3>>& rVoxelSizes,const array_1d<double,3>& rGlobalVoxelSize);

    /**
     * @brief Computes partitions and voxel sizes by direction based on refinement areas.
     * @details This method calculates the partitions and their respective voxel sizes in each direction, considering the provided refinement areas and voxel sizes. It's useful for adaptive mesh refinement.
     * @param rRefinementAreas Reference to a vector of BoundingBoxType indicating the areas of refinement.
     * @param rVoxelSizes Reference to a vector of voxel sizes for each refinement area.
     * @param rGlobalBoundingBox The global bounding box within which the calculations are performed.
     * @param rGlobalVoxelSize The global voxel size, serving as a baseline for calculations.
     * @return A pair of array_1d objects, each containing vectors of doubles representing the partition limits and the theoretical voxel sizes in each direction.
     */
    std::pair<array_1d<std::vector<double>,3>,array_1d<std::vector<double>,3>> ComputePartitionsAndVoxelSizeByDirection(std::vector<BoundingBoxType>& rRefinementAreas,
                                                    std::vector<array_1d<double,3>>& rVoxelSizes,
                                                    BoundingBoxType& rGlobalBoundingBox,
                                                    const array_1d<double,3>& rGlobalVoxelSize);

   /**
    * @brief Returns the theoretical voxel size for a specific position and direction within the global bounding box.
    * @details This method calculates the theoretical voxel size at a given position and direction, taking into account the refinement areas and their voxel sizes.
    * @param rRefinementAreas Reference to a vector of BoundingBoxType indicating the areas of refinement.
    * @param rVoxelSizes Reference to a vector of voxel sizes for each refinement area.
    * @param rGlobalBoundingBox The global bounding box within which the calculations are performed.
    * @param rGlobalVoxelSize The baseline global voxel size.
    * @param Direction The direction in which the voxel size is calculated.
    * @param Position The position at which the voxel size is calculated.
    * @return double The calculated theoretical voxel
    */
    double ReturnLocalTheoreticalVoxelSize(std::vector<BoundingBoxType>& rRefinementAreas,
                                                    std::vector<array_1d<double,3>>& rVoxelSizes,
                                                    BoundingBoxType& rGlobalBoundingBox,
                                                    const array_1d<double,3>& rGlobalVoxelSize,
                                                    int Direction,
                                                    double Position);

    /**
    * @brief Generates keyplanes based on partition limits and theoretical voxel sizes.
    * @details This function is designed to process the partitioning of a computational domain into keyplanes based on the provided partition limits and the theoretical voxel sizes. Keyplanes are essential for optimizing spatial queries and data access in computational geometry and volume processing tasks.
    * @param PartitionLimits An array of three std::vector<double> elements representing the partition limits in each of the three dimensions (X, Y, Z). Each vector contains the boundary values that define the partitions in its respective dimension.
    * @param TheoreticalVoxelSize An array of three std::vector<double> elements representing the theoretical voxel sizes in each of the three dimensions (X, Y, Z). Each vector contains values that define the desired resolution of the voxels in its respective dimension.
    */
    void GenerateKeyplanes(array_1d<std::vector<double>,3>& rPartitionLimits, array_1d<std::vector<double>,3>& rTheoreticalVoxelSize);

    /**
     * @brief Computes the effective voxel size given the minimum and maximum coordinates of a region and the input voxel size.
     * The effective voxel size is calculated to fit within the specified region defined by minimum and maximum coordinates while adhering to the input voxel size as closely as possible. This calculation is crucial for adjusting the resolution of a discretized space to match specific computational requirements.
     * @param rMinCoordinate An array_1d<double,3> representing the minimum coordinates (Xmin, Ymin, Zmin) of the region for which the effective voxel size is being computed.
     * @param rMaxCoordinate An array_id<double,3> representing the maximum coordinates (Xmax, Ymax, Zmax) of the region for which the effective voxel size is being computed.
     * @param rInputVoxelSize An array_id<double,3> representing the initial voxel size (Xsize, Ysize, Zsize) as input for the computation.
     * @return array_1d<double,3> representing the effective voxel size (Xeff, Yeff, Zeff) that best fits within the specified region while considering the input voxel size.
     */
    array_1d<double,3> ComputEffectiveVoxelSize(const array_1d<double,3>& rMinCoordinate,
                                                const array_1d<double,3>& rMaxCoordinate,
                                                const array_1d<double,3>& rInputVoxelSize);

};

}
