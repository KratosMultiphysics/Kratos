//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:     Máté Kelemen
//

#ifndef KRATOS_HDF5APPLICATION_POINT_SET_OUTPUT_PROCESS_H
#define KRATOS_HDF5APPLICATION_POINT_SET_OUTPUT_PROCESS_H

// Application includes
#include "custom_io/hdf5_file.h"
#include "custom_utilities/vertex.h"
#include "includes/ublas_interface.h"

// Core inlcudes
#include "processes/process.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "utilities/interval_utility.h"
#include "utilities/brute_force_point_locator.h"


namespace Kratos
{
namespace HDF5
{


class PointSetOutputProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PointSetOutputProcess);

    PointSetOutputProcess(Parameters parameters, Model& rModel);

    virtual void ExecuteFinalizeSolutionStep() override;

    virtual const Parameters GetDefaultParameters() const override;

protected:
    /// Constructor for derived @ref{LineOutputProcess}
    PointSetOutputProcess(Parameters parameters,
                          Model& rModel,
                          Detail::VertexContainerType&& rVertices);

private:
    /// BruteForcePointLocator with configuration and tolerance persistence
    class Locator
    {
    public:
        Locator(ModelPart& rModelPart, const Globals::Configuration configuration, const double tolerance);

        int FindElement(const Point& rPoint, Kratos::Vector& rShapeFunctionValues) const;

    private:
        BruteForcePointLocator mLocator;

        const Globals::Configuration mConfiguration;

        const double mTolerance;
    };

    /** Initialize common members
     *  (except Points and members that must be initialized before the constructor body gets executed)
     */
    void InitializeFromParameters(Parameters parameters);

    File::Pointer MakeHDF5File(const std::size_t stepID) const;

    ModelPart& mrModelPart;

    IntervalUtility mInterval;

    Detail::VertexContainerType mVertices;

    std::unique_ptr<Locator> mpLocator;

    std::vector<std::string> mVariableNames;

    std::string mOutputPath;
};


} // namespace HDF5
} // namespace Kratos


#endif