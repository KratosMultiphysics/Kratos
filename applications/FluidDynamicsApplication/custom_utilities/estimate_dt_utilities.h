//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela, Ruben Zorrilla
//
//

#ifndef KRATOS_ESTIMATE_DT_UTILITIES_H
#define	KRATOS_ESTIMATE_DT_UTILITIES_H

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "includes/serializer.h"
#include "utilities/openmp_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/// Estimate the time step in a fluid problem to obtain a given Courant number.
template< unsigned int TDim >
class EstimateDtUtility
{
public:

    ///@name Life Cycle
    ///@{

    /// Constructor
    /**
     * @param ModelPart The model part containing the problem mesh
     * @param CFL The user-defined Courant-Friedrichs-Lewy number
     * @param DtMin user-defined minimum time increment allowed
     * @param DtMax user-defined maximum time increment allowed
     */
    EstimateDtUtility(ModelPart &ModelPart, const double CFL, const double DtMin, const double DtMax):
      mrModelPart(ModelPart)
    {
        mCFL = CFL;
        mDtMin = DtMin;
        mDtMax = DtMax;
    }

    /// Constructor with Kratos parameters
    /**
     * @param ModelPart The model part containing the problem mesh
     * @param rParameters Kratos parameters containing the CFL number and max time step
     */
    EstimateDtUtility(ModelPart& ModelPart, Parameters& rParameters):
      mrModelPart(ModelPart)
    {
        Parameters defaultParameters(R"({
            "automatic_time_step"   : true,
            "CFL_number"            : 1.0,
            "minimum_delta_time"    : 1e-4,
            "maximum_delta_time"    : 0.1
        })");

        rParameters.ValidateAndAssignDefaults(defaultParameters);

        mCFL = rParameters["CFL_number"].GetDouble();
        mDtMin = rParameters["minimum_delta_time"].GetDouble();
        mDtMax = rParameters["maximum_delta_time"].GetDouble();
    }

    /// Destructor
    ~EstimateDtUtility()
    {}

    ///@}
    ///@name Operations
    ///@{

    /// Set the CFL value.
    /**
    * @param CFL the user-defined CFL number used in the automatic time step computation
    */
    void SetCFL(const double CFL)
    {
        mCFL = CFL;
    }

    /// Set the maximum time step allowed value.
    /**
    * @param CFL the user-defined CFL number used in the automatic time step computation
    */
    void SetDtMin(const double DtMin)
    {
        mDtMin = DtMin;
    }

    /// Set the maximum time step allowed value.
    /**
    * @param CFL the user-defined CFL number used in the automatic time step computation
    */
    void SetDtMax(const double DtMax)
    {
        mDtMax = DtMax;
    }

    /// Calculate the maximum time step that satisfies the Courant-Friedrichs-Lewy (CFL) condition.
    /**
     * @return A time step value that satisfies the CFL condition for the current mesh and velocity field
     */
    double EstimateDt()
    {
        KRATOS_TRY;

        unsigned int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ElementPartition;
        OpenMPUtils::DivideInPartitions(mrModelPart.NumberOfElements(),NumThreads,ElementPartition);

        double CurrentDt = mrModelPart.GetProcessInfo().GetValue(DELTA_TIME);

        std::vector<double> MaxCFL(NumThreads,0.0);

        #pragma omp parallel shared(MaxCFL)
        {
            int k = OpenMPUtils::ThisThread();
            ModelPart::ElementIterator ElemBegin = mrModelPart.ElementsBegin() + ElementPartition[k];
            ModelPart::ElementIterator ElemEnd = mrModelPart.ElementsBegin() + ElementPartition[k+1];

            GeometryDataContainer GeometryInfo;

            double MaxLocalCFL = 0.0;

            for( ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
            {
                double ElementCFL = CalculateElementCFL(*itElem,GeometryInfo,CurrentDt);
                if (ElementCFL > MaxLocalCFL)
                {
                   MaxLocalCFL = ElementCFL;
                }
            }

            MaxCFL[k] = MaxLocalCFL;
        }

        // Reduce to maximum the thread results
        // Note that MSVC14 does not support max reductions, which are part of OpenMP 3.1
        double CurrentCFL = MaxCFL[0];
        for (unsigned int k = 1; k < NumThreads; k++)
        {
            if (CurrentCFL > MaxCFL[k]) CurrentCFL = MaxCFL[k];
        }

        double NewDt = 0.0;

        // Avoid division by 0 when the maximum CFL number is close to 0 (e.g. problem initialization)
        if (CurrentCFL < 1e-10)
        {
            KRATOS_INFO("EstimateDtUtility") << "Setting minimum delta time " << mDtMin << " as current time step." << std::endl;
            NewDt = mDtMin;
        }
        else
        {
            // Compute new Dt
            NewDt = mCFL * CurrentDt / CurrentCFL;
            // Limit max and min Dt
            if (NewDt > mDtMax)
            {
                NewDt = mDtMax;
            }
            else if (NewDt < mDtMin)
            {
                NewDt = mDtMin;
            }
        }

        // Perform MPI sync if needed
        mrModelPart.GetCommunicator().MinAll(NewDt);

        return NewDt;

        KRATOS_CATCH("")
    }

    /// Calculate each element's CFL for the current time step.
    /**
     * The elemental CFL is stored in the CFL_NUMBER elemental variable.
     * To view it in the post-process file, remember to print CFL_NUMBER as a Gauss Point result.
     */
    void CalculateLocalCFL()
    {
        KRATOS_TRY;

        unsigned int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ElementPartition;
        OpenMPUtils::DivideInPartitions(mrModelPart.NumberOfElements(),NumThreads,ElementPartition);

        const double CurrentDt = mrModelPart.GetProcessInfo().GetValue(DELTA_TIME);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();
            ModelPart::ElementIterator ElemBegin = mrModelPart.ElementsBegin() + ElementPartition[k];
            ModelPart::ElementIterator ElemEnd = mrModelPart.ElementsBegin() + ElementPartition[k+1];

            GeometryDataContainer GeometryInfo;

            for( ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
            {
                double ElementCFL = CalculateElementCFL(*itElem,GeometryInfo,CurrentDt);
                itElem->SetValue(CFL_NUMBER,ElementCFL);
            }
        }

        KRATOS_CATCH("")
    }

    ///@} // Operators

private:

    ///@name Auxiliary Data types
    ///@{

    struct GeometryDataContainer {
        double Area;
        array_1d<double, TDim+1> N;
        BoundedMatrix<double, TDim+1, TDim> DN_DX;
    };

    ///@}
    ///@name Member Variables
    ///@{

    double    mCFL;         // User-defined CFL number
    double    mDtMax;       // User-defined maximum time increment allowed
    double    mDtMin;       // User-defined minimum time increment allowed
    ModelPart &mrModelPart; // The problem's model part

    ///@} // Member variables
    ///@name Private Operations
    ///@{

    double CalculateElementCFL(Element &rElement, GeometryDataContainer& rGeometryInfo, double Dt)
    {
        double Proj = 0.0;

        // Get the element's geometric parameters
        Geometry< Node<3> >& rGeom = rElement.GetGeometry();
        GeometryUtils::CalculateGeometryData(rGeom, rGeometryInfo.DN_DX, rGeometryInfo.N, rGeometryInfo.Area);

        // Elemental Velocity
        array_1d<double,3> ElementVel = rGeometryInfo.N[0]*rGeom[0].FastGetSolutionStepValue(VELOCITY);
        for (unsigned int i = 1; i < TDim+1; ++i)
            ElementVel += rGeometryInfo.N[i]*rGeom[i].FastGetSolutionStepValue(VELOCITY);

        // Calculate u/h as the maximum projection of the velocity along element heights
        for (unsigned int i = 0; i < TDim+1; ++i)
        {
            for (unsigned int d = 0; d < TDim; ++d)
                Proj += ElementVel[d]*rGeometryInfo.DN_DX(i,d);
            Proj = fabs(Proj);
        }

        return Proj*Dt;
    }


    ///@} // Private Operations

};

///@} // Kratos classes

///@}

} // namespace Kratos.


#endif	/* KRATOS_ESTIMATE_DT_UTILITIES_H */
