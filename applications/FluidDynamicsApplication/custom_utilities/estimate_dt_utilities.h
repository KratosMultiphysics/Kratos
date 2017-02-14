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

///@addtogroup IncompressibleFluidApplication
///@{

///@name Kratos Classes
///@{

/// Estimate the time step in a fluid problem to obtain a given Courant number.
template< unsigned int TDim, unsigned int TNumNodes=TDim+1 >
class EstimateDtUtility
{
public:

    ///@name Life Cycle
    ///@{

    /// Constructor
    /**
     * @param rModelPart The model part containing the problem mesh
     * @param CFL The user-defined Courant-Friedrichs-Lewy number
     * @param DtMax user-defined maximum time increment allowed
     */
    EstimateDtUtility(ModelPart& rModelPart, const double CFL, const double DtMax)
    {
        mCFL = CFL;
        mDtMax = DtMax;
        mrModelPart = rModelPart;
    }

    /// Constructor with Kratos parameters
    /**
     * @param rModelPart The model part containing the problem mesh
     * @param rParameters Kratos parameters containing the CFL number and max time step
     */
    EstimateDtUtility(ModelPart& rModelPart, Parameters& rParameters)
    {
        Parameters defaultParameters(R"({
            "automatic_time_step"   : true,
            "CFL_number"            : 1.0,
            "maximum_delta_time"    : 0.01
        })");

        rParameters.ValidateAndAssignDefaults(defaultParameters);

        mCFL = rParameters["CFL_number"].GetDouble();
        mDtMax = rParameters["maximum_delta_time"].GetDouble();
        mrModelPart = rModelPart;
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

        int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ElementPartition;
        OpenMPUtils::DivideInPartitions(mrModelPart.NumberOfElements(),NumThreads,ElementPartition);

        std::vector<double> MaxProj(NumThreads,0.0);

        #pragma omp parallel shared(MaxProj)
        {
            int k = OpenMPUtils::ThisThread();
            ModelPart::ElementIterator ElemBegin = mrModelPart.ElementsBegin() + ElementPartition[k];
            ModelPart::ElementIterator ElemEnd = mrModelPart.ElementsBegin() + ElementPartition[k+1];

            double& rMaxProj = MaxProj[k];

            double Area;
            array_1d<double, TNumNodes> N;
            boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;

            for( ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
            {
                // Get the element's geometric parameters
                Geometry< Node<3> >& rGeom = itElem->GetGeometry();
                GeometryUtils::CalculateGeometryData(rGeom, DN_DX, N, Area);

                // Elemental Velocity
                array_1d<double,3> ElementVel = N[0]*itElem->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
                for (unsigned int i = 1; i < TNumNodes; ++i)
                    ElementVel += N[i]*rGeom[i].FastGetSolutionStepValue(VELOCITY);

                // Maximum element size along the direction of velocity
                for (unsigned int i = 0; i < TNumNodes; ++i)
                {
                    double Proj = 0.0;
                    for (unsigned int d = 0; d < TDim; ++d)
                        Proj += ElementVel[d]*DN_DX(i,d);
                    Proj = fabs(Proj);
                    if (Proj > rMaxProj) rMaxProj = Proj;
                }
            }
        }

        // Obtain the maximum projected element size (compare thread results)
        double Max = 0.0;
        for (int k = 0; k < NumThreads; ++k)
            if (Max < MaxProj[k]) Max = MaxProj[k];

        // Dt to obtain user-defined CFL
        double dt = mCFL / Max;
        if(dt > mDtMax)
            dt = mDtMax;

        // Perform MPI sync if needed
        double global_dt = dt;
        mrModelPart.GetCommunicator().MinAll(global_dt);
        dt = global_dt;

        return dt;

        KRATOS_CATCH("")
    }

    /// Calculate each element's CFL for a given time step value.
    /**
     * This function is mainly intended for test and debug purposes, but can
     * be sometimes useful to detect where a mesh is inadequate. CFL values
     * are stored as the elemental value of DIVPROJ, so be careful with elements
     * that use it, such as Fluid and VMS (in particular, avoid printing both
     * nodal and elemental values of the variable in GiD)
     * @param Dt The time step used by the fluid solver
     */
    void CalculateLocalCFL(const double Dt)
    {
        KRATOS_TRY;

        int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ElementPartition;
        OpenMPUtils::DivideInPartitions(mrModelPart.NumberOfElements(),NumThreads,ElementPartition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();
            ModelPart::ElementIterator ElemBegin = mrModelPart.ElementsBegin() + ElementPartition[k];
            ModelPart::ElementIterator ElemEnd = mrModelPart.ElementsBegin() + ElementPartition[k+1];

            double Area;
            array_1d<double, TNumNodes> N;
            boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;

            for( ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
            {
                // Get the element's geometric parameters
                Geometry< Node<3> >& rGeom = itElem->GetGeometry();
                GeometryUtils::CalculateGeometryData(rGeom, DN_DX, N, Area);

                // Elemental Velocity
                array_1d<double,3> ElementVel = N[0]*itElem->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
                for (unsigned int i = 1; i < TNumNodes; ++i)
                    ElementVel += N[i]*rGeom[i].FastGetSolutionStepValue(VELOCITY);

                // Maximum element size along the direction of velocity
                double ElemProj = 0.0;
                for (unsigned int i = 0; i < TNumNodes; ++i)
                {
                    double Proj = 0.0;
                    for (unsigned int d = 0; d < TDim; ++d)
                        Proj += ElementVel[d]*DN_DX(i,d);
                    Proj = fabs(Proj);
                    if (Proj > ElemProj) ElemProj = Proj;
                }
                itElem->SetValue(DIVPROJ,ElemProj*Dt);
            }
        }

        KRATOS_CATCH("")
    }

    ///@} // Operators

private:

    ///@name Member Variables
    ///@{

    double              mCFL; // A reference to the user-defined CFL number
    double              mDtMax; // A reference to the user-defined maximum time increment allowed
    ModelPart           mrModelPart; // A reference to the problem's model part

    ///@} // Member variables
    ///@name Serialization
    ///@{

    //friend class Serializer;

    //virtual void save(Serializer& rSerializer) const
    //{
    //    rSerializer.save("mrModelPart",mrModelPart);
    //}

    //virtual void load(Serializer& rSerializer)
    //{
    //    rSerializer.load("mrModelPart",mrModelPart);
    //}

    ///@}

};

///@} // Kratos classes

///@}

} // namespace Kratos.


#endif	/* KRATOS_ESTIMATE_DT_UTILITIES_H */
