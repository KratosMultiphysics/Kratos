/*
==============================================================================
KratosIncompressibleFluidApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */
/*
 * File:   estimate_dt_utilities.h
 * Author: jcotela
 *
 * Created on January 18, 2011, 6:45 PM
 */

#ifndef KRATOS_ESTIMATE_DT_UTILITIES_H
#define	KRATOS_ESTIMATE_DT_UTILITIES_H

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "utilities/geometry_utilities.h"
#include "utilities/openmp_utils.h"


namespace Kratos {

    ///@addtogroup IncompressibleFluidApplication
    ///@{

    ///@name Kratos Classes
    ///@{

    /// Estimate the time step in a fluid problem to obtain a given Courant number.
    template< unsigned int TDim >
    class EstimateDtUtil
    {
    public:

        ///@name Life Cycle
        ///@{

        /// Constructor
        /**
         * @param rModelPart The model part containing the problem mesh
         */
        EstimateDtUtil(ModelPart& rModelPart):
            mrModelPart(rModelPart)
        {}

        /// Destructor
        ~EstimateDtUtil()
        {}
        
        ///@}
        ///@name Operations
        ///@{

        /// Calculate the maximum time step that satisfies the Courant-Friedrichs-Lewy (CFL) condition.
        /**
         * @param CFL The upper limit for the Courant number
         * @param dt_max Maximum admissible time step (upper bound to be used for situations with very low velocity fields)
         * @return A time step value that satisfies the CFL condition for the current mesh and velocity field
         */
        double EstimateDt(double CFL, double dt_max)
        {
            KRATOS_TRY;

            const unsigned int NumNodes = TDim +1;

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
                array_1d<double, NumNodes> N;
                boost::numeric::ublas::bounded_matrix<double, NumNodes, TDim> DN_DX;

                for( ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
                {
                    // Get the element's geometric parameters
                    Geometry< Node<3> >& rGeom = itElem->GetGeometry();
                    GeometryUtils::CalculateGeometryData(rGeom, DN_DX, N, Area);
                    
                    // Elemental Velocity
                    array_1d<double,3> ElementVel = N[0]*itElem->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
                    for (unsigned int i = 1; i < NumNodes; ++i)
                        ElementVel += N[i]*rGeom[i].FastGetSolutionStepValue(VELOCITY);
                    
                    // Velocity norm
                    double VelNorm = ElementVel[0]*ElementVel[0];
                    for (unsigned int d = 1; d < TDim; ++d)
                        VelNorm += ElementVel[d]*ElementVel[d];
                    VelNorm = sqrt(VelNorm);

                    // Maximum element size along the direction of velocity
                    for (unsigned int i = 0; i < NumNodes; ++i)
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

            // Dt to obtain desired CFL
            double dt = CFL / Max;
            if(dt > dt_max)
                dt = dt_max;

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
        void CalculateLocalCFL(double Dt)
        {
            KRATOS_TRY;

            const unsigned int NumNodes = TDim +1;
            
            int NumThreads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::PartitionVector ElementPartition;
            OpenMPUtils::DivideInPartitions(mrModelPart.NumberOfElements(),NumThreads,ElementPartition);

            #pragma omp parallel
            {
                int k = OpenMPUtils::ThisThread();
                ModelPart::ElementIterator ElemBegin = mrModelPart.ElementsBegin() + ElementPartition[k];
                ModelPart::ElementIterator ElemEnd = mrModelPart.ElementsBegin() + ElementPartition[k+1];
                
                double Area;
                array_1d<double, NumNodes> N;
                boost::numeric::ublas::bounded_matrix<double, NumNodes, TDim> DN_DX;

                for( ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
                {
                    // Get the element's geometric parameters
                    Geometry< Node<3> >& rGeom = itElem->GetGeometry();
                    GeometryUtils::CalculateGeometryData(rGeom, DN_DX, N, Area);
                    
                    // Elemental Velocity
                    array_1d<double,3> ElementVel = N[0]*itElem->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
                    for (unsigned int i = 1; i < NumNodes; ++i)
                        ElementVel += N[i]*rGeom[i].FastGetSolutionStepValue(VELOCITY);
                    
                    // Velocity norm
                    double VelNorm = ElementVel[0]*ElementVel[0];
                    for (unsigned int d = 1; d < TDim; ++d)
                        VelNorm += ElementVel[d]*ElementVel[d];
                    VelNorm = sqrt(VelNorm);

                    double ElemProj = 0.0;
                    // Maximum element size along the direction of velocity
                    for (unsigned int i = 0; i < NumNodes; ++i)
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

        /// A reference to the problem's model part
        ModelPart& mrModelPart;
        
        ///@} // Member variables
        ///@name Serialization
        ///@{

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const
        {
            rSerializer.save("mrModelPart",mrModelPart);
        }

        virtual void load(Serializer& rSerializer)
        {
            rSerializer.load("mrModelPart",mrModelPart);
        }

        ///@}

    };

    ///@} // Kratos classes

    ///@}

} // namespace Kratos.


#endif	/* KRATOS_ESTIMATE_DT_UTILITIES_H */

