/*
==============================================================================
KratosFluidDynamicsApplication
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
 * File:   periodic_condition_utilities.h
 * Author: jcotela
 *
 * Created on September 22, 2011, 3:22 PM
 */

#ifndef KRATOS_PERIODIC_CONDITION_UTILITIES_H
#define	KRATOS_PERIODIC_CONDITION_UTILITIES_H

// System includes
#include <string>
#include <iostream>
#include <algorithm>
#include <math.h>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/model_part.h"
#include "includes/properties.h"
#include "geometries/point_3d.h"
#include "spatial_containers/spatial_containers.h"
#include "containers/variables_list.h"

#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/periodic_variables_container.h"

namespace Kratos
{
    ///@addtogroup FluidDynamicsApplication
    ///@{

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

    /// Auxiliary utilitiy to define periodic boundary conditions for flow problems.
    /**
     * @todo Write new detailed description
     Create a condition linking node pairs, one in each side of the periodic boundary.
     This class assumes conformant meshes at both sides of the boundary (that is,
     pairs of nodes with concordant coordinates on both sides of the periodic
     boundary). The class works by identifying both "sides" of the periodic boundary,
     marked by the user with IS_INTERFACE 1.0 or 2.0 respectively, and creating a
     condition that matches each node on one side to the nearest node to its image
     on the other side.

     In addition to identifying the boundary, the user must also provide the vector
     that provides the translation from the IS_INTERFACE 1.0 side to the other with
     SetShiftTranslation and provide a weight to the penalty function that imposes
     the periodic condition with SetWeight.
     */
    class PeriodicConditionUtilities
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of PeriodicConditionUtilities
        KRATOS_CLASS_POINTER_DEFINITION(PeriodicConditionUtilities);

        typedef std::size_t IndexType;

        typedef std::size_t SizeType;

        //defintions for spatial search
        typedef Node<3>                                     PointType;
        typedef Node<3>::Pointer                            PointTypePointer;
        typedef std::vector<PointType::Pointer>             PointVector;
        typedef std::vector<PointType::Pointer>::iterator   PointIterator;
        typedef std::vector<double>                         DistanceVector;
        typedef std::vector<double>::iterator               DistanceIterator;

        // Bucket types
//        typedef Bucket< 3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator >       BucketType;
        typedef Bins< 3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator >         StaticBins;
//        typedef BinsDynamic< 3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator >  DynamicBins;

        // DynamicBins
//        typedef Tree< KDTreePartition<BucketType> > tree; 		//Kdtree;
//        typedef Tree< OCTreePartition<BucketType> > tree; 		//Octree;
        typedef Tree< StaticBins > tree; 		     		//Binstree;
//        typedef Tree< KDTreePartition<StaticBins> > tree; 		//KdtreeBins;

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        PeriodicConditionUtilities(ModelPart& ThisModelPart,SizeType ThisDomainSize):
            mrModelPart(ThisModelPart),
            mDomainSize(ThisDomainSize),
            mpSearchStrategy()
        {}

        /// Destructor.
        virtual ~PeriodicConditionUtilities()
        {}


        ///@}
        ///@name Operators
        ///@{



        ///@}
        ///@name Operations
        ///@{

        void SetUpSearchStructure(Variable<double> const& rFlagVar,
                                  const double FlagValue)
        {
            KRATOS_TRY

            // Find candidate nodes for periodic conditions
            mCandidateNodes = PointVector();

            for(ModelPart::NodesContainerType::ptr_iterator itNode = mrModelPart.Nodes().ptr_begin();
                itNode != mrModelPart.Nodes().ptr_end(); itNode++)
            {
                if( (**itNode).FastGetSolutionStepValue(rFlagVar,0) == FlagValue)
                    mCandidateNodes.push_back(*itNode);
            }

            if(mCandidateNodes.size() == 0) KRATOS_ERROR(std::invalid_argument,"No nodes found for periodic conditions generation","");

            // Initialize Search Strategy
            SizeType BucketSize = 20;
            mpSearchStrategy = boost::shared_ptr<tree>(new tree(mCandidateNodes.begin(),mCandidateNodes.end(),BucketSize));

//            std::cout << *mpSearchStrategy << std::endl;

            KRATOS_CATCH("")
        }

        /// Generate a set of conditions linking each node in the periodic boundary to its image on the other side.
        template< class TReference >
        void GenerateConditions(const TReference& MovementRef,
                                Properties::Pointer pProperties,
                                const double Tolerance = 1e-4)
        {
            KRATOS_TRY

            // check that the spatial seach structure was initialized
            if(mpSearchStrategy == 0)
                KRATOS_ERROR(std::logic_error,"PeriodicConditionUtilities error: GenerateConditions() called without a spatial search structure. Please call SetUpSearchStructure() first.","")

            // Get reference condition
            const Condition& rCondition = KratosComponents<Condition>::Get("PeriodicCondition");

            // Create a geometry for new conditions
            ModelPart::ConditionsContainerType& rConditions = mrModelPart.Conditions();
            Geometry< Node<3> >::PointsArrayType ConditionNodes(2);

            // Get Id for new conditions
            IndexType Id;
            if( mrModelPart.ConditionsArray().size() == 0)
                Id = 0;
            else
                Id = (mrModelPart.ConditionsEnd()-1)->Id();

//            // Spatial search setup (for SearchInRadius)
//            const double Distance = sqrt(Translation[0]*Translation[0] + Translation[1]*Translation[1] + Translation[2]*Translation[2]);
//            const double SearchRadius = Tolerance * Distance;
//
//            PointType ImageNode;
//            SizeType NodesFound;
//            SizeType MaxResults = 20;
//
//            PointIterator Results;
//            DistanceIterator Distances;

            // Spatial search setup (for SearchNearestPoint)
            double ResultDistance = 0.0;
            PointTypePointer SecondNode;
            PointType ImageNode;
            
            // Spatial search, create new conditions with results
            unsigned int ConditionCount = 0;
            for(PointVector::iterator itNode = mCandidateNodes.begin(); itNode != mCandidateNodes.end(); itNode++)
            {
                this->MoveNode(**itNode,ImageNode,MovementRef); // * for iterator + * for pointer
//                NodesFound = mpSearchStrategy->SearchInRadius(ImageNode,SearchRadius,Results,MaxResults);
//                if (NodesFound > 0)
//                {
//                    ConditionNodes.GetContainer()[0] = *itNode;
//
//                    // Find nearest node in radius
//                    double MinDistance = Distance;
//                    for(SizeType i = 0; i < NodesFound; ++i)
//                    {
//                        if(*Distances < MinDistance)
//                        {
//                            MinDistance = *Distances;
//                            ConditionNodes.GetContainer()[1] = *Results;
//                        }
//                        Distances++;
//                        Results++;
//                    }
//
//                    // Add condition to model part
//                    Condition::Pointer pNewCondition = rCondition.Create(++Id,ConditionNodes,pProperties);
//                    rConditions.push_back(pNewCondition);
//                }

                SecondNode = mpSearchStrategy->SearchNearestPoint(ImageNode,ResultDistance);
                if(ResultDistance < Tolerance)
                {
                    ConditionNodes.GetContainer()[0] = *itNode;
                    ConditionNodes.GetContainer()[1] = SecondNode;
                    // Add condition to model part
                    Condition::Pointer pNewCondition = rCondition.Create(++Id,ConditionNodes,pProperties);
                    rConditions.push_back(pNewCondition);
                    ++ConditionCount;
                }

            }

            std::cout << "Found " << ConditionCount << " node pairs in periodic boundary." << std::endl;

            KRATOS_CATCH("")
        }


        void DefinePeriodicBoundary(const double PenaltyWeight,
                                    const double TranslationX,
                                    const double TranslationY,
                                    const double TranslationZ = 0.0)
        {
            KRATOS_TRY

             // check that the spatial seach structure was initialized
            if(mpSearchStrategy == 0)
                KRATOS_ERROR(std::logic_error,"PeriodicConditionUtilities error: DefinePeriodicBoundary() called without a spatial search structure. Please call SetUpSearchStructure() first.","")


            const double Tolerance = 1e-4; // Relative tolerance when searching for node pairs

            array_1d<double,3> Translation;
            Translation[0] = TranslationX;
            Translation[1] = TranslationY;
            Translation[2] = TranslationZ;

            Properties::Pointer pNewProperties = boost::shared_ptr<Properties>( new Properties() );
            SetPropertiesForVelocity(pNewProperties);
            SetSymmetry(pNewProperties,PenaltyWeight);

            GenerateConditions(Translation,pNewProperties,Tolerance);

            KRATOS_CATCH("")
        }

        void DefineCentralSymmetry(const double PenaltyWeight,
                                   const double CentreX,
                                   const double CentreY,
                                   const double CentreZ = 0.0)
        {
            KRATOS_TRY

             // check that the spatial seach structure was initialized
            if(mpSearchStrategy == 0)
                KRATOS_ERROR(std::logic_error,"PeriodicConditionUtilities error: DefinePeriodicBoundary() called without a spatial search structure. Please call SetUpSearchStructure() first.","")


            const double Tolerance = 1e-4; // Relative tolerance when searching for node pairs

            Node<3> Centre;
            Centre.Coordinate(0) = CentreX;
            Centre.Coordinate(1) = CentreY;
            Centre.Coordinate(2) = CentreZ;

            Properties::Pointer pNewProperties = boost::shared_ptr<Properties>( new Properties() );
            SetPropertiesForVelocity(pNewProperties);
            SetSymmetry(pNewProperties,PenaltyWeight);

            GenerateConditions(Centre,pNewProperties,Tolerance);

            KRATOS_CATCH("")
        }
        void DefineCentralAntimetry(const double PenaltyWeight,
                                    const double CentreX,
                                    const double CentreY,
                                    const double CentreZ = 0.0)
        {
            KRATOS_TRY

             // check that the spatial seach structure was initialized
            if(mpSearchStrategy == 0)
                KRATOS_ERROR(std::logic_error,"PeriodicConditionUtilities error: DefinePeriodicBoundary() called without a spatial search structure. Please call SetUpSearchStructure() first.","")


            const double Tolerance = 1e-4; // Relative tolerance when searching for node pairs

            Node<3> Centre;
            Centre.Coordinate(0) = CentreX;
            Centre.Coordinate(1) = CentreY;
            Centre.Coordinate(2) = CentreZ;

            Properties::Pointer pNewProperties = boost::shared_ptr<Properties>( new Properties() );
            SetPropertiesForVelocity(pNewProperties);
            SetAntimetry(pNewProperties,PenaltyWeight);

            GenerateConditions(Centre,pNewProperties,Tolerance);

            KRATOS_CATCH("")
        }


        void DefinePeriodicBoundaryViscosity(const double PenaltyWeight,
                                             const double TranslationX,
                                             const double TranslationY,
                                             const double TranslationZ = 0.0)
        {
            KRATOS_TRY

            if(mpSearchStrategy == 0)
                KRATOS_ERROR(std::logic_error,"PeriodicConditionUtilities error: DefinePeriodicBoundaryViscosity() called without a spatial search structure. Please call SetUpSearchStructure() first.","")


            const double Tolerance = 1e-4; // Relative tolerance when searching for node pairs

            array_1d<double,3> Translation;
            Translation[0] = TranslationX;
            Translation[1] = TranslationY;
            Translation[2] = TranslationZ;

            Properties::Pointer pNewProperties = boost::shared_ptr<Properties>( new Properties() );
            SetPropertiesForViscosity(pNewProperties);
            SetSymmetry(pNewProperties,PenaltyWeight);

            GenerateConditions(Translation,pNewProperties,Tolerance);

            KRATOS_CATCH("")
        }

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
        virtual std::string Info() const
        {
            return "PeriodicConditionUtilities";
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "PeriodicConditionUtilities";
        }

        /// Print object's data.
        virtual void PrintData(std::ostream& rOStream) const
        {}


        ///@}
        ///@name Friends
        ///@{


        ///@}

    protected:
        ///@name Protected static Member Variables
        ///@{


        ///@}
        ///@name Protected member Variables
        ///@{


        ///@}
        ///@name Protected Operators
        ///@{


        ///@}
        ///@name Protected Operations
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

        /// ModelPart containing the problem domain
        ModelPart& mrModelPart;

        /// Number of spatial dimensions
        SizeType mDomainSize;

        /// Nodes in mrModelPart identified as candidates for the periodic condition
        PointVector mCandidateNodes;

        /// Spatial search strategy
        typename tree::Pointer mpSearchStrategy;

        ///@}
        ///@name Private Operators
        ///@{


        ///@}
        ///@name Private Operations
        ///@{

        void SetSymmetry(Properties::Pointer pProperties,
                         const double PenaltyWeight) const
        {
            pProperties->GetValue(PERIODIC_VARIABLES).SetSymmetricCondition(PenaltyWeight);
        }

        void SetAntimetry(Properties::Pointer pProperties,
                          const double PenaltyWeight) const
        {
            pProperties->GetValue(PERIODIC_VARIABLES).SetAntimetricCondition(PenaltyWeight);
        }


        void SetPropertiesForVelocity(Properties::Pointer pProperties) const
        {
            pProperties->GetValue(PERIODIC_VARIABLES) = PeriodicVariablesContainer();
            PeriodicVariablesContainer& rPeriodicVariables = pProperties->GetValue(PERIODIC_VARIABLES);

            rPeriodicVariables.Add(VELOCITY_X);
            rPeriodicVariables.Add(VELOCITY_Y);
            if(mDomainSize == 3)
                rPeriodicVariables.Add(VELOCITY_Z);
//            rPeriodicVariables.Add(PRESSURE);
        }

        void SetPropertiesForViscosity(Properties::Pointer pProperties) const
        {
            pProperties->GetValue(PERIODIC_VARIABLES) = PeriodicVariablesContainer();
            PeriodicVariablesContainer& rPeriodicVariables = pProperties->GetValue(PERIODIC_VARIABLES);

            rPeriodicVariables.Add(TURBULENT_VISCOSITY);
        }

        void MoveNode( const Node<3>& rInNode, Node<3>& rOutNode, const array_1d<double,3>& rTranslation) const
        {
            rOutNode.Coordinates() = rInNode.Coordinates() + rTranslation;
        }

        void MoveNode( const Node<3>& rInNode, Node<3>& rOutNode, const Node<3>& rCentreNode )
        {
            rOutNode.Coordinates() = rInNode.Coordinates() + 2.0 * ( rCentreNode.Coordinates() - rInNode.Coordinates() );
        }

        ///@}
        ///@name Private  Access
        ///@{


        ///@}
        ///@name Private Inquiry
        ///@{


        ///@}
        ///@name Un accessible methods
        ///@{



        ///@}

    }; // Class PeriodicConditionUtilities

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function
    inline std::istream& operator >> (std::istream& rIStream,
                                      PeriodicConditionUtilities& rThis)
    {
        return rIStream;
    }

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream,
                                      const PeriodicConditionUtilities& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}

    ///@} //Group

}  // namespace Kratos.

#endif	/* KRATOS_PERIODIC_CONDITION_UTILITIES_H */

