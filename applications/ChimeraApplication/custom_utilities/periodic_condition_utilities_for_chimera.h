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
#include "containers/periodic_variables_container.h"

#include "includes/variables.h"

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
 * This utility will try to find node pairs where one of the nodes is on one side
 * of the periodic boundary and the other is its image on the other side. For each
 * pair, a PeriodicCondition object linking them will be created and appended to the
 * ModelPart's Conditions.
 *
 * This class is used as follows:
 * - Initialize a PeriodicConditionUtilitiesForChimera, passing the ModelPart containing the nodes
 * and the domain size.
 * - Define the nodal unknowns of the problem with AddPeriodicVariable
 * - Define a spatial search strategy with SetUpSearchStrategy
 * - Define a periodic boundary with GenerateConditions
 * - Additional periodic boundaries can be defined by calling DefinePeriodicBoundary
 * again, if the new boundary is identified by the same variable and value (Otherwise,
 * first set a new search structure with the new variable and value).
 * @see PeriodicCondition
 * @note If the PeriodicCondition are already defined in the ModelPart (for example,
 * if they are already defined in the mdpa file) one can simply call AddPeriodicVariable
 * to define the system unknowns to be linked and DefinePeriodicBoundary to initialize
 * the conditions.
 */
class PeriodicConditionUtilitiesForChimera
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PeriodicConditionUtilitiesForChimera
    KRATOS_CLASS_POINTER_DEFINITION(PeriodicConditionUtilitiesForChimera);

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
    /**
      * @param ThisModelPart The problem's ModelPart
      * @param ThisDomainSize The domain size
      */
    PeriodicConditionUtilitiesForChimera(ModelPart& ThisModelPart,
                               SizeType ThisDomainSize):
        mrModelPart(ThisModelPart),
	  /*mDomainSize(ThisDomainSize),*/
        mpSearchStrategy()
    {}

    /// Destructor.
    virtual ~PeriodicConditionUtilitiesForChimera()
    {}


    ///@}
    ///@name Operators
    ///@{



    ///@}
    ///@name Operations
    ///@{

    /// Set a spatial search structure that will be used to find the periodic boundary node pairs.
    /**
      * This function generates a spatial search structure containing all nodes in the periodic
      * boundary, which will be used to find node pairs. Note that both sides of the periodic
      * boundary have to be identified with the same value of rFlagVar.
      * @param rFlagVar The Kratos (double) Variable used to identify the periodic boundary
      * @param FlagValue The value of Variable rFlagVer used to identify the periodic boundary
      */
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

        if(mCandidateNodes.size() == 0) KRATOS_THROW_ERROR(std::invalid_argument,"No nodes found for periodic conditions generation","");

        // Initialize Search Strategy
        SizeType BucketSize = 20;
        mpSearchStrategy = Kratos::shared_ptr<tree>(new tree(mCandidateNodes.begin(),mCandidateNodes.end(),BucketSize));

//            std::cout << *mpSearchStrategy << std::endl;

        KRATOS_CATCH("")
    }

    /// Generate a set of conditions linking each node in the periodic boundary to its image on the other side.
    /**
      * @param MovementRef If TReference == array_1d<double,3>, MovementRef is assumed to be the transaltion
      * vector between the two sides of the periodic boundary. If TReference == Node<3>, MovementRef is assumed
      * to be the center of symmetry for node pairs.
      * @param pProperties Pointer to the properties that will be assigned to new conditions. Note that PeriodicConditon
      * objects need to have a value for PERIODIC_VARIABLES in their properties.
      * @param rConditionLabel Label of the periodic condition to be generated.
      * @param Tolerance Spatial search tolerance. Two nodes will be considered each other's image if the distance
      * between one and the image of the other is less than this value.
      */
    template< class TReference >
    void GenerateConditions(const TReference& MovementRef,
                            Properties::Pointer pProperties,
			    const std::string& rConditionLabel,
                            const double Tolerance = 1e-4)
    {
        KRATOS_TRY

        // check that the spatial seach structure was initialized
        if(mpSearchStrategy == 0)
            KRATOS_THROW_ERROR(std::logic_error,"PeriodicConditionUtilitiesForChimera error: GenerateConditions() called without a spatial search structure. Please call SetUpSearchStructure() first.","")

        // Get reference condition
	const Condition& rCondition = KratosComponents<Condition>::Get(rConditionLabel);

        // Create a geometry for new conditions
        ModelPart::ConditionsContainerType& rConditions = mrModelPart.Conditions();
        Geometry< Node<3> >::PointsArrayType ConditionNodes(2);

        // Get Id for new conditions
        IndexType Id;
        if( mrModelPart.ConditionsArray().size() == 0)
            Id = 0;
        else
            Id = (mrModelPart.ConditionsEnd()-1)->Id();

        // Spatial search setup (for SearchNearestPoint)
        double ResultDistance = 0.0;
        PointTypePointer SecondNode;
        PointType ImageNode;

        // Spatial search, create new conditions with results
        std::size_t ConditionCount = 0;
        for(PointVector::iterator itNode = mCandidateNodes.begin(); itNode != mCandidateNodes.end(); itNode++)
        {
            this->MoveNode(**itNode,ImageNode,MovementRef); // * for iterator + * for pointer

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

        KRATOS_INFO("Found ") << ConditionCount << " node pairs in periodic boundary." << std::endl;

        KRATOS_CATCH("")
    }


    /// Find node pairs to define periodic boundary conditions.
    /** This function uses GenerateConditions to find node pairs
      * where one is the image of the other by the translation
      * defined by the arguments. The resulting conditions will enforce
      * equal values of velocity for each node pair.
      * @param rConditionLabel Label of the periodic condition to be generated.
      * @param TranslationX X component of the vector that transforms each
      * node in one side of the periodic boundary to its image in the other.
      * @param TranslationY Y component of the vector that transforms each
      * node in one side of the periodic boundary to its image in the other.
      * @param TranslationZ Z component of the vector that transforms each
      * node in one side of the periodic boundary to its image in the other.
      */
    void DefinePeriodicBoundary(Properties::Pointer pNewProperties,
				const std::string& rConditionLabel,
                                const double TranslationX,
                                const double TranslationY,
                                const double TranslationZ = 0.0)
    {
        KRATOS_TRY

        // check that the spatial seach structure was initialized
        if(mpSearchStrategy == 0)
            KRATOS_THROW_ERROR(std::logic_error,"PeriodicConditionUtilitiesForChimera error: DefinePeriodicBoundary() called without a spatial search structure. Please call SetUpSearchStructure() first.","")


            const double Tolerance = 1e-4; // Relative tolerance when searching for node pairs

        array_1d<double,3> Translation;
        Translation[0] = TranslationX;
        Translation[1] = TranslationY;
        Translation[2] = TranslationZ;

	GenerateConditions(Translation,pNewProperties,rConditionLabel,Tolerance);

        KRATOS_CATCH("")
    }

    void AddPeriodicVariable(Properties& rProperties,
                             Variable<double>& rVariable)
    {
        rProperties.GetValue(PERIODIC_VARIABLES).Add(rVariable);
    }

    void AddPeriodicVariable(Properties &rProperties,
                             VariableComponent< VectorComponentAdaptor< array_1d<double, 3> > >&rVariable)
    {
        rProperties.GetValue(PERIODIC_VARIABLES).Add(rVariable);
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
        return "PeriodicConditionUtilitiesForChimera";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "PeriodicConditionUtilitiesForChimera";
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
    //SizeType mDomainSize;

    /// Nodes in mrModelPart identified as candidates for the periodic condition
    PointVector mCandidateNodes;

    /// Spatial search strategy
    tree::Pointer mpSearchStrategy;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

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

}; // Class PeriodicConditionUtilitiesForChimera

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  PeriodicConditionUtilitiesForChimera& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const PeriodicConditionUtilitiesForChimera& rThis)
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

