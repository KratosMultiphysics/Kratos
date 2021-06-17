//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Masó Sotomayor
//


#ifndef KRATOS_CONVECTION_OPERATOR_H_INCLUDED
#define KRATOS_CONVECTION_OPERATOR_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/binbased_fast_point_locator.h"


namespace Kratos
{
///@addtogroup ShallowWaterOperator
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

/// Short class definition.
/** Detail class definition.
*/
template<size_t TDim>
class ConvectionOperator
{
public:
    ///@name Type Definitions
    ///@{

    typedef Point PointType;
    typedef Geometry<Node<3>> GeometryType;

    /// Pointer definition of ConvectionOperator
    KRATOS_CLASS_POINTER_DEFINITION(ConvectionOperator);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    ConvectionOperator(ModelPart& rModelPart)
     : mrModelPart(rModelPart)
     , mSearchStructure(rModelPart)
    {
        UpdateSearchDatabase();
    }

    /// Destructor
    virtual ~ConvectionOperator(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual void SetMaxResults(size_t MaxResults)
    {
        mMaxResults = MaxResults;
    }

    /**
     * @brief This method updates the search structure if the mesh has been modified
     */
    virtual void UpdateSearchDatabase()
    {
        mSearchStructure.UpdateSearchDatabase();
    }

    virtual void InitializeSearch()
    {
        mResults = typename BinBasedFastPointLocator<TDim>::ResultContainerType(mMaxResults);
        miResultBegin = mResults.begin();
    }

    /**
     * @brief This method convects a point with RK2 and substepping
     */
    virtual bool Convect(
        const double Dt,
        PointType& rPosition,
        Element::Pointer& pElement,
        Vector& rN,
        const int Substeps = 1)
    {
        bool is_found = false;
        is_found = mSearchStructure.FindPointOnMesh(rPosition, rN, pElement, miResultBegin, mMaxResults);
        if (is_found)
        {
            GeometryType geom = pElement->GetGeometry();
            array_1d<double,3> vel = 0.5 * rN[0] * (geom[0].FastGetSolutionStepValue(VELOCITY) + geom[0].FastGetSolutionStepValue(VELOCITY,1));
            Convect(Dt, rPosition, pElement, rN, vel, Substeps);
        }
        return is_found;
    }

    /**
     * @brief This method convects a point starting from a known location with RK2 and substepping
     */
    virtual bool Convect(
        const double Dt,
        PointType& rPosition,
        Element::Pointer& pElement,
        Vector& rN,
        const array_1d<double,3>& rInitialVelocity,
        const int Substeps = 1)
    {
        bool is_found = false;
        int sub_steps = std::abs(Substeps);
        const double sub_dt = Dt / static_cast<double>(Substeps);
        for (int i = 0; i < sub_steps; ++i)
        {
            is_found = ConvectSingleStep(sub_dt, rPosition, pElement, rN, rInitialVelocity);
            if (!is_found) break;
        }
        return is_found;
    }

    /**
     * @brief This method convects a point starting from a known location with RK2
     */
    virtual bool ConvectSingleStep(
        const double Dt,
        PointType& rPosition,
        Element::Pointer& pElement,
        Vector& rN,
        const array_1d<double,3>& rInitialVelocity)
    {
        bool is_found = false;
        array_1d<double,3> pos_step1 = rPosition + 0.5 * Dt * rInitialVelocity;
        is_found = mSearchStructure.FindPointOnMesh(pos_step1, rN, pElement, miResultBegin, mMaxResults);
        if (is_found)
        {
            Geometry<Node<3>>& geom = pElement->GetGeometry();
            array_1d<double,3> vel_step1 = 0.5 * rN[0] * (geom[0].FastGetSolutionStepValue(VELOCITY) + geom[0].FastGetSolutionStepValue(VELOCITY,1));
            for (std::size_t i = 1; i < geom.size(); ++i) {
                vel_step1 += 0.5 * rN[i] * (geom[i].FastGetSolutionStepValue(VELOCITY) + geom[i].FastGetSolutionStepValue(VELOCITY,1));
            }
            rPosition += Dt * vel_step1;
            is_found = mSearchStructure.FindPointOnMesh(rPosition, rN, pElement, miResultBegin, mMaxResults);
        }
        return is_found;
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
        std::stringstream buffer;
        buffer << "ConvectionOperator" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "ConvectionOperator";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

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

    ModelPart& mrModelPart;
    BinBasedFastPointLocator<TDim> mSearchStructure;
    typename BinBasedFastPointLocator<TDim>::ResultContainerType mResults;
    typename BinBasedFastPointLocator<TDim>::ResultIteratorType miResultBegin;
    size_t mMaxResults;

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


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ConvectionOperator& operator=(ConvectionOperator const& rOther){}

    /// Copy constructor.
    ConvectionOperator(ConvectionOperator const& rOther){}

    ///@}

}; // Class ConvectionOperator

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                 ConvectionOperator& rThis){}

/// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                 const ConvectionOperator& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CONVECTION_OPERATOR_H_INCLUDED  defined
