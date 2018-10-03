//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//


#ifndef KRATOS_POTENTIAL_WALL_CONDITION_H
#define KRATOS_POTENTIAL_WALL_CONDITION_H

// System includes
#include <string>
#include <iostream>

#include "includes/kratos_flags.h"
#include "includes/deprecated_variables.h"

// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/process_info.h"
#include "includes/cfd_variables.h"

// Application includes
#include "compressible_potential_flow_application_variables.h"

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

/// Implements a wall condition for the potential flow formulation
template< unsigned int TDim, unsigned int TNumNodes = TDim >
class PotentialWallCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PotentialWallCondition
    KRATOS_CLASS_POINTER_DEFINITION(PotentialWallCondition);

    typedef Node < 3 > NodeType;

    typedef Properties PropertiesType;

    typedef Geometry<NodeType> GeometryType;

    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    typedef VectorMap<IndexType, DataValueContainer> SolutionStepsConditionalDataContainerType;

    typedef Element::WeakPointer ElementWeakPointerType;
    
    typedef Element::Pointer ElementPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
      @param NewId Index for the new condition
      */
    PotentialWallCondition(IndexType NewId = 0):
        Condition(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    PotentialWallCondition(IndexType NewId,
                           const NodesArrayType& ThisNodes):
        Condition(NewId,ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    PotentialWallCondition(IndexType NewId,
                           GeometryType::Pointer pGeometry):
        Condition(NewId,pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    PotentialWallCondition(IndexType NewId,
                           GeometryType::Pointer pGeometry,
                           PropertiesType::Pointer pProperties):
        Condition(NewId,pGeometry,pProperties)
    {
    }

    /// Copy constructor.
    PotentialWallCondition(PotentialWallCondition const& rOther):
        Condition(rOther)
    {
    }

    /// Destructor.
    ~PotentialWallCondition() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Copy constructor
    PotentialWallCondition & operator=(PotentialWallCondition const& rOther)
    {
        Condition::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create a new PotentialWallCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Condition::Pointer(new PotentialWallCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }


    Condition::Pointer Create(IndexType NewId, Condition::GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Condition::Pointer(new PotentialWallCondition(NewId, pGeom, pProperties));
    }

    /**
     * Clones the selected element variables, creating a new one
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */

    Condition::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const override
    {
        Condition::Pointer pNewCondition = Create(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

        pNewCondition->SetData(this->GetData());
        pNewCondition->SetFlags(this->GetFlags());

        return pNewCondition;
    }

    //Find the condition's parent element.
    void Initialize() override
    {
        KRATOS_TRY;
     
        // const array_1d<double, 3> &rNormal = this->GetValue(NORMAL);

        // KRATOS_ERROR_IF(norm_2(rNormal) < std::numeric_limits<double>::epsilon())
        //     << "Error on condition -> " << this->Id() << "\n"
        //     << "NORMAL must be calculated before using this condition." << std::endl;
       
        if (mInitializeWasPerformed)
            return;

        mInitializeWasPerformed = true;

        GeometryType &rGeom = this->GetGeometry();
        WeakPointerVector<Element> ElementCandidates;
        GetElementCandidates(ElementCandidates, rGeom);

        std::vector<IndexType> NodeIds, ElementNodeIds;
        GetSortedIds(NodeIds, rGeom);
        FindParentElement(NodeIds, ElementNodeIds, ElementCandidates);

        KRATOS_ERROR_IF(!mpElement.lock()) << "error in condition # " << this->Id() << "\n"
                                           << "Condition cannot find parent element" << std::endl;
        KRATOS_CATCH("");
    }

    void CalculateLeftHandSide(MatrixType &rLeftHandSideMatrix,
                               ProcessInfo &rCurrentProcessInfo) override
    {
        VectorType RHS;
        this->CalculateLocalSystem(rLeftHandSideMatrix, RHS, rCurrentProcessInfo);
    }

    /// Calculate wall stress term for all nodes with IS_STRUCTURE != 0.0
    /**
      @param rDampingMatrix Left-hand side matrix
      @param rRightHandSideVector Right-hand side vector
      @param rCurrentProcessInfo ProcessInfo instance (unused)
      */
    void CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
                              VectorType &rRightHandSideVector,
                              ProcessInfo &rCurrentProcessInfo) override
    {
        if (rLeftHandSideMatrix.size1() != TNumNodes)
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
        if (rRightHandSideVector.size() != TNumNodes)
            rRightHandSideVector.resize(TNumNodes, false);
        rLeftHandSideMatrix.clear();
        
        array_1d<double, 3> An;
        if (TDim == 2)
            CalculateNormal2D(An);
        else
            CalculateNormal3D(An);

        const PotentialWallCondition &r_this = *this;
        const array_1d<double, 3> &v = r_this.GetValue(VELOCITY_INFINITY);
        const double value = -inner_prod(v, An) / static_cast<double>(TNumNodes);

        for (unsigned int i = 0; i < TNumNodes; ++i)
            rRightHandSideVector[i] = value;
    }

    /// Check that all data required by this condition is available and reasonable
    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        int Check = Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area > 0

        if (Check != 0)
        {
            return Check;
        }
        else
        {
            // Check that all required variables have been registered
            if(POSITIVE_FACE_PRESSURE.Key() == 0)
                KRATOS_ERROR << "POSITIVE_FACE_PRESSURE Key is 0. Check if the application was correctly registered.";
            if(NEGATIVE_FACE_PRESSURE.Key() == 0)
                KRATOS_ERROR << "NEGATIVE_FACE_PRESSURE Key is 0. Check if the application was correctly registered.";

            // Checks on nodes

            // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
            for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
            {

                if(this->GetGeometry()[i].SolutionStepsDataHas(POSITIVE_FACE_PRESSURE) == false)
                    KRATOS_ERROR << "missing POSITIVE_FACE_PRESSURE variable on solution step data for node " << this->GetGeometry()[i].Id();
                if(this->GetGeometry()[i].SolutionStepsDataHas(NEGATIVE_FACE_PRESSURE) == false)
                    KRATOS_ERROR << "missing NEGATIVE_FACE_PRESSURE variable on solution step data for node " << this->GetGeometry()[i].Id();


                return Check;
            }
        }
        return 0;

            KRATOS_CATCH("");
        }

        /// Provides the global indices for each one of this element's local rows.
        /** This determines the elemental equation ID vector for all elemental DOFs
         * @param rResult A vector containing the global Id of each row
         * @param rCurrentProcessInfo the current process info object (unused)
         */
        void EquationIdVector(EquationIdVectorType& rResult,
                                      ProcessInfo& rCurrentProcessInfo) override
        {
            if (rResult.size() != TNumNodes)
                rResult.resize(TNumNodes, false);

            for (unsigned int i = 0; i < TNumNodes; i++)
                rResult[i] = GetGeometry()[i].GetDof(POSITIVE_FACE_PRESSURE).EquationId();
        }


        /// Returns a list of the element's Dofs
        /**
         * @param ElementalDofList the list of DOFs
         * @param rCurrentProcessInfo the current process info instance
         */
        void GetDofList(DofsVectorType& ConditionDofList,
                                ProcessInfo& CurrentProcessInfo) override
        {
            if (ConditionDofList.size() != TNumNodes)
                ConditionDofList.resize(TNumNodes);

            for (unsigned int i = 0; i < TNumNodes; i++)
                ConditionDofList[i] = GetGeometry()[i].pGetDof(POSITIVE_FACE_PRESSURE);

        }

        void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override
        {
            std::vector<double> rValues;
            ElementPointerType pElem = pGetElement();
            pElem->GetValueOnIntegrationPoints(PRESSURE, rValues, rCurrentProcessInfo);
            this->SetValue(PRESSURE,rValues[0]);
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
        std::string Info() const override
        {
            std::stringstream buffer;
            this->PrintInfo(buffer);
            return buffer.str();
        }

        /// Print information about this object.
        void PrintInfo(std::ostream& rOStream) const override
        {
            rOStream << "PotentialWallCondition" << TDim << "D #" << this->Id();
        }

        /// Print object's data.
        void PrintData(std::ostream& rOStream) const override
        {
            this->pGetGeometry()->PrintData(rOStream);
        }


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

        inline ElementPointerType pGetElement()
        {
            KRATOS_ERROR_IF_NOT(mpElement.lock() != 0)
                << "No element found for condition #" << this->Id() << std::endl;
            return mpElement.lock();
        }

        void GetElementCandidates(WeakPointerVector<Element> &ElementCandidates, GeometryType &rGeom)
        {
            for (SizeType i = 0; i < TDim; i++)
            {
                WeakPointerVector<Element> &rNodeElementCandidates = rGeom[i].GetValue(NEIGHBOUR_ELEMENTS);
                for (SizeType j = 0; j < rNodeElementCandidates.size(); j++)
                    ElementCandidates.push_back(rNodeElementCandidates(j));
            }
        }

        void GetSortedIds(std::vector<IndexType> &Ids,
                          const GeometryType &rGeom)
        {
            Ids.resize(rGeom.PointsNumber());
            for (SizeType i = 0; i < Ids.size(); i++)
                Ids[i] = rGeom[i].Id();
            std::sort(Ids.begin(), Ids.end());
        }

        void FindParentElement(std::vector<IndexType> &NodeIds,
                               std::vector<IndexType> &ElementNodeIds,
                               WeakPointerVector<Element> ElementCandidates)
        {
            for (SizeType i = 0; i < ElementCandidates.size(); i++)
            {
                GeometryType &rElemGeom = ElementCandidates[i].GetGeometry();
                GetSortedIds(ElementNodeIds, rElemGeom);

                if (std::includes(ElementNodeIds.begin(), ElementNodeIds.end(), NodeIds.begin(), NodeIds.end()))
                {
                    mpElement = ElementCandidates(i);
                    return;
                }
            }
        }

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

        bool mInitializeWasPerformed = false;
        double mMinEdgeLength;
        ElementWeakPointerType mpElement;

        void CalculateNormal2D(array_1d<double, 3> &An)
        {
            Geometry<Node<3>> &pGeometry = this->GetGeometry();

            An[0] = pGeometry[1].Y() - pGeometry[0].Y();
            An[1] = -(pGeometry[1].X() - pGeometry[0].X());
            An[2] = 0.00;
        }

        void CalculateNormal3D(array_1d<double, 3> &An)
        {
            Geometry<Node<3>> &pGeometry = this->GetGeometry();

            array_1d<double, 3> v1, v2;
            v1[0] = pGeometry[1].X() - pGeometry[0].X();
            v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
            v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

            v2[0] = pGeometry[2].X() - pGeometry[0].X();
            v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
            v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

            MathUtils<double>::CrossProduct(An, v1, v2);
            An *= 0.5;
        }

        ///@}
        ///@name Serialization
        ///@{

        friend class Serializer;

        void save(Serializer& rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
        }

        void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
        }

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


        ///@}

    }; // Class PotentialWallCondition


    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function
    template< unsigned int TDim, unsigned int TNumNodes >
    inline std::istream& operator >> (std::istream& rIStream,
                                      PotentialWallCondition<TDim,TNumNodes>& rThis)
    {
        return rIStream;
    }

    /// output stream function
    template< unsigned int TDim, unsigned int TNumNodes >
    inline std::ostream& operator << (std::ostream& rOStream,
                                      const PotentialWallCondition<TDim,TNumNodes>& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }

    ///@}

    ///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_POTENTIAL_WALL_CONDITION_H
