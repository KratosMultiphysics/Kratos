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


#ifndef KRATOS_INCOMPRESSIBLE_ADJOINT_POTENTIAL_WALL_CONDITION_H
#define KRATOS_INCOMPRESSIBLE_ADJOINT_POTENTIAL_WALL_CONDITION_H

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
#include "utilities/geometry_utilities.h"

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
class IncompressibleAdjointPotentialWallCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of IncompressibleAdjointPotentialWallCondition
    KRATOS_CLASS_POINTER_DEFINITION(IncompressibleAdjointPotentialWallCondition);

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

    typedef Element::WeakPointer ElementWeakPointerType;
    
    typedef Element::Pointer ElementPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
      @param NewId Index for the new condition
      */
    IncompressibleAdjointPotentialWallCondition(IndexType NewId = 0):
        Condition(NewId)
    {
    }

    IncompressibleAdjointPotentialWallCondition(Condition::Pointer pPrimalCondition)
                    : Condition(pPrimalCondition->Id(), pPrimalCondition->pGetGeometry(), pPrimalCondition->pGetProperties())
                    , mpPrimalCondition(pPrimalCondition)
    {
    };


    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    IncompressibleAdjointPotentialWallCondition(IndexType NewId,
                           const NodesArrayType& ThisNodes):
        Condition(NewId,ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    IncompressibleAdjointPotentialWallCondition(IndexType NewId,
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
    IncompressibleAdjointPotentialWallCondition(IndexType NewId,
                           GeometryType::Pointer pGeometry,
                           PropertiesType::Pointer pProperties):
        Condition(NewId,pGeometry,pProperties)
    {
    }

    /// Copy constructor.
    IncompressibleAdjointPotentialWallCondition(IncompressibleAdjointPotentialWallCondition const& rOther):
        Condition(rOther)
    {
    }

    /// Destructor.
    ~IncompressibleAdjointPotentialWallCondition() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Copy constructor
    IncompressibleAdjointPotentialWallCondition & operator=(IncompressibleAdjointPotentialWallCondition const& rOther)
    {
        Condition::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create a new IncompressibleAdjointPotentialWallCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Condition::Pointer(new IncompressibleAdjointPotentialWallCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }


    Condition::Pointer Create(IndexType NewId, Condition::GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Condition::Pointer(new IncompressibleAdjointPotentialWallCondition(NewId, pGeom, pProperties));
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

    void Initialize() override
    {   
        mpPrimalCondition->Initialize();
    }

    void ResetConstitutiveLaw() override
    {
        mpPrimalCondition->ResetConstitutiveLaw();
    }

    void CleanMemory() override
    {
        mpPrimalCondition->CleanMemory();
    }

    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->Data() = this->Data();
        mpPrimalCondition->Set(Flags(*this));
        mpPrimalCondition->InitializeSolutionStep(rCurrentProcessInfo);
    }
    //Find the condition's parent element.
    void GetValuesVector(Vector& rValues, int Step=0) override
    {
    
        KRATOS_TRY
  
        if(rValues.size() != TNumNodes)
            rValues.resize(TNumNodes, false);
        for (unsigned int i = 0; i < TNumNodes; i++)
            rValues[i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_POSITIVE_POTENTIAL);
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
            if(ADJOINT_POSITIVE_POTENTIAL.Key() == 0)
                KRATOS_ERROR << "ADJOINT_POSITIVE_POTENTIAL Key is 0. Check if the application was correctly registered.";
            if(ADJOINT_NEGATIVE_POTENTIAL.Key() == 0)
                KRATOS_ERROR << "ADJOINT_NEGATIVE_POTENTIAL Key is 0. Check if the application was correctly registered.";

            // Checks on nodes

            // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
            for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
            {

                if(this->GetGeometry()[i].SolutionStepsDataHas(ADJOINT_POSITIVE_POTENTIAL) == false)
                    KRATOS_ERROR << "missing ADJOINT_POSITIVE_POTENTIAL variable on solution step data for node " << this->GetGeometry()[i].Id();
                if(this->GetGeometry()[i].SolutionStepsDataHas(ADJOINT_NEGATIVE_POTENTIAL) == false)
                    KRATOS_ERROR << "missing ADJOINT_NEGATIVE_POTENTIAL variable on solution step data for node " << this->GetGeometry()[i].Id();


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
                rResult[i] = GetGeometry()[i].GetDof(ADJOINT_POSITIVE_POTENTIAL).EquationId();
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
                ConditionDofList[i] = GetGeometry()[i].pGetDof(ADJOINT_POSITIVE_POTENTIAL);
        }

        void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override
        {
            mpPrimalCondition -> FinalizeSolutionStep(rCurrentProcessInfo);
        }

        void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                         std::vector<double>& rValues,
                         const ProcessInfo& rCurrentProcessInfo) override
        {
            mpPrimalCondition ->GetValueOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
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
            rOStream << "IncompressibleAdjointPotentialWallCondition" << TDim << "D #" << this->Id();
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
        Condition::Pointer mpPrimalCondition;
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

    }; // Class IncompressibleAdjointPotentialWallCondition


    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function
    template< unsigned int TDim, unsigned int TNumNodes >
    inline std::istream& operator >> (std::istream& rIStream,
                                      IncompressibleAdjointPotentialWallCondition<TDim,TNumNodes>& rThis)
    {
        return rIStream;
    }

    /// output stream function
    template< unsigned int TDim, unsigned int TNumNodes >
    inline std::ostream& operator << (std::ostream& rOStream,
                                      const IncompressibleAdjointPotentialWallCondition<TDim,TNumNodes>& rThis)
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
