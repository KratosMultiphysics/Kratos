// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_BASE_LOAD_CONDITION_3D_H_INCLUDED )
#define  KRATOS_BASE_LOAD_CONDITION_3D_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

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
    
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION)  BaseLoadCondition
    : public Condition
{
public:

    ///@name Type Definitions
    ///@{
    
    // Counted pointer of BaseLoadCondition
    KRATOS_CLASS_POINTER_DEFINITION( BaseLoadCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    BaseLoadCondition()
    {};

    // Constructor using an array of nodes
    BaseLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry ):Condition(NewId,pGeometry)
    {};

    // Constructor using an array of nodes with properties
    BaseLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ):Condition(NewId,pGeometry,pProperties)
    {};

    // Destructor
    virtual ~BaseLoadCondition()
    {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        if (rResult.size() != dim*NumberOfNodes)
            rResult.resize(dim*NumberOfNodes,false);

        const unsigned int pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

        if(dim == 2)
        {
            for (unsigned int i = 0; i < NumberOfNodes; ++i)
            {
                const unsigned int index = i * 2;
                rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos).EquationId();
                rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos+1).EquationId();
            }
        }
        else
        {
            for (unsigned int i = 0; i < NumberOfNodes; ++i)
            {
                const unsigned int index = i * 3;
                rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos).EquationId();
                rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos+1).EquationId();
                rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z,pos+2).EquationId();
            }
        }
        KRATOS_CATCH("")
    };


    virtual void GetDofList(
        DofsVectorType& ElementalDofList,
        ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int dim =  GetGeometry().WorkingSpaceDimension();
        ElementalDofList.resize(0);
        ElementalDofList.reserve(dim*NumberOfNodes);

        if(dim == 2)
        {
            for (unsigned int i = 0; i < NumberOfNodes; ++i)
            {
                ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            }
        }
        else
        {
            for (unsigned int i = 0; i < NumberOfNodes; ++i)
            {
                ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
            }
        }
        KRATOS_CATCH("")
    };

    virtual void GetValuesVector(
        Vector& values,
        int Step = 0 )
    {
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int MatSize = NumberOfNodes * dim;
        if (values.size() != MatSize)
            values.resize(MatSize, false);
        for (unsigned int i = 0; i < NumberOfNodes; i++)
        {
            const array_1d<double, 3 > & Displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            unsigned int index = i * dim;
            for(unsigned int k=0; k<dim; ++k)
            {
                values[index+k] = Displacement[k];
            }
        }
    }


    virtual void GetFirstDerivativesVector(
        Vector& values,
        int Step = 0 )
    {
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int MatSize = NumberOfNodes * dim;
        if (values.size() != MatSize)
        {
            values.resize(MatSize, false);
        }
        
        for (unsigned int i = 0; i < NumberOfNodes; i++)
        {
            const array_1d<double, 3 > & Velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
            unsigned int index = i * dim;
            for(unsigned int k=0; k<dim; ++k)
            {
                values[index+k] = Velocity[k];
            }
        }
    }

    virtual void GetSecondDerivativesVector(
        Vector& values,
        int Step = 0 )
    {
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int MatSize = NumberOfNodes * dim;
        if (values.size() != MatSize)
            values.resize(MatSize, false);
        for (unsigned int i = 0; i < NumberOfNodes; i++)
        {
            const array_1d<double, 3 > & Acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
            unsigned int index = i * dim;
            for(unsigned int k=0; k<dim; ++k)
            {
                values[index+k] = Acceleration[k];
            }
        }
    }

    virtual void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo )
    {
        if(rMassMatrix.size1() != 0)
        {
            rMassMatrix.resize(0, 0, false);
        }
    }

    virtual void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        ProcessInfo& rCurrentProcessInfo )
    {
        if(rDampingMatrix.size1() != 0)
        {
            rDampingMatrix.resize(0, 0, false);
        }
    }

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    virtual int Check( const ProcessInfo& rCurrentProcessInfo )
    {
        if ( DISPLACEMENT.Key() == 0 )
        {
            KRATOS_ERROR <<  "DISPLACEMENT has Key zero! (check if the application is correctly registered" << std::endl;
        }

        //verify that the dofs exist
        for ( unsigned int i = 0; i < this->GetGeometry().size(); i++ )
        {
            if ( this->GetGeometry()[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
            {
                KRATOS_ERROR << "missing variable DISPLACEMENT on node " << this->GetGeometry()[i].Id() << std::endl;
            }

            if ( this->GetGeometry()[i].HasDofFor( DISPLACEMENT_X ) == false ||
                 this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Y ) == false ||
                 this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Z ) == false )
            {
                KRATOS_ERROR << "missing one of the dofs for the variable DISPLACEMENT on node " << GetGeometry()[i].Id() << " of condition " << Id() << std::endl;
            }
        }
        return 0;
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

    ///@}
    ///@name Friends
    ///@{
    
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
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
    }

}; // class BaseLoadCondition.

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // KRATOS_BASE_LOAD_CONDITION_3D_H_INCLUDED  defined 
