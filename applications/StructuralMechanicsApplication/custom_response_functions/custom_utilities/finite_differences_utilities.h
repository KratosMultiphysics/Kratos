// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder 
//   


#ifndef FINITE_DIFFERENCES_UTILITIES_H
#define FINITE_DIFFERENCES_UTILITIES_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include <boost/python.hpp>
#include <boost/numeric/ublas/io.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "processes/process.h"
#include "includes/kratos_flags.h"

#include "includes/element.h"
#include "includes/condition.h"
#include "includes/process_info.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "structural_mechanics_application_variables.h"

// ==============================================================================

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

/// Short class definition.
/** Detail class definition.

*/

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) FiniteDifferencesUtilities
{
public:
    ///@name Type Definitions
    ///@{

    typedef array_1d<double,3> array_3d;

    /// Pointer definition of FiniteDifferencesUtilities
    KRATOS_CLASS_POINTER_DEFINITION(FiniteDifferencesUtilities);



    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FiniteDifferencesUtilities()
    {
    }

    /// Destructor.
    virtual ~FiniteDifferencesUtilities()
    {
    }

    virtual FiniteDifferencesUtilities::Pointer Clone() const
    {
      KRATOS_ERROR << "Called the virtual function for Clone" << std::endl;
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    void SetDesignVariable(std::string _DesignVariable) { mDesignVariable = _DesignVariable; }
    std::string GetDesignVariable() { return mDesignVariable; }

    void SetDerivedObject(std::string _DerivedObject) { mDerivedObject = _DerivedObject; }
    std::string GetDerivedObject() { return mDerivedObject; }

    void DisturbElementDesignVariable(Element& rTracedElement,std::string VariableLabel , double DisturbanceMeasure )
    {
        const Variable<double> & rDesignVariable =
                 KratosComponents<Variable<double> >::Get(VariableLabel);

        if ( rTracedElement.GetProperties().Has(rDesignVariable) )
        {
            mElementID = rTracedElement.Id();
            // Save properties and its pointer
            Properties& r_global_property= rTracedElement.GetProperties();
            mpGlobalProperties = rTracedElement.pGetProperties();

            // Create new property and assign it to the element
            Properties::Pointer p_local_property(new Properties(r_global_property));
            rTracedElement.SetProperties(p_local_property);


            // Disturb the design variable
            const double current_property_value = rTracedElement.GetProperties()[rDesignVariable];
            p_local_property->SetValue(rDesignVariable, (current_property_value + DisturbanceMeasure));
        }
        else
            KRATOS_ERROR << "Chosen Design Variable not availible!!" << std::endl;

    }

    void UndisturbElementDesignVariable(Element& rTracedElement)
    {
        if ( rTracedElement.Id() != mElementID  )
        {
            KRATOS_ERROR << "Undisturbing failed!" << std::endl;
        }
        else
        {
            rTracedElement.SetProperties(mpGlobalProperties);
        }

    }

    Vector GetStressResultantBeam(Element& rTracedElement, std::string location, std::string stress_label,
                                    const ProcessInfo& rCurrentProcessInfo )
    {
        Vector output_vector;
        
        if(location == "STRESS_ON_GP" || location == "STRESS_ON_NODE")
        {
            std::string traced_stress_type = stress_label;

            const char item_1 = traced_stress_type.at(0);
            const char item_2 = traced_stress_type.at(1);

            int direction_1 = 0;
            std::vector< array_1d<double, 3 > > stress_vector;

            if(item_1 == 'M')
                rTracedElement.GetValueOnIntegrationPoints(MOMENT, stress_vector, rCurrentProcessInfo);
            else if(item_1 == 'F')
                rTracedElement.GetValueOnIntegrationPoints(FORCE, stress_vector, rCurrentProcessInfo);
            else
                KRATOS_ERROR << "Invalid stress type! " << traced_stress_type << (" is not supported!")  << std::endl;

            if(item_2 == 'X')
                direction_1 = 0;
            else if(item_2 == 'Y')
                direction_1 = 1;
            else if(item_2 == 'Z')
                direction_1 = 2;
            else
                KRATOS_ERROR << "Invalid stress type! " << traced_stress_type << (" is not supported!")  << std::endl;

            if(location == "STRESS_ON_GP")
            {
                const unsigned int&  GP_num = rTracedElement.GetGeometry().IntegrationPointsNumber(Kratos::GeometryData::GI_GAUSS_3);

                output_vector.resize(GP_num);
                output_vector.clear();
                for(unsigned int i = 0; i < GP_num ; i++)
                {
                    output_vector(i) = stress_vector[i][direction_1];
                }
            }
            else if(location == "STRESS_ON_NODE")
            {
                output_vector.resize(2);
                output_vector(0) = 2 * stress_vector[0][direction_1] - stress_vector[1][direction_1];
                output_vector(1) = 2 * stress_vector[2][direction_1] - stress_vector[1][direction_1];
            }
        }
        return output_vector;

    }

    Vector GetStressResultantShell(Element& rTracedElement, std::string location, std::string stress_label,
                                    const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY;

        Vector output_vector;
        if(location == "STRESS_ON_GP")
        {
            std::string traced_stress_type = stress_label;

            const char item_1 = traced_stress_type.at(0);
            const char item_2 = traced_stress_type.at(1);
            const char item_3 = traced_stress_type.at(2);
            int direction_1 = 0;
            int direction_2 = 0;
            std::vector<Matrix> stress_vector;

            if(item_1 == 'M')
                rTracedElement.GetValueOnIntegrationPoints(SHELL_MOMENT_GLOBAL, stress_vector, rCurrentProcessInfo);
            else if(item_1 == 'F')
                rTracedElement.GetValueOnIntegrationPoints(SHELL_FORCE_GLOBAL, stress_vector, rCurrentProcessInfo);
            else
                KRATOS_ERROR << "Invalid stress type! " << traced_stress_type << (" is not supported!")  << std::endl;

            if(item_2 == 'X')
                direction_1 = 0;
            else if(item_2 == 'Y')
                direction_1 = 1;
            else if(item_2 == 'Z')
                direction_1 = 2;
            else
                KRATOS_ERROR << "Invalid stress type! " << traced_stress_type << (" is not supported!")  << std::endl;

            if(item_3 == 'X')
                direction_2 = 0;
            else if(item_3 == 'Y')
                direction_2 = 1;
            else if(item_3 == 'Z')
                direction_2 = 2;
            else
                KRATOS_ERROR << "Invalid stress type! " << traced_stress_type << (" is not supported!")  << std::endl;

            unsigned int num_GP = stress_vector.size();
            output_vector.resize(num_GP);
            for(size_t i = 0; i < num_GP; i++)
            {
                output_vector(i) = stress_vector[i](direction_1, direction_2);
            }
        }
        else
        {
            output_vector.resize(1);
            output_vector.clear();
        }

        return output_vector;

        KRATOS_CATCH("")
    }

    double GetNodalDisplacement(const int NodeId, std::string TracedDofLabel, ModelPart& rModelPart )
    {
        KRATOS_TRY;

        typedef Node<3>::Pointer PointTypePointer;
        PointTypePointer traced_pNode = rModelPart.pGetNode(NodeId);

        typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> VariableComponentType;
        const VariableComponentType& rTRACED_DOF =
            KratosComponents<VariableComponentType>::Get(TracedDofLabel);

        double displacement_value = traced_pNode->FastGetSolutionStepValue(rTRACED_DOF);

        return displacement_value;

        KRATOS_CATCH("");
    }

    double GetStrainEnergy( ModelPart& rModelPart )
    {
        KRATOS_TRY;

        ProcessInfo &r_current_process_info = rModelPart.GetProcessInfo();
        double strain_energy = 0.0;

        // Sum all elemental strain energy values calculated as: W_e = u_e^T K_e u_e
        for (ModelPart::ElementIterator elem_i = rModelPart.ElementsBegin(); elem_i != rModelPart.ElementsEnd(); ++elem_i)
        {
            Matrix LHS;
            Vector RHS;
            Vector u;

            // Get state solution relevant for energy calculation
            elem_i->GetValuesVector(u,0);

            elem_i->CalculateLocalSystem(LHS,RHS,r_current_process_info);

            // Compute strain energy
            strain_energy += 0.5 * inner_prod(u,prod(LHS,u));
        }

        return strain_energy;


        KRATOS_CATCH("");
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "FiniteDifferencesUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "FiniteDifferencesUtilities";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
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

    std::string mDesignVariable;
    std::string mDerivedObject;
    unsigned int mElementID;
    Properties::Pointer mpGlobalProperties;

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
//      FiniteDifferencesUtilities& operator=(FiniteDifferencesUtilities const& rOther);

    /// Copy constructor.
//      FiniteDifferencesUtilities(FiniteDifferencesUtilities const& rOther);


    ///@}


    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        rSerializer.save("mDesignVariable",mDesignVariable);
        rSerializer.save("mDerivedObject",mDerivedObject);
    }

    virtual void load( Serializer& rSerializer )
    {
        rSerializer.load("mDesignVariable",mDesignVariable);
        rSerializer.load("mDerivedObject",mDerivedObject);
    }


}; // Class FiniteDifferencesUtilities

///@}

///@name Type Definitions
///@{

/**
* Definition of FiniteDifferencesUtilities variable
*/
//KRATOS_DEFINE_VARIABLE_IMPLEMENTATION( STRUCTURAL_MECHANICS_APPLICATION, FiniteDifferencesUtilities::Pointer, FINITE_DIFFERENCE_INFORMATION )


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // FINITE_DIFFERENCES_UTILITIES_H
