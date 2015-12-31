// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



#if !defined(KRATOS_KRATOS_APPLICATION_H_INCLUDED )
#define  KRATOS_KRATOS_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "includes/kratos_components.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/periodic_condition.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// This class defines the interface with kernel for all applications in Kratos.
/** The application class defines the interface necessary for providing the information
    needed by Kernel in order to configure the whole sistem correctly.

*/

class KRATOS_API(KRATOS_CORE) KratosApplication
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosApplication);


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosApplication();


    /// Copy constructor.
    KratosApplication(KratosApplication const& rOther) :
        mpVariableData(rOther.mpVariableData),
        mpIntVariables(rOther.mpIntVariables),
        mpUnsignedIntVariables(rOther.mpUnsignedIntVariables),
        mpDoubleVariables(rOther.mpDoubleVariables),
        mpArray1DVariables(rOther.mpArray1DVariables),
        mpVectorVariables(rOther.mpVectorVariables),
        mpMatrixVariables(rOther.mpMatrixVariables),
        mpArray1DVariableComponents(rOther.mpArray1DVariableComponents),
        mpElements(rOther.mpElements),
        mpConditions(rOther.mpConditions) {}



    /// Destructor.
    virtual ~KratosApplication() {}



    ///@}
    ///@name Operations
    ///@{

    virtual void Register()

    {

        RegisterVariables();

    }


    void RegisterVariables();
    
    ///////////////////////////////////////////////////////////////////
    void RegisterDeprecatedVariables(); //TODO: remove, this variables should not be there
    void RegisterC2CVariables(); //TODO: move to application
    void RegisterCFDVariables(); //TODO: move to application
    void RegisterDEMVariables(); //TODO: move to application
    void RegisterLegacyStructuralAppVariables(); //TODO: move to application

    ///@}
    ///@name Access
    ///@{



//	template<class TComponentType>
//		typename KratosComponents<TComponentType>::ComponentsContainerType& GetComponents(TComponentType const& rComponentType)
//	{
//		return KratosComponents<TComponentType>::GetComponents();
//	}







    // I have to see why the above version is not working for multi thread ...
    // Anyway its working with these functions.Pooyan.
    KratosComponents<Variable<int> >::ComponentsContainerType& GetComponents(Variable<int> const& rComponentType)
    {
        return *mpIntVariables;
    }

    KratosComponents<Variable<unsigned int> >::ComponentsContainerType& GetComponents(Variable<unsigned int> const& rComponentType)
    {
        return *mpUnsignedIntVariables;
    }

    KratosComponents<Variable<double> >::ComponentsContainerType& GetComponents(Variable<double> const& rComponentType)
    {
        return *mpDoubleVariables;
    }

    KratosComponents<Variable<array_1d<double, 3> > >::ComponentsContainerType& GetComponents(Variable<array_1d<double, 3> >  const& rComponentType)
    {
        return *mpArray1DVariables;
    }

    KratosComponents<Variable<Vector> >::ComponentsContainerType& GetComponents(Variable<Vector> const& rComponentType)
    {
        return *mpVectorVariables;
    }

    KratosComponents<Variable<Matrix> >::ComponentsContainerType& GetComponents(Variable<Matrix>  const& rComponentType)
    {
        return *mpMatrixVariables;
    }

    KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >::ComponentsContainerType& GetComponents(VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > const& rComponentType)
    {
        return *mpArray1DVariableComponents;
    }

    KratosComponents<VariableData>::ComponentsContainerType& GetVariables()
    {
        return *mpVariableData;

    }

    KratosComponents<Element>::ComponentsContainerType& GetElements()
    {
        return *mpElements;
    }

    KratosComponents<Condition>::ComponentsContainerType& GetConditions()
    {
        return *mpConditions;
    }

    void SetComponents(KratosComponents<VariableData>::ComponentsContainerType const& VariableDataComponents)

    {
        for(KratosComponents<VariableData>::ComponentsContainerType::iterator i = mpVariableData->begin() ;

                i != mpVariableData->end() ; i++)

        {
            std::string const& variable_name = i->second->Name();
            KratosComponents<VariableData>::ComponentsContainerType::const_iterator i_variable = VariableDataComponents.find(variable_name);

            if(i_variable == VariableDataComponents.end())

                KRATOS_THROW_ERROR(std::logic_error, "This variable is not registered in Kernel : ",   *(i_variable->second));

            unsigned int variable_key = i_variable->second->Key();

            if(variable_key == 0)

                KRATOS_THROW_ERROR(std::logic_error, "This variable is not initialized in Kernel : ",   *(i_variable->second));



            //			KRATOS_WATCH(i_variable->second.get());

            //			KRATOS_WATCH(i->second.get().Key());

            //			KRATOS_WATCH(variable_key);

            i->second->SetKey(variable_key);

        }

        //			KRATOS_WATCH("!!!!!!!!!!!!!!!!!!!!! END SETTING COMPONENETS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

    }




    void SetComponents(KratosComponents<Element>::ComponentsContainerType const& ElementComponents)

    {
        // It's better to make a loop over new components and add them if they are NOT already exist in application. Or make an ERROR for incompatibility between applications.

        mpElements->insert(ElementComponents.begin(), ElementComponents.end());

    }



    void SetComponents(KratosComponents<Condition>::ComponentsContainerType const& ConditionComponents)

    {

        mpConditions->insert(ConditionComponents.begin(), ConditionComponents.end());

    }


    Serializer::RegisteredObjectsContainerType& GetRegisteredObjects()
    {
        return *mpRegisteredObjects;
    }

    Serializer::RegisteredObjectsNameContainerType& GetRegisteredObjectsName()
    {
        return *mpRegisteredObjectsName;
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

        return "KratosApplication";

    }



    /// Print information about this object.

    virtual void PrintInfo(std::ostream& rOStream) const

    {

        rOStream << Info();

    }



    /// Print object's data.

    virtual void PrintData(std::ostream& rOStream) const

    {

        rOStream << "Variables:" << std::endl;

        KratosComponents<VariableData>().PrintData(rOStream);

        rOStream << std::endl;

        rOStream << "Elements:" << std::endl;

        KratosComponents<Element>().PrintData(rOStream);

        rOStream << std::endl;

        rOStream << "Conditions:" << std::endl;

        KratosComponents<Condition>().PrintData(rOStream);

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

    //general conditions must be defined
    const Condition  mCondition;
    const Condition  mCondition3D;
    const Condition  mCondition2D;

    // Periodic Condition 
    const PeriodicCondition mPeriodicCondition;
    const PeriodicCondition mPeriodicConditionEdge;
    const PeriodicCondition mPeriodicConditionCorner;


    //general elements must be defined
    const Element  mElement;
    const Element  mElement3D4N;
    const Element  mElement2D3N;

    KratosComponents<VariableData>::ComponentsContainerType* mpVariableData;

    KratosComponents<Variable<int> >::ComponentsContainerType* mpIntVariables;

    KratosComponents<Variable<unsigned int> >::ComponentsContainerType* mpUnsignedIntVariables;

    KratosComponents<Variable<double> >::ComponentsContainerType* mpDoubleVariables;

    KratosComponents<Variable<array_1d<double, 3> > >::ComponentsContainerType* mpArray1DVariables;

    KratosComponents<Variable<Vector> >::ComponentsContainerType* mpVectorVariables;

    KratosComponents<Variable<Matrix> >::ComponentsContainerType* mpMatrixVariables;

    KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >::ComponentsContainerType* mpArray1DVariableComponents;

    KratosComponents<Element>::ComponentsContainerType* mpElements;

    KratosComponents<Condition>::ComponentsContainerType* mpConditions;

    Serializer::RegisteredObjectsContainerType* mpRegisteredObjects;

    Serializer::RegisteredObjectsNameContainerType* mpRegisteredObjectsName;





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

    KratosApplication& operator=(KratosApplication const& rOther);





    ///@}



}; // Class KratosApplication



///@}



///@name Type Definitions

///@{





///@}

///@name Input and output

///@{





/// input stream function

inline std::istream& operator >> (std::istream& rIStream,

                                  KratosApplication& rThis);



/// output stream function

inline std::ostream& operator << (std::ostream& rOStream,

                                  const KratosApplication& rThis)

{

    rThis.PrintInfo(rOStream);

    rOStream << std::endl;

    rThis.PrintData(rOStream);



    return rOStream;

}

///@}





}  // namespace Kratos.



#endif // KRATOS_KRATOS_APPLICATION_H_INCLUDED  defined 





