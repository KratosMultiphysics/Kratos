//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

#if !defined(KRATOS_KRATOS_APPLICATION_H_INCLUDED)
#define KRATOS_KRATOS_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "includes/kratos_components.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/periodic_condition.h"
#include "utilities/quaternion.h"
#include "includes/master_slave_constraint.h"
#include "includes/linear_master_slave_constraint.h"

namespace Kratos {
///@name Kratos Classes
///@{

/// This class defines the interface with kernel for all applications in Kratos.
/** The application class defines the interface necessary for providing the information
    needed by Kernel in order to configure the whole sistem correctly.

*/

class KRATOS_API(KRATOS_CORE) KratosApplication {
   public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosApplication(const std::string ApplicationName);

    KratosApplication() = delete;

    /// Copy constructor.
    KratosApplication(KratosApplication const& rOther)
        : mpVariableData(rOther.mpVariableData),
          mpIntVariables(rOther.mpIntVariables),
          mpUnsignedIntVariables(rOther.mpUnsignedIntVariables),
          mpDoubleVariables(rOther.mpDoubleVariables),
          mpArray1DVariables(rOther.mpArray1DVariables),
          mpArray1D4Variables(rOther.mpArray1D4Variables),
          mpArray1D6Variables(rOther.mpArray1D6Variables),
          mpArray1D9Variables(rOther.mpArray1D9Variables),
          mpVectorVariables(rOther.mpVectorVariables),
          mpMatrixVariables(rOther.mpMatrixVariables),
          mpArray1DVariableComponents(rOther.mpArray1DVariableComponents),
          mpArray1D4VariableComponents(rOther.mpArray1D4VariableComponents),
          mpArray1D6VariableComponents(rOther.mpArray1D6VariableComponents),
          mpArray1D9VariableComponents(rOther.mpArray1D9VariableComponents),
          mpElements(rOther.mpElements),
          mpConditions(rOther.mpConditions),
          mpMasterSlaveConstraints(rOther.mpMasterSlaveConstraints) {}

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
    void
    RegisterDeprecatedVariables();  //TODO: remove, this variables should not be there
    void RegisterC2CVariables();                  //TODO: move to application
    void RegisterCFDVariables();                  //TODO: move to application
    void RegisterALEVariables();                  //TODO: move to application
    void RegisterMappingVariables();              //TODO: move to application
    void RegisterDEMVariables();                  //TODO: move to application
    void RegisterFSIVariables();                  //TODO: move to application
    void RegisterMATVariables();                  //TODO: move to application
    void RegisterLegacyStructuralAppVariables();  //TODO: move to application

    const std::string& Name() const { return mApplicationName; }

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
    KratosComponents<Variable<int> >::ComponentsContainerType& GetComponents(
        Variable<int> const& rComponentType) {
        return *mpIntVariables;
    }

    KratosComponents<Variable<unsigned int> >::ComponentsContainerType&
    GetComponents(Variable<unsigned int> const& rComponentType) {
        return *mpUnsignedIntVariables;
    }

    KratosComponents<Variable<double> >::ComponentsContainerType& GetComponents(
        Variable<double> const& rComponentType) {
        return *mpDoubleVariables;
    }

    KratosComponents<Variable<array_1d<double, 3> > >::ComponentsContainerType&
    GetComponents(Variable<array_1d<double, 3> > const& rComponentType) {
        return *mpArray1DVariables;
    }

    KratosComponents<Variable<array_1d<double, 4> > >::ComponentsContainerType&
    GetComponents(Variable<array_1d<double, 4> > const& rComponentType) {
        return *mpArray1D4Variables;
    }

    KratosComponents<Variable<array_1d<double, 6> > >::ComponentsContainerType&
    GetComponents(Variable<array_1d<double, 6> > const& rComponentType) {
        return *mpArray1D6Variables;
    }

    KratosComponents<Variable<array_1d<double, 9> > >::ComponentsContainerType&
    GetComponents(Variable<array_1d<double, 9> > const& rComponentType) {
        return *mpArray1D9Variables;
    }

    KratosComponents<Variable<Quaternion<double> > >::ComponentsContainerType&
    GetComponents(Variable<Quaternion<double> > const& rComponentType) {
        return *mpQuaternionVariables;
    }

    KratosComponents<Variable<Vector> >::ComponentsContainerType& GetComponents(
        Variable<Vector> const& rComponentType) {
        return *mpVectorVariables;
    }

    KratosComponents<Variable<Matrix> >::ComponentsContainerType& GetComponents(
        Variable<Matrix> const& rComponentType) {
        return *mpMatrixVariables;
    }

    KratosComponents<VariableComponent<VectorComponentAdaptor< array_1d<double, 3> > > >::ComponentsContainerType& GetComponents(
        VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > const& rComponentType) {
        return *mpArray1DVariableComponents;
    }

    KratosComponents<VariableComponent<VectorComponentAdaptor< array_1d<double, 4> > > >::ComponentsContainerType& GetComponents(
        VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > > const& rComponentType) {
        return *mpArray1D4VariableComponents;
    }

    KratosComponents<VariableComponent<VectorComponentAdaptor< array_1d<double, 6> > > >::ComponentsContainerType& GetComponents(
        VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > > const& rComponentType) {
        return *mpArray1D6VariableComponents;
    }

    KratosComponents<VariableComponent<VectorComponentAdaptor< array_1d<double, 9> > > >::ComponentsContainerType& GetComponents(
        VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > > const& rComponentType) {
        return *mpArray1D9VariableComponents;
    }

    KratosComponents<VariableData>::ComponentsContainerType& GetVariables() {
        return *mpVariableData;
    }

    KratosComponents<Element>::ComponentsContainerType& GetElements() {
        return *mpElements;
    }

    KratosComponents<Condition>::ComponentsContainerType& GetConditions() {
        return *mpConditions;
    }

    KratosComponents<MasterSlaveConstraint>::ComponentsContainerType& GetMasterSlaveConstraints() {
        return *mpMasterSlaveConstraints;
    }

    void SetComponents(
        KratosComponents<VariableData>::ComponentsContainerType const&
            VariableDataComponents)

    {
        for (KratosComponents<VariableData>::ComponentsContainerType::iterator
                 i = mpVariableData->begin();

             i != mpVariableData->end(); i++)

        {
            std::string const& variable_name = i->second->Name();
            KratosComponents<VariableData>::ComponentsContainerType::
                const_iterator i_variable =
                    VariableDataComponents.find(variable_name);

            if (i_variable == VariableDataComponents.end())

                KRATOS_THROW_ERROR(std::logic_error,
                    "This variable is not registered in Kernel : ",
                    *(i_variable->second));

            unsigned int variable_key = i_variable->second->Key();

            if (variable_key == 0)

                KRATOS_THROW_ERROR(std::logic_error,
                    "This variable is not initialized in Kernel : ",
                    *(i_variable->second));

            //			KRATOS_WATCH(i_variable->second.get());

            //			KRATOS_WATCH(i->second.get().Key());

            //			KRATOS_WATCH(variable_key);

            i->second->SetKey(variable_key);
        }

        //			KRATOS_WATCH("!!!!!!!!!!!!!!!!!!!!! END SETTING COMPONENETS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    }

    void SetComponents(KratosComponents<Element>::ComponentsContainerType const&
            ElementComponents)

    {
        // It's better to make a loop over new components and add them if they are NOT already exist in application. Or make an ERROR for incompatibility between applications.

        mpElements->insert(ElementComponents.begin(), ElementComponents.end());
    }

    void SetComponents(KratosComponents<MasterSlaveConstraint>::ComponentsContainerType const&
            MasterSlaveConstraintComponents)

    {
        mpMasterSlaveConstraints->insert(MasterSlaveConstraintComponents.begin(), MasterSlaveConstraintComponents.end());
    }

    void SetComponents(
        KratosComponents<Condition>::ComponentsContainerType const&
            ConditionComponents)

    {
        mpConditions->insert(
            ConditionComponents.begin(), ConditionComponents.end());
    }

    Serializer::RegisteredObjectsContainerType& GetRegisteredObjects() {
        return *mpRegisteredObjects;
    }

    Serializer::RegisteredObjectsNameContainerType& GetRegisteredObjectsName() {
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

        rOStream << std::endl;

        rOStream << "MasterSlaveConstraints:" << std::endl;

        KratosComponents<MasterSlaveConstraint>().PrintData(rOStream);
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
    std::string mApplicationName;

    // General conditions must be defined

    // Point conditions
    const Condition mPointCondition2D1N;
    const Condition mPointCondition3D1N;
    // Line conditions
    const Condition mLineCondition2D2N;
    const Condition mLineCondition2D3N;
    const Condition mLineCondition3D2N;
    const Condition mLineCondition3D3N;
    // Surface conditions
    const Condition mSurfaceCondition3D3N;
    const Condition mSurfaceCondition3D6N;
    const Condition mSurfaceCondition3D4N;
    const Condition mSurfaceCondition3D8N;
    const Condition mSurfaceCondition3D9N;


    // Master-Slave base constraint
    const MasterSlaveConstraint mMasterSlaveConstraint;
    const LinearMasterSlaveConstraint mLinearMasterSlaveConstraint;

    // Deprecated conditions start
    const Condition mCondition;
    const Condition mCondition2D;
    const Condition mCondition2D2N;
    const Condition mCondition2D3N;
    const Condition mCondition3D;
    const Condition mCondition3D2N;
    const Condition mCondition3D3N;
    const Condition mCondition3D6N;
    const Condition mCondition3D4N;
    const Condition mCondition3D8N;
    const Condition mCondition3D9N;
    // Deprecated conditions end

    // Periodic Condition
    const PeriodicCondition mPeriodicCondition;
    const PeriodicCondition mPeriodicConditionEdge;
    const PeriodicCondition mPeriodicConditionCorner;

    // General elements must be defined
    const Element mElement;
    const Element mElement2D2N;
    const Element mElement2D3N;
    const Element mElement2D4N;

    const Element mElement3D2N;
    const Element mElement3D3N;
    const Element mElement3D4N;
    const Element mElement3D6N;
    const Element mElement3D8N;
    const Element mElement3D10N;



    const ConstitutiveLaw mConstitutiveLaw;

    KratosComponents<VariableData>::ComponentsContainerType* mpVariableData;

    KratosComponents<Variable<int> >::ComponentsContainerType* mpIntVariables;

    KratosComponents<Variable<unsigned int> >::ComponentsContainerType* mpUnsignedIntVariables;

    KratosComponents<Variable<double> >::ComponentsContainerType* mpDoubleVariables;

    KratosComponents<Variable<array_1d<double, 3> > >::ComponentsContainerType* mpArray1DVariables;

    KratosComponents<Variable<array_1d<double, 4> > >::ComponentsContainerType* mpArray1D4Variables;

    KratosComponents<Variable<array_1d<double, 6> > >::ComponentsContainerType* mpArray1D6Variables;

    KratosComponents<Variable<array_1d<double, 9> > >::ComponentsContainerType* mpArray1D9Variables;

    KratosComponents<Variable<Quaternion<double> > >::ComponentsContainerType* mpQuaternionVariables;

    KratosComponents<Variable<Vector> >::ComponentsContainerType* mpVectorVariables;

    KratosComponents<Variable<Matrix> >::ComponentsContainerType* mpMatrixVariables;

    KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >::ComponentsContainerType* mpArray1DVariableComponents;

    KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > > >::ComponentsContainerType* mpArray1D4VariableComponents;

    KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > > >::ComponentsContainerType* mpArray1D6VariableComponents;

    KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > > >::ComponentsContainerType* mpArray1D9VariableComponents;

    KratosComponents<Element>::ComponentsContainerType* mpElements;

    KratosComponents<Condition>::ComponentsContainerType* mpConditions;

    KratosComponents<MasterSlaveConstraint>::ComponentsContainerType* mpMasterSlaveConstraints;

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

};  // Class KratosApplication

///@}

///@name Type Definitions

///@{

///@}

///@name Input and output

///@{

/// input stream function

inline std::istream& operator>>(std::istream& rIStream,

    KratosApplication& rThis);

/// output stream function

inline std::ostream& operator<<(std::ostream& rOStream,

    const KratosApplication& rThis)

{
    rThis.PrintInfo(rOStream);

    rOStream << std::endl;

    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.

#endif  // KRATOS_KRATOS_APPLICATION_H_INCLUDED  defined
