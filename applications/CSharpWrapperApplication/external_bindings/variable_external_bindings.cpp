#define EXPORT __declspec(dllexport)

#include "includes/variables.h"

extern "C" {
#if defined(KRATOS_COMPILED_IN_WINDOWS)

EXPORT const Kratos::Variable<double> *__stdcall Variable_GetVar1d(char *variableName) {
    return &Kratos::KratosComponents<Kratos::Variable<double >>::Get(variableName);
}

EXPORT const Kratos::Variable<Kratos::array_1d<double, 3>> *__stdcall Variable_GetVar3d(char *variableName) {
    return &Kratos::KratosComponents<Kratos::Variable<Kratos::array_1d<double, 3>>>::Get(variableName);
}

EXPORT const Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > *
Variable_GetVarComponent(char *variableName) {
    return &Kratos::KratosComponents<Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > >>::Get(
            variableName);
}

EXPORT bool __stdcall Variable_HasVar1d(char *variableName) {
    return Kratos::KratosComponents<Kratos::Variable<double>>::Has(variableName);
}

EXPORT bool __stdcall Variable_HasVar3d(char *variableName) {
    return Kratos::KratosComponents<Kratos::Variable<Kratos::array_1d<double, 3>>>::Has(variableName);
}

EXPORT bool __stdcall Variable_HasVariableComponent(char *variableName) {
    return Kratos::KratosComponents<Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > >>::Has(
            variableName);
}
#endif
}
