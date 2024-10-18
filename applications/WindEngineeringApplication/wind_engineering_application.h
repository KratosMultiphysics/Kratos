//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

#ifndef KRATOS_WIND_ENGINEERING_APPLICATION_H
#define KRATOS_WIND_ENGINEERING_APPLICATION_H

// Project includes
#include "includes/kratos_application.h"


namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(WIND_ENGINEERING_APPLICATION) KratosWindEngineeringApplication : public KratosApplication {
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(KratosWindEngineeringApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosWindEngineeringApplication();

    /// Destructor.
    ~KratosWindEngineeringApplication() override {}

    ///@}
    ///@name Operations
    ///@{

    void Register() override;

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override;

    void PrintInfo(std::ostream& rOStream) const override;

    void PrintData(std::ostream& rOStream) const override;

    ///@}
private:
    ///@name Un accessible methods
    ///@{

    KratosWindEngineeringApplication& operator=(const KratosWindEngineeringApplication& rOther) = delete;

    KratosWindEngineeringApplication(const KratosWindEngineeringApplication& rOther) = delete;

    ///@}

}; // class KratosWindEngineeringApplication

///@}


} // namespace Kratos

#endif // KRATOS_WIND_ENGINEERING_APPLICATION_H
