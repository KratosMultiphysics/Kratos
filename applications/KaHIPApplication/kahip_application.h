//     __ __      __  __________  ___                ___            __  _           
//    / //_/___ _/ / / /  _/ __ \/   |  ____  ____  / (_)________ _/ /_(_)___  ____ 
//   / ,< / __ `/ /_/ // // /_/ / /| | / __ \/ __ \/ / / ___/ __ `/ __/ / __ \/ __ \ 
//  / /| / /_/ / __  // // ____/ ___ |/ /_/ / /_/ / / / /__/ /_/ / /_/ / / /_/ / / /
// /_/ |_\__,_/_/ /_/___/_/   /_/  |_/ .___/ .___/_/_/\___/\__,_/\__/_/\____/_/ /_/ 
//                                  /_/   /_/                                       
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "custom_modeler/kahip_partitioning_modeler.h"

namespace Kratos
{

///@addtogroup KaHIPApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class KratosKaHIPApplication
 * @ingroup KaHIPApplication
 * @brief KaHIP-based mesh partitioning application for Kratos Multiphysics.
 * @details This application provides graph-based mesh partitioning using KaHIP
 *          (Karlsruhe High Quality Partitioning) as an alternative to MetisApplication.
 *          KaHIP typically produces lower edge cuts than METIS, especially with the
 *          ECO and STRONG preconfigurations.
 *
 *          Supported workflows:
 *          - Serial partitioning via KaHIPDivideHeterogeneousInputProcess
 *          - Parallel (MPI) distributed partitioning via ParHIP (when built with MPI)
 *
 *          Python API mirrors MetisApplication for drop-in compatibility.
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KAHIP_APPLICATION) KratosKaHIPApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosKaHIPApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosKaHIPApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosKaHIPApplication();

    /// Destructor.
    ~KratosKaHIPApplication() override {}

    ///@}
    ///@name Operations
    ///@{

    void Register() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "KratosKaHIPApplication";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        KRATOS_WATCH("in KratosKaHIPApplication");
        KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size());
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

private:
    ///@name Private Operators
    ///@{

    /// Assignment operator.
    KratosKaHIPApplication& operator=(KratosKaHIPApplication const& rOther) = delete;

    /// Copy constructor.
    KratosKaHIPApplication(KratosKaHIPApplication const& rOther) = delete;

    /// Prototype modeler instance used for factory registration.
    const KaHIPPartitioningModeler mKaHIPPartitioningModeler;

    ///@}

}; // Class KratosKaHIPApplication

///@}

} // namespace Kratos
