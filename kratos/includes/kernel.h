//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

#pragma once

// System includes
#include <string>
#include <iostream>
#include <unordered_set>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

namespace Kratos {

///@name Kratos Classes
///@{

/// Kernel is in charge of synchronization the whole system of Kratos itself and its application.
/** Kernel is the first component of the Kratos to be created and then used to plug the application into Kratos.
    Kernel takes the list of variables defined in the Variables file in Kratos and then by adding each application
    synchronizes the variables of this application with its variables and add the new ones to the Kratos.
    After adding all applications its time to initialize the Kratos to assign variables key to the list of all variables
    in Kratos and all added applications. Finally the initialized variables with keys are synchronized in each
    application in time of calling InitializeApplication method for each of them.
    The sequence of using Kernel is as follow:
1. Creating the Kernel using its default constructor
2. Adding applications to Kernel using ImportApplication method
3. Initializing the Kernel using Initialize method
4. Initializing the applications using InitializeApplication method

    It is very important to perform all this step exactly in the same order as described above.

    @see ImportApplication
    @see Initialize
    @see InitializeApplication
    @see KratosApplication
*/
class KRATOS_API(KRATOS_CORE) Kernel {
   public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Kernel
    KRATOS_CLASS_POINTER_DEFINITION(Kernel);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** The default constructor creates a list of registered variables and classes in kratos_core.cpp
        by calling the RegisterKratosCore method of application class.

        @see KratosApplication

    */
    Kernel();

    Kernel(bool IsDistributedRun);

    /// Copy constructor.
    /** This constructor is empty
    */
    Kernel(Kernel const& rOther);

    /// Destructor.
    virtual ~Kernel();

    ///@}
    ///@name Operations
    ///@{

    /// Plugging an application into Kratos.
    /** This method first call the register method of the new application in order to create the
        components list of the application and then synchronizes the lists of its components with Kratos ones.
        The synchronized lists are
      - Variables
      - Elements
      - Conditions

    @param NewApplication The application to be added and synchronized
    */
    void ImportApplication(KratosApplication::Pointer pNewApplication);

    /// This method is used to print the information about the kernel
    /** This method is used to print the information about the kernel which contains the main information
        about the compiler usnsed, the release type and also number of threads and processes used to run in this execution.
    */
    void PrintInfo();

    /// To be deprecated because variables have their own hash key.
    /** The keys of Variables are not sequential anymore, so this method will be deprecated
    */
    void Initialize();

    /// Initializes and synchronizes the list of variables, elements and conditions in each application.
    /** This method gives the application the list of all variables, elements and condition which is registered
        by kratos and all other added applications.
        @see ImportApplication
        @see Initialize
    */
    void InitializeApplication(KratosApplication& NewApplication) {}

    bool IsImported(const std::string& rApplicationName) const;

    bool IsLibraryAvailable(const std::string& rLibraryName) const;

    static bool IsDistributedRun();

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;

    static std::unordered_set<std::string>& GetApplicationsList();

    /**
     * @brief Get the list of libraries available in the system
     * @return The list of libraries available in the system (stored in a set)
     */
    static std::unordered_set<std::string> GetLibraryList();

    static std::string Version();

    static std::string BuildType();

    static std::string OSName();

    static std::string PythonVersion();

    static std::string Compiler();

    static void SetPythonVersion(std::string);

    void PrintParallelismSupportInfo() const;

    /**
     * @brief This method prints all the system hardware information
     */
    void PrintSystemInfo() const;

    ///@}
   private:
    ///@name Static Member Variables
    ///@{

    KratosApplication::Pointer mpKratosCoreApplication;

    static bool         mIsDistributedRun;
    static std::string  mPyVersion;

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void RegisterKratosCore();

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    Kernel& operator=(Kernel const& rOther);

    ///@}

};  // Class Kernel

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator>>(std::istream& rIStream, Kernel& rThis);

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream, const Kernel& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.
