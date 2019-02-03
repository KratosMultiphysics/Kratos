//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_BASE_FACTORY_H_INCLUDED )
#define  KRATOS_BASE_FACTORY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/kratos_components.h"
#include "includes/shared_pointers.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "spaces/ublas_space.h"

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

/**
 * @class BaseFactory
 * @ingroup KratosCore
 * @brief Here we add the functions needed for the registration of class
 * @details Defines the base class factory
 * @author Vicente Mataix Ferrandiz
 * @tparam TClass The class to create the factory
 * @tparam TAuxiliarClass The auxiliar class to create the factory
 */
template<typename TClass, typename TAuxiliarClass = TClass>
class BaseFactory
{
public:
    ///@name Type Definitions
    ///@{

    /// The definition of the custom class
    typedef BaseFactory<TClass, TAuxiliarClass> FactoryType;

    /// The definition of the class
    typedef TClass ClassType;

    /// The definition of the auxiliar class
    typedef TAuxiliarClass AuxiliarClassType;

    /// Pointer definition of BaseFactory
    KRATOS_CLASS_POINTER_DEFINITION(BaseFactory);

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method checks if the linear class is registered
     * @return True if registered, false otherwise
     */
    virtual bool Has(const std::string& rSolverType)
    {
        return KratosComponents< FactoryType >::Has( rSolverType );
    }

    /**
     * @brief This method creates a new class
     * @param Settings The settings of the factory
     * @return The pointer to the class of interest
     */
    virtual typename ClassType::Pointer Create(Kratos::Parameters Settings)
    {
        const std::string& name = Settings["name"].GetString();
        KRATOS_ERROR_IF_NOT(Has( name )) << "Trying to construct a class with type name= " << name << std::endl <<
                                            "Which does not exist. The list of available options (for currently loaded applications) is: " << std::endl <<
                                            KratosComponents< FactoryType >() << std::endl;
        const auto& aux = KratosComponents< FactoryType >::Get( name );
        return aux.CreateClass(Settings);
    }

    /**
     * @brief This method creates a new class
     * @param rModel The model containing the problem
     * @param Settings The settings of the factory
     * @return The pointer to the class of interest
     */
    virtual typename ClassType::Pointer Create(Model& rModel, Kratos::Parameters Settings)
    {
        const std::string& name = Settings["name"].GetString();
        KRATOS_ERROR_IF_NOT(Has( name )) << "Trying to construct a class with type name= " << name << std::endl <<
                                            "Which does not exist. The list of available options (for currently loaded applications) is: " << std::endl <<
                                            KratosComponents< FactoryType >() << std::endl;
        const auto& aux = KratosComponents< FactoryType >::Get( name );
        return aux.CreateClass(rModel, Settings);
    }

    /**
     * @brief This method creates a new class
     * @param rModelPart The model part containing the problem
     * @param Settings The settings of the factory
     * @return The pointer to the class of interest
     */
    virtual typename ClassType::Pointer Create(ModelPart& rModelPart, Kratos::Parameters Settings)
    {
        const std::string& name = Settings["name"].GetString();
        KRATOS_ERROR_IF_NOT(Has( name )) << "Trying to construct a class with type name= " << name << std::endl <<
                                            "Which does not exist. The list of available options (for currently loaded applications) is: " << std::endl <<
                                            KratosComponents< FactoryType >() << std::endl;
        const auto& aux = KratosComponents< FactoryType >::Get( name );
        return aux.CreateClass(rModelPart, Settings);
    }

    /**
     * @brief This method creates a new class
     * @param pAuxiliarClass The pointer to the auxiliar class
     * @param Settings The settings of the factory
     * @return The pointer to the class of interest
     */
    virtual typename ClassType::Pointer Create(typename AuxiliarClassType::Pointer pAuxiliarClass, Kratos::Parameters Settings)
    {
        const std::string& name = Settings["name"].GetString();
        KRATOS_ERROR_IF_NOT(Has( name )) << "Trying to construct a class with type name= " << name << std::endl <<
                                            "Which does not exist. The list of available options (for currently loaded applications) is: " << std::endl <<
                                            KratosComponents< FactoryType >() << std::endl;
        const auto& aux = KratosComponents< FactoryType >::Get( name );
        return aux.CreateClass(pAuxiliarClass, Settings);
    }

    ///@}
protected:
    ///@name Protected Operators
    ///@{

    /**
     * @brief This method is an auxiliar method to create a new class with settings
     * @param Settings The settings of the factory
     * @return The pointer to the class of interest
     */
    virtual typename ClassType::Pointer CreateClass(Kratos::Parameters Settings)  const
    {
        KRATOS_ERROR << "Calling the base class BaseFactory" << std::endl;
    }

    /**
     * @brief This method is an auxiliar method to create a new class with Model and settings
     * @param rModel The model containing the problem
     * @param Settings The settings of the factory
     * @return The pointer to the class of interest
     */
    virtual typename ClassType::Pointer CreateClass(Model& rModel, Kratos::Parameters Settings)  const
    {
        KRATOS_ERROR << "Calling the base class BaseFactory" << std::endl;
    }

    /**
     * @brief This method is an auxiliar method to create a new class with model part and settings
     * @param rModelPart The model part containing the problem
     * @param Settings The settings of the factory
     * @return The pointer to the class of interest
     */
    virtual typename ClassType::Pointer CreateClass(ModelPart& rModelPart, Kratos::Parameters Settings)  const
    {
        KRATOS_ERROR << "Calling the base class BaseFactory" << std::endl;
    }

    /**
     * @brief This method is an auxiliar method to create a new class with auxiliar class and settings
     * @param pAuxiliarClass The pointer to the auxiliar class
     * @param Settings The settings of the factory
     * @return The pointer to the class of interest
     */
    virtual typename ClassType::Pointer CreateClass(typename AuxiliarClassType::Pointer pAuxiliarClass, Kratos::Parameters Settings)  const
    {
        KRATOS_ERROR << "Calling the base class BaseFactory" << std::endl;
    }
};

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
template<class TClass>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const BaseFactory<TClass>& rThis)
{
    rOStream << "Factory" << std::endl;

    return rOStream;
}

///@}

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSparseSpaceType;
typedef LinearSolver<SparseSpaceType,LocalSparseSpaceType> LinearSolverType;

typedef SolvingStrategy<SparseSpaceType, LocalSparseSpaceType, LinearSolverType> SolvingStrategyType;
typedef BuilderAndSolver<SparseSpaceType, LocalSparseSpaceType, LinearSolverType> BuilderAndSolverType;
typedef Scheme<SparseSpaceType,LocalSparseSpaceType> SchemeType;
typedef ConvergenceCriteria<SparseSpaceType,LocalSparseSpaceType> ConvergenceCriteriaType;

typedef BaseFactory<SolvingStrategyType> StrategyFactoryType;
typedef BaseFactory<BuilderAndSolverType, LinearSolverType> BuilderAndSolverFactoryType;
typedef BaseFactory<SchemeType> SchemeFactoryType;
typedef BaseFactory<ConvergenceCriteriaType> ConvergenceCriteriaFactoryType;


#ifdef KRATOS_REGISTER_STRATEGY
#undef KRATOS_REGISTER_STRATEGY
#endif
#define KRATOS_REGISTER_STRATEGY(name, reference) \
    KratosComponents<StrategyFactoryType>::Add(name, reference);

#ifdef KRATOS_REGISTER_BUILDER_AND_SOLVER
#undef KRATOS_REGISTER_BUILDER_AND_SOLVER
#endif
#define KRATOS_REGISTER_BUILDER_AND_SOLVER(name, reference) \
    KratosComponents<BuilderAndSolverFactoryType>::Add(name, reference);

#ifdef KRATOS_REGISTER_SCHEME
#undef KRATOS_REGISTER_SCHEME
#endif
#define KRATOS_REGISTER_SCHEME(name, reference) \
    KratosComponents<SchemeFactoryType>::Add(name, reference);

#ifdef KRATOS_REGISTER_CONVERGENCE_CRITERIA
#undef KRATOS_REGISTER_CONVERGENCE_CRITERIA
#endif
#define KRATOS_REGISTER_CONVERGENCE_CRITERIA(name, reference) \
    KratosComponents<ConvergenceCriteriaFactoryType>::Add(name, reference);

}  // namespace Kratos.

#endif // KRATOS_BASE_FACTORY_H_INCLUDED  defined
