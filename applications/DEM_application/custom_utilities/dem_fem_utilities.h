#ifndef DEM_FEM_UTILITIES_H
#define DEM_FEM_UTILITIES_H

// System includes

#include <limits>
#include <iostream>
#include <iomanip>

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "utilities/openmp_utils.h"
#include "GeometryFunctions.h"
#include "utilities/timer.h"

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

namespace Kratos
{

class KRATOS_API(DEM_APPLICATION) DEMFEMUtilities {



    public:

    typedef ModelPart::ConditionType                                  ConditionType;
    typedef ModelPart::ElementsContainerType                          ElementsArrayType;
    typedef ModelPart::NodesContainerType                             NodesArrayType;
    typedef ModelPart::PropertiesType                                 PropertiesType;
    typedef WeakPointerVector<Element>                                ParticleWeakVectorType;
    typedef WeakPointerVector<Element >::iterator                     ParticleWeakIteratorType;

    KRATOS_CLASS_POINTER_DEFINITION(DEMFEMUtilities);

    /// Default constructor

    DEMFEMUtilities() {}

    /// Destructor

    virtual ~DEMFEMUtilities();

    void MoveAllMeshes(ModelPart& r_model_part, double time, double dt);
    void MoveAllMeshesUsingATable(ModelPart& r_model_part, double time, double dt);

    void CreateRigidFacesFromAllElements(ModelPart& r_model_part, PropertiesType::Pointer pProps);

    /// Turn back information as a string
    virtual std::string Info() const;

    /// Print information about this object
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data
    virtual void PrintData(std::ostream& rOStream) const;

    protected:

        DenseVector<unsigned int> mElementPartition;

    private:

        array_1d<double, 3> mInitialCenterOfMassAndMass;
        double mInitialMass;

    /// Assignment operator
    DEMFEMUtilities & operator=(DEMFEMUtilities const& rOther);

    }; // Class DEMFEMUtilities

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
// 	template<std::size_t TDim>
// 	inline std::ostream& operator << (std::ostream& rOStream)
// 	{
// 		rThis.PrintInfo(rOStream);
// 		rOStream << std::endl;
// 		rThis.PrintData(rOStream);
//
// 		return rOStream;
// 	}
///@}

} // namespace Kratos.

#endif // DEM_FEM_UTILITIES_H
