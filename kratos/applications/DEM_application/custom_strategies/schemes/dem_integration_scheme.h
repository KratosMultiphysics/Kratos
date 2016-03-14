//        
// Author: Miguel AngelCeligueta, maceli@cimne.upc.edu
//


#if !defined(KRATOS_DEM_INTEGRATION_SCHEME_H_INCLUDED)
#define KRATOS_DEM_INTEGRATION_SCHEME_H_INCLUDED


// System includes
#include <float.h>
#include <string>
#include <iostream> 

// External includes 

// Project includes
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) DEMIntegrationScheme {
    public:

        typedef ModelPart::NodesContainerType NodesArrayType;
        KRATOS_CLASS_POINTER_DEFINITION(DEMIntegrationScheme);

        DEMIntegrationScheme();

        virtual ~DEMIntegrationScheme();

        virtual void AddSpheresVariables(ModelPart & r_model_part);
        virtual void AddClustersVariables(ModelPart & r_model_part);
        
        void SetRotationOption(const int rotation_option);
        
        virtual void UpdateLinearDisplacementAndVelocityOfSpheres(ModelPart & rcluster_model_part);         
        virtual void Calculate(ModelPart& model_part, int StepFlag = -1);
                
        virtual void UpdateTranslationalVariables(
                            int StepFlag, 
                            const Node < 3 > & i,
                            array_1d<double, 3 >& coor,
                            array_1d<double, 3 >& displ,
                            array_1d<double, 3 >& delta_displ,
                            array_1d<double, 3 >& vel,
                            const array_1d<double, 3 >& initial_coor,
                            const array_1d<double, 3 >& force,
                            const double force_reduction_factor,
                            const double mass,
                            const double delta_t,
                            const bool Fix_vel[3]);
        virtual void CalculateTranslationalMotion(ModelPart& model_part, NodesArrayType& pNodes, int StepFlag );
        
        virtual void CalculateLocalAngularAcceleration(
                                const Node < 3 > & i,
                                const double moment_of_inertia,
                                const array_1d<double, 3 >& torque, 
                                const double moment_reduction_factor,
                                array_1d<double, 3 >& angular_acceleration);

        virtual void UpdateRotationalVariables(
                int StepFlag,
                const Node < 3 > & i,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& angular_acceleration,
                const double delta_t,
                const bool Fix_Ang_vel[3]);

        virtual void CalculateRotationalMotion(ModelPart& model_part, NodesArrayType& pNodes, int StepFlag);
        
        virtual void CalculateLocalAngularAccelerationByEulerEquations(
                                    const Node < 3 > & i,
                                    const array_1d<double, 3 >& local_angular_velocity,
                                    const array_1d<double, 3 >& moments_of_inertia,
                                    const array_1d<double, 3 >& local_torque, 
                                    const double moment_reduction_factor,
                                    array_1d<double, 3 >& local_angular_acceleration);

        virtual void CalculateRotationalMotionOfClusters(ModelPart& rcluster_model_part, int StepFlag);

        virtual std::string Info() const {
            std::stringstream buffer;
            buffer << "DEMIntegrationScheme";
            return buffer.str();
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream & rOStream) const {
            rOStream << "DEMIntegrationScheme";
        }

        /// Print object's data.

        virtual void PrintData(std::ostream & rOStream) const {
        }

        protected:

        private:
            
        bool mRotationOption;

        DEMIntegrationScheme& operator=(DEMIntegrationScheme const& rOther) {
            return *this;
        }

        /// Copy constructor.

        DEMIntegrationScheme(DEMIntegrationScheme const& rOther) {
            *this = rOther;
        }
    }; // Class DEMIntegrationScheme

    /// input stream function

    inline std::istream& operator>>(std::istream& rIStream, DEMIntegrationScheme& rThis) {
        return rIStream;
    }

    /// output stream function

    inline std::ostream& operator<<(std::ostream& rOStream, const DEMIntegrationScheme& rThis) {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);
        return rOStream;
    }

} // namespace Kratos.

#endif // KRATOS_DEM_INTEGRATION_SCHEME_H_INCLUDED  defined 

