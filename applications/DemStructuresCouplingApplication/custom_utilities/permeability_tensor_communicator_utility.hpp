/*
 * Author: Miguel Angel Celigueta
 *
 *  maceli@cimne.upc.edu
 */

#ifndef PERMEABILITY_TENSOR_COMMUNICATOR_UTILITY_H
#define PERMEABILITY_TENSOR_COMMUNICATOR_UTILITY_H

#include "includes/variables.h"
#include <limits>
#include <iostream>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "includes/define.h"
#include "includes/condition.h"
#include "includes/model_part.h"
#include "dem_structures_coupling_application_variables.h"
#include "../../PoromechanicsApplication/poromechanics_application_variables.h"
#include "../../DEMApplication/DEM_application_variables.h"
#include "utilities/binbased_fast_point_locator.h"

// Project includes
#include "spatial_containers/dem_search.h"
#include "utilities/openmp_utils.h"

// Configures
#include "../../DEMApplication/custom_utilities/discrete_particle_configure.h"
#include "../../DEMApplication/custom_utilities/geometrical_object_configure.h"
#include "../../DEMApplication/custom_utilities/node_configure.h"
#include "../../DEMApplication/custom_utilities/omp_dem_search.h"

// Search
#include "spatial_containers/bins_dynamic_objects.h"
#include "spatial_containers/bins_dynamic.h"
#include "custom_search/bins_dynamic_objects_periodic.h"

namespace Kratos {

    class PermeabilityTensorCommunicatorUtility {

        public:

        typedef ModelPart::NodesContainerType::ContainerType::iterator  NodesIteratorType;
        typedef SpatialSearch                                           SearchType;
        typedef SearchType::ElementsContainerType                       ElementsContainerType;
        typedef SearchType::ElementsContainerType::ContainerType        ContainerType;
        typedef SearchType::ElementsContainerType                       ElementsContainerType;
        typedef SearchType::NodesContainerType                          NodesContainerType;
        typedef SearchType::ElementType                                 ElementType;
        typedef ContainerType::value_type                               PointerType;
        typedef ContainerType::iterator                                 IteratorType;
        typedef ElementsContainerType::iterator                         ElementIteratorType;
        typedef SearchType::ElementsContainerType::ContainerType        ResultContainerType;
        typedef DiscreteParticleConfigure<3>                  ElementConfigureType;
        typedef BinsObjectDynamic<ElementConfigureType>               BinsType;
        typedef std::unique_ptr<BinsType>                             BinsUniquePointerType;
        typedef ElementsContainerType::ContainerType              ResultElementsContainerType;
        typedef SpatialSearch::DistanceType DistanceType;

        KRATOS_CLASS_POINTER_DEFINITION(PermeabilityTensorCommunicatorUtility);

        PermeabilityTensorCommunicatorUtility(ModelPart& r_source_model_part, ModelPart& destination_model_part):mrDEMModelPart(r_source_model_part), mrFEMModelPart(destination_model_part) {
            Check();
        }

        virtual ~PermeabilityTensorCommunicatorUtility() {
            delete mpSearchStructure;
        }

        void Initialize() {
            mpSearchStructure = new BinBasedFastPointLocator<2>(mrDEMModelPart);
            mpSearchStructure->UpdateSearchDatabase();
        }

        void Diagonalize(const BoundedMatrix<double, 3, 3> A, BoundedMatrix<double, 3, 3>& Q, BoundedMatrix<double, 3, 3>& D)
        {
            // Slightly modified version of  Stan Melax's code for 3x3 matrix diagonalization
            // source: http://www.melax.com/diag.html?attredirects=0
            // A must be a symmetric matrix.
            // returns Q and D such that 
            // Diagonal matrix D = QT * A * Q;  and  A = Q * D * QT
            const int maxsteps = 24;  // certainly wont need that many.
            int k0, k1, k2;
            double o[3], m[3];
            double q [4] = {0.0, 0.0, 0.0, 1.0};
            double jr[4];
            double sqw, sqx, sqy, sqz;
            double tmp1, tmp2, mq;
            BoundedMatrix<double, 3, 3> AQ = ZeroMatrix(3,3);
            double thet, sgn, t, c;
            for (int i = 0; i < maxsteps; ++i) {
                // quat to matrix
                sqx      = q[0]*q[0];
                sqy      = q[1]*q[1];
                sqz      = q[2]*q[2];
                sqw      = q[3]*q[3];
                Q(0,0) = ( sqx - sqy - sqz + sqw);
                Q(1,1) = (-sqx + sqy - sqz + sqw);
                Q(2,2) = (-sqx - sqy + sqz + sqw);
                tmp1     = q[0]*q[1];
                tmp2     = q[2]*q[3];
                Q(1,0) = 2.0 * (tmp1 + tmp2);
                Q(0,1) = 2.0 * (tmp1 - tmp2);
                tmp1     = q[0]*q[2];
                tmp2     = q[1]*q[3];
                Q(2,0) = 2.0 * (tmp1 - tmp2);
                Q(0,2) = 2.0 * (tmp1 + tmp2);
                tmp1     = q[1]*q[2];
                tmp2     = q[0]*q[3];
                Q(2,1)  = 2.0 * (tmp1 + tmp2);
                Q(1,2)  = 2.0 * (tmp1 - tmp2);

                AQ = prod(A, Q);
                BoundedMatrix<double, 3, 3> QT = Transpose(Q);
                D = prod(QT, AQ);

                o[0]    = D(1,2);
                o[1]    = D(0,2);
                o[2]    = D(0,1);
                m[0]    = fabs(o[0]);
                m[1]    = fabs(o[1]);
                m[2]    = fabs(o[2]);
                k0      = (m[0] > m[1] && m[0] > m[2])? 0: (m[1] > m[2])? 1 : 2; // index of largest element of offdiag
                k1      = (k0+1)%3;
                k2      = (k0+2)%3;

                if (!o[k0])
                {
                    break;  // diagonal already
                }
                thet    = (D(k2,k2) - D(k1,k1)) / (2.0 * o[k0]);
                sgn     = (thet > 0.0)? 1.0: -1.0;
                thet   *= sgn; // make it positive
                t       = sgn /(thet +((thet < 1.E6)? sqrt(thet * thet + 1.0): thet)) ; // sign(T)/(|T|+sqrt(T^2+1))
                c       = 1.0 / sqrt(t * t + 1.0); //  c= 1/(t^2+1) , t=s/c
                
                if (c == 1.0)
                {
                    break;  // no room for improvement - reached machine precision.
                }
                jr[0 ]  = jr[1] = jr[2] = jr[3] = 0.0;
                jr[k0]  = sgn * sqrt((1.0 - c) / 2.0);  // using 1/2 angle identity sin(a/2) = sqrt((1-cos(a))/2)  
                jr[k0] *= -1.0; // since our quat-to-matrix convention was for v*M instead of M*v
                jr[3 ]  = sqrt(1.0f - jr[k0] * jr[k0]);

                if (jr[3] == 1.0)
                {
                    break; // reached limits of floating point precision
                }
                q[0]  = (q[3]*jr[0] + q[0]*jr[3] + q[1]*jr[2] - q[2]*jr[1]);
                q[1]  = (q[3]*jr[1] - q[0]*jr[2] + q[1]*jr[3] + q[2]*jr[0]);
                q[2]  = (q[3]*jr[2] + q[0]*jr[1] - q[1]*jr[0] + q[2]*jr[3]);
                q[3]  = (q[3]*jr[3] - q[0]*jr[0] - q[1]*jr[1] - q[2]*jr[2]);
                mq    = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
                q[0] /= mq;
                q[1] /= mq;
                q[2] /= mq;
                q[3] /= mq;
            }
        }

        BoundedMatrix<double, 3, 3> Transpose(BoundedMatrix<double, 3, 3> A) {
            BoundedMatrix<double, 3, 3> AT = ZeroMatrix(3,3);
            for (int i = 0; i < 3; i++ ) {
                for (int j = 0; j < 3; j++ ) {
                    AT(j,i) = A(i,j);
                }
            }
            return AT;
        }

        void TrasferUpdatedPermeabilityTensor() {

            KRATOS_TRY

            // This should take the diagonalised strain tensor coming from DEM and its diagonalising matrices
            // and compute the updated permeability tensor using them and the initial porosity, the initial permeability
            // and the updated porosities

            int property_id = 0;
            for (int i = 0; i < (int)mrDEMModelPart.Elements().size(); i++) {
                const auto elem_it = mrDEMModelPart.ElementsBegin() + i;
                property_id = elem_it->GetProperties().Id();
                break;
            }
            const double initial_porosity = mrDEMModelPart.GetProperties(property_id)[POROSITY];

            property_id = 0;
            for (int i = 0; i < (int)mrFEMModelPart.Elements().size(); i++) {
                const auto elem_it = mrFEMModelPart.ElementsBegin() + i;
                property_id = elem_it->GetProperties().Id();
                break;
            }
            const double initial_permeability_xx = mrFEMModelPart.GetProperties(property_id)[PERMEABILITY_XX];
            const double initial_permeability_yy = mrFEMModelPart.GetProperties(property_id)[PERMEABILITY_YY];
            const double initial_permeability_zz = mrFEMModelPart.GetProperties(property_id)[PERMEABILITY_ZZ];

            const ElementsContainerType& r_dem_elements = mrDEMModelPart.Elements();
            ElementsContainerType::ContainerType& elements_dem_model_part = const_cast<ElementsContainerType::ContainerType&>(r_dem_elements.GetContainer());
            BinsUniquePointerType p_bins = std::unique_ptr<BinsType>(new BinsType(elements_dem_model_part.begin(), elements_dem_model_part.end()));
            int max_number_of_results = r_dem_elements.size();
            ResultElementsContainerType   localResults(max_number_of_results);
            DistanceType                  localResultsDistances(max_number_of_results);

            for (int k = 0; k < (int) mrFEMModelPart.Elements().size(); k++) {

                const auto elem_it = mrFEMModelPart.ElementsBegin() + k;
                std::size_t number_of_results = 0;
                Element::GeometryType& rGeom = elem_it->GetGeometry();
                Point centroid = rGeom.Center();
                Geometry<Node<3>>::PointsArrayType node_centroid;
                node_centroid.push_back(Node<3>::Pointer(new Node<3>(1, centroid[0], centroid[1], centroid[2])));
                Element::Pointer p_particle = Kratos::make_intrusive<SphericParticle>(1, node_centroid);
                const double radius = rGeom.Length();
                number_of_results = p_bins->SearchObjectsInRadiusExclusive(p_particle, radius, localResults.begin(), localResultsDistances.begin(), max_number_of_results);
                BoundedMatrix<double, 3, 3> K_full_total = ZeroMatrix(3,3);

                for (int i = 0; i <  number_of_results; ++i) {

                    SphericParticle* particle = dynamic_cast<SphericParticle*>(&*localResults[i]);
                    const BoundedMatrix<double, 3, 3> strain_tensor = *particle->mStrainTensor;

                    // Diagonalise matrix strain_tensor, diag(strain_tensor) = D = QT * strain_tensor * Q
                    BoundedMatrix<double, 3, 3> Q = ZeroMatrix(3,3);
                    BoundedMatrix<double, 3, 3> D = ZeroMatrix(3,3);
                    Diagonalize(strain_tensor, Q, D);
                    BoundedMatrix<double, 3, 3> QT = Transpose(Q);

                    // Compute diagonal permeability matrix
                    BoundedMatrix<double, 3, 3> K_diag = ZeroMatrix(3,3);
                    K_diag(0,0) = initial_permeability_xx * (initial_porosity + D(1,1)) * (initial_porosity + D(2,2)) / (initial_porosity * initial_porosity * (1.0 + D(1,1)) * (1.0 + D(2,2)));
                    K_diag(1,1) = initial_permeability_yy * (initial_porosity + D(0,0)) * (initial_porosity + D(2,2)) / (initial_porosity * initial_porosity * (1.0 + D(0,0)) * (1.0 + D(2,2)));
                    K_diag(2,2) = initial_permeability_zz * (initial_porosity + D(0,0)) * (initial_porosity + D(1,1)) / (initial_porosity * initial_porosity * (1.0 + D(0,0)) * (1.0 + D(1,1)));
                    K_diag(1,0) = K_diag(0,1) = K_diag(2,0) = K_diag(0,2) = 0.0;

                    // De-diagonalise K_diag so K_full = Q * K_diag * QT
                    BoundedMatrix<double, 3, 3> Temp = prod(Q, K_diag);
                    BoundedMatrix<double, 3, 3> K_full = prod(Temp, QT);

                    // Sum up each DEM full permeability matrices
                    K_full_total += K_full;
                }

                // Divide obtained full permeability matrix by the number_of_results
                if (number_of_results) {
                    K_full_total /= number_of_results;
                } else {
                    continue;
                }

                // Transfer full permeability matrix to FEM
                unsigned int dim = rGeom.WorkingSpaceDimension();
                Matrix FEMPermeabilityTensor;
                FEMPermeabilityTensor.resize(dim, dim, false);
                for (unsigned int i = 0; i < dim; i++) {
                    for (unsigned int j = 0; j < dim; j++){
                        FEMPermeabilityTensor(i,j) = K_full_total(i,j);
                    }
                }
                unsigned int NumGPoints = rGeom.IntegrationPointsNumber(elem_it->GetIntegrationMethod());
                std::vector<Matrix> FEMPermeabilityTensorVector(NumGPoints);
                for (unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++) {
                    FEMPermeabilityTensorVector[GPoint] = FEMPermeabilityTensor;
                }
                const ProcessInfo& CurrentProcessInfo = mrFEMModelPart.GetProcessInfo();
                elem_it->SetValuesOnIntegrationPoints(PERMEABILITY_MATRIX, FEMPermeabilityTensorVector, CurrentProcessInfo);
            }

            KRATOS_CATCH("")
        }

        virtual std::string Info() const { return "";}
        virtual void PrintInfo(std::ostream& rOStream) const {}
        virtual void PrintData(std::ostream& rOStream) const {}

        private:

        ModelPart& mrDEMModelPart;
        ModelPart& mrFEMModelPart;
        BinBasedFastPointLocator<2>* mpSearchStructure;

        PermeabilityTensorCommunicatorUtility& operator= (PermeabilityTensorCommunicatorUtility const& rOther);

        void Check() {}

    }; // class PermeabilityTensorCommunicatorUtility

} // namespace Kratos

#endif // PERMEABILITY_TENSOR_COMMUNICATOR_UTILITY_H
