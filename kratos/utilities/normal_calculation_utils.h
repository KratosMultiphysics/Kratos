//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//
//

#if !defined(KRATOS_NORMAL_CALCULATION_UTILS )
#define  KRATOS_NORMAL_CALCULATION_UTILS

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "utilities/math_utils.h"
#include "utilities/variable_utils.h"
#include "includes/deprecated_variables.h"

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
 * @namespace NormalCalculationUtils
 * @ingroup KratosCore
 * @brief Tool to evaluate the normals on nodes based on the normals of a set of surface conditions
 * @author Pooyan Dadvand
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_CORE) NormalCalculationUtils
{
public:
    ///@name Type Definitions
    ///@{

    /// The index type definition
    typedef std::size_t IndexType;

    /// The size type definition
    typedef std::size_t SizeType;

    // Node definitions
    typedef Node<3> NodeType;

    /// Definition of geometries
    typedef Geometry<NodeType> GeometryType;

    /// Condition type definition
    typedef ModelPart::ConditionType ConditionType;

    /// Conditions array definition
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

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
     * @brief It computes the normal in the conditions
     * @param rModelPart The model part to compute
     */
    template<class TContainerType>
    void CalculateNormalsInContainer(ModelPart& rModelPart);

    /**
     * @brief It computes the mean of the normal in the entities and in all the nodes
     * @param rModelPart The model part to compute
     * @param EnforceGenericGeometryAlgorithm If enforce the generic algorithm for any kind of geometry
     * @tparam TEntity The entity type considered
     */
    template<class TEntity>
    void CalculateNormals(
        ModelPart& rModelPart,
        const bool EnforceGenericGeometryAlgorithm = false
        );

    /**
     * @brief It computes the mean of the normal in the entities and in all the nodes (unit normal version)
     * @param rModelPart The model part to compute
     * @param EnforceGenericGeometryAlgorithm If enforce the generic algorithm for any kind of geometry
     * @tparam TEntity The entity type considered
     */
    template<class TEntity>
    void CalculateUnitNormals(
        ModelPart& rModelPart,
        const bool EnforceGenericGeometryAlgorithm = false
        );

    /**
     * @brief Calculates the "area normal" (vector oriented as the normal with a dimension proportional to the area).
     * @details This is done on the base of the Conditions provided which should be understood as the surface elements of the area of interest.
     * @param rConditions A set of conditions defining the "skin" of a model
     * @param Dimension Spatial dimension (2 or 3)
     * @note This function is not recommended for distributed (MPI) runs, as the user has to ensure that the calculated normals are assembled between processes. The overload of this function that takes a ModelPart is preferable in ths case, as it performs the required communication.
     */
    void CalculateOnSimplex(
        ConditionsArrayType& rConditions,
        const std::size_t Dimension
        );

    /**
     * @brief Calculates nodal area normal shape sensitivities w.r.t. nodal coordinates of the condition.
     *
     * @param rConditions   List of conditions where shape sensitivities need to be calculated.
     * @param Dimension     Dimensionality of the conditions
     */
    void CalculateNormalShapeDerivativesOnSimplex(
        ConditionsArrayType& rConditions,
        const std::size_t Dimension
    );

    /**
     * @brief Calculates the area normal (vector oriented as the normal with a dimension proportional to the area).
     * @details This is done on the base of the Conditions provided which should be  understood as the surface elements of the area of interest.
     * @param rModelPart ModelPart of the problem. Must have a set of conditions defining the "skin" of the domain
     * @param Dimension Spatial dimension (2 or 3)
     * @note Use this fuction instead of its overload taking a Conditions array for MPI applications, as it will take care of communication between partitions.
     */
    void CalculateOnSimplex(
        ModelPart& rModelPart,
        const std::size_t Dimension
        );

    /**
     * @brief Calculates the area normal (vector oriented as the normal with a dimension proportional to the area).
     * @details This is done on the base of the Conditions provided which should be  understood as the surface elements of the area of interest.
     * @param rModelPart ModelPart of the problem. Must have a set of conditions defining the "skin" of the domain
     * @note Use this fuction instead of its overload taking a Conditions array for MPI applications, as it will take care of communication between partitions.
     */
    void CalculateOnSimplex(
        ModelPart& rModelPart
        );

    /**
     * @brief This function swaps the normal of all of the conditions in a model part
     * @details This is done by swapping the two first nodes in the geometry and is thus appropriate for simplicial elements
     * @param rModelPart ModelPart of the problem. Must have a set of conditions defining the "skin" of the domain
     */
    void SwapNormals(ModelPart& rModelPart);

    /**
     * @brief Calculates the area normal (vector oriented as the normal with a dimension proportional to the area) using only nodes marked with a flag variable.
     * @details This function is equivalent to other implementations of CalculateOnSimplex, but instead of using all conditions in the array, it only uses those that contain a value of rVariable != Zero. This is useful in problems where a part of the boundary is a slip condition, as it provides more reasonable values for the normals on the border between this area and other parts of the boundary. This function is safe to use in MPI.
     * @param rModelPart ModelPart of the problem. Must have a set of conditions defining the "skin" of the domain.
     * @param Dimension Spatial dimension (2 or 3).
     * @param rVariable The Kratos::Variable used to indicate which parts of the boundary will be used to calculate the normals.
     * @param Zero The 'off' value for the flag. Conditions where rVariable == Zero will be skipped for normal calculation.
     */
    template<class TValueType>
    void CalculateOnSimplex(
        ModelPart& rModelPart,
        const std::size_t Dimension,
        Variable<TValueType>& rVariable,
        const TValueType Zero
        )
    {
        KRATOS_TRY;

        // Reset normals
        const array_1d<double,3> ZeroNormal(3,0.0);

        for(ModelPart::NodesContainerType::iterator it = rModelPart.NodesBegin(); it !=rModelPart.NodesEnd(); it++) {
            noalias(it->FastGetSolutionStepValue(NORMAL)) = ZeroNormal;
        }

        // Calculate new condition normals, using only conditions with rVariable == rValue
        array_1d<double,3> An(3,0.0);

        if ( Dimension == 2 ) {
            for ( ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); ++itCond ) {
                if ( itCond->GetValue(rVariable) != Zero )
                    CalculateNormal2D(*itCond,An);
            }
        } else if ( Dimension == 3 ) {
            array_1d<double,3> v1(3,0.0);
            array_1d<double,3> v2(3,0.0);

            for ( ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); ++itCond ) {
                if ( itCond->GetValue(rVariable) != Zero )
                    CalculateNormal3D(*itCond,An,v1,v2);
            }
        }

        // Transfer normals to nodes
        for ( ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); ++itCond ) {
            Condition::GeometryType& rGeom = itCond->GetGeometry();
            const double Coef = 1.0 / rGeom.PointsNumber();
            const array_1d<double,3>& rNormal = itCond->GetValue(NORMAL);
            for ( Condition::GeometryType::iterator itNode = rGeom.begin(); itNode != rGeom.end(); ++itNode)
                noalias( itNode->FastGetSolutionStepValue(NORMAL) ) += rNormal * Coef;
        }

        // For MPI: correct values on partition boundaries
        rModelPart.GetCommunicator().AssembleCurrentData(NORMAL);

        KRATOS_CATCH("");
    }

    ///
    /**
     * @brief Calculates the area normal (vector oriented as the normal with a dimension proportional to the area) using only nodes marked with a flag variable.
     * @details This function is equivalent to other implementations of CalculateOnSimplex, but instead of using all conditions in the array, it only uses those that contain a value of rVariable != Zero. This is useful in problems where a part of the boundary is a slip condition, as it provides more reasonable values for the normals on the border between this area and other parts of the boundary. This function is safe to use in MPI.
     * @param rModelPart ModelPart of the problem. Must have a set of conditions defining the "skin" of the domain.
     * @param Dimension Spatial dimension (2 or 3).
     * @param rVariable The Kratos::Variable used to indicate which parts of the boundary will be used to calculate the normals. Conditions where rVariable == Zero will be skipped.
     */
    template<class TValueType>
    void CalculateOnSimplex(
        ModelPart& rModelPart,
        const std::size_t Dimension,
        Variable<TValueType>& rVariable
        )
    {
        CalculateOnSimplex(rModelPart,Dimension,rVariable,TValueType());
    }

    /**
     *  @brief Calculates the area normal (vector oriented as the normal with a dimension proportional to the area) using only nodes marked with a flag variable and detecting corners. Corners are defined as nodes that recieves more than 2 normals from their neighbor conditions with a difference in angle greater than Alpha .
     * @details This function is equivalent to other implementations of CalculateOnSimplex, but instead of using all conditions in the array, it only uses those that contain a value of rVariable != Zero. This is useful in problems where a part of the boundary is a slip condition, as it provides more reasonable values for the normals on the border between this area and other parts of the boundary. This function is safe to use in MPI.
     * @param rModelPart ModelPart of the problem. Must have a set of conditions defining the "skin" of the domain.
     * @param Dimension Spatial dimension (2 or 3).
     * @param rVariable The Kratos::Variable used to indicate which parts of the boundary will be used to calculate the normals. Conditions where rVariable == Zero will be skipped.
     * @param rAlpha the maximum angle to distinguish normals.
     */
    template<class TValueType>
    void CalculateOnSimplex(
        ModelPart& rModelPart,
        const std::size_t Dimension,
        Variable<TValueType>& rVariable,
        const TValueType Zero,
        const double rAlpha
        )
    {
        KRATOS_TRY;

        // Reset normals
        const array_1d<double,3> ZeroNormal(3,0.0);

        for(ModelPart::NodesContainerType::iterator it =  rModelPart.NodesBegin(); it !=rModelPart.NodesEnd(); it++) {
            noalias(it->FastGetSolutionStepValue(NORMAL)) = ZeroNormal;
            it->FastGetSolutionStepValue(NODAL_PAUX) = 0.0;
        }

        // Calculate new condition normals, using only conditions with rVariable == rValue
        array_1d<double,3> An(3,0.0);

        if ( Dimension == 2 ) {
            for ( ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); ++itCond ) {
                if ( itCond->GetValue(rVariable) != Zero )
                    CalculateNormal2D(*itCond,An);
            }
        } else if ( Dimension == 3 ) {
            array_1d<double,3> v1(3,0.0);
            array_1d<double,3> v2(3,0.0);

            for ( ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); ++itCond ) {
                if ( itCond->GetValue(rVariable) != Zero )
                    CalculateNormal3D(*itCond,An,v1,v2);
            }
        }

        // Loop over nodes to set normals
        for(ModelPart::NodesContainerType::iterator it =  rModelPart.NodesBegin(); it !=rModelPart.NodesEnd(); it++) {
            std::vector< array_1d<double,3> > N_Mat;
            N_Mat.reserve(10);
            double nodal_area = 0.0;

            GlobalPointersVector<Condition >& ng_cond = it->GetValue(NEIGHBOUR_CONDITIONS);

            if(ng_cond.size() != 0) {
                for(GlobalPointersVector<Condition >::iterator ic = ng_cond.begin(); ic!=ng_cond.end(); ic++) {
                    Condition::GeometryType& pGeom = ic->GetGeometry();
                    const array_1d<double,3>&  rNormal = ic->GetValue(NORMAL);
                    const double Coef = 1.0 / pGeom.PointsNumber();
                    double norm_normal = norm_2( rNormal );

                    if(norm_normal != 0.0) {
                        nodal_area += Coef * norm_normal;

                        if(N_Mat.size() == 0.0) {
                            N_Mat.push_back( rNormal * Coef );
                        } else {
                            int added = 0;
                            for(unsigned int ii=0; ii<N_Mat.size();++ii) {
                                const array_1d<double,3>& temp_normal = N_Mat[ii];
                                double norm_temp = norm_2( temp_normal );

                                double cos_alpha=temp_normal[0]*rNormal[0] + temp_normal[1]*rNormal[1] +temp_normal[2]*rNormal[2];
                                cos_alpha /= (norm_temp*norm_normal);

                                if( cos_alpha > std::cos(0.017453293*rAlpha)) {
                                    N_Mat[ii] += rNormal * Coef;
                                    added = 1;
                                }
                            }

                            if(!added)
                                N_Mat.push_back( rNormal*Coef );
                        }
                    }
                }
            }

            // Compute NORMAL and mark
            array_1d<double,3> sum_Normal(3,0.0);

            for(unsigned int ii=0; ii<N_Mat.size(); ++ii) {
                sum_Normal += N_Mat[ii];
            }

            noalias( it->FastGetSolutionStepValue(NORMAL) ) = sum_Normal;
            it->FastGetSolutionStepValue(NODAL_PAUX) = nodal_area;
            //assign IS_SLIP = 0 for vertices
            if(N_Mat.size() == 2) {
//                 it->SetValue(IS_SLIP,0);
                it->FastGetSolutionStepValue(IS_SLIP)=20.0;
            } else if(N_Mat.size() == 3) {
                it->FastGetSolutionStepValue(IS_SLIP)=30.0;
            } else if(N_Mat.size() == 1) {
                it->FastGetSolutionStepValue(IS_SLIP)=10.0;
            }
        }

        // For MPI: correct values on partition boundaries
        rModelPart.GetCommunicator().AssembleCurrentData(NORMAL);
        rModelPart.GetCommunicator().AssembleCurrentData(NODAL_PAUX);

        KRATOS_CATCH("");
    }

    /**
     *  @brief Calculates the area normal (vector oriented as the normal with a dimension proportional to the area) using only nodes marked with a flag variable and detecting corners. Corners are defined as nodes that recieves more than 2 normals from their neighbor conditions with a difference in angle greater than Alpha . (Low memory version)
     * @details This function is equivalent to other implementations of CalculateOnSimplex, but instead of using all conditions in the array, it only uses those that contain a value of rVariable != Zero. This is useful in problems where a part of the boundary is a slip condition, as it provides more reasonable values for the normals on the border between this area and other parts of the boundary. This function is safe to use in MPI.
     * @param rModelPart ModelPart of the problem. Must have a set of conditions defining the "skin" of the domain.
     * @param Dimension Spatial dimension (2 or 3).
     * @param rVariable The Kratos::Variable used to indicate which parts of the boundary will be used to calculate the normals. Conditions where rVariable == Zero will be skipped.
     * @param rAlpha the maximum angle to distinguish normals.
     */
    template< class TValueType >
    void CalculateOnSimplexLowMemory(
        ModelPart& rModelPart,
        int Dimension,
        Variable<TValueType>& rVariable,
        const TValueType Zero,const double rAlpha)
    {
        KRATOS_TRY;

        // Reset normals
        const array_1d<double,3> ZeroNormal(3,0.0);

        for(ModelPart::NodesContainerType::iterator it =  rModelPart.NodesBegin(); it !=rModelPart.NodesEnd(); it++) {
            noalias(it->GetValue(NORMAL)) = ZeroNormal;
            it->GetValue(NODAL_PAUX) = 0.0;
        }

        // Calculate new condition normals, using only conditions with rVariable == rValue
        array_1d<double,3> An(3,0.0);

        if ( Dimension == 2 ) {
            for ( ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); ++itCond ) {
                if ( itCond->GetValue(rVariable) != Zero )
                    CalculateNormal2D(*itCond,An);
            }
        } else if ( Dimension == 3 ) {
            array_1d<double,3> v1(3,0.0);
            array_1d<double,3> v2(3,0.0);

            for ( ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); ++itCond ) {
                if ( itCond->GetValue(rVariable) != Zero )
                    CalculateNormal3D(*itCond,An,v1,v2);
            }
        }

        // Loop over nodes to set normals
        for(ModelPart::NodesContainerType::iterator it =  rModelPart.NodesBegin(); it !=rModelPart.NodesEnd(); it++) {
            std::vector< array_1d<double,3> > N_Mat;
            N_Mat.reserve(10);
            double nodal_area = 0.0;

            GlobalPointersVector<Condition >& ng_cond = it->GetValue(NEIGHBOUR_CONDITIONS);

            if(ng_cond.size() != 0){
                for(GlobalPointersVector<Condition >::iterator ic = ng_cond.begin(); ic!=ng_cond.end(); ic++) {
                Condition::GeometryType& pGeom = ic->GetGeometry();
                const array_1d<double,3>&  rNormal = ic->GetValue(NORMAL);
                const double Coef = 1.0 / pGeom.PointsNumber();
                double norm_normal = norm_2( rNormal );

                if(norm_normal != 0.0) {
                    nodal_area += Coef * norm_normal;

                    if(N_Mat.size() == 0.0) {
                        N_Mat.push_back( rNormal * Coef );
                    } else{
                        int added = 0;
                        for(unsigned int ii=0; ii<N_Mat.size();++ii) {
                            const array_1d<double,3>& temp_normal = N_Mat[ii];
                            double norm_temp = norm_2( temp_normal );

                            double cos_alpha=temp_normal[0]*rNormal[0] + temp_normal[1]*rNormal[1] +temp_normal[2]*rNormal[2];
                            cos_alpha /= (norm_temp*norm_normal);

                            if( cos_alpha > cos(0.017453293*rAlpha) ){
                                N_Mat[ii] += rNormal * Coef;
                                added = 1;}
                            }
                            if(!added)
                                N_Mat.push_back( rNormal*Coef );

                        }
                    }
                }
            }

            // Compute NORMAL and mark
            array_1d<double,3> sum_Normal(3,0.0);

            for(unsigned int ii=0; ii<N_Mat.size(); ++ii){
                sum_Normal += N_Mat[ii];
            }

            noalias( it->FastGetSolutionStepValue(NORMAL) ) = sum_Normal;
            it->FastGetSolutionStepValue(NODAL_PAUX) = nodal_area;
            //assign IS_SLIP = 0 for vertices
            if(N_Mat.size() == 2){
//                 it->SetValue(IS_SLIP,0);
                it->FastGetSolutionStepValue(IS_SLIP)=20.0;
            } else if(N_Mat.size() == 3) {
                it->FastGetSolutionStepValue(IS_SLIP)=30.0;
            } else if(N_Mat.size() == 1) {
                it->FastGetSolutionStepValue(IS_SLIP)=10.0;

            }
        }

        // For MPI: correct values on partition boundaries
        rModelPart.GetCommunicator().AssembleCurrentData(NORMAL);
        rModelPart.GetCommunicator().AssembleCurrentData(NODAL_PAUX);

        KRATOS_CATCH("");
    }

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

    /**
     * @brief It initializes the normal in the entites and in all the nodes
     * @param rModelPart The model part to compute
     * @tparam TEntity The entity type considered
     */
    template<class TEntity>
    void InitializeNormals(ModelPart& rModelPart);

    /**
     * @brief It computes the unit normals from the area normals
     * @param rModelPart The model part to compute
     */
    void ComputeUnitNormalsFromAreaNormals(ModelPart& rModelPart);

    /**
     * @brief This function adds the Contribution of one of the geometries to the corresponding nodes
     * @param rCondition Reference to the target condition
     * @param rAn Area normal
     */
    static void CalculateNormal2D(
        Condition& rCondition,
        array_1d<double,3>& rAn
        );

    /**
     * @brief Calculates 2D condition area normals shape sensitivity
     *
     * @param rCondition    Reference to the targe condition
     */
    static void CalculateNormalShapeDerivative2D(
        ConditionType& rCondition
        );

    /**
     * @brief This function adds the Contribution of one of the geometries to the corresponding nodes
     * @param rCondition Reference to the target condition
     * @param rAn Area normal
     * @param rv1 First tangent vector
     * @param rv2 Second tangent vector
     */
    static void CalculateNormal3D(
        Condition& rCondition,
        array_1d<double,3>& rAn,
        array_1d<double,3>& rv1,
        array_1d<double,3>& rv2
        );

    /**
     * @brief Calculates 3D condition area normals shape sensitivity
     *
     * @param rCondition    Reference to the targe condition
     */
    static void CalculateNormalShapeDerivative3D(
        ConditionType& rCondition
        );

    template<class TContainerType>
    TContainerType& GetContainer(ModelPart& rModelPart);

    template<class TContainerType>
    void CalculateNormalsUsingGenericAlgorithm(ModelPart& rModelPart);

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}

}; // Class NormalCalculationUtils

} // namespace Kratos

#endif /* KRATOS_NORMAL_CALCULATION_UTILS  defined */
