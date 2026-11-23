---
type: C++ Function
title: MakePRestrictionOperator
resource: solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_utilities.hpp#L305-L628
generated:
  by: okf-rs/0.7.0
---

# Signature

`void MakePRestrictionOperator(ModelPart& rModelPart, const std::size_t FineSystemSize, const PointerVectorSet<Dof<double>>& rParentIndirectDofSet, typename TUblasSparseSpace<TValue>::MatrixType& rRestrictionOperator, const VariablesList::Pointer& rpVariableList, std::vector<detail::DofData>& rDofSet, PointerVectorSet<Dof<double>>& rIndirectDofSet, std::vector<std::size_t>& rDofMap)`