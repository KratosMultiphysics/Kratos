---
type: C++ Method
title: AugmentedLagrangeConstraintAssembler
resource: solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler.cpp#L167-L196
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`AugmentedLagrangeConstraintAssembler<TSparse,TDense>::AugmentedLagrangeConstraintAssembler(Parameters Settings, std::string&& rInstanceName) : Base(ConstraintImposition::AugmentedLagrange, std::move(rInstanceName)), mpImpl(new Impl{/*mConstraintIdToIndexMap: */ CSRHashMap<std::size_t,std::pair<std::size_t,std::size_t>>(), /*mMaybeTransposeRelationMatrix: */std::optional<typename TSparse::MatrixType>(), /*mpScaling*/ nullptr, /*mVerbosity: */1})`