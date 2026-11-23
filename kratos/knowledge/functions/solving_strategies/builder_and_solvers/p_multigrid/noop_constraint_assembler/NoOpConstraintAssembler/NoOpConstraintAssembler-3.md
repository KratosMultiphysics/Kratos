---
type: C++ Method
title: NoOpConstraintAssembler
resource: solving_strategies/builder_and_solvers/p_multigrid/noop_constraint_assembler.hpp#L39-L42
generated:
  by: okf-rs/0.7.0
---

# Signature

`NoOpConstraintAssembler(Parameters Settings, std::string&& rInstanceName) noexcept : Base(ConstraintImposition::None, std::move(rInstanceName)), mVerbosity(Settings.Has("verbosity") ? Settings["verbosity"].Get<int>() : 1)`