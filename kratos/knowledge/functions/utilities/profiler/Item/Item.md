---
type: C++ Method
title: Item
resource: utilities/profiler.h#L81-L88
generated:
  by: okf-rs/0.7.0
---

# Signature

`explicit Item(CodeLocation&& rLocation) noexcept : Item(0, // <== .mCallCount Duration(0), // <== .mCumulative Duration(std::numeric_limits<typename Duration::rep>::max()), // <== .mMin Duration(0), // <== .mMax std::move(rLocation)) // <== .mLocation`