---
type: C++ Method
title: Item
resource: utilities/profiler.h#L91-L103
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`Item(std::size_t CallCount, Duration CumulativeDuration, Duration MinDuration, Duration MaxDuration, CodeLocation&& rLocation) noexcept : mRecursionLevel(0) , mCallCount(CallCount) , mCumulative(CumulativeDuration) , mMin(MinDuration) , mMax(MaxDuration) , mLocation(std::move(rLocation))`