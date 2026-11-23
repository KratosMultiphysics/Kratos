---
type: C++ Method
title: Exception
resource: sources/exception.cpp#L43-L46
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`Exception::Exception(const Exception& Other) : std::exception(Other), mMessage(Other.mMessage), mWhat(Other.mWhat), mCallStack(Other.mCallStack)`