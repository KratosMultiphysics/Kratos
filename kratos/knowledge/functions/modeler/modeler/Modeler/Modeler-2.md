---
type: C++ Method
title: Modeler
resource: modeler/modeler.h#L80-L88
generated:
  by: okf-rs/0.7.0
---

# Signature

`explicit Modeler( Model& rModel, Parameters ModelerParameters = Parameters()) : mParameters(ModelerParameters) , mEchoLevel( ModelerParameters.Has("echo_level") ? ModelerParameters["echo_level"].GetInt() : 0)`