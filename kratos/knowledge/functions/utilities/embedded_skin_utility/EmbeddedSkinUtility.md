---
type: C++ Function
title: EmbeddedSkinUtility
resource: utilities/embedded_skin_utility.h#L95-L104
generated:
  by: okf-rs/0.7.0
---

# Signature

`EmbeddedSkinUtility( ModelPart &rModelPart, ModelPart &rSkinModelPart, const std::string LevelSetType = "continuous", const std::vector<std::string>& InterpolatedSkinVariables = {}) : mrModelPart(rModelPart), mrSkinModelPart(rSkinModelPart), mLevelSetType(LevelSetType == "continuous" ? Continuous : Discontinuous), mrConditionPrototype(KratosComponents<Condition>::Get(this->GetConditionType())), mInterpolatedSkinVariables(InterpolatedSkinVariables)`