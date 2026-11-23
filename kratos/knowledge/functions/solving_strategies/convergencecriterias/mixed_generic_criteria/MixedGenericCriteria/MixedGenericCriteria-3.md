---
type: C++ Method
title: MixedGenericCriteria
resource: solving_strategies/convergencecriterias/mixed_generic_criteria.h#L92-L133
generated:
  by: okf-rs/0.7.0
---

# Signature

`MixedGenericCriteria(const ConvergenceVariableListType& rConvergenceVariablesList) : BaseType() , mVariableSize([&] (const ConvergenceVariableListType& rList) -> int {return rList.size();} (rConvergenceVariablesList)) , mVariableDataVector([&] (const ConvergenceVariableListType& rList) -> std::vector<const VariableData*> { int i = 0; std::vector<const VariableData*> aux_vect(mVariableSize); for (const auto &r_tup : rList) { aux_vect[i++] = std::get<0>(r_tup); } return aux_vect; } (rConvergenceVariablesList)) , mRatioToleranceVector([&] (const ConvergenceVariableListType& rList) -> std::vector<TDataType> { int i = 0; std::vector<TDataType> aux_vect(mVariableSize); for (const auto &r_tup : rList) { aux_vect[i++] = std::get<1>(r_tup); } return aux_vect; } (rConvergenceVariablesList)) , mAbsToleranceVector([&] (const ConvergenceVariableListType& rList) -> std::vector<TDataType> { int i = 0; std::vector<TDataType> aux_vect(mVariableSize); for (const auto &r_tup : rList) { aux_vect[i++] = std::get<2>(r_tup); } return aux_vect; } (rConvergenceVariablesList)) , mLocalKeyMap([&] (const ConvergenceVariableListType& rList) -> std::unordered_map<KeyType, KeyType> { KeyType local_key = 0; std::unordered_map<KeyType, KeyType> aux_map; for (const auto &r_tup : rList) { const auto *p_var_data = std::get<0>(r_tup); if (aux_map.find(p_var_data->Key()) != aux_map.end()) { KRATOS_ERROR << "Convergence variable " << p_var_data->Name() << " is repeated. Check the input convergence variable list." << std::endl; } else { KRATOS_ERROR_IF(p_var_data->IsComponent()) << "Trying to check convergence with the " << p_var_data->Name() << " component variable. Use the corresponding vector one." << std::endl; aux_map[p_var_data->Key()] = local_key++; } } return aux_map; } (rConvergenceVariablesList))`