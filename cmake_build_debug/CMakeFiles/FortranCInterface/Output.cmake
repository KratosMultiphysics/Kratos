# Global symbol without underscore.
set(FortranCInterface_GLOBAL_SYMBOL  "mysub_")
set(FortranCInterface_GLOBAL_PREFIX  "")
set(FortranCInterface_GLOBAL_SUFFIX  "_")
set(FortranCInterface_GLOBAL_CASE    "LOWER")
set(FortranCInterface_GLOBAL_MACRO   "(name,NAME) name##_")

# Global symbol with underscore.
set(FortranCInterface_GLOBAL__SYMBOL "my_sub_")
set(FortranCInterface_GLOBAL__PREFIX "")
set(FortranCInterface_GLOBAL__SUFFIX "_")
set(FortranCInterface_GLOBAL__CASE   "LOWER")
set(FortranCInterface_GLOBAL__MACRO  "(name,NAME) name##_")

# Module symbol without underscore.
set(FortranCInterface_MODULE_SYMBOL  "__mymodule_MOD_mysub")
set(FortranCInterface_MODULE_PREFIX  "__")
set(FortranCInterface_MODULE_MIDDLE  "_MOD_")
set(FortranCInterface_MODULE_SUFFIX  "")
set(FortranCInterface_MODULE_CASE    "LOWER")
set(FortranCInterface_MODULE_MACRO   "(mod_name,name, mod_NAME,NAME) __##mod_name##_MOD_##name")

# Module symbol with underscore.
set(FortranCInterface_MODULE__SYMBOL "__my_module_MOD_my_sub")
set(FortranCInterface_MODULE__PREFIX "__")
set(FortranCInterface_MODULE__MIDDLE "_MOD_")
set(FortranCInterface_MODULE__SUFFIX "")
set(FortranCInterface_MODULE__CASE   "LOWER")
set(FortranCInterface_MODULE__MACRO  "(mod_name,name, mod_NAME,NAME) __##mod_name##_MOD_##name")

# Summarize what was found.
set(FortranCInterface_GLOBAL_FOUND 1)
set(FortranCInterface_MODULE_FOUND 1)

