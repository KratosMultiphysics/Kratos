# Remove fortran modules provided by this target.
FILE(REMOVE
  "my_module.mod"
  "MY_MODULE.mod"
  "CMakeFiles/myfort.dir/my_module.mod.stamp"

  "mymodule.mod"
  "MYMODULE.mod"
  "CMakeFiles/myfort.dir/mymodule.mod.stamp"
  )
