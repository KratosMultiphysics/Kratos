# Developers guide
## I/ Documenting your code using Doxygen
### 1) How
#### General case
Our project use **Doxygen** to automatically generate the developer documentation. If you implement a new function in **Mmg**, please, comment it and give at least its interface description (function arguments and return values).

For example a minimal documentation for the function that saves the mesh may be this one:
<!-- do not mark this as C code or Doxygen will remove the interesting part -->
```
/**
 * \param mesh  pointer to the mesh structure.
 * \param filename  pointer to the name of the file.
 * \return 0 if failed, 1 otherwise.
 *
 * Save mesh data.
 *
 * \warning you must call the \a _MMG3D_packMesh function before saving your
 * mesh.
 */
int MMG3D_saveMesh(MMG5_pMesh mesh, char *filename);
```
Additionaly, it is good practice to include text inside the routine to explain the work carried out.

You can refer to the [Doxygen documentation](http://www.stack.nl/~dimitri/doxygen/) for a description of the **Doxygen** commands.

#### API functions
Because the library header for Fortran users is automatically generated from the C header, you must add the interface of the fortran function to your documentation. Each line of this interface must begin with the `>` symbol and end with `\n` (a backslash and the letter n).

For example, if the previous function is an API function, its documentation becomes the following:
<!-- do not mark this as C code or Doxygen will remove the interesting part -->
```
/**
 * \param mesh pointer to the mesh structure.
 * \param filename pointer to the name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Save mesh data.
 *
 * \warning you must call the \a _MMG3D_packMesh function before saving your
 * mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_saveMesh(mesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 int MMG3D_saveMesh(MMG5_pMesh mesh, char *filename);
```

### 2) Where
Please, comments only your functions in the `.c` file, except for the **API** functions; these must be documentated in the appropriate `libmmg<X>.h` file (and only there).

## II/ Memory management: dynamic allocations and deallocations
We need to control the memory consumption in our applications so the memory used by dynamic allocations is counted and updated at each allocation and deallocation.

Note that with a high verbosity (at least 6), you can check that at the end of the process the memory count is 0.

To make the update of memory consumption easier, we have wrapped the `malloc`, `calloc`, `realloc` and `free` functions into macros that must be called in place of these functions.

| `C function`  | `Mmg macro`   |
|----------------------|-----------------------------------|
| `ptr = (type *) malloc(size*sizeof(type));`  |  `MMG5_SAFE_MALLOC(ptr,size,type,law);` |
| `ptr = (type *) calloc(size,sizeof(type));`  |  `MMG5_SAFE_CALLOC(ptr,size,type,law);` |
| `ptr = (type *) realloc(ptr,size*sizeof(type));`<br>`if ( high_verbosity )`<br>&nbsp;&nbsp;&nbsp;&nbsp;`printf("  ## Warning:%s:%d: %s reallocation.\n",__FILE__,__LINE__,tab_info);`| `MMG5_SAFE_REALLOC(ptr,prevSize,newSize,type,tab_name,law);` |
| `Decrease_memory_count(size); `<br>`free(ptr); ptr = NULL; `| `MMG5_DEL_MEM(mesh,ptr)`|

Note that other macros which aim to help to manage the memory have been implemented.

### 1) Allocations
To check that we have enough memory to allocate a pointer of size `siz` and to increment the memory counter, you must precede your allocation by a call to the `MMG5_ADD_MEM(mesh, siz, "tab_name", law)` macro.

For example, the following allocation:
```c
 Increase_memory_count(5*sizeof(double));
 ptr = (double *) malloc (5*sizeof(double));
 if ( !ptr ) {
    fprintf(stdout,"  ## Error: unable to allocate my array.");
    exit(EXIT_FAILURE);
 }
```

must be replaced by the following one in the **Mmg** code:
```c
MMG5_ADD_MEM(mesh,5*sizeof(double),"my array",exit(EXIT_FAILURE));
MMG5_SAFE_MALLOC(ptr,5,double,exit(EXIT_FAILURE));
```

That said, calling `exit` from a library function is not polite. You may wish to do something else in the `law` argument, for example setting a flag that tells your function to free any memory it allocated and  return a value indicating failure.

### 2) Deallocations
To decrement the memory counter, to deallocate your pointer and to leave it pointing toward `NULL`, you just need to call the `MMG5_DEL_MEM(mesh,ptr)` macro.

To deallocate the memory allocated in the previous example, instead of the following code:
```c
Decrease_memory_count(5*sizeof(double));
free(ptr);
ptr = NULL;
```
just write:
```c
MMG5_DEL_MEM(mesh,ptr);
```

## III/ Coding style
Please, use the following configuration in your editor:
  * no tabs;
  * 1 indent = 2 spaces;
  * no trailing whitespaces;
  * limit the size of your lines to 80 characters;

Besides, try to respect the following rules:
  * declaration of variables in the top of the function;
  * do not use exit(), use a return value instead;
  * do not implement void API functions;
  * the main library functions returns `MMG<X>_SUCCESS` if successful, `MMG<X>_LOWFAILURE` if the function fails but we can save the mesh, and `MMG<X>_STRONGFAILURE` if fails and we can't save a conforming mesh.
