# Developers guide
## I/ Documents your code using Doxygen
### 1) How
#### General case
Our project use **Doxygen** to automatically generate the developer documentation. If you implement a new function in **Mmg**, please, comment it and give at least its interface description (function's arguments and return values).  

For example a minimal documentation for the function that save the mesh may be this one:
```c
/**
 * \param mesh pointer toward the mesh structure.
 * \param filename pointer toward the name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Save mesh data.
 *
 * \warning you must call the \a _MMG3D_packMesh function before saving your
 * mesh.
 */
int MMG3D_saveMesh(MMG5_pMesh mesh, char *filename);
```
Additionaly, it is a good practice to include text inside the routine to explain the work carried out.  

You can refer to the [Doxygen documentation](http://www.stack.nl/~dimitri/doxygen/) for a description of the **Doxygen** commands.

#### API's functions
Because the library header for Fortran users is automatically generated from the C header, you must add to your documentation the interface of the fortran function. Each line of this interface must begin with the **>** symbol and end with the **\n** one.

For example, if the previous function is an API function, its documentation becames the following:
```c
/**
 * \param mesh pointer toward the mesh structure.
 * \param filename pointer toward the name of file.
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
Please, comments only your functions in the **.c** file, except for the **API**'s functions that must be documentated in the suitable **_libmmg*.h_** file (and only here).

## II/ Memory management: dynamic allocations and deallocations
We need to control the memory consumption in our applications so the memory used by dynamic allocations is counted and updated at each allocation and deallocation.  

Note that with a high verbosity (at least 6), you can check that at the end of the process the memory count is 0.

To make the update of memory consumption easier, we have wrapped the **malloc**, **calloc**, **realloc** and **free** functions into macros that must be called in place of the matching function.

| **_C function_**  | **_Mmg macro_**   |
|----------------------|-----------------------------------|
| ``ptr = (type *) malloc(size*sizeof(type));``  |  ``_MMG5_SAFE_MALLOC(ptr,size,type);`` |
| ``ptr = (type *) calloc(size,sizeof(type));``  |  ``_MMG5_SAFE_CALLOC(ptr,size,type);`` |
| ``ptr = (type *) realloc(ptr,size*sizeof(type));``<br>``if ( high_verbosity )``<br>&nbsp;&nbsp;&nbsp;&nbsp;``printf("  ## Warning:%s:%d: %s reallocation.\n",__FILE__,__LINE__,tab_info);``| ``_MMG5_SAFE_REALLOC(ptr,size,type,tab_info);`` |
| ``Decrease_memory_count(size); ``<br>``free(ptr); ptr = NULL; ``| ``_MMG5_DEL_MEM(mesh,ptr,size)``|

### 1) Allocations
To check that we have enough memory to allocate a pointer of size **_siz_** and to increment the memory counter, you must precede your allocation by a call to the **_MMG5_ADD_MEM(mesh, siz, tab_info, law)** macro.

For example, the following allocation:
```c
 Increase_memory_count(5*sizeof(double));
 ptr = (double *) malloc (5*sizeof(double));
 if ( !ptr ) {
    fprintf(stdout,"  ## Error: unable to allocate ""table of double"".");  
    exit(EXIT_FAILURE);
 }
```

must be replaced by the following one in the **Mmg** code:
```c 
_MMG5_ADD_MEM(mesh,5*sizeof(double),"table of double",exit(EXIT_FAILURE));
_MMG5_SAFE_MALLOC(ptr,5,double); 
```

### 2) Deallocations
To decrement the memory counter, to deallocate your pointer and to leave it pointing toward ``NULL``, you just need to call the ``MMG5_DEL_MEM(mesh,ptr,siz)`` macro. 

To deallocate the memory allocated in the previous example, instead of the following code:
```c
Decrease_memory_count(5*sizeof(double));
free(ptr);
ptr = NULL;
```
just write:
```c
MMG5_DEL_MEM(mesh,ptr,5*sizeof(double));
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
  * do not implement void API function;
  * the main library functions returns MMG*_SUCCESS if success, MMG*_LOWFAILURE if fails but we can save the mesh, MMG*_STRONGFAILURE if fails and we can't save a conform mesh; 
