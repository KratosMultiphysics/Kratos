
#-------------------------------------------------------------------

macro( array2d_get_item out_value offset )
  math( EXPR _finalindex "${_array2d_index}+${offset}" )
  list( GET _array2d_array ${_finalindex} _item )
  set( ${out_value} "${_item}" )
endmacro()

#-------------------------------------------------------------------

macro( array2d_begin_loop out_advanced array width var_names )
  set( _array2d_out_advanced ${out_advanced} )
  set( _array2d_index 0 )
  set( _array2d_array ${array} )
  set( _array2d_width ${width} )
  set( _array2d_var_names ${var_names} )
  array2d_advance()
endmacro()

#-------------------------------------------------------------------

macro( array2d_advance )
  if( NOT _array2d_array )
    set( ${_array2d_out_advanced} false )
  else()	
    list( LENGTH _array2d_array _size )
    math( EXPR _remaining "${_size}-${_array2d_index}" )
    
    if( (_array2d_width LESS 1) OR (_size LESS _array2d_width) OR (_remaining LESS _array2d_width) )
      set( ${_array2d_out_advanced} false )
    else()
      math( EXPR _adjusted_width "${_array2d_width}-1" )
      foreach( offset RANGE ${_adjusted_width} )
	list( GET _array2d_var_names ${offset} _var_name )
	array2d_get_item( ${_var_name} ${offset} )
      endforeach()
      
      math( EXPR _index "${_array2d_index}+${_array2d_width}" )
      set( _array2d_index ${_index} )
      set( ${_array2d_out_advanced} true )
    endif()
  endif()
endmacro()

#-------------------------------------------------------------------