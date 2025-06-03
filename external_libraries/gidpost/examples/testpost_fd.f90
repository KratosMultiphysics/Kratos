PROGRAM EXAMPLE
  use gidpost
  implicit none
  REAL(8) :: nodes(1:9,1:3), rx, ry, rz, value
  INTEGER :: elems(1:4,1:4)
  INTEGER :: idx, x, y, nid(1:4)
  CHARACTER(LEN=1024) :: name
  CHARACTER(LEN=1024) :: component_names(2)
  TYPE(GiD_File) :: fdm, fdr
  REAL(8) :: gpa
  INTEGER :: do_write_array, num_nodes, idx_coords, num_elements, idx_con
  REAL(8) :: coords_xyz(1:27), scalar_value(1:9)
  INTEGER :: connectivities(1:16), ids(1:9)
  num_nodes = 9
  num_elements = 4

  idx = 1
  do x=0,2
    do y=0,2
      nodes(idx,1) = x
      nodes(idx,2) = y
      nodes(idx,3) = 0
      idx = idx + 1
    end do
  end do
  idx = 1
  do x=1,2
    do y=1,2
      elems(idx,1) = x + (y-1)*3
      elems(idx,2) = elems(idx,1) + 1
      elems(idx,3) = elems(idx,1) + 4
      elems(idx,4) = elems(idx,1) + 3
      idx = idx + 1
    end do
  end do

  CALL GiD_PostInit

  !CALL GID_OPENPOSTMESHFILE('testfortran90.post.msh',GiD_PostAscii)      
  !CALL GiD_OpenPostResultFile('testfortran90.post.res',GiD_PostBinary)      
  !CALL GiD_OpenPostResultFile('testfortran90.post.res',GiD_PostHDF5)  
  fdm=GiD_fOpenPostMeshFile('test_fdfortran90.post.msh',GiD_PostAscii)  
  fdr=GiD_fOpenPostResultFile('test_fdfortran90.post.res',GiD_PostAscii)
  
  !fdm=fdr
  
  CALL GiD_fBeginMeshColor(fdm,'quadmesh',GiD_2D,GiD_Quadrilateral,4,0.7d0,0.7d0,0.4d0)
  CALL GiD_fMeshUnit(fdm,'mm')

  do_write_array = 1

  if ( do_write_array == 1) then
    ! write coordinates
    idx_coords = 1
    num_nodes = 9
    do idx=1,num_nodes
        coords_xyz( idx_coords + 0) = nodes( idx, 1)
        coords_xyz( idx_coords + 1) = nodes( idx, 2)
        coords_xyz( idx_coords + 2) = nodes( idx, 3)
        idx_coords = idx_coords + 3
    end do
    CALL GiD_fWriteCoordinatesBlock( fdm, num_nodes, coords_xyz)
    ! write connectivities
    idx_con = 1
    num_elements = 4
    do idx=1,num_elements
        connectivities( idx_con + 0) = elems(idx,1)
        connectivities( idx_con + 1) = elems(idx,2)
        connectivities( idx_con + 2) = elems(idx,3)
        connectivities( idx_con + 3) = elems(idx,4)
        idx_con = idx_con + 4
    end do
    CALL GiD_fWriteElementsBlock( fdm, num_elements, connectivities)
  else
    ! write coordinates
    CALL GiD_fBeginCoordinates(fdm)
    do idx=1,9
      CALL GiD_fWriteCoordinates(fdm,idx,nodes(idx,1),nodes(idx,2),nodes(idx,3))      
    end do
    CALL GiD_fEndCoordinates(fdm)
    ! write connectivities
    CALL GiD_fBeginElements(fdm)  
    do idx=1,4
      nid(1) = elems(idx,1)
      nid(2) = elems(idx,2)
      nid(3) = elems(idx,3)
      nid(4) = elems(idx,4)
      CALL GiD_fWriteElement(fdm,idx,nid)    
    end do      
    CALL GiD_fEndElements(fdm)
  endif
  CALL GiD_fEndMesh(fdm)
  
  !Gauss points, declared but not used by any result  
  CALL GiD_fBeginGaussPoint(fdr,'element_gp',GiD_Quadrilateral,GID_STRING_NULL,4,0,0)
  gpa = 0.57735027  
  CALL GiD_fWriteGaussPoint2D(fdr,-gpa,-gpa)
  CALL GiD_fWriteGaussPoint2D(fdr,gpa,-gpa)
  CALL GiD_fWriteGaussPoint2D(fdr,0.57735027d0,0.57735027d0)
  CALL GiD_fWriteGaussPoint2D(fdr,-gpa,gpa)
  CALL GiD_fEndGaussPoint(fdr)
  
  !Scalar result
  if ( do_write_array == 1) then
    do idx=1,9
      scalar_value( idx) = idx * 1.5
      ids(idx)=idx
    end do
    CALL GiD_fWriteResultBlock( fdr, 'Pressure', 'Analysis', 1.0d0, GiD_Scalar, GiD_onNodes, GID_STRING_NULL,    &
              GID_STRING_NULL, 0, GID_STRING_NULL, 'Pa', num_nodes, ids, 1, scalar_value)
  else
    CALL GiD_fBeginResultHeader(fdr,'Pressure','Analysis',1.0d0,GiD_Scalar,GiD_onNodes,GID_STRING_NULL)  
    !CALL GiD_fScalarComp(fdr,'Comp')
    CALL GiD_fResultUnit(fdr,'Pa')
    CALL GiD_fResultValues(fdr)
    !CALL GiD_fBeginScalarResult(fdr,'Pressure','Analysis',1.0d0,GiD_onNodes,GID_STRING_NULL,GID_STRING_NULL,GID_STRING_NULL) 
    do idx=1,9
      value=idx*1.5
      CALL GiD_fWriteScalar(fdr,idx,value)
    end do
    CALL GiD_fEndResult(fdr)
  endif
  CALL GiD_fFlushPostFile(fdr)
  !Vectorial result
  CALL GiD_fBeginResultHeader(fdr,'Velocity','Analysis',1.0d0,GiD_Vector,GiD_onNodes,GID_STRING_NULL)
  CALL GiD_fVectorComp(fdr,'Comp X','Comp Y','Comp Z','Module')
  CALL GiD_fResultUnit(fdr,'m/s')
  CALL GiD_fResultValues(fdr)
  do idx=1,9
    value=idx*1.5
    CALL GiD_fWrite2DVector(fdr,idx,value,value*3)
  end do
  CALL GiD_fEndResult(fdr)

  !Complex Scalar result
  CALL GiD_fBeginResultHeader(fdr,'Potential','Analysis',1.0d0,GiD_ComplexScalar,GiD_onNodes,GID_STRING_NULL)    
  component_names(1)='real_part'
  component_names(2)='imaginary_part'
  CALL GiD_fResultComponents(fdr,2,component_names)
  CALL GiD_fResultUnit(fdr,'Volt')
  CALL GiD_fResultValues(fdr)
  !CALL GiD_fBeginComplexScalarResult(fdr,'Complex scalar','Analysis',1.0d0,GiD_onNodes,GID_STRING_NULL,GID_STRING_NULL,'real_part','imaginary_part')
  do idx=1,9
    value=idx*1.5
    CALL GiD_fWriteComplexScalar(fdr,idx,value,value*1.2)
  end do
  CALL GiD_fEndResult(fdr)

  !CALL GiD_fClosePostMeshFile(fdm)
  CALL GiD_fClosePostResultFile(fdr)
  
  CALL GiD_PostDone
END PROGRAM
