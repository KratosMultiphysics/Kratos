      PROGRAM EXAMPLE

      REAL*8 nodes(1:9,1:3), rx, ry, rz, value
      INTEGER elems(1:4,1:4)
      INTEGER*4 idx, x, y, nid(1:4)
      CHARACTER*4 NULL
      CHARACTER*100 nombre

      NULL = CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)
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
            elems(idx,4) = elems(idx,1) + 3
            elems(idx,3) = elems(idx,4) + 1
            idx = idx + 1
         end do
      end do

C      CALL GID_OPENPOSTMESHFILE('testfortran.post.msh',0)
      CALL GID_OPENPOSTRESULTFILE('testfortran.post.res',2)
      CALL GID_BEGINMESHCOLOR('quadmesh',2,4,4,0.7d0,0.7d0,0.4d0)
      CALL GID_BEGINCOORDINATES
      idx = 1
      do x=0,2
         do y=0,2
            rx = x
            ry = y	
		  rz = 0.0;      
            CALL GID_WRITECOORDINATES(idx, rx, ry, rz)
            idx = idx + 1
         end do
      end do
      CALL GID_ENDCOORDINATES
      CALL GID_BEGINELEMENTS
      idx = 1
      
      do x=1,2
         do y=1,2
            nid(1) = x + (y-1)*3
            nid(2) = nid(1) + 1
            nid(4) = nid(1) + 3
            nid(3) = nid(4) + 1
            CALL GID_WRITEELEMENT(idx, nid)
            idx = idx + 1
         end do
      end do      
      CALL GID_ENDELEMENTS
      CALL GID_ENDMESH
C     Scalar result
      CALL GID_BEGINSCALARRESULT('Scalar'//char(0),'An.'//char(0),1.0,0,
     .NULL,NULL,NULL)      
      do idx=1,9
         value=idx*1.5;
         CALL GID_WRITESCALAR(idx,value)
      end do
      CALL GID_ENDRESULT
      CALL GID_FLUSHPOSTFILE
C     Vectorial result
      CALL GID_BEGINRESULTHEADER('Vectorial'//char(0),
     .'An.'//char(0), 1.0, 1, 0, NULL)
C     CALL GID_VECTORCOMP('Comp X', 'Comp Y', 'Comp Z', "Length")
      CALL GID_RESULTVALUES
      do idx=1,9
         value=idx*1.5;
         CALL GID_WRITEVECTOR(idx,value,value*3,value*5)
      end do
      CALL GID_ENDRESULT

C      CALL GID_CLOSEPOSTMESHFILE
      CALL GID_CLOSEPOSTRESULTFILE
      END
