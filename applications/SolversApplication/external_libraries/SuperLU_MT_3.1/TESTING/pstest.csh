#!/bin/csh

set ofile = pstest.out			# output file
if ( -e $ofile ) then
    rm -f $ofile
endif
echo 'Single-precision testing output' > $ofile

set NVAL        = (10 19)
set NRHS        = (2)
set LWORK       = (0 100000000)
set PANELSIZE   = (2)
set RELAX       = (2)
set NPROCS	= (1 4)

setenv OMP_NUM_THREADS 4

#
# Loop through all matrices ...
#
foreach m (LAPACK g10)
echo $m
    #--------------------------------------------
    # Test matrix types generated in LAPACK-style
    #--------------------------------------------
    if ($m == 'LAPACK') then
      	echo '' >> $ofile
      	echo '** LAPACK test matrices' >> $ofile
      	foreach n ($NVAL)
            foreach s ($NRHS)
              	foreach l ($LWORK)
		    foreach p ($NPROCS)
	    	      	echo '' >> $ofile
            	      	echo 'n='$n 'nrhs='$s 'lwork='$l 'nprocs='$p >> $ofile
            	      	./pstest -t "LA" -l $l -n $n -s $s -p $p >> $ofile
		    end
              	end
            end
        end
    #--------------------------------------------
    # Test a specified sparse matrix
    #--------------------------------------------
    else
      	echo '' >> $ofile
      	echo '** sparse matrix:' $m >> $ofile
      	foreach  w ($PANELSIZE)
            foreach r ($RELAX)
                foreach s ($NRHS)
                    foreach l ($LWORK)
			foreach p ($NPROCS)
	                    echo '' >> $ofile
                      	    echo 'w='$w 'relax='$r 'nrhs='$s 'lwork='$l \
					'nprocs='$p >> $ofile
            	      	    ./pstest -t "SP" -w $w -r $r -s $s -l $l \
                             		< ../EXAMPLE/$m >> $ofile
	   		end
        	  end
      	      end
    	  end
      end
  endif

end


