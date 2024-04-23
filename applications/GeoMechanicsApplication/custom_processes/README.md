# Custom processes

## c-$\phi$ reduction process
For the assesment of a safety factor to characterize slope stability, a Mohr-Coulomb material based c-$\phi$ reduction 
scheme is implemented. The apex of the Mohr_Coulomb cone shaped failure surface is kept in the same position, 
therefore both c and $\tan \phi$ will diminish at the same rate.

### Incrementation scheme
The c-$\phi$ reduction process requires the existence of a stress state in your model and the use of a Mohr=Coulomb material 
(in a UDSM or UMAT formulation). Preferably this stress state is an equilibrium state, such that no stresses in integration 
points violate the given Mohr-Coulomb failure surface. During the stage with the active c-$\phi$ reduction process, 
c and $\tan \phi$ will be incrementally reduced in steps with an initial size of 10%. For each reduction step stresses are 
mapped back onto the reduced Mohr-Coulomb yield surface and equilibrium is found if possible. When equilibrium is no longer 
found in the given number of iterations for the Newton-Raphson scheme, the step is retried with a halved reduction increment. This is repeated until the allowed number of cycles for a step is reached.   

### Safety factor
The safety factor is computed as the inverse of the reduction factor:

$$SF = \frac{1}{\alpha}$$

where the reduction factor $\alpha$ is the ratio of the found critical values for cohesion and friction angle $c_c$ or $\phi_c$ and original material parameters.

$$\alpha = \frac{c_c}{c} = \frac{\tan \phi_c}{\tan \phi}$$ 

## Bibliography
Brinkgreve, R.B.J., Bakker, H.L., 1991. Non-linear finite element analysis of safety factors, Computer Methods and Advances in Geomechanics, Beer, Booker & Carterr (eds), Balkema, Rotterdam.
