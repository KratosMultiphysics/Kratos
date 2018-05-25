sessile_droplet_test:


This example is an illustration of the water sessile droplet behavior in term of the surface tension, curvature, pressure distribution, and contact angle evolution.


The initial contact angle can be modified as required from lagrangian_sessile_droplet_test.py (line 130). 

The initial condition for the advancing and receding contact angles can be further modified using "SurfaceTension_monolithic_solver.py" which is inside the python_script solver (lines 328-330):
def cont_angle_cond(self):
        theta_adv = self.contact_angle + 1.0
        theta_rec = self.contact_angle - 1.0


The user can use Gid to have this code working with different initial sessile droplet geometry so the droplet can be either hydrophilic or hydrophobic.


To run this test, run the following command in the terminal:
python3 lagrangian_sessile_droplet_test.py