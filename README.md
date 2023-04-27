# Hydra_General

General code for simulating hydra organism. Performed on a triangulated discretisation of the spheroid, characterised by spatial position of each node, and a chemistry value $\Psi$ assigned to each triangular face. 

A set of IC's are defined at the bottom of 'hydra_3d_main', a simulation can be run using an initial condition using e.g. Sim1 = RunSimulation(SphereMonopole, 50, 0.02), which runs for T=50 timesteps with dt=0.02. Optional argumenst V0grad, V0thresh, V0rest are used to characterise the sawtooth value of the contained fluid. By default the gradient is V0grad=0.0, make this number positive for sawtoothing. 

After completion this simulation can be visualised with AnimateSimulation(Sim1), which saves the animation into the Animations folder. The function Show(Sim1, t) can be used to view an interactable plot of a particular timestep. 
