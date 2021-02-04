# MRMPCHalo

# Overview
 - [:orange_book: The general idea](#orange_book-some-theory-behind-the-code)
 - [:page_facing_up: Dependencies](#page_facing_up-dependencies)

# :orange_book: The general idea
This repo tests the recently proposed multirate planning and nmpc control scheme on the problem of translunar station-keeping for a _quasi_ halo 
orbit.

## The mode

For modelling the equations of motion of a satellite/space-craft in the Earth moon system, we opt for the following representation

![model](https://github.com/mebbaid/MRMPCHalo/blob/main/MR%20MPC%20simulink%20implm/image/the_model-1.png)




Additionally, a simplified model of the solar radiation pressure effect on the position of the space-craft is model, so getting a more realistic model. To this end, the
following section highlights the simulation mode coupled with the proposed planning and control scheme.

## The simulation set-up

For the simulation purposes, the following scheme serves to illustrate the general diagram, in a very rough way.

![simulation](https://github.com/mebbaid/MRMPCHalo/blob/main/MR%20MPC%20simulink%20implm/image/simulation_scheme.png)


# :page_facing_up: Dependencies
1. **CASADI and IPOPT:** Install the CASADI and IPOPT if you intend to use different solver, by default we use Matlab's ```fmincon``` [CASADI](https://github.com/casadi/casadi/wiki/InstallationInstructions)


# :hammer: Run the simulation
## Windows

You can simply run ```initMRMPCHalo.m```. At the start of this init file, most relevant parameters to the simulation are present, adjust for different
scenarios e.g. primaries eccentricity, thrust saturation limits, different weights ...etc here. 
