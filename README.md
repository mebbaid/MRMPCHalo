# MRMPCHalo

# Overview
 - [:orange_book: The general idea](#orange_book-some-theory-behind-the-code)
 - [:page_facing_up: Dependencies](#page_facing_up-dependencies)

# :orange_book: The general idea
This repo tests the recently proposed multirate planning and nmpc control scheme on the problem of translunar station-keeping for a _quasi_ halo 
orbit.
![alt text](https://github.com/mebbaid/MRMPCHalo/blob/main/image/Header.JPG)

# :page_facing_up: Dependencies
1. **CASADI and IPOPT:** Install the CASADI and IPOPT if you intend to use different solver, by default we use Matlab's ```fmincon``` [CASADI](https://github.com/casadi/casadi/wiki/InstallationInstructions)


# :hammer: Run the simulation
## Windows

You can simply run ```initMRMPCHalo.m```. At the start of this init file, most relevant parameters to the simulation are present, adjust for different
scenarios e.g. primaries eccentricity, thrust saturation limits, different weights ...etc here. 
