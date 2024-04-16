# Beam propagation software

Example code used to simulate beam propagation in a system with real optical components (lenses, mirrors, etc.). 

The example presented here propagates the simulated near-field emission output of a terahertz quantum cascade laser (QCL) through two optical lenses, which first collimate and then focus the beam into the imaging plane. Further simulation regarding reflection from the imaging plane and backward propagation through both optical lenses onto the QCL is also shown. Simulation is presented for the following planes: the QCL emission plane, collimating lens L1, focusing lens L2, imaging plane, backward propagation to L2, then to L1, and final reinjection at the QCL's emission plane.  


## Installation
You will need two custom modules 'QCL_methods.py' and 'import_module_from_anywhere.ipynb'. The input QCL near-field profile data is imported from 'QCL_laser_output.mat'. In case there are some nested dependencies, please pip install those. I hope I did not forget nested dependencies of my own. If there is a module name that you can't find online that is required, please let me know and I will fix it.

## Support
In case of any encountered issues (very likely at this stage), let me know via issue tracker or at my email: m.ploschner@uq.edu.au

## Roadmap
The code is in a pretty raw state. It contains a lot of debugging flags and I made no attempts to optimise, beutify or streamline it. I am happy to make the code much more accessible if there is interest in the community. 

## Contributing
I am happy to implement contributions and suggestions from others and make the code of higher-quality.

## Authors and acknowledgment
Author: Martin Ploschner. 
Acknowledgement: Thanks to all those amazing people on https://stackoverflow.com/ finding solutions to all my questions. 

## License
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

## Project status
Active
