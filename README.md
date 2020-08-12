# Magnetic-Fields-RHIC
In Relativistic Heavy Ion Collisions(RHIC) the magnetic fields created by spectators protons are extremely strong. Therefore, magnetic fields are extensively study in order to characterize and analyze the effects on particle observables and transports effects. The model used to calculate the electromagnetic field, we used the Lienard-Wiechert model obtained by considering the retardation time in classical electrodynamic equations. 

<img src="https://latex.codecogs.com/svg.latex?\centering&space;e\mathbf{B}(\mathbf{r},t)=\alpha&space;\frac{\mathbf{v}\times&space;\mathbf{R}&space;(&space;1-v^{2})}&space;{&space;R^{3}(1-\frac{(|\mathbf{R}\times\mathbf{v}|)^{2}}{R^{2}})^{3/2}&space;}&space;\qquad&space;e\mathbf{E}(\mathbf{r},t)=\alpha&space;\frac{\mathbf{R}&space;(&space;1-v^{2})}&space;{&space;R^{3}(1-\frac{(|\mathbf{R}\times\mathbf{v}|)^{2}}{R^{2}})^{3/2}&space;}" title="\centering e\mathbf{B}(\mathbf{r},t)=\alpha \frac{\mathbf{v}\times \mathbf{R} ( 1-v^{2})} { R^{3}(1-\frac{(|\mathbf{R}\times\mathbf{v}|)^{2}}{R^{2}})^{3/2} } \qquad e\mathbf{E}(\mathbf{r},t)=\alpha \frac{\mathbf{R} ( 1-v^{2})} { R^{3}(1-\frac{(|\mathbf{R}\times\mathbf{v}|)^{2}}{R^{2}})^{3/2} }" /></a>

To present the results that can be easily compared to recent publications, the equations above are scaled to the Pion´s mass and since we are working in natural units the Lienard-Wiechert electrmoagnetic fields are:




## Monte Carlo Generator
The Monte Carlo generator is used to simulate the initial condition of nuclei collisions, the Glauber Monte Carlo is very accurate with experiments results and is the one used in this work. After the initial condition is setted up, in order to take into account the temporal evolution of the system, several models are available. In UrQMD(arXiv:hep-ph/9909407), the nucleus propagation is considered without any potential for all particles. However is possible to change this configuration, if we want to consider Coulomb, Pauli, among other potentials. 
