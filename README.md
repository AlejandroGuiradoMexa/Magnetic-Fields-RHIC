# Magnetic-Fields-RHIC
In Relativistic Heavy Ion Collisions(RHIC) the magnetic fields created by spectators protons are extremely strong. Therefore, magnetic fields are extensively study in order to characterize and analyze the effects on particle observables and transports effects. The model used to calculate the electromagnetic field, we used the Lienard-Wiechert model obtained by considering the retardation time in classical electrodynamic equations. 

<img src="https://latex.codecogs.com/svg.latex?\centering&space;e\mathbf{B}(\mathbf{r},t)=\alpha&space;\frac{\mathbf{v}\times&space;\mathbf{R}&space;(&space;1-v^{2})}&space;{&space;R^{3}(1-\frac{(|\mathbf{R}\times\mathbf{v}|)^{2}}{R^{2}})^{3/2}&space;}&space;\qquad&space;e\mathbf{E}(\mathbf{r},t)=\alpha&space;\frac{\mathbf{R}&space;(&space;1-v^{2})}&space;{&space;R^{3}(1-\frac{(|\mathbf{R}\times\mathbf{v}|)^{2}}{R^{2}})^{3/2}&space;}" title="\centering e\mathbf{B}(\mathbf{r},t)=\alpha \frac{\mathbf{v}\times \mathbf{R} ( 1-v^{2})} { R^{3}(1-\frac{(|\mathbf{R}\times\mathbf{v}|)^{2}}{R^{2}})^{3/2} } \qquad e\mathbf{E}(\mathbf{r},t)=\alpha \frac{\mathbf{R} ( 1-v^{2})} { R^{3}(1-\frac{(|\mathbf{R}\times\mathbf{v}|)^{2}}{R^{2}})^{3/2} }" /></a>

To present the results that can be easily compared to recent publications, the equations above are scaled to the Pion´s mass and since we are working in natural units the Lienard-Wiechert electrmoagnetic fields are:

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{e\mathbf{B}(\mathbf{r},t)}{m_{\pi}^{2}}=\alpha\left(\frac{197}{135}&space;\right)^{2}&space;\frac{\mathbf{p}\times&space;\mathbf{R}&space;(&space;E^{2}-|\mathbf{p}|^{2}&space;)&space;}&space;{&space;((R&space;E)^{2}-|\mathbf{R}\times\mathbf{p}|^{2})^{3/2}&space;}&space;\qquad&space;\frac{e\mathbf{E}(\mathbf{r},t)}{m_{\pi}^{2}}=\alpha\left(\frac{197}{135}&space;\right)^{2}&space;\frac{\mathbf{R}&space;(&space;E^{2}-|\mathbf{p}|^{2}&space;)&space;}&space;{&space;((R&space;E)^{2}-|\mathbf{R}\times\mathbf{p}|^{2})^{3/2}&space;}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\frac{e\mathbf{B}(\mathbf{r},t)}{m_{\pi}^{2}}=\alpha\left(\frac{197}{135}&space;\right)^{2}&space;\frac{\mathbf{p}\times&space;\mathbf{R}&space;(&space;E^{2}-|\mathbf{p}|^{2}&space;)&space;}&space;{&space;((R&space;E)^{2}-|\mathbf{R}\times\mathbf{p}|^{2})^{3/2}&space;}&space;\qquad&space;\frac{e\mathbf{E}(\mathbf{r},t)}{m_{\pi}^{2}}=\alpha\left(\frac{197}{135}&space;\right)^{2}&space;\frac{\mathbf{R}&space;(&space;E^{2}-|\mathbf{p}|^{2}&space;)&space;}&space;{&space;((R&space;E)^{2}-|\mathbf{R}\times\mathbf{p}|^{2})^{3/2}&space;}" title="\frac{e\mathbf{B}(\mathbf{r},t)}{m_{\pi}^{2}}=\alpha\left(\frac{197}{135} \right)^{2} \frac{\mathbf{p}\times \mathbf{R} ( E^{2}-|\mathbf{p}|^{2} ) } { ((R E)^{2}-|\mathbf{R}\times\mathbf{p}|^{2})^{3/2} } \qquad \frac{e\mathbf{E}(\mathbf{r},t)}{m_{\pi}^{2}}=\alpha\left(\frac{197}{135} \right)^{2} \frac{\mathbf{R} ( E^{2}-|\mathbf{p}|^{2} ) } { ((R E)^{2}-|\mathbf{R}\times\mathbf{p}|^{2})^{3/2} }" /></a>

The last equations are implemented in all the codes, some conditions are added in order to consider all charged particles, spectators or participants. A brief description of the codes:

* OutReacionPlane. -> Magnetic field produced by spectator protons at different times of collisions.
* TemporalEvol_OutReactionPlaneGIF. -> Produce a GIF of the OutReactionPlane results.
* TemporalEvolution. -> Temporal evolution of magnetic fields created by protons at x,y,z directions.
* CentralityBin_y. -> Temporal evolution of magnetic fields for different centrality ranges.

To use this codes from scratch, I mean generate the simulated data is needed to run urqmd and use the macro urqmdtoroot.C to create a C++ struct in order to save the data in ROOT tree files.




## Monte Carlo Generator
The Monte Carlo generator is used to simulate the initial condition of nuclei collisions, the Glauber Monte Carlo is very accurate with experiments results and is the one used in this work. After the initial condition is setted up, in order to take into account the temporal evolution of the system, several models are available. In UrQMD(<a href="arXiv:hep-ph/9909407"> arXiv:hep-ph/9909407</a>  ), the nucleus propagation is considered without any potential for all particles. However it is possible to add the Coulomb, Pauli, among other potentials.  
