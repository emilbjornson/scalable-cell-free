Scalable Cell-Free Massive MIMO Systems
==================

This is a code package is related to the follow scientific article:

Emil Björnson and Luca Sanguinetti, “[Scalable Cell-Free Massive MIMO Systems](https://arxiv.org/pdf/1908.03119),” IEEE Transactions on Communications, to appear.

The package contains a simulation environment, based on Matlab, that reproduces some of the numerical results and figures in the article. *We encourage you to also perform reproducible research!*


## Abstract of Article

Imagine a coverage area with many wireless access points that cooperate to jointly serve the users, instead of creating autonomous cells. Such a cell-free network operation can potentially resolve many of the interference issues that appear in current cellular networks. This ambition was previously called Network MIMO (multiple-input multiple-output) and has recently reappeared under the name Cell-Free Massive MIMO. The main challenge is to achieve the benefits of cell-free operation in a practically feasible way, with computational complexity and fronthaul requirements that are scalable to large networks with many users. We propose a new framework for scalable Cell-Free Massive MIMO systems by exploiting the dynamic cooperation cluster concept from the Network MIMO literature. We provide a novel algorithm for joint initial access, pilot assignment, and cluster formation that is proved to be scalable. Moreover, we adapt the standard channel estimation, precoding, and combining methods to become scalable. A new uplink and downlink duality is proved and used to heuristically design the precoding vectors on the basis of the combining vectors. Interestingly, the proposed scalable precoding and combining outperform conventional maximum ratio processing and also performs closely to the best unscalable alternatives.


## Content of Code Package

The article contains 5 simulation figures, numbered 4, 5(a), 5(b), 6(a) and 6(b). simulationFigure4 generates Figure 4, simulationFigure5 generates Figures 5(a) and 5(b), and simulationFigure6 generates Figures 6(a) and 6(b). The package also contains 5 Matlab functions that are used by some of the scripts.

See each file for further documentation.


## Acknowledgements

E. Björnson was supported by ELLIIT and the Wallenberg AI, Autonomous Systems and Software Program (WASP). L. Sanguinetti was supported by the University of Pisa under the PRA 2018-2019 Research Project CONCEPT.


## License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
