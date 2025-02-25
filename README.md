# Mapped Adjoint Control Transformation for Low Thrust Trajectory Design

Spacecraft trajectory optimization is a critical task in space mission design. The propulsion system of a spacecraft significantly influences the types of trajectories it can achieve. Over the past few decades, electric propulsion systems characterized by high specific impulse but low thrust magnitudes—have revolutionized space trajectory design.

Low-thrust trajectory design can be formulated as boundary-value problems, which are often challenging to solve due to:
- A small domain of convergence.
- The lack of knowledge about initial costates when using the indirect formalism of optimal control.

Estimating the missing values of non-intuitive costates is a crucial step in solving these boundary-value problems. By leveraging the costate vector mapping theorem, the Adjoint Control Transformation (ACT) method is extended to alternative sets of coordinates/elements for solving low-thrust trajectory optimization problems. This extension is called Mapped Adjoint Control Transformation (MACT).

This repository contains the MATLAB codes for the interplanetary rendezvous maneuver from Earth to Dionysus, as presented in my MSc thesis. The thesis can be accessed at https://www.researchgate.net/publication/373044382_Mapped_Adjoint_Control_Transformation_for_Low-Thrust_Space_Mission_Design

