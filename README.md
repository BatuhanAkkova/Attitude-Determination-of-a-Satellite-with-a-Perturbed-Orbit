# Attitude-Determination-of-a-Satellite-with-a-Perturbed-Orbit
This repository contains my MATLAB codes for my Bachelor Thesis. It calculates the perturbations of a satellite orbiting Earth using the formulation of Encke. Then it calculates the attitude determination via QUEST method. Lastly it calculates the RMS error for the attitude determination.

Initial Values contains all the initial parameters used in the code and thesis.
Main is the (obviously) the main code to run the complete simulation.

The Abstract of my thesis is given below:

The Wahba's optimization problem is one of the most important problems in attitude determination and control theory. In this thesis, we implement the Quaternion Estimator algorithm to solve the optimization problem. An elliptical Keplerian orbit around Earth is chosen for demonstrating the orbit. In order to create a measurement model, a hypothetical magnetometer and a hypothetical sun sensor which are onboard our hypothetical spacecraft are modeled. Following the sensor modeling, various orbital perturbations are introduced to emulate real-world scenarios. These perturbations include Earth's aspherical gravitational effect, atmospheric drag, the gravitational influence of the Moon, and the effects of solar radiation pressure. The perturbed orbit is reconstructed to reflect these influences. The perturbation effects on the attitude determination algorithm are investigated. By incorporating all these perturbations, the algorithm's robustness and accuracy in real-world conditions are assessed. Root mean square errors (RMSE) are computed and analyzed to quantify the attitude determination's performance and reliability. The results of the error analysis indicate that when all perturbation effects are accounted for, the algorithm yields the minimum error, thereby demonstrating the desired accuracy and reliability. This comprehensive approach underscores the importance of considering multiple perturbation effects in enhancing the precision of attitude determination algorithms in spacecraft operations.
