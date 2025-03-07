# **Robust Tube Linear MPC with Artificial Reference**  

This repository contains a simulation demonstration for the paper:  

**Limón, Daniel, et al.**  
*"MPC for tracking piecewise constant references for constrained linear systems."*  
_Automatica, vol. 44, no. 9, 2008, pp. 2382-2387._  

## **Recommended References**  
For a deeper understanding of robust tube-based MPC and related concepts, the following papers are recommended:  

1. **Limon, D., et al.**  
   *"On the design of robust tube-based MPC for tracking."*  
   _IFAC Proceedings Volumes, vol. 41, no. 2, 2008, pp. 15333-15338._  

2. **Mayne, David Q., María M. Seron, and Saša V. Raković.**  
   *"Robust model predictive control of constrained linear systems with bounded disturbances."*  
   _Automatica, vol. 41, no. 2, 2005, pp. 219-224._  

3. **Gilbert, Elmer G., and K. Tin Tan.**  
   *"Linear systems with state and control constraints: The theory and application of maximal output admissible sets."*  
   _IEEE Transactions on Automatic Control, vol. 36, no. 9, 1991, pp. 1008-1020._  

4. **Rakovic, Sasa V., et al.**
   *"Invariant approximations of the minimal robust positively invariant set."*
   _IEEE Transactions on automatic control 50.3 (2005): 406-410._
   
## **Required Toolboxes**  
The following MATLAB toolboxes are required to run the simulations:  

1. **[MPT3](https://www.mpt3.org/)** – for set operations.  
2. **[MATLAB Optimization Toolbox](https://www.mathworks.com/help/optim/)** – for `fmincon`.  
3. **[YALMIP](https://yalmip.github.io/)** – for solving LMIs (included in MPT3 if installed properly).  

---
