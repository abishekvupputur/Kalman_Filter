System Identification Assignment

To run Kalman Filter run the following file: simdata2019/EKF.m
Within EKF:
Kal_betterNames: Is used to generate the Jacobian for state transition matrix (F matrix)
KalHnames: Is used to generate Jacobian for observation matrix. (H matrix)
T2: As in Task 2 of question 2 is used to prepare dataset for preprocessing
rl4: Is used to perform one step ahead prediction
kf_calc_F: Is used to evaluate the jacobian of state transition.
kf_calc_H: Is used to evaluate the jacobian of observation matrix.
kf_calc_f: Is used to evaluate the state transition.
kf_calc_h: Is used to evaluate the observation.

LSE: Used to run System Identification
LSE_red: Used to run System Identification alternate model.