%% README

%Following CFD simulation, run the following functions to extract some useful information

%NOTE: CFD simulations should output: examplespline.txt, and Eel3d/
%%
%% Extract Energy expenditure from CFD simulations
[E_tot] =  energetics_calculator();
%% Validate fish movement
[sim_eye_vel] = velocity_calculator();
[turn_angle] = turning_calculator();