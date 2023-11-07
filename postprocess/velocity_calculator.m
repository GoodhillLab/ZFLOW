function [sim_eye_vel] = velocity_calculator()
%% Import the CFD body tracking data
IBAMR_tracking = table2array(readtable("example_simulation/examplespline.txt"));
%% Extract eye midpoint positions
sim_eye = IBAMR_tracking(IBAMR_tracking(:,4)==1,1:3); 
%% Velocity calculation
dt = (1/49)/500; %Simulation time-step. This is derived from our preprocessing pipeline (500 fps recording, with interpolation between frames)
sim_eye_vel = diff(sim_eye)/dt;
end