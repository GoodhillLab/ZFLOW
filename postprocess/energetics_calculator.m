function [E_tot] =  energetics_calculator()
%% This function loads the CFD power expenditure outputs, computes the total power expenditure, and computes the total energy usage
%% Load the simulation power data
ibamr_power=load('/ZFISH3d/Eel3d_Power_spent_struct_no_0');
%% Extract power
timeColumn = 1;
noiseOFF = 10; %Small offset to exclude spurious power noise near end of simulation
power_spent_x = -1.*(ibamr_power(1:end-noiseOFF,2) - ibamr_power(1:end-noiseOFF,5)); %NOTE1: Eel3d_Power_spent_struct_no_0 represents power exerted by the fish into the environment. As such, it is negative. For ease of interpretation, we flip the axis. 
power_spent_y = -1.*(ibamr_power(1:end-noiseOFF,3) - ibamr_power(1:end-noiseOFF,6)); %NOTE2: An explanation of why these columns are subtracted from each other subtractions can be found in the https://groups.google.com/g/ibamr-users, or by reading the IBAMR paper, or code base
power_spent_z = -1.*(ibamr_power(1:end-noiseOFF,4) - ibamr_power(1:end-noiseOFF,7)); %NOTE3: Z-axis power is very small in our simulations, and is excluded from this point onwards.
%% Compute energy
E_x = trapz(ibamr_power(1:end-noiseOFF,timeColumn), power_spent_x);
E_y = trapz(ibamr_power(1:end-noiseOFF,timeColumn), power_spent_y);
E_tot = E_x+E_y;
end