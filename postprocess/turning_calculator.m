function [turn_angle] = turning_calculator()
%% Import the CFD body tracking data
IBAMR_tracking = table2array(readtable("example_simulation/examplespline.txt"));
%% Extract eye midpoint & swim bladder positions
sim_eye = IBAMR_tracking(IBAMR_tracking(:,4)==1,1:3); 
sim_bladder = IBAMR_tracking(IBAMR_tracking(:,4)==2,1:3);
%% Compute total turn angle
noiseOFF = 10; %Small offset to exclude spurious power noise near end of simulation
%Starting heading angle
start_eye = sim_eye(1,1:2); %column 1&2 are the x&y positions
start_bladder = sim_bladder(1,1:2);
finish_eye = sim_eye(end-noiseOFF,1:2);
finish_bladder = sim_bladder(end-noiseOFF,1:2);
a = start_eye(1);
aa = start_bladder(1);
b = start_eye(2);
bb = start_bladder(2);
heading_start = rad2deg(atan2(b-bb,a-aa));  
%Finishing heading angle 
a = finish_eye(1);
aa = finish_bladder(1);
b = finish_eye(2);
bb = finish_bladder(2);
heading_finish = rad2deg(atan2(b-bb,a-aa)); 
%Total turn angle
turn_angle = heading_start - heading_finish; %Depending on the fish-system orientation, you will need to do further steps on this turn angle
%% Sanity-check
%plot(a,b,'.r');hold on;plot(aa,bb,'.r')
end