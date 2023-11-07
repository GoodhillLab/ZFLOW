%%
% This script visualizes the 2D tracking data for our example swim bout, following preprocessing. 
% 
%%
clear
clc
close all
load('zebrafish_2D_coords.txt')
i = 1;
ii = 125;
for k = 1:length(zebrafish_2D_coords)/125
    if mod(k,10)==0
        plot(zebrafish_2D_coords(i:ii,1),zebrafish_2D_coords(i:ii,2))
        xlim([0 1.5])
        ylim([0 1.5])
        pause(0.000001)
    end
    i = ii+1;
    ii = ii + 125;
end