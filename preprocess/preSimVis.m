%% This script visualizes the 2D tracking data for our example swim bout, following preprocessing. 

load('../simulation/zebrafish_2D_coords.txt')

figure;

assert(mod(length(zebrafish_2D_coords), 125) == 0, "Must be 125 tracked points")
T = uint32(length(zebrafish_2D_coords) / 125);

for t = 1:T
    if mod(t, 10)==0
        title(['Frame ', num2str(t)])
        i = (t - 1) * 125 + 1;
        ii = t * 125;

        plot(zebrafish_2D_coords(i:ii,1),zebrafish_2D_coords(i:ii,2))
        xlim([0 1.5])
        ylim([0 1.5])
        pause(0.000001)
    end    
end
