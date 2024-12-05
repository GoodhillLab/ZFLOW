%% load data
fname = "Eel3d/Eel3d_COM_coordinates_struct_no_0";
data_table = readtable(fname, "ReadRowNames", false, "ReadVariableNames", false);
data_table.Properties.VariableNames = {'t', 'x', 'y', 'z'}; % names of columns
data_table = data_table(5000: end - 10, :); % remove simulation noise at the end

%%

dx = diff(data_table.x);
dy = diff(data_table.y);
distance = sqrt(dx.^2 + dy.^2);
t = data_table.t(2: end);
speed = distance ./ t;
plot(t, speed)
xlabel("Times (s)")
ylabel("")
