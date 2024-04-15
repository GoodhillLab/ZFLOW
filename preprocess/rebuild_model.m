% This script is used to rebuild the 3D model from the point clouds and the
% 2D coordinates of the zebrafish.
% The model is scaled to the desired length and translated to the desired
% position. The 2D coordinates of the zebrafish are used to determine the
% position of the model in the x-direction. The model is then shifted to
% the appropriate position in the x-direction.
% Original author: Thomas Darveniza from
% https://github.com/ThomasDarveniza/ZFHYDRO_CFD/blob/main/CFD_setup/preprocess/code/rebuild_model.m
% Modified by: Robert Wong

%% load the model

source = "\\wsl.localhost\Ubuntu\data\test\CFD_test";
dest = "..\simulation";
NUM_KEYPOINTS = 125;
sizes = load_file(fullfile(source, 'sizes.txt'), [1, Inf], 3);
eel3d = load_file(fullfile(source, 'eel3d.vertex'), [2, Inf], 3);
zebrafish2Dcoords = load_file(fullfile(dest, 'zebrafish_2D_coords.txt'), [1, NUM_KEYPOINTS], 2);

%% Rebuild the model

%Shift to x=0
eel3d(:,1) = eel3d(:,1) - min(eel3d(:,1));

%Change length (large-scale)
%Current length
curr_len = abs(max(eel3d(:,1))-min(eel3d(:,1)));
%Desired length
des_len = 0.45; % hard-coded for now
final_shape = eel3d*(des_len/curr_len);
%Translate x-direction (Again)
desired_x = 0.850185;
actual_x = mean(final_shape(:,2));
shift_x = desired_x - actual_x;
final_shape(:,2) = final_shape(:,2) + shift_x;

%Change length (small-scale)
%Shift the eel3d shape to appropriate position
mps = linspace(0, des_len, NUM_KEYPOINTS)';
xval = 0.850185;

%Calculate the initial shifts
reference_shape = zebrafish2Dcoords;
shiftxx = reference_shape(:,1)-mps;
shiftyy = reference_shape(:,2)-xval;
%Do the initial shifts
mps = mps+shiftxx;
xval = xval+shiftyy;

start=1;
finish=sizes(1,1);
p = nan(size(eel3d));

for ii = 1: NUM_KEYPOINTS
    for j=start:finish
        p(j,1)=final_shape(j,1)+shiftxx(ii);
        p(j,2)=final_shape(j,2)+shiftyy(ii);
        p(j,3)=final_shape(j,3);
    end
    if ii < NUM_KEYPOINTS
        start=finish+1;
        finish=start+sizes(ii+1)-1;
    end
end

%% output new eel3d.vertex and stacked.txt

% Write to eel3d.vertex
fid = fopen(fullfile(dest, 'eel3d.vertex'), 'wt');
fprintf(fid, '%i\n', size(p, 1));
fprintf(fid, '%f\t%f\t%f\n', p.');
fclose(fid);

fid = fopen(fullfile(dest, 'stacked.txt'), 'wt');
fprintf(fid, '%f\t%f\t%f\n', p.');
fclose(fid);


%% helper function
function data = load_file(file_path, data_lines, num_vars)
    opts = delimitedTextImportOptions("NumVariables", num_vars);
    opts.DataLines = data_lines;
    opts.Delimiter = "\t";
    opts.VariableNames = arrayfun(@(n) sprintf('VarName%d', n), 1:num_vars, 'UniformOutput', false);
    opts.VariableTypes = repmat("double", 1, num_vars);
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    data = readtable(file_path, opts);
    data = table2array(data);
end
