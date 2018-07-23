clear ; close all; clc;

dirs = {'optimizers','auxiliary'};
disp('Added to path:')    
for i = 1:length(dirs)    
    addpath(genpath(strcat(pwd,filesep,dirs{i})));
    disp(dirs{i});
end
