%%
% required files
% opt1_3c1v.m : 
%           -  input files : connect_region1.csv, node_region1.csv, GT2.mat
% opt2_3c4v.m ~ opt4
%           -  input files: node_fullmodel.csv, connect_fullmodel.csv, GT2.mat


%%
clear; close all;
addpath('~/matlabtool/spm12/')
addpath('../lib/')
addpath('../lib/cf/')

%parpool(2)
step_names = { 'opt1_3c1v', 'opt2_3c4v', 'opt3_3c4v', 'opt4_3c4v' };

for ii=1:4 %length(step_names)

    sname = step_names{ii};

    % coarse/fine tuning with Bayesian optimization
    for run_type = { 'coarse', 'fine' }
        run_options = set_run_options(run_type);
        run_options.do_parallel = false;
        run_options.preDCM = [];
        if ii > 1             
            run_options.preDCM = sprintf('%s_DCM.mat', step_names{ii-1});            
        end
        cmd = sprintf("%s(sname, run_options);", sname); eval(cmd)
    end

end

function run_options = set_run_options(run_type)
if strcmp(run_type,'coarse')
    run_options.type = 'coarse'; 
    run_options.Bayes_NMax = 2; % fast run
    %run_options.Bayes_NMax = 1000; 
    run_options.DCM_NMax = 1; 
else
    run_options.type = 'fine'; 
    run_options.Bayes_NMax = 2; % fast run
    %run_options.Bayes_NMax = 32;    
    run_options.DCM_NMax = 1; % fast run
    %run_options.DCM_NMax = 128;
end
end



