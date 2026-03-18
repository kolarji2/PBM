%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Load data to memory if ran for the first time.
%   (If varibale "readdata" is set to true or if does not exist)
%
%   data_nd = normal dissolution data 
%   data_rd = reactive dissolution data 
%   data_kin = kinetics exp data 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fname_kinetics='paliperidone_kinetics.xlsx'; % file with kinetics data
fname_normaldiss='nd_allexps.xlsx'; % file with normal dissolution data
fname_reactivediss='rd_allexps.xlsx';% file with reactive diss data
basedir='SonntagEtAl/expdata'; % default path to files
%filename_distribution='paliperidone_phiCamSizer2.xlsx'; % path to particle size distribution
filename_distribution='paliperidone_phi0.xlsx'; % horiba ps only
% 'paliperidone_phiCamSizer2.xlsx'; % combined horiba for milled, camsizer for crude


if ~exist('readdata')
    readdata=true;
end

if readdata
    fname_dist=fullfile(basedir, filename_distribution);
    data_nd=load_exp_fromfile( fullfile(basedir, fname_normaldiss),fname_dist);
    data_rd=load_exp_fromfile( fullfile(basedir, fname_reactivediss),fname_dist);
    data_kin = load_exp_fromfile( fullfile(basedir, fname_kinetics),fname_dist);
    fprintf('Experimental data loaded to memory\n as variables ''data_nd; data_rd; data_kin''.\n')
    warning('Use {datavar}.all.dexp{i} to access diameter of each experiment. In pars.dexp is only last distribution!')
    readdata=false;
else
    warning('Data loading skipped. Variable ''readdata'' is set to false. Set variable readdata to true to load new data.')
end