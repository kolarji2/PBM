%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Load data to memory if ran for the first time.
%   (If varibale "readdataKCl" is set to true or if does not exist)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('readdataKCl')
    readdataKCl=true;
end
calib=load('DjukajEtAl/expdata/conductivity.mat');
filenames={'monomodalArticle.xlsx';'bimodalArticle.xlsx';'trimodalArticle.xlsx'};

Nexp=length(filenames);
if readdataKCl      
    dataKCl=cell(Nexp,1);
    infoKCl=cell(Nexp,1);    
    for iexp=1:Nexp;
        [dataKCl{iexp},infoKCl{iexp}] = load_exp_dataKCl( ['DjukajEtAl/expdata/' filenames{iexp}],calib.cond2conc);
        dataKCl{iexp}.expid=iexp;
    end
    fprintf('Experimental data set loaded to memory as variable ''data''.\n')
    readdataKCl=false;
else
    warning('Data loading skipped. Variable ''readdataKCl'' is set to false. Set variable readdataKCl to true to load new data.');
end