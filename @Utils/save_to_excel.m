function save_to_excel(obj,data,fields,fpath,sheet)
%SAVE_TO_EXCELL(data,fields,fpath) Save data to excel
%   Input:
%       data: cell array of structures or PBM classes
%               {PBM,expdata}
%
%       fields: cell array of cell arrays of corresponding fields to save, column wise
%               {{tm,cm},{texp,cexp,cstd}};
k=1;
utils=Utils;
for i=1:length(data)
    datai=data{i};
    fieldsi=fields{i};
    for j=1:length(fieldsi)        
        coldata=datai.(fieldsi{j});
        [Nrow,Ncol]=size(coldata);
        cols=utils.xlsColNum2Str([k,k+Ncol-1]);        
        range_0=sprintf('%s1:%s1',cols{1},cols{2});
        writecell(repmat({fieldsi{j}},1,Ncol),fpath,'Range',range_0,'Sheet',sheet); 
        range_=sprintf('%s2:%s%d',cols{1},cols{2},Nrow+1);
        writematrix(coldata,fpath,'Range',range_,'Sheet',sheet);        
        k=k+Ncol;
    end    
end

end

