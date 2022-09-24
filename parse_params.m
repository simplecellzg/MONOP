% clear
% load read_ff_mat
% file=file;
% ffdata=reax;
% ffid_pointer=ffid_pointer;
% 
% 
% paramfile='params';
function params=parse_params(ffid_pointer,file,ffdata)

% file.param=paramfile;
fo=fopen(file.log,'a+');
fi=fopen(file.params,'r');
fprintf(fo, '2, READING PARAMETERS FILE:  %s ...\n',  file.params);
s=fgetl(fi);
i=1;
while ~feof(fi)
    
    s=fgetl(fi);
    strtmp=strsplit(strtrim(s));
    section=str2num(strtmp{1});
    entry=str2num(strtmp{2});
    parameter=str2num(strtmp{3});
    offset=str2num(strtmp{4});
    if offset == 0.0
    params.Lower(i)=str2num(strtmp{5});
    params.Upper(i)=str2num(strtmp{6});
    end
    
    switch section
        case 1  % general parameters
            params.params_ptr(i) = parameter;
            section_label='GENERAL';
        case 2 % single body parameters
            if ffdata.lg_reax == 1
                params.params_ptr(i) = ffdata.gp_idx + (entry -1) * 34 + parameter;
            else
                params.params_ptr(i) = ffdata.gp_idx + (entry -1) * 32 + parameter;
            end
            section_label='ATOMS';
        case 3 % two body parameters
            params.params_ptr(i) = ffdata.gp_idx + ffdata.sbp_idx + (entry -1) * 16 + parameter;
            section_label='BONDS';
        case 4 % two body off-diagonal parameters
            if ffdata.lg_reax == 1
                params.params_ptr(i) = ffdata.gp_idx + ffdata.sbp_idx + ffdata.tbp_idx + (entry -1) * 7 + parameter;
            else
                params.params_ptr(i) = ffdata.gp_idx + ffdata.sbp_idx + ffdata.tbp_idx + (entry -1) * 6 + parameter;
            end
            section_label='OFF-DIAGONALS';
        case 5 % three body parameters
            params.params_ptr(i) = ffdata.gp_idx + ffdata.sbp_idx + ffdata.tbp_idx + ffdata.tbodp_idx + ...
                (entry -1) * 7 + parameter;
            section_label='ANGLES';
        case 6 % four body parameters
            params.params_ptr(i) = ffdata.gp_idx + ffdata.sbp_idx + ffdata.tbp_idx + ffdata.tbodp_idx + ...
            ffdata.thbp_idx + (entry -1) * 7 + parameter;
            section_label='TORSIONS';
        case 7 % hydrogen body parameters
            params.params_ptr(i) = ffdata.gp_idx + ffdata.sbp_idx + ffdata.tbp_idx + ffdata.tbodp_idx + ...
            ffdata.thbp_idx + ffdata.fbp_idx + (entry -1) * 4 + parameter;
            section_label='H-BONDS';
    end
    if offset ~= 0.0
        if  ffid_pointer(params.params_ptr(i)) >= 0
    	params.Lower(i) =  ffid_pointer(params.params_ptr(i))* (1.0 - offset);
    	params.Upper(i) =  ffid_pointer(params.params_ptr(i))* (1.0 + offset);
        else
        params.Lower(i) =  ffid_pointer(params.params_ptr(i))* (1.0 + offset);
    	params.Upper(i) =  ffid_pointer(params.params_ptr(i))* (1.0 - offset);
        end
    end
    if section == 1
        fprintf (fo,'--> %s: change parameter %d [with global index %d] within [%f:%f]\n', ...
            section_label, parameter, params.params_ptr(i), params.Lower(i), params.Upper(i));
    else
        fprintf (fo,'--> %s: change entry %d, parameter %d [with global index %d] within [%f:%f]\n', ...
            section_label, entry, parameter, params.params_ptr(i), params.Lower(i), params.Upper(i));
    end
    i = i + 1;
end
 params.param_num=i-1;

fclose(fo);
fclose(fi);

 % save params_data file params
 % clear

end