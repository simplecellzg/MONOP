
function [ffid_pointer,file,reax]=read_reaxc_ff_file(ff,log,lgvdw)

%ff='ffield.reax.weak2017m';
%log='logfile_m';
% lgvdw=1; or lgvdw=0; use or not use lg
reax.lg_reax=lgvdw;
file.ff=ff;
file.log=log;
% reax.lg_reax=0;
%  initialization
reax.sbp=[]; reax.tbp=[]; reax.tbodp=[]; reax.thbp=[]; reax.fbp=[]; reax.hbp=[];
% create log file
fid=fopen(log,'w+');
fclose(fid);
% start to write log file
fid=fopen(log,'a+');
fprintf (fid,'1, READING REAXFF FORCE FIELD FILE: ffield ...\n');
% global parameters
% reading parameters from ffield file
ffid=fopen(ff,'r');
% reading first header comment
s=fgetl(ffid);
% line 2 is number of global parameters
% strstim remove blankspace of row header and tail;strsplit split string in a cell
s=fgetl(ffid);
strtmp=strsplit(strtrim(s));
n=str2num(strtmp{1});
reax.gp.n_global=n;
% get general types, and skip
for i=1:n
    s=fgetl(ffid);
end
% get number of atom types, and skip
s=fgetl(ffid);
strtmp=strsplit(strtrim(s));
reax.num_atom_types=str2num(strtmp{1});
% if low-gradient used, then atom entries have 5 lines, else 4
% and off-diagonal parameters are 7, vs 6 if no lg
if reax.lg_reax==1
    sblines = 5;
    sb_params = 34;
    off_diag_params = 7;
else
    sblines = 4;
    sb_params = 32;
    off_diag_params = 6;
end
for i = 1 : reax.num_atom_types*sblines+3
    s=fgetl(ffid);
end
%  get number of bond types, and skip
s=fgetl(ffid);
strtmp = strsplit(strtrim(s));
reax.num_bond_types = str2num(strtmp{1});
for i = 1 : reax.num_bond_types * 2 + 1;
    s=fgetl(ffid);
end
%  get number of off-diagonals types, and skip
s=fgetl(ffid);
strtmp=strsplit(strtrim(s));
reax.num_off_diag_types = str2num(strtmp{1});
for i = 1 : reax.num_off_diag_types;
    s=fgetl(ffid);
end
%  get number of angle types, and skip
s=fgetl(ffid);
strtmp=strsplit(strtrim(s));
reax.num_angle_types = str2num(strtmp{1});
for i = 1 :  reax.num_angle_types;
    s=fgetl(ffid);
end
%  get number of torsion types, and skip
s=fgetl(ffid);
strtmp=strsplit(strtrim(s));
reax.num_torsion_types = str2num(strtmp{1});
for i = 1 :  reax.num_torsion_types;
    s=fgetl(ffid);
end
%  get number of hbond types, and skip
s=fgetl(ffid);
strtmp=strsplit(strtrim(s));
reax.num_hbond_types = str2num(strtmp{1});
for i = 1 : reax.num_hbond_types;
    s=fgetl(ffid);
end
% output ffield info
fprintf (fid, 'FFIELD: %d GENERAL, %d ATOM, %d BOND, %d OFF-DIAG, %d ANGLES, %d TORSION, and %d H-BONDS\n', ...
    reax.gp.n_global, reax.num_atom_types, reax.num_bond_types, reax.num_off_diag_types, ...
    reax.num_angle_types, reax.num_torsion_types, reax.num_hbond_types);
% return to start, skip globals and comments for general and reax.sbp
frewind(ffid);

%  Define number of entries for each ffield section
reax.gp_idx = reax.gp.n_global;
reax.sbp_idx = reax.num_atom_types * sb_params;
reax.tbp_idx = reax.num_bond_types * 16;
reax.tbodp_idx = reax.num_off_diag_types * off_diag_params;
reax.thbp_idx = reax.num_angle_types * 7;
reax.fbp_idx = reax.num_torsion_types * 7;
reax.hbp_idx = reax.num_hbond_types * 4;

%  vdWaals type: 1: Shielded Morse, no inner-wall
%                2: inner wall, no shielding
%                3: inner wall+shielding
reax.gp.vdw_type = 0;
a = reax.gp_idx + reax.sbp_idx + reax.tbp_idx + reax.tbodp_idx + reax.thbp_idx + reax.fbp_idx + reax.hbp_idx;
reax.num_allp=a;
fprintf (fid, 'FFIELD: %d parameters\n', a);
fprintf (fid, '1: GENERAL PARAMETERS [%d]\n', reax.gp_idx);
a = 0;
%  reading general parameters
%  reading first header comment
s=fgetl(ffid);
%  line 2 is number of global parameters
s=fgetl(ffid);
%  see reax_types.h for mapping between l(i) and the lambdas used in ff
for i = 1 : reax.gp.n_global;
    s=fgetl(ffid);
    strtmp=strsplit(strtrim(s));
    val = str2num((strtmp{1}));
    reax.gp.l(i) = val;
    a=a+1; ffid_pointer(a) = reax.gp.l(i);
    fprintf (fid, 'GENERAL: [%-2d]%8.4f\n', a, reax.gp.l(i));
end

%  skip atom header lines
for i = 1 : 4
    s=fgetl(ffid);
end

%  reading single atom parameters
%  there are 4 lines of each single atom parameters in ff files. these
%  parameters later determine some of the pair and triplet parameters using
%  combination rules.

fprintf (fid, '2: ATOM PARAMETERS [%d]\n', reax.sbp_idx);

for i = 1 : reax.num_atom_types;
    % line one
    s=fgetl(ffid);
    strtmp=strsplit(strtrim(s));
    reax.sbp(i).name = strtmp{1};
    
    val = str2num (strtmp{2});
    reax.sbp(i).r_s = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).r_s;
    val = str2num (strtmp{3});
    reax.sbp(i).valency = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).valency;
    val = str2num (strtmp{4});
    reax.sbp(i).mass = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).mass;
    val = str2num (strtmp{5});
    reax.sbp(i).r_vdw = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).r_vdw;
    val = str2num (strtmp{6});
    reax.sbp(i).epsilon = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).epsilon;
    val = str2num (strtmp{7});
    reax.sbp(i).gamma = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).gamma;
    val = str2num (strtmp{8});
    reax.sbp(i).r_pi = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).r_pi;
    val = str2num (strtmp{9});
    reax.sbp(i).valency_e = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).valency_e;
    reax.sbp(i).nlp_opt = 0.5 * (reax.sbp(i).valency_e - reax.sbp(i).valency);
    
    % line two
    s=fgetl(ffid);
    strtmp=strsplit(strtrim(s));
    
    val = str2num (strtmp{1});
    reax.sbp(i).alpha = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).alpha;
    val = str2num (strtmp{2});
    reax.sbp(i).gamma_w = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).gamma_w;
    val = str2num (strtmp{3});
    reax.sbp(i).valency_boc = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).valency_boc;
    val = str2num (strtmp{4});
    reax.sbp(i).p_ovun5 = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).p_ovun5;
    val = str2num (strtmp{5});	% not used
    reax.sbp(i).p_ovun6 = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).p_ovun6;
    val = str2num (strtmp{6});
    reax.sbp(i).chi = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).chi;
    val = str2num (strtmp{7});
    reax.sbp(i).eta = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).eta;
    val = str2num (strtmp{8});
    reax.sbp(i).p_hbond = val;
    a=a+1; ffid_pointer(a) = double(reax.sbp(i).p_hbond);
    %    reax.sbp(i).p_hbond = val;
    %    a=a+1; ffid_pointer(a) = reax.sbp(i).p_hbond;
    
    % line 3
    s=fgetl(ffid);
    strtmp=strsplit(strtrim(s));
    
    val = str2num (strtmp{1});
    reax.sbp(i).r_pi_pi = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).r_pi_pi;
    val = str2num (strtmp{2});
    reax.sbp(i).p_lp2 = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).p_lp2;
    val = str2num (strtmp{3});	% not used
    reax.sbp(i).heat = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).heat;
    val = str2num (strtmp{4});
    reax.sbp(i).p_boc3 = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).p_boc3;
    val = str2num (strtmp{5});
    reax.sbp(i).p_boc4 = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).p_boc4;
    val = str2num (strtmp{6});
    reax.sbp(i).p_boc5 = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).p_boc5;
    val = str2num (strtmp{7});	% not used
    reax.sbp(i).reax.sbp_23 = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).reax.sbp_23;
    val = str2num (strtmp{8});	% not used
    reax.sbp(i).reax.sbp_24 = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).reax.sbp_24;
    
    % line 4
    s=fgetl(ffid);
    strtmp=strsplit(strtrim(s));
    
    val = str2num (strtmp{1});
    reax.sbp(i).p_ovun2 = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).p_ovun2;
    val = str2num (strtmp{2});
    reax.sbp(i).p_val3 = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).p_val3;
    val = str2num (strtmp{3});	% not used
    reax.sbp(i).reax.sbp_27 = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).reax.sbp_27;
    val = str2num (strtmp{4});
    reax.sbp(i).valency_val = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).valency_val;
    val = str2num (strtmp{5});
    reax.sbp(i).p_val5 = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).p_val5;
    val = str2num (strtmp{6});
    reax.sbp(i).rcore2 = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).rcore2;
    val = str2num (strtmp{7});
    reax.sbp(i).ecore2 = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).ecore2;
    val = str2num (strtmp{8});
    reax.sbp(i).acore2 = val;
    a=a+1; ffid_pointer(a) = reax.sbp(i).acore2;
    
    % line 5 - iff sb_params
    if reax.lg_reax == 1
        s=fgetl(ffid);
        strtmp=strsplit(strtrim(s));
        
        val = str2num (strtmp{1});
        reax.sbp(i).lg1 = val;
        a=a+1; ffid_pointer(a) = reax.sbp(i).lg1;
        val = str2num (strtmp{2});
        reax.sbp(i).lg2 = val;
        a=a+1; ffid_pointer(a) = reax.sbp(i).lg2;
    end
    
    for idx = 1 : sb_params
        if idx == 1
            fprintf (fid, '--%2d--%3s\n', i, reax.sbp(i).name);
        end
        fprintf (fid, '[%4d]%8.4f ', idx + (i -1)* sb_params + reax.gp_idx, ...
            ffid_pointer( idx + (i-1) * sb_params + reax.gp_idx));
        if (idx == 8 || idx == 16 || idx == 24 || idx == 32 || idx == 34)
            fprintf (fid, '\n');
        end
    end
    
    if (reax.sbp(i).rcore2 > 0.01 && reax.sbp(i).acore2 > 0.01) 	% Inner-wall
        if (reax.sbp(i).gamma_w > 0.5)  % Shielding vdWaals
            if (reax.gp.vdw_type ~= 0 && reax.gp.vdw_type ~= 3)
                fprintf (['WARNING: inconsistent vdWaals-parameters\n' ...
                    'Force field parameters for element %s\n' ...
                    'indicate inner wall+shielding, but earlier\n' ...
                    'atoms indicate different vdWaals-method.\n' ...
                    'This may cause division-by-zero errors.\n' ...
                    'Keeping vdWaals-setting for earlier atoms.\n'], reax.sbp(i).name)
            else
                reax.gp.vdw_type = 3;
                fprintf ('----> vdWaals type for element %s: Shielding+inner-wall\n', ...
                    reax.sbp(i).name);
            end
        else 			% No shielding vdWaals parameters present
            if (reax.gp.vdw_type ~= 0 && reax.gp.vdw_type ~= 2)
                fprintf (['WARNING: inconsistent vdWaals-parameters\n' ...
                    'Force field parameters for element %s\n' ...
                    'indicate inner wall without shielding, but earlier\n' ...
                    'atoms indicate different vdWaals-method.\n' ...
                    'This may cause division-by-zero errors.\n' ...
                    'Keeping vdWaals-setting for earlier atoms.\n'], reax.sbp(i).name);
            else
                reax.gp.vdw_type = 2;
                fprintf (stderr, '----> vdWaals type for element %s: No Shielding,inner-wall\n', ...
                    reax.sbp(i).name);
            end
        end
        
    else 			% No Inner wall parameters present
        if (reax.sbp(i).gamma_w > 0.5) 	% Shielding vdWaals
            if (reax.gp.vdw_type ~= 0 && reax.gp.vdw_type ~= 1)
                fprintf (['WARNING: inconsistent vdWaals-parameters\n' ...
                    'Force field parameters for element %s\n' ...
                    'indicate  shielding without inner wall, but earlier\n' ...
                    'atoms indicate different vdWaals-method.\n' ...
                    'This may cause division-by-zero errors.\n' ...
                    'Keeping vdWaals-setting for earlier atoms.\n'], reax.sbp(i).name);
            else
                reax.gp.vdw_type = 1;
                fprintf ('----> vdWaals type for element %s: Shielding,no inner-wall\n', ...
                    reax.sbp(i).name);
            end
        else
            fprintf (['ERROR: inconsistent vdWaals-parameters\n' ...
                'No shielding or inner-wall set for element %s\n'], reax.sbp(i).name);
        end
    end
end



fprintf (fid, '>> Done with ATOM interaction (%i types and %i parameters)\n', ...
    reax.num_atom_types,a);
fprintf ('----> vdWaals type: %d\n', reax.gp.vdw_type);

% Equate vval3 to valf for first-row elements (25/10/2004)
for i = 1 : reax.num_atom_types
    if (reax.sbp(i).mass < 21 && reax.sbp(i).valency_val ~= reax.sbp(i).valency_boc)
        fprintf ('\nWARNING: changed valency_val to valency_boc for %s\n', reax.sbp(i).name);
        reax.sbp(i).valency_val = reax.sbp(i).valency_boc;
    end
end

% next line is number of two body combination and some comments
s=fgetl(ffid);
% a line of comments
s=fgetl(ffid);
fprintf (fid, '3: 2-BODY PARAMETERS [%d]\n', reax.tbp_idx);


for (i = 1 : reax.num_bond_types)
    %    line 1
    s=fgetl(ffid);
    strtmp=strsplit(strtrim(s));
    
    j = str2num (strtmp{1});
    k = str2num (strtmp{2});
    
    %  if (j <  reax.num_atom_types && k <  reax.num_atom_types)
    reax.tbp(j,k).used = 1;
    reax.tbp(j,k).ffid_ps = a;
    reax.tbp(j,k).i = j;
    reax.tbp(j,k).j = k;
    
    val = str2num (strtmp{3});
    reax.tbp(j,k).De_s = val;
    reax.tbp(k,j).De_s = val;
    a=a+1; ffid_pointer(a) =  reax.tbp(j,k).De_s;
    val = str2num (strtmp{4});
    reax.tbp(j,k).De_p = val;
    reax.tbp(k,j).De_p = val;
    a=a+1; ffid_pointer(a) =  reax.tbp(j,k).De_p;
    val = str2num (strtmp{5});
    reax.tbp(j,k).De_pp = val;
    reax.tbp(k,j).De_pp = val;
    a=a+1; ffid_pointer(a) =  reax.tbp(j,k).De_pp;
    val = str2num (strtmp{6});
    reax.tbp(j,k).p_be1 = val;
    reax.tbp(k,j).p_be1 = val;
    a=a+1; ffid_pointer(a) =  reax.tbp(j,k).p_be1;
    val = str2num (strtmp{7});
    reax.tbp(j,k).p_bo5 = val;
    reax.tbp(k,j).p_bo5 = val;
    a=a+1; ffid_pointer(a) =  reax.tbp(j,k).p_bo5;
    val = str2num (strtmp{8});
    reax.tbp(j,k).v13cor = val;
    reax.tbp(k,j).v13cor = val;
    a=a+1; ffid_pointer(a) =  reax.tbp(j,k).v13cor;
    val = str2num (strtmp{9});
    reax.tbp(j,k).p_bo6 = val;
    reax.tbp(k,j).p_bo6 = val;
    a=a+1; ffid_pointer(a) =  reax.tbp(j,k).p_bo6;
    val = str2num (strtmp{10});
    reax.tbp(j,k).p_ovun1 = val;
    reax.tbp(k,j).p_ovun1 = val;
    a=a+1; ffid_pointer(a) =  reax.tbp(j,k).p_ovun1;
    
    %    line 2
    s=fgetl(ffid);
    strtmp=strsplit(strtrim(s));
    
    val = str2num (strtmp{1});
    reax.tbp(j,k).p_be2 = val;
    reax.tbp(k,j).p_be2 = val;
    a=a+1; ffid_pointer(a) =  reax.tbp(j,k).p_be2;
    val = str2num (strtmp{2});
    reax.tbp(j,k).p_bo3 = val;
    reax.tbp(k,j).p_bo3 = val;
    a=a+1; ffid_pointer(a) =  reax.tbp(j,k).p_bo3;
    val = str2num (strtmp{3});
    reax.tbp(j,k).p_bo4 = val;
    reax.tbp(k,j).p_bo4 = val;
    a=a+1; ffid_pointer(a) =  reax.tbp(j,k).p_bo4;
    val = str2num (strtmp{4});
    reax.tbp(j,k).reax.tbp_12 = val;
    reax.tbp(k,j).reax.tbp_12 = val;
    a=a+1; ffid_pointer(a) =  reax.tbp(j,k).reax.tbp_12;
    val = str2num (strtmp{5});
    reax.tbp(j,k).p_bo1 = val;
    reax.tbp(k,j).p_bo1 = val;
    a=a+1; ffid_pointer(a) =  reax.tbp(j,k).p_bo1;
    val = str2num (strtmp{6});
    reax.tbp(j,k).p_bo2 = val;
    reax.tbp(k,j).p_bo2 = val;
    a=a+1; ffid_pointer(a) =  reax.tbp(j,k).p_bo2;
    val = str2num (strtmp{7});
    reax.tbp(j,k).reax.tbp_15 = val;
    reax.tbp(k,j).reax.tbp_15 = val;
    a=a+1; ffid_pointer(a) =  reax.tbp(j,k).reax.tbp_15;
    val = str2num (strtmp{8});
    reax.tbp(j,k).reax.tbp_16 = val;
    reax.tbp(k,j).reax.tbp_16 = val;
    a=a+1; ffid_pointer(a) =  reax.tbp(j,k).reax.tbp_16;
    
    for (idx = 1 : 16)
        if (idx == 1)
            fprintf (fid, '%3d: %3d-%-2d\n', i ,  reax.tbp(j,k).i,  reax.tbp(j,k).j);
        end
        fprintf (fid, '[%4d]%8.4f ', idx + (i-1) * 16 + reax.gp_idx + reax.sbp_idx, ...
            ffid_pointer(idx + (i-1) * 16 + reax.gp_idx + reax.sbp_idx));
        if (idx == 8 || idx == 16)
            fprintf (fid, '\n');
        end
    end
end

fprintf (fid,'>> Done with 2-BODY interactions (%i types and %i parameters)\n', ...
    reax.num_bond_types, a);

%   %    calculating combination rules and filling up remaining fields.
%   %
%      for (i=0; i <  reax.num_atom_types; i++)
%      for (j=i; j <  reax.num_atom_types; j++)
%
%       reax.tbp(i,j).r_s = 0.5 * ( reax.sbp(i).r_s +  reax.sbp(j).r_s);
%       reax.tbp(j,i).r_s = 0.5 * ( reax.sbp(j).r_s +  reax.sbp(i).r_s);
%
%       reax.tbp(i,j).r_p = 0.5 * ( reax.sbp(i).r_pi +  reax.sbp(j).r_pi);
%       reax.tbp(j,i).r_p = 0.5 * ( reax.sbp(j).r_pi +  reax.sbp(i).r_pi);
%
%       reax.tbp(i,j).r_pp = 0.5 * ( reax.sbp(i).r_pi_pi +  reax.sbp(j).r_pi_pi);
%       reax.tbp(j,i).r_pp = 0.5 * ( reax.sbp(j).r_pi_pi +  reax.sbp(i).r_pi_pi);
%
%       reax.tbp(i,j).p_boc3 = sqrt( reax.sbp(i).b_o_132 *  reax.sbp(j).b_o_132);
%       reax.tbp(j,i).p_boc3 = sqrt( reax.sbp(j).b_o_132 *  reax.sbp(i).b_o_132);
%
%       reax.tbp(i,j).p_boc4 = sqrt( reax.sbp(i).b_o_131 *  reax.sbp(j).b_o_131);
%       reax.tbp(j,i).p_boc4 = sqrt( reax.sbp(j).b_o_131 *  reax.sbp(i).b_o_131);
%
%       reax.tbp(i,j).p_boc5 = sqrt( reax.sbp(i).b_o_133 *  reax.sbp(j).b_o_133);
%       reax.tbp(j,i).p_boc5 = sqrt( reax.sbp(j).b_o_133 *  reax.sbp(i).b_o_133);
%
%       reax.tbp(i,j).D = sqrt( reax.sbp(i).epsilon *  reax.sbp(j).epsilon);
%       reax.tbp(j,i).D = sqrt( reax.sbp(j).epsilon *  reax.sbp(i).epsilon);
%
%       reax.tbp(i,j).alpha = sqrt( reax.sbp(i).alpha *  reax.sbp(j).alpha);
%       reax.tbp(j,i).alpha = sqrt( reax.sbp(j).alpha *  reax.sbp(i).alpha);
%
%       reax.tbp(i,j).r_vdW = 2.0 * sqrt( reax.sbp(i).r_vdw *  reax.sbp(j).r_vdw);
%       reax.tbp(j,i).r_vdW = 2.0 * sqrt( reax.sbp(j).r_vdw *  reax.sbp(i).r_vdw);
%
%       reax.tbp(i,j).gamma_w = sqrt( reax.sbp(i).gamma_w *  reax.sbp(j).gamma_w);
%       reax.tbp(j,i).gamma_w = sqrt( reax.sbp(j).gamma_w *  reax.sbp(i).gamma_w);
%
%       reax.tbp(i,j).gamma = pow( reax.sbp(i).gamma *  reax.sbp(j).gamma,-1.5);
%       reax.tbp(j,i).gamma = pow( reax.sbp(j).gamma *  reax.sbp(i).gamma,-1.5);
%
%      %   additions for additional vdWaals interaction types - inner core
%       reax.tbp(i,j).rcore =  reax.tbp(j,i).rcore = sqrt(  reax.sbp(i).rcore2 *  reax.sbp(j).rcore2 );
%       reax.tbp(i,j).ecore =  reax.tbp(j,i).ecore = sqrt(  reax.sbp(i).ecore2 *  reax.sbp(j).ecore2 );
%       reax.tbp(i,j).acore =  reax.tbp(j,i).acore = sqrt(  reax.sbp(i).acore2 *  reax.sbp(j).acore2 );



%    next line is number of two body offdiagonal combinations and comments
%    these are two body offdiagonal terms that are different from the
%      combination rules defined above

s=fgetl(ffid);
fprintf (fid, '4: 2-BODY OFF-DIAGONAL PARAMETERS [%d]\n', reax.tbodp_idx);

for (i = 1 : reax.num_off_diag_types)
    s=fgetl(ffid);
    strtmp=strsplit(strtrim(s));
    j = str2num (strtmp{1});
    k = str2num (strtmp{2});
    
  %  if (j <=  reax.num_atom_types && k <=  reax.num_atom_types)
        reax.tbp(j,k).used_offdiag = 1;
        reax.tbp(j,k).od_ffid_ps = a;
        reax.tbp(j,k).i = j;
        reax.tbp(j,k).j = k;
        
        val = str2num (strtmp{3});
        %  if (val > 0.0)
        reax.tbp(j,k).D = val;
        reax.tbp(k,j).D = val;
        a=a+1; ffid_pointer(a) =  reax.tbp(j,k).D;
        %
        
        val = str2num (strtmp{4});
        %  if (val > 0.0)
        reax.tbp(j,k).r_vdW = val;
        reax.tbp(k,j).r_vdW = val;
        a=a+1; ffid_pointer(a) =  reax.tbp(j,k).r_vdW;
        %
        
        val = str2num (strtmp{5});
        %  if (val > 0.0)
        reax.tbp(j,k).alpha = val;
        reax.tbp(k,j).alpha = val;
        a=a+1; ffid_pointer(a) =  reax.tbp(j,k).alpha;
        %
        
        val = str2num (strtmp{6});
        %  if (val > 0.0)
        reax.tbp(j,k).r_s = val;
        reax.tbp(k,j).r_s = val;
        a=a+1; ffid_pointer(a) =  reax.tbp(j,k).r_s;
        %
        
        val = str2num (strtmp{7});
        %  if (val > 0.0)
        reax.tbp(j,k).r_p = val;
        reax.tbp(k,j).r_p = val;
        a=a+1; ffid_pointer(a) =  reax.tbp(j,k).r_p;
        %
        
        val = str2num (strtmp{8});
        %  if (val > 0.0)
        reax.tbp(j,k).r_pp = val;
        reax.tbp(k,j).r_pp = val;
        a=a+1; ffid_pointer(a) =  reax.tbp(j,k).r_pp;
        
        if (reax.lg_reax == 1)
            val = str2num (strtmp{9});
            %  if (val > 0.0)
            reax.tbp(j,k).lg_cij = val;
            reax.tbp(k,j).lg_cij = val;
            a=a+1; ffid_pointer(a) =  reax.tbp(j,k).lg_cij;
        end
  %  end
    
    for (idx = 1 : off_diag_params)
        if (idx == 1)
            fprintf (fid, '%3d: %3d-%-2d\n', i,  reax.tbp(j,k).i,  reax.tbp(j,k).j);
        end
        fprintf (fid, '[%4d]%8.4f ', idx + (i-1) * off_diag_params + reax.gp_idx + reax.sbp_idx + reax.tbp_idx, ...
            ffid_pointer(idx + (i-1) * off_diag_params + reax.gp_idx + reax.sbp_idx + reax.tbp_idx));
        if (idx == off_diag_params)
            fprintf (fid, '\n');
        end
    end
end

fprintf (fid, '>> Done with Off-diag interactions (%i types and %i parameters)\n', ...
    reax.num_off_diag_types, a);

%    3-body parameters -
%   supports multi-well potentials (upto MAX_3BODY_PARAM in mytypes.h)
%    clear entries first
for (i = 1 : reax.num_atom_types)
    for (j = 1 : reax.num_atom_types)
        for (k = 1 : reax.num_atom_types)
            reax.thbp(i,j,k).cnt = 1;
        end
    end
end

%    next line is number of 3-body params and some comments
s=fgetl(ffid);
fprintf (fid, '5: 3-BODY PARAMETERS [%d]\n', reax.thbp_idx);

for (i = 1 :  reax.num_angle_types)
    s=fgetl(ffid);
    strtmp=strsplit(strtrim(s));
    
    j = str2num (strtmp{1});
    k = str2num (strtmp{2});
    m = str2num (strtmp{3});
    
    %  if (j <  reax.num_atom_types && k <  reax.num_atom_types &&
    %  m <  reax.num_atom_types)
    reax.thbp(j,k,m).used = 1;
    reax.thbp(j,k,m).ffid_ps = a;
    reax.thbp(j,k,m).i = j;
    reax.thbp(j,k,m).j = k;
    reax.thbp(j,k,m).k = m;
    cnt =  reax.thbp(j,k,m).cnt;
    reax.thbp(j,k,m).cnt=reax.thbp(j,k,m).cnt+1;
    reax.thbp(m,k,j).cnt=reax.thbp(m,k,j).cnt+1;
    
    val = str2num (strtmp{4});
    reax.thbp(j,k,m).prm(cnt).theta_00 = val;
    reax.thbp(m,k,j).prm(cnt).theta_00 = val;
    a=a+1; ffid_pointer(a) =  reax.thbp(j,k,m).prm(cnt).theta_00;
    
    val = str2num (strtmp{5});
    reax.thbp(j,k,m).prm(cnt).p_val1 = val;
    reax.thbp(m,k,j).prm(cnt).p_val1 = val;
    a=a+1; ffid_pointer(a) =  reax.thbp(j,k,m).prm(cnt).p_val1;
    
    val = str2num (strtmp{6});
    reax.thbp(j,k,m).prm(cnt).p_val2 = val;
    reax.thbp(m,k,j).prm(cnt).p_val2 = val;
    a=a+1; ffid_pointer(a) =  reax.thbp(j,k,m).prm(cnt).p_val2;
    
    val = str2num (strtmp{7});
    reax.thbp(j,k,m).prm(cnt).p_coa1 = val;
    reax.thbp(m,k,j).prm(cnt).p_coa1 = val;
    a=a+1; ffid_pointer(a) =  reax.thbp(j,k,m).prm(cnt).p_coa1;
    
    val = str2num (strtmp{8});
    reax.thbp(j,k,m).prm(cnt).p_val7 = val;
    reax.thbp(m,k,j).prm(cnt).p_val7 = val;
    a=a+1; ffid_pointer(a) =  reax.thbp(j,k,m).prm(cnt).p_val7;
    
    val = str2num (strtmp{9});
    reax.thbp(j,k,m).prm(cnt).p_pen1 = val;
    reax.thbp(m,k,j).prm(cnt).p_pen1 = val;
    a=a+1; ffid_pointer(a) =  reax.thbp(j,k,m).prm(cnt).p_pen1;
    
    val = str2num (strtmp{10});
    reax.thbp(j,k,m).prm(cnt).p_val4 = val;
    reax.thbp(m,k,j).prm(cnt).p_val4 = val;
    a=a+1; ffid_pointer(a) =  reax.thbp(j,k,m).prm(cnt).p_val4;
    
    for (idx = 1 : 7)
        if (idx == 1)
            fprintf (fid, '%3d: %3d %2d %2d\n', i,  reax.thbp(j,k,m).i, ...
                reax.thbp(j,k,m).j,  reax.thbp(j,k,m).k);
        end
        fprintf (fid, '[%4d]%8.4f ', idx + (i-1) * 7 + reax.gp_idx + reax.sbp_idx + reax.tbp_idx + reax.tbodp_idx, ...
            ffid_pointer(idx + (i-1) * 7 + reax.gp_idx + reax.sbp_idx + reax.tbp_idx + reax.tbodp_idx));
        if (idx == 7)
            fprintf (fid, '\n');
        end
    end
end
fprintf (fid, '>> Done with three body interaction (%i types and %i parameters)\n', ...
    reax.num_angle_types, a);

%   %    4-body parameters are entered in compact form. i.e. 0-X-Y-0
%      correspond to any type of pair of atoms in 1 and 4
%      position. However, explicit X-Y-Z-W takes precedence over the
%      default description.
%      supports multi-well potentials (upto MAX_4BODY_PARAM in mytypes.h)
%      IMPORTANT: for now, directions on how to read multi-entries from ffield
%      is not clear

%    next line is number of 4-body params and some comments
s=fgetl(ffid);
fprintf (fid, '6: 4-BODY PARAMETERS [%d]\n', reax.fbp_idx);

for (i = 1 : reax.num_torsion_types)
    s=fgetl(ffid);
    strtmp=strsplit(strtrim(s));
    
    j = str2num (strtmp{1});
    k = str2num (strtmp{2});
    m = str2num (strtmp{3});
    n = str2num (strtmp{4});
    
    % simple style of 4-body parameters
    % if j and n equal 0 , give max atom type plus 1
    if (j == 0 && n == 0)
        j=reax.num_atom_types+1; n=reax.num_atom_types+1;
    end
    
    reax.fbp(j,k,m,n).used = 1;
    reax.fbp(j,k,m,n).ffid_ps = a;
    
    val = str2num (strtmp{5});
    reax.fbp(j,k,m,n).prm(1).V1 = val;
    reax.fbp(n,m,k,j).prm(1).V1 = val;
    a=a+1; ffid_pointer(a) =  reax.fbp(j,k,m,n).prm(1).V1;
    
    val = str2num (strtmp{6});
    reax.fbp(j,k,m,n).prm(1).V2 = val;
    reax.fbp(n,m,k,j).prm(1).V2 = val;
    a=a+1; ffid_pointer(a) =  reax.fbp(j,k,m,n).prm(1).V2;
    
    val = str2num (strtmp{7});
    reax.fbp(j,k,m,n).prm(1).V3 = val;
    reax.fbp(n,m,k,j).prm(1).V3 = val;
    a=a+1; ffid_pointer(a) =  reax.fbp(j,k,m,n).prm(1).V3;
    
    val = str2num (strtmp{8});
    reax.fbp(j,k,m,n).prm(1).p_tor1 = val;
    reax.fbp(n,m,k,j).prm(1).p_tor1 = val;
    a=a+1; ffid_pointer(a) =  reax.fbp(j,k,m,n).prm(1).p_tor1;
    
    val = str2num (strtmp{9});
    reax.fbp(j,k,m,n).prm(1).p_cot1 = val;
    reax.fbp(n,m,k,j).prm(1).p_cot1 = val;
    a=a+1; ffid_pointer(a) =  reax.fbp(j,k,m,n).prm(1).p_cot1;
    
    val = str2num (strtmp{10});
    reax.fbp(j,k,m,n).prm(1).reax.fbp_6 = val;
    reax.fbp(n,m,k,j).prm(1).reax.fbp_6 = val;
    a=a+1; ffid_pointer(a) =  reax.fbp(j,k,m,n).prm(1).reax.fbp_6;
    
    val = str2num (strtmp{11});
    reax.fbp(j,k,m,n).prm(1).reax.fbp_7 = val;
    reax.fbp(n,m,k,j).prm(1).reax.fbp_7 = val;
    a=a+1; ffid_pointer(a) =  reax.fbp(j,k,m,n).prm(1).reax.fbp_7;
    
    for (idx = 1 : 7)
        if (idx == 1)
            if (j == reax.num_atom_types+1 && n == reax.num_atom_types+1)
                j = 0; n = 0;
            end
            fprintf (fid, '%3d: %3d %2d %2d %2d: ', i, j, k, m, n);
        end
        fprintf (fid, '[%4d]%8.4f ', ...
            idx + (i-1) * 7 + reax.gp_idx + reax.sbp_idx + reax.tbp_idx + reax.tbodp_idx + reax.thbp_idx, ...
            ffid_pointer(idx + (i-1) * 7 + reax.gp_idx + reax.sbp_idx + reax.tbp_idx + reax.tbodp_idx + reax.thbp_idx));
        if (idx == 7)
            fprintf (fid, '\n');
        end
    end
end

fprintf (fid, '>> Done with 4-BODY interactions (%i types and %i parameters)\n', ...
    reax.num_torsion_types, a);
%    next line is number of hydrogen bond params and some comments
s=fgetl(ffid);
fprintf (fid, '7: H-BOND PARAMETERS [%d]\n', reax.hbp_idx);

reax.hbp=[];
for (i = 1 : reax.num_hbond_types)
    s=fgetl(ffid);
    strtmp=strsplit(strtrim(s));
    
    j = str2num (strtmp{1});
    k = str2num (strtmp{2});
    m = str2num (strtmp{3});
    
    %  if( j <  reax.num_atom_types && m <  reax.num_atom_types )
    reax.hbp(j,k,m).used = 1;
    reax.hbp(j,k,m).ffid_ps = a;
    reax.hbp(j,k,m).i = j;
    reax.hbp(j,k,m).j = k;
    reax.hbp(j,k,m).k = m;
    val = str2num (strtmp{4});
    reax.hbp(j,k,m).r0_hb = val;
    a=a+1; ffid_pointer(a) =  reax.hbp(j,k,m).r0_hb;
    
    val = str2num (strtmp{5});
    reax.hbp(j,k,m).p_hb1 = val;
    a=a+1; ffid_pointer(a) =  reax.hbp(j,k,m).p_hb1;
    
    val = str2num (strtmp{6});
    reax.hbp(j,k,m).p_hb2 = val;
    a=a+1; ffid_pointer(a) =  reax.hbp(j,k,m).p_hb2;
    
    val = str2num (strtmp{7});
    reax.hbp(j,k,m).p_hb3 = val;
    a=a+1; ffid_pointer(a) =  reax.hbp(j,k,m).p_hb3;
    
    for (idx = 1 : 4)
        if (idx == 1)
            fprintf (fid, '%3d: %3d %2d %2d: ', i,  reax.hbp(j,k,m).i,  reax.hbp(j,k,m).j, ...
                reax.hbp(j,k,m).k);
        end
        fprintf (fid, '[%4d]%8.4f ', ...
            idx + (i-1) * 4 + reax.gp_idx + reax.sbp_idx + reax.tbp_idx + reax.tbodp_idx + reax.thbp_idx + reax.fbp_idx, ...
            ffid_pointer(idx + (i-1) * 4 + reax.gp_idx + reax.sbp_idx + reax.tbp_idx + reax.tbodp_idx + reax.thbp_idx + ...
            reax.fbp_idx));
        if (idx == 4)
            fprintf (fid, '\n');
        end
    end
end
fprintf (fid, '>> Done with H-BOND interaction (%i types and %i parameters)\n', ...
    reax.num_hbond_types, a);
clear s strtmp;

fprintf (fid, '>> Read in %i parameters from ffield\n', reax.num_allp);

fclose(ffid);
fclose(fid);

save read_ff_mat file reax ffid_pointer
 %clear
% load read_ff_mat

end
