% read parameters data from ascii file
% clear
% path='ffield_cl';
% num_allp=738;
% fid=fopen(path, 'r');
%  for  i=1:num_allp
%       ss=fgetl(fid);
%       dt=sscanf(ss,'%f');
%       ffid_pointer(i,1) = dt; % ffid_pointer is all parameters in single column
%  end
%  fclose(fid);
% % read atom types from ascii file
% path='atom_types_cl';
% reax.num_atom_types=5;
% fid=fopen(path, 'r');
%   for  i=1:reax.num_atom_types
%       ss=fgetl(fid);
%      %  dt=sscanf(ss,'%3s');
%       reax.sbp_name(i,:)=ss; % reax.sbp_name is atom type in single column
%  end
%  fclose(fid);

%  gp_n_global=39;
%  lg_flag='nlg';

function write_reaxc_ff_file(ffid_pointer_new,fflmp)

load read_ff_mat
file=file;
ffdata=reax;
ffid_pointer=ffid_pointer_new;
token=0;
% write ff file in ascii
% fflmp=ffname;
% fflmp='ffield.reax.opted';
file.fflmp=fflmp;
% create a ff file
fid=fopen(fflmp,'w+');
fclose(fid);
fid=fopen(fflmp,'a+');
fprintf (fid, 'Reactive MD-fidrce field\n');
fprintf (fid, '%3d       ! Number of general parameters\n', reax.gp.n_global);
for i=1:reax.gp.n_global
    token=token+1; fprintf (fid, '%10.4f !Comment here\n', ffid_pointer(i));
end
% 1-body data
fprintf (fid, '%3d    !Nr of atoms; cov.r; valency;a.m;Rvdw;Evdw;gammaEEM;cov.r2;\n', reax.num_atom_types);
fprintf (fid, '            alfa;gammavdW;valency;Eunder;Eover;chiEEM;etaEEM;n.u.\n');
fprintf (fid, '            cov r3;Elp;Heat inc.;n.u.;n.u.;n.u.;n.u.\n');
fprintf (fid, '            ov/un;val1;n.u.;val3,vval4\n');
for i =1:reax.num_atom_types
    fprintf (fid, '%3s', reax.sbp(i).name);
    
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
    token=token+1;  fprintf (fid, '%9.4f\n', ffid_pointer(token));
    
    token=token+1;  fprintf (fid, '   %9.4f', ffid_pointer(token));
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));	% not used
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
    token=token+1;  fprintf (fid, '%9.4f\n', ffid_pointer(token));
    
    token=token+1;  fprintf (fid, '   %9.4f', ffid_pointer(token));
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));	% not used
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));	% not used
    token=token+1;  fprintf (fid, '%9.4f\n', ffid_pointer(token));	% not used
    
    token=token+1;  fprintf (fid, '   %9.4f', ffid_pointer(token));
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));	% not used
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
    token=token+1;  fprintf (fid, '%9.4f\n', ffid_pointer(token));
    
    if reax.lg_reax == 1
        token=token+1;    fprintf (fid, '   %9.4f', ffid_pointer(token));
        token=token+1;    fprintf (fid, '%9.4f\n', ffid_pointer(token));
    end
end


% 2-body data
fprintf (fid, '%3d      ! Nr of bonds; Edis1;LPpen;n.u.;pbe1;pbo5;13corr;pbo6\n', reax.num_bond_types);
fprintf (fid, '                         pbe2;pbo3;pbo4;n.u.;pbo1;pbo2;ovcorr\n');
for (i = 1 : reax.num_atom_types)
    for (j = 1 :  reax.num_atom_types)
        %  Avoid subscript overflow
        [a,b]=size(reax.tbp);
        if (i > a || j > b)
            reax.tbp(i,j).used=[];
        end
        if (reax.tbp(i,j).used)
            token = reax.tbp(i,j).ffid_ps;
            fprintf (fid, '%3d%3d', reax.tbp(i,j).i, reax.tbp(i,j).j);
            
            token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
            token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
            token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
            token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
            token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
            token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
            token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
            token=token+1;  fprintf (fid, '%9.4f\n', ffid_pointer(token));
            
            token=token+1;  fprintf (fid, '      %9.4f', ffid_pointer(token));
            token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
            token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
            token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
            token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
            token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
            token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
            token=token+1;  fprintf (fid, '%9.4f\n', ffid_pointer(token));
        end
    end
end




% off-diagonal
fprintf (fid, '%3d    ! Nr of off-diagonal terms; Ediss;Ro;gamma;rsigma;rpi;rpi2\n', ...
    reax.num_off_diag_types);
for (i = 1 : reax.num_atom_types)
    for (j = 1 : reax.num_atom_types)
        [a,b]=size(reax.tbp);
        if (i > a || j > b)
            reax.tbp(i,j).used=[];
        end
        if (reax.tbp(i,j).used_offdiag)
            token = reax.tbp(i,j).od_ffid_ps;
            
            fprintf (fid, '%3d%3d', reax.tbp(i,j).i, reax.tbp(i,j).j);
            token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
            token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
            token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
            token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
            token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));  
            if reax.lg_reax == 1
                token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
                token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
            end
            fprintf (fid, '\n');
        end
    end
end




%  3-body
fprintf (fid, '%3d    ! Nr of angles;at1;at2;at3;Thetao,o;ka;kb;pv1;pv2\n', ...
    reax.num_angle_types);
for (i = 1 : reax.num_atom_types)
    for (j = 1 : reax.num_atom_types)
        for (k = 1 :  reax.num_atom_types)
            [a,b,c]=size(reax.thbp);
            if (i > a || j > b || k >c)
                reax.thbp(i,j,k).used=[];
            end
            if (reax.thbp(i,j,k).used)
                token = reax.thbp(i,j,k).ffid_ps;
                
                fprintf (fid, '%3d%3d%3d', reax.thbp(i,j,k).i, reax.thbp(i,j,k).j, ...
                    reax.thbp(i,j,k).k);
                
                token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
                token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
                token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
                token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
                token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
                token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
                token=token+1;  fprintf (fid, '%9.4f\n', ffid_pointer(token));
            end
        end
    end
end


%  4-body
%   4/15/2015: thanks to Byung-Hyun Kim (Uppsala University) for reporting a bug when
%   the number of atoms is equal to the number of torsions in the force field file.
%   The bug appeared due to the support for compact torsion representation, and it was fixed
%   by looping over reax.num_atom_types+1 for i,j,k, and l
fprintf (fid, '%3d    ! Nr of torsions;at1;at2;at3;at4;;V1;V2;V3;V2(BO);vconj;n.u;n\n', ...
    reax.num_torsion_types);
for (i = 1 : reax.num_atom_types+1)
    for (j = 1 : reax.num_atom_types+1 )
        for (k = 1 : reax.num_atom_types+1)
            for (l = 1 : reax.num_atom_types+1)
                [a,b,c,d]=size(reax.fbp);
                if (i > a || j > b || k > c || l > d)
                    reax.fbp(i,j,k,l).used=[];
                end
                if (reax.fbp(i,j,k,l).used)
                    token = reax.fbp(i,j,k,l).ffid_ps;
                    % ir lr for i l replace
                    ir=i; lr=l;
                    if (i == reax.num_atom_types+1 && l == reax.num_atom_types+1)
                        ir = 0; lr = 0;
                    end
                    fprintf (fid, '%3d%3d%3d%3d', ir, j, k, lr);
                    
                    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
                    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
                    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
                    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
                    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
                    token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
                    token=token+1;  fprintf (fid, '%9.4f\n', ffid_pointer(token));
                end
            end
        end
    end
end






%  h-bond
fprintf (fid, '%3d    ! Nr of hydrogen bonds;at1;at2;at3;Rhb;Dehb;vhb1\n', ...
    reax.num_hbond_types);
for (i = 1 : reax.num_atom_types)
    for (j = 1 :  reax.num_atom_types)
        for (k = 1 :  reax.num_atom_types)
            [a,b,c]=size(reax.hbp);
            if (i > a || j > b || k >c)
                reax.hbp(i,j,k).used=[];
            end
            if (reax.hbp(i,j,k).used)
                token = reax.hbp(i,j,k).ffid_ps;
                fprintf (fid, '%3d%3d%3d', reax.hbp(i,j,k).i, reax.hbp(i,j,k).j, ...
                    reax.hbp(i,j,k).k);
                
                token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
                token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
                token=token+1;  fprintf (fid, '%9.4f', ffid_pointer(token));
                token=token+1;  fprintf (fid, '%9.4f\n', ffid_pointer(token));
            end
        end
    end
end
fclose(fid);
save write_ff_mat file reax ffid_pointer
clear



