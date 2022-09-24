% get trainset data and make in file

        % construct ff file and start lammps initialization

  function f = evaluate_objective(x, M, V)

%% function f = evaluate_objective(x, M, V)
% Function to evaluate the objective functions for the given input vector
% x. x is an array of decision variables and f(1), f(2), etc are the
% objective functions. The algorithm always minimizes the objective
% function hence if you would like to maximize the function then multiply
% the function by negative one. M is the numebr of objective functions and
% V is the number of decision variables. 
%
% This functions is basically written by the user who defines his/her own
% objective function. Make sure that the M and V matches your initial user
% input. Make sure that the 
%
% An example objective function is given below. It has two six decision
% variables are two objective functions.
% change ffid_pointer to new chromosome
ffid_pointer_new=x(1:V);
fflmp='ffield.reax.opted';
write_reaxc_ff_file(ffid_pointer_new,fflmp);

% ! cp ffield ./lg_trainset/
traindir=dir();
[num,l]=size(traindir);
for i = 3: num
    if traindir(i).isdir == 1
        dir_s=traindir(i).name;
        path=[dir_s,'\',fflmp];
        copyfile(fflmp,path);
        %  s=s+1;
    end
end

% start lammps to caculate
% lmp_startup;
% extract lammps result
traindir=dir();
[num,l]=size(traindir);
for i = 3: num
    if traindir(i).isdir == 1
        dir_s=traindir(i).name;
       cd(dir_s);
       
       % Replace in file potential function file
      uu=dir(['in.reaxc','.*']);
        %  s=s+1;
        len=length(uu);
        tmp_in='tmp_in';
         % a=0;
         for k=1:len
           ssid  = fopen(uu(k).name,'r+');
           aaid  = fopen(tmp_in,'w+');
           while ~feof(ssid)
             tl = fgetl(ssid);
             id = findstr(tl,'pair_coeff');
              if ~isempty(id)
                  fprintf(aaid,'pair_coeff      * * %s ${elem}\n',fflmp);
              else
                  fprintf(aaid,'%s\n',tl);
              end
           %   a=a+1;
             % fprintf('a=%f\n',a);             
          end;
            fclose(ssid);
             fclose(aaid); 
            xx= [uu(k).name];
            delete(xx);
            command = ['rename',32,tmp_in,32,xx];
            system(command);
         end
%        system(['for input in `ls in.reaxc.*`;' ...
%                 'do INPUT_FILE=${input};' ... 
%                 'OUTPUT_FILE=${input:9}.lammps.log;' ...
%                 ' mpirun -np 20 lmp_mpi < $INPUT_FILE > ${OUTPUT_FILE};' ...
%                 ' done'])
%          fprintf( '%s\n', dirs(i).name);
        ss=dir(['*.','lammps.log']);
        %  s=s+1;
        len=length(ss);
for k=1:len
logdata=readlog(ss(k).name);
logave(k).name=ss(k).name(1:end-11);
logdata.relax=str2num(logdata.data{1});
logdata.data_num=str2num(logdata.data{2});
idx.step=find(strcmp(logdata.Chead(2,:),'Step'));
idx.density=find(strcmp(logdata.Chead(2,:),'Density'));
idx.cella=find(strcmp(logdata.Chead(2,:),'Cella'));
idx.cellb=find(strcmp(logdata.Chead(2,:),'Cellb'));
idx.cellc=find(strcmp(logdata.Chead(2,:),'Cellc'));

ave(k).endstep=length(logdata.data_num(:,1));
% extract last 100 rows and average
ave(k).num=100;
ave(k).density=mean(logdata.data_num(ave(k).endstep-ave(k).num+1:1:ave(k).endstep,idx.density));
ave(k).cella=mean(logdata.data_num(ave(k).endstep-ave(k).num+1:1:ave(k).endstep,idx.cella));
ave(k).cellb=mean(logdata.data_num(ave(k).endstep-ave(k).num+1:1:ave(k).endstep,idx.cellb));
ave(k).cellc=mean(logdata.data_num(ave(k).endstep-ave(k).num+1:1:ave(k).endstep,idx.cellc));
end

lgfile = 'lg_train';
fid = fopen(lgfile,'r');
loop = 1;
% first row
s = fgetl(fid);
strtmp=strsplit(strtrim(s));
lgdata.title=strtmp;
% name supercella supercellb supercellc density cella cellb cellc weight
while feof(fid) == 0
    s = fgetl(fid);
   strtmp=strsplit(strtrim(s));
   [row,cl]=size(strtmp);
   lgdata(loop).name=strtmp{1};
   data=strrep(s,lgdata(loop).name,'');
   lgdata(loop).data=data;
  loop=loop+1;
end
fclose(fid);

for j=1:loop-1

tmp=str2num(lgdata(j).data);
spcella(j)=tmp(1);
spcellb(j)=tmp(2);
spcellc(j)=tmp(3);
density(j)=tmp(4);
cella(j)=tmp(5);
cellb(j)=tmp(6);
cellc(j)=tmp(7);

% if weight is given, use the weight, if not ,weight is equal.
weight(j)=tmp(8);
if weight(j) ~=0
    weight(j)=weight(j);
else
    weight(j)=1.0/len;
end

%lammps data
ave(j).density=ave(j).density;
ave(j).cella=ave(j).cella/spcella(j);
ave(j).cellb=ave(j).cellb/spcellb(j);
ave(j).cellc=ave(j).cellc/spcellc(j);
end
f(1)=0;
f(2)=0;
for z=1:len
    filename=logave(z).name;
    filename_postion=find(strcmp({lgdata.name},filename));
    g= filename_postion;
    % logave num is z; lg num is g
    % density error function
    f(1)=f(1)+((density(g)-ave(z).density)/weight(g))^2;
    % cella cellb cellc
     f(2)=f(2)+(((cella(g)-ave(z).cella)+ (cellb(g)-ave(z).cellb) +(cellc(g)-ave(z).cellc))/weight(g))^2;
    
end
cd ..
    end
end

% f = [];
% %% Objective function one
% % Decision variables are used to form the objective function.
% f(1) = 1 - exp(-4*x(1))*(sin(6*pi*x(1)))^6;
% sum = 0;
% for i = 2 : 6
%     sum = sum + x(i)/4;
% end
% %% Intermediate function
% g_x = 1 + 9*(sum)^(0.25);
% 
% %% Objective function two
% f(2) = g_x*(1 - ((f(1))/(g_x))^2);

%% Kursawe proposed by Frank Kursawe.
% Take a look at the following reference
% A variant of evolution strategies for vector optimization.
% In H. P. Schwefel and R. Männer, editors, Parallel Problem Solving from
% Nature. 1st Workshop, PPSN I, volume 496 of Lecture Notes in Computer 
% Science, pages 193-197, Berlin, Germany, oct 1991. Springer-Verlag. 
%
% Number of objective is two, while it can have arbirtarly many decision
% variables within the range -5 and 5. Common number of variables is 3.
% Objective function one
%{
sum = 0;
for i = 1 : V - 1
    sum = sum - 10*exp(-0.2*sqrt((x(i))^2 + (x(i + 1))^2));
end
% Decision variables are used to form the objective function.
%}


% Objective function two
%{
sum = 0;
for i = 1 : V
    sum = sum + (abs(x(i))^0.8 + 5*(sin(x(i)))^3);
end
% Decision variables are used to form the objective function.
%}

%% Check for error
if length(f) ~= M
    error('The number of decision variables does not match you previous input. Kindly check your objective function');
end
