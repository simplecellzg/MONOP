% main
clear
%% reading ff_file from given potential file; use lgvdw =1,not use lgvdw=0;
% set log file for watch;
ff='ffield.reax.weak2017m';
lgvdw=1;
log='logfile_m';
[ffid_pointer,file,reax]=read_reaxc_ff_file(ff,log,lgvdw);
fprintf ('>> Finished reading force field\n');



%% reading parameters that needs to be changed;
file=file;
ffdata=reax;
ffid_pointer=ffid_pointer;
paramfile='params';
file.params=paramfile;
params=parse_params(ffid_pointer,file,ffdata);
fprintf ('>> Finished reading params\n');


%% extract trainset lammps in file
fo=fopen(file.log,'a+');
fprintf(fo, '3, READING TRAIN FILE ...\n');
dirs=dir();
[num,l]=size(dirs);
for i = 3: num
    if dirs(i).isdir == 1
        dir_s=dirs(i).name;
        cd(dir_s)
        ss=dir(['*.','lammps.log']);
        len=length(ss);
         cd ..
         fprintf(fo, '-->in %s dir, %3d train file\n', dirs(i).name,len);
        %  s=s+1;
    end
end
fprintf('>> Finished reading train file\n');


%% reading is done start to use NSGA2
% pop - Population size
% gen - Total number of generations
pop=21;
gen=10;
flag.density=1;
flag.cell=1;
nsga_2(pop,gen,flag,ffid_pointer,params,file);
fprintf('>> Finished NSGA2 optimization\n');







