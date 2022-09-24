
%clear
%diary('ssss.log')
%diary on
dirs=dir();
[num,l]=size(dirs);
% s=0;
for i = 3: num
    if dirs(i).isdir == 1
        dir_s=dirs(i).name;
        cd(dir_s)
        system(['for input in `ls in.reaxc.*`;' ...
                'do INPUT_FILE=${input};' ... 
                'OUTPUT_FILE=${input:9}.lammps.log;' ...
                ' mpirun -np 20 lmp_mpi < $INPUT_FILE > ${OUTPUT_FILE};' ...
                ' done'])
        cd ..
         fprintf( '%s\n', dirs(i).name);
        %  s=s+1;
    end
end
% diary off;
%clear
