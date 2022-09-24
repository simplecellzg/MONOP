uu=dir(['in.reaxc','.*']);
fflmp='ffield.reax.opted';
        %  s=s+1;
        len=length(uu);
        tmp_in='tmp_in';
        a=0;
         for k=1:len
           ssid  = fopen(uu(k).name,'r+');
           aaid  = fopen(tmp_in,'w+');
             R = 'aaa';
           while ~feof(ssid)
             tl = fgetl(ssid);
             id = findstr(tl,'pair_coeff');
              if ~isempty(id)
                  fprintf(aaid,'pair_coeff      * * %s ${elem}\n',fflmp);
              else
                  fprintf(aaid,'%s\n',tl);
              end
             a=a+1;
             fprintf('a=%f\n',a);             
          end;
            fclose(ssid);
             fclose(aaid);
             
            xx= [uu(k).name];
            delete(xx);
            command = ['rename',32,tmp_in,32,xx];
            system(command);
         end
        