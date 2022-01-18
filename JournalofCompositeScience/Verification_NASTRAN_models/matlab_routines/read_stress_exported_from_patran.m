function Stress=read_stress_exported_from_patran(rpt_fname,minElmid,maxElmid,Layerlimit)

% ************************************************************************

% Stress output format
% -Entity ID--El. Pos. ID--X Component---Y Component---XY Component-


fid111=fopen(rpt_fname);
status=fseek(fid111,0,'eof');
EOF=ftell(fid111);
currentFPI=fseek(fid111,0,'bof');
% modeshape=zeros(maxmode,6);
modeflag=1;

while currentFPI<EOF && modeflag
    
    linef06=fgetl(fid111);currentFPI=ftell(fid111);
    str=findstr(linef06,'Layer Layer');
    
    if length(linef06)>65 && isempty(str)==0
        layer_num=str2num(linef06(81:end));
        
        linef06=fgetl(fid111);currentFPI=ftell(fid111);
        linef06=fgetl(fid111);currentFPI=ftell(fid111);
        linef06=fgetl(fid111);currentFPI=ftell(fid111);
        linef06=fgetl(fid111);currentFPI=ftell(fid111);
        
        jj=1; % The first line, jj=1, is the frequency line, set jj=2 if include frequency in the first line
        ElmID=minElmid;
        flag=1;
        
        if layer_num<=Layerlimit
            
            disp(['----------------- Reading the stress of Layer #' num2str(layer_num) '----------------']);
            
            % uncomment this code if include frequency
            %             eval(['modeshape' num2str(mode) '(1,1:7)=[NaturalFrequency(' num2str(mode) '),0,0,0,0,0,0];']);
            
            while currentFPI<EOF && flag && ElmID<maxElmid
                
                linef06=fgetl(fid111);currentFPI=ftell(fid111);
%               linef06
                ElmID=str2num(linef06(1:12));
                
                
                if isempty(ElmID)==1
                    flag=0;
                else
                    %                 eval(['modeshape' num2str(jj) '=mode' num2str(jj) '(5)']);
                    
%                     str2num(linef06(25:end))
                   Stress(jj,1:4,layer_num) = [ElmID str2num(linef06(25:end))];
                    
                    
                    %                     eval(['modeshape' num2str(mode_num) '(' num2str(jj) ',1:4)=[nodeid,str2num(linef06(20:72))];']);
                    
%                     linef06=fgetl(fid111);currentFPI=ftell(fid111);
%                     ModeShape.modeshape(jj,5:7,layer_num) = str2num(linef06(20:72));
                    %                     eval(['modeshape' num2str(mode_num) '(' num2str(jj) ',5:7)=str2num(linef06(20:72));']);
                    
                    jj=jj+1;
                end
            end
            %             eval(['ModeShape.modeshape' num2str(mode_num) '=modeshape' num2str(mode_num) ';']);
            
        else
            modeflag=0;
        end
    end
end
fclose(fid111);
