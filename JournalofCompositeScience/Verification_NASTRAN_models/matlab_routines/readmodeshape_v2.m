function ModeShape=readmodeshape_v2(pchfname,minodeid,maxnodeid,Modelimit)

% ************************************************************************
% Read mode shape from *pch file
% Delete the last blank lines

% ModeShape            --------   structure data
% ModeShape.modeshape1 --------   mode shape of the mode 1
% ...
% ModeShape.modeshapeN --------   mode shape of the mode n

fid111=fopen(pchfname);
status=fseek(fid111,0,'eof');
EOF=ftell(fid111);
currentFPI=fseek(fid111,0,'bof');
% modeshape=zeros(maxmode,6);
modeflag=1;

while currentFPI<EOF && modeflag
    
    linef06=fgetl(fid111);currentFPI=ftell(fid111);
    str=findstr(linef06,'MODE');
    
    if length(linef06)>70 && isempty(str)==0
        mode_num=str2num(linef06(40:42));
        
        jj=1; % The first line, jj=1, is the frequency line, set jj=2 if include frequency in the first line
        nodeid=minodeid;
        flag=1;
        
        if mode_num<=Modelimit
            
            disp(['----------------- Reading the mode shape of Mode #' num2str(mode_num) '----------------']);
            
            % uncomment this code if include frequency
            %             eval(['modeshape' num2str(mode) '(1,1:7)=[NaturalFrequency(' num2str(mode) '),0,0,0,0,0,0];']);
            
            while currentFPI<EOF && flag && nodeid<=maxnodeid
                linef06=fgetl(fid111);currentFPI=ftell(fid111);
                nodeid=str2num(linef06(1:14));
                
                if isempty(nodeid)==1
                    flag=0;
                else
                    %                 eval(['modeshape' num2str(jj) '=mode' num2str(jj) '(5)']);
                    
                    
                    ModeShape.modeshape(jj,1:4,mode_num) = [nodeid,str2num(linef06(20:72))];
                    
                    
                    %                     eval(['modeshape' num2str(mode_num) '(' num2str(jj) ',1:4)=[nodeid,str2num(linef06(20:72))];']);
                    
                    linef06=fgetl(fid111);currentFPI=ftell(fid111);
                    ModeShape.modeshape(jj,5:7,mode_num) = str2num(linef06(20:72));
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
