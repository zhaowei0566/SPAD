function Beam_stress =read_cbeam_stress(pch_filename,minElemid,maxElemid)

% ************************************************************************
% Read mode shape from *pch file
% Delete the last blank lines

% ModeShape            --------   structure data
% ModeShape.modeshape1 --------   mode shape of the mode 1
% ...
% ModeShape.modeshapeN --------   mode shape of the mode n

fid111=fopen(pch_filename);
status=fseek(fid111,0,'eof');
EOF=ftell(fid111);
currentFPI=fseek(fid111,0,'bof');
% modeshape=zeros(maxmode,6);



jj=1;
elemid = minElemid;
flag = 1;
while currentFPI<EOF
    
    linef06=fgetl(fid111);currentFPI=ftell(fid111);
    
    str=findstr(linef06,'$ELEMENT TYPE =           2  BEAM');
    
    if  isempty(str)==0
        
        
        while currentFPI<EOF && flag
            
            linef06=fgetl(fid111);currentFPI=ftell(fid111);
            
            elemid =str2num(linef06(1:16));
            
            if isempty(elemid) == 0
                
                if elemid<=maxElemid && elemid>=minElemid
                    
                    temp = str2num(linef06)
                    
                    
                    if isempty(temp) == 0
                        Beam_stress(jj,:) = str2num(linef06);
                        jj = jj+1;
                    else
                        
                        flag = 0;
                    end
%                 else
%                     flag = 0;
                end
                
            end
        end
        
        
        %          while currentFPI<EOF && flag && elemid<=maxElemid
        %              linef06=fgetl(fid111);currentFPI=ftell(fid111);
        %
        %              con_flag = strcmp(linef06(1:6),'-CONT-');
        %
        %
        %
        %          end
        
        
    end
    
    
    
    
    %
    %     if length(linef06)>=81 && isempty(str)==0
    %
    %         linef06=fgetl(fid111);currentFPI=ftell(fid111);
    %
    %         elemid =str2num(linef06(1:16));
    %
    %
    %  flag = 1;
    %         % uncomment this code if include frequency
    %         %             eval(['modeshape' num2str(mode) '(1,1:7)=[NaturalFrequency(' num2str(mode) '),0,0,0,0,0,0];']);
    %
    %         while currentFPI<EOF && flag && elemid<=maxElemid
    %             linef06=fgetl(fid111);currentFPI=ftell(fid111);
    %             con_flag = strcmp(linef06(1:6),'-CONT-');
    %
    %             if isempty(con_flag )==1
    %                 flag=0;
    %             else
    %                 %                 eval(['modeshape' num2str(jj) '=mode' num2str(jj) '(5)']);
    %
    %
    %                 linef06=fgetl(fid111);currentFPI=ftell(fid111);
    %
    %                 linef06=fgetl(fid111);currentFPI=ftell(fid111);
    %                 str2num(linef06)
    %                 Beam_stress(jj,:) = str2num(linef06);
    %                 %                     eval(['modeshape' num2str(mode_num) '(' num2str(jj) ',5:7)=str2num(linef06(20:72));']);
    %
    %                 jj=jj+1;
    %             end
    %         end
    %         %             eval(['ModeShape.modeshape' num2str(mode_num) '=modeshape' num2str(mode_num) ';']);
    %
    %     else
    %         modeflag=0;
    %     end
end

fclose(fid111);
