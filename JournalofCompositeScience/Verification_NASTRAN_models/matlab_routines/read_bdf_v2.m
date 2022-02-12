function [QUAD4_Elem,TRIA3_Elem,BEAM3_Elem,BEAM2_Elem,bdf_node,bdf_cord,...
    QUAD8_Elem]=read_bdf_v2(filename)

% this function is to read nodes and element connectivity and save in one
% file bdf_node_elem
% Revise the read of scientific data; March-19-2016
bdf_elem=[];
QUAD4_Elem=[];
TRIA3_Elem=[];
BEAM3_Elem=[];
BEAM2_Elem=[];
bdf_node=[];
bdf_cord=[];
QUAD8_Elem=[];
fid100=fopen(filename);

if fid100>0
    % string1='CQUAD4';string2='GRID';string3='CBAR';
    
    status=fseek(fid100,0,'eof');
    EOF=ftell(fid100);
    currentFPI=fseek(fid100,0,'bof');
    
    
    quad4_elem_num = 1;
    quad8_elem_num = 1;
    tria3_elem_num = 1;
    beam3_elem_num = 1;
    beam2_elem_num = 1;
    
    node_num=1;
    cord_num=1;
    
    bdf_cord=zeros(1,10);
    
    while currentFPI<EOF
        
        linebdf=fgetl(fid100);currentFPI=ftell(fid100);
        
        %         linebdf
        
        if length(linebdf)>=7
            %
            if strcmp(linebdf(1:6),'CQUAD4')
                
                bdf_elem_bdf{quad4_elem_num}=linebdf;
                bdf_elem_temp=linebdf;
                QUAD4_Elem(quad4_elem_num,:)=str2num(bdf_elem_temp(9:end));
                
                quad4_elem_num=quad4_elem_num+1;
                
            elseif strcmp(linebdf(1:6),'CTRIA3')
                
                bdf_elem_bdf{tria3_elem_num}=linebdf;
                bdf_elem_temp=linebdf;
                TRIA3_Elem(tria3_elem_num,:)=str2num(bdf_elem_temp(9:end));
                
                tria3_elem_num = tria3_elem_num+1;
                
                
            elseif strcmp(linebdf(1:6),'CBEAM3')
                
                bdf_elem_bdf{beam3_elem_num}=linebdf;
                bdf_elem_temp=linebdf;
                BEAM3_Elem(beam3_elem_num,:)=str2num(bdf_elem_temp(9:54));
                
                beam3_elem_num = beam3_elem_num+1;
                
                
            elseif strcmp(linebdf(1:5),'CBEAM')
                
                bdf_elem_bdf{beam2_elem_num}=linebdf;
                bdf_elem_temp=linebdf;
                BEAM2_Elem(beam2_elem_num,:)=str2num(bdf_elem_temp(9:40));
                
                beam2_elem_num = beam2_elem_num+1;
                
                
            elseif strcmp(linebdf(1:6),'CQUAD8')
                %                 bdf_elem_bdf{quad8_elem_num}=linebdf;
                %                 bdf_elem_temp=linebdf;
                node1 = str2num(linebdf(9:end));
                
                 linebdf=fgetl(fid100);currentFPI=ftell(fid100);
                  node2 = str2num(linebdf(9:end));
                 
                
                QUAD8_Elem(quad8_elem_num,:)= [node1(3:end) node2];
                
                quad8_elem_num=quad8_elem_num+1;
                
            end
            
            
            %
            if strcmp(linebdf(1:4),'GRID')
                
                bdf_node(node_num,1)=str2num(linebdf(9:16)); % node ID
                
                % replace e with blank
                
              e_places=  strfind(linebdf,'e');
                
              if ~isempty(e_places)
                linebdf(e_places) = blanks(length(e_places));
              end
                %% X
%                 if strfind(linebdf(25:32),'-')% THE FIRST negative sign could be a negative value
%                     
%                     %                 disp('============= yes =============');
%                     negative_loc=strfind(linebdf(25:32),'-');
%                     
%                     bdf_node(node_num,2)=str2num(linebdf(25:24+ negative_loc-1))*10^(-str2num(linebdf(24+ negative_loc+1:32)));
%                     
%                 else
%                     
%                     bdf_node(node_num,2)=str2num(linebdf(25:32)); % node Y cord
%                     
%                 end
                
                                if strfind(linebdf(25:32),'-')% THE FIRST negative sign could be a negative value
                    %                     linebdf(33:40)
                    
                    % ~isspace(linebdf(33:40));
                    %
                    negative_sign_id = strfind(linebdf(25:32),'-');
                    %
                    %
                    %                     %                 disp('============= yes =============');
                    
                    
                    temp_string = strtrim(linebdf(25:32));
                    
                    negative_loc=strfind(temp_string ,'-');
                    
                    if length(negative_loc)>1
                        
                        bdf_node(node_num,2)=str2num(linebdf(24+1:24+ negative_sign_id(2)-1))*10^(-str2num(linebdf(24+ negative_sign_id(2)+1:32)));
                        
                    elseif length(negative_loc)==1 && negative_loc>1
                        %
                        bdf_node(node_num,2)=str2num(linebdf(24+1:24+ negative_loc-1))*10^(-str2num(linebdf(24+ negative_loc+1:32)));
                    else
                        
                          bdf_node(node_num,2)=str2num(linebdf(25:32)); % node Y cord
                        
                    end
                    
                else
%                     linebdf
                    bdf_node(node_num,2)=str2num(linebdf(25:32)); % node Y cord
                    
                end
                %% Y
                if strfind(linebdf(33:40),'-')% THE FIRST negative sign could be a negative value
                    %                     linebdf(33:40)
                    
                    % ~isspace(linebdf(33:40));
                    %
                    negative_sign_id = strfind(linebdf(33:40),'-');
                    %
                    %
                    %                     %                 disp('============= yes =============');
                    
                    
                    temp_string = strtrim(linebdf(33:40));
                    
                    negative_loc=strfind(temp_string ,'-');
                    
                    if length(negative_loc)>1
                        

                        bdf_node(node_num,3)=str2num(linebdf(32+1:32+ negative_sign_id(2)-1))*10^(-str2num(linebdf(32+ negative_sign_id(2)+1:40)));
                        
                    elseif length(negative_loc)==1 && negative_loc>1
                        %
                        bdf_node(node_num,3)=str2num(linebdf(32+1:32+ negative_loc-1))*10^(-str2num(linebdf(32+ negative_loc+1:40)));
                    else
                        
                          bdf_node(node_num,3)=str2num(linebdf(33:40)); % node Y cord
                        
                    end
                    
                else
                    
                    bdf_node(node_num,3)=str2num(linebdf(33:40)); % node Y cord
                    
                end
                
                %                 linebdf
                %                 str2num(linebdf(41:end))
                
                %%  Z
                
                string_total_length = length(linebdf);
                
                if string_total_length>48
                    z_end_string_length = 48;
                else
                    z_end_string_length = string_total_length;
                end
                
                
                if strfind(linebdf(41: z_end_string_length ),'-')
                    
                    negative_sign_id = strfind(linebdf(41: z_end_string_length ),'-');
                    
                    temp_string = strtrim(linebdf(41: z_end_string_length ));
                    negative_loc=strfind(temp_string ,'-');
                    
                    
                    if length(negative_loc)>1
                        bdf_node(node_num,4)=str2num(linebdf(41:40+negative_sign_id(2)-1))*10^(-str2num(linebdf(40+ negative_sign_id(2)+1: z_end_string_length )));
                        
                    elseif length(negative_loc)==1 && negative_loc>1
%                      
%                         linebdf
%                         str2num(linebdf(41:40+negative_loc-1))*10^(-str2num(linebdf(40+ negative_loc+1: z_end_string_length )))
                        bdf_node(node_num,4)=str2num(linebdf(41:40+negative_loc-1))*10^(-str2num(linebdf(40+ negative_loc+1: z_end_string_length )));
                    else
                        
                          bdf_node(node_num,4)=str2num(linebdf(41: z_end_string_length ));
                    end
                    
                else
                    
                    bdf_node(node_num,4)=str2num(linebdf(41: z_end_string_length )); % node Z cord
                    
   
                end
                
                
                
                
                
                
                node_num=node_num+1;
            end
            %
            
            % read coordinate frames 302 (302 is updated as the spars and ribs locations change)
            if strcmp(linebdf(1:6),'CORD2R') && str2double(linebdf(9:16))==302
                
                bdf_cord(cord_num,1)=str2double(linebdf(9:16));
                
                bdf_cord(cord_num,2)=str2double(linebdf(25:32));
                bdf_cord(cord_num,3)=str2double(linebdf(33:40));
                bdf_cord(cord_num,4)=str2double(linebdf(41:48));
                
                bdf_cord(cord_num,5)=str2double(linebdf(49:56));
                bdf_cord(cord_num,6)=str2double(linebdf(57:64));
                bdf_cord(cord_num,7)=str2double(linebdf(65:end));
                
                % next line
                linebdf=fgetl(fid100);currentFPI=ftell(fid100);
                bdf_cord(cord_num,8)=str2double(linebdf(9:16));
                bdf_cord(cord_num,9)=str2double(linebdf(17:24));
                bdf_cord(cord_num,10)=str2double(linebdf(25:end));
                
                cord_num=cord_num+1;
                
            end
        end
        
    end
    
else
    bdf_elem=NaN;
    bdf_node=NaN;
    disp('2D planform mesh failed!!!!!!!!!!!!!!!!!!')
end
