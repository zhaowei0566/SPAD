%% Geometry Coordinate and thickness of plate
% Read coordinates from MSC Patran directly;
% CQUAD8
function FEM=geometryinfoplate(Stru)
CBAR2=[];CQUAD4=[];CQUAD8=[];
%% This is used to read CQUAD4 element
inputfile=Stru.bdfname;
fid100=fopen([Stru.pathfile filesep inputfile]);
string1='CQUAD4';string2='GRID';string3='CBAR';
status=fseek(fid100,0,'eof');
EOF=ftell(fid100);
currentFPI=fseek(fid100,0,'bof');
elem=1;grid=1;cbarNo=1;
while currentFPI<EOF
    linebdf=fgetl(fid100);currentFPI=ftell(fid100);
    str1=findstr(linebdf,'CQUAD4');%%element node connection
    str2=findstr(linebdf,'GRID'); %% node coordinates
    str3=findstr(linebdf,'CBAR');
    
    if length(linebdf)>=50 && min(linebdf(1:6)=='CQUAD4')
        CQUAD4(elem,:)=str2num(linebdf(7:end));
        elem=elem+1;
    end
    if length(linebdf)>=58 && min(linebdf(1:4)=='CBAR')
        CBAR2(cbarNo,:)= str2num(linebdf(7:end));
        cbarNo=cbarNo+1;
    end
    if length(linebdf)>=43 && min(linebdf(1:4)=='GRID')
        NodeCoord(grid,:)=[str2num(linebdf(9:16)) str2num(linebdf(25:32)) str2num(linebdf(33:40)) ...
            str2num( linebdf(41:end))];
        grid=grid+1;
    end
    if length(linebdf)>=67 && min(linebdf(1:6)=='CQUAD8')
        CQUAD8(elem,1:8)=str2num(linebdf(7:end));
        linebdf2=fgetl(fid100);currentFPI=ftell(fid100);
        CQUAD8(elem,9:10)=str2num(linebdf2);
        elem=elem+1;
    end
end

%% Information of FEMFEM.nodeCoordinates=NodeCoord(:,2:3);
if isempty(CQUAD4)==0
    FEM.elementNodes=CQUAD4(:,3:6);
end
if isempty(CQUAD8)==0
    Number=[1,3:10];
    FEM.elementNodes=CQUAD8(:,3:10);
end
FEM.numberElements=size(FEM.elementNodes,1);
FEM.cbar=CBAR2;
FEM.nodesCord=NodeCoord;
fclose(fid100);
% FEM.GDof;
% FEM.numberElements
% FEM.elementNodes
% FEM.numberNodes
% FEM.nodeCoordinates