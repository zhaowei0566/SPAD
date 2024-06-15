function [Mass,Melem]=LinearMassLaminatedStiffener(FEM,Mat,ebar,Stiffener,Gausspointbend)
% ------------------------------------
% subroutine for stiffener mass matrix
% ------------------------------------

disp('-------------Assembling Stiffener Mass matrix-------------')
gDOF=Stiffener.nodedof*Stiffener.nodenum;

rho=Mat.density;

Mass=zeros(gDOF,gDOF); % initialize mass matrix for stiffener mass matrix

numberElements=size(FEM.Stiffener.element,1); % Stiffener elements number

As=Stiffener.width*Stiffener.height; % Cross sectional area of stiffeners

In=1/12*Stiffener.width*Stiffener.height^3+ebar^2*As;

Ib=1/12*Stiffener.width^3*Stiffener.height;

num_elem_nodes=3;

% Cycle for each element, element stiffness
for e=1:numberElements
%     disp(['-------------stiffener No#' num2str(e) '-------------']);
    
    stiffenerNODE=FEM.Stiffener.element(e,2:end); %% Node NO. for one element
    XY= Stiffener.PointsCoord(stiffenerNODE,2:3);
    
    elementDof2=stiffenerNODE;
            
    elementDof=[elementDof2 elementDof2+Stiffener.nodenum  elementDof2+Stiffener.nodenum*2 ...
                elementDof2+Stiffener.nodenum*3 elementDof2+Stiffener.nodenum*4];
        
    Melem=zeros(size(elementDof2,2)*5,size(elementDof2,2)*5);
    %
    [GaussWeights,GaussLocations]=gaussQuadrature1D(Gausspointbend);
    
    % Loop for Gauss Points
    for ee=1:size(GaussWeights,1)
        
        GaussPoint=GaussLocations(ee);
        xi=GaussPoint;
        [shape,naturalderivatives,d2Nds2]=shapefunctionbeam(xi,FEM);

        detJ=sqrt((naturalderivatives*XY(:,1))^2+(naturalderivatives*XY(:,2))^2);
        
        %% --------Membrane strain energy-------------
        
        NmassShape=zeros(5,num_elem_nodes*5);
        
        NmassShape(1,1:num_elem_nodes)=shape';
        NmassShape(2,1+num_elem_nodes:num_elem_nodes*2)=shape';
        NmassShape(3,1+num_elem_nodes*2:num_elem_nodes*3)=shape';
        NmassShape(4,1+num_elem_nodes*3:num_elem_nodes*4)=shape';
        NmassShape(5,1+num_elem_nodes*4:num_elem_nodes*5)=shape';
        
        mass=rho*[As,      0,     0,      As*ebar,    0;
                0,       As,    0,         0,    As*ebar;
                0,       0,     As,        0,       0;
             As*ebar,    0,     0,         In,      0;
                0,    As*ebar,  0,         0,     In+Ib];
            
       TransfMatrix=[0 0 0 1 0;
           0 0 0 0 1;
           1 0 0 0 0;
           0 1 0 0 0;
           0 0 1 0 0];
        
       mass=TransfMatrix'*mass*TransfMatrix;
       
        Mass(elementDof,elementDof)=Mass(elementDof,elementDof)+...
            NmassShape'*mass*NmassShape*GaussWeights(ee)*detJ;
        
        Melem=Melem+rho*NmassShape'*mass*NmassShape*GaussWeights(ee)*detJ;
        
    end
end


