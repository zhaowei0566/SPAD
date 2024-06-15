function [shape,naturalderivatives,d2Nds2]=shapefunctionbeam(xi,FEM)
type=FEM.typestiff;
switch type
    
    case 'CBAR3'
        shape=[xi*(xi-1)/2;
            xi*(xi+1)/2;
            (1-xi^2)];
        naturalderivatives=[xi-1/2,xi+1/2,-2*xi];
        
        d2Nds2=[1,1,-2];
        
       
    case 'CBAR2'
        %C0
        shape=[1/2*(1-xi); 1/2*(1+xi)];
        naturalderivatives=[-1/2,1/2];
        
        d2Nds2=[0,0];

end
