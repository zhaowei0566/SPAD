function [Ds,C] = Orthotropic_beam_Ds(Mat,Stiffener,ebar)

kappa = Mat.kappa; % shear correction factor for paroblic stress distribuion

bs = Stiffener.width;
hs = Stiffener.height;
As = Stiffener.CrossArea;
In = Stiffener.Istiffener;

Ds = zeros(5,5);


[C11bar,C55bar,C66bar,C] = Coeff_stress_strain_4_beam(Mat);

% % % % C11bar = C(1,1);
% C(1,1)



Ds(1,1) = C11bar*As; Ds(1,4) = ebar*As*C11bar; 

Ds(2,2) = kappa*As*C66bar; Ds(2,5) = ebar*kappa*As*C66bar;

Ds(3,3) = kappa*As*C55bar;

Ds(4,1) = Ds(1,4); Ds(4,4) = In*C11bar;

Ds(5,2) = Ds(2,5);



%% Effective GJ 

Gs = sqrt(Mat.G12*Mat.G13);

if hs/bs>=10
    
    temp_sum = 0;
    
    for p = 1:1000
        
        
        temp1 = (1-(-1)^p)/p^5 *tanh(p*pi*hs/2/bs*sqrt(Mat.G13/Mat.G12));
        
        temp_sum = temp_sum + temp1;
        
    end
    
    Js = sqrt(Mat.G13/Mat.G12)* bs^3*hs/3 *(1-96/pi^5*bs/hs*sqrt(Mat.G13/Mat.G12)*temp_sum);

    GJ = Gs*Js;
    
else
    
    GJ = Gs*bs^3*hs/3;
    
end

%% 
Ds(5,5) = GJ;

