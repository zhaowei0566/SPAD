function [Ds,Qbar,Stiffene]=LaminatedBeamDs(Stiffener,FEM,ebar,Mat,GJflag)
% Calculate Constitutive Matrix for laminated beam
switch Stiffener.orientation
    
    case 'parallel' % parallel
        
        disp('---------------Stiffener PARALLEL -------------');
        
        Stiffener.thickness=Stiffener.height/Stiffener.layer*ones(1,Stiffener.layer);%% lamina along thickness
        Amatrix=zeros(2,2);
        Bmatrix=zeros(2,2);
        Dmatrix=zeros(2,2);
        Ashear=zeros(1,1);
        for layer=1:length(Stiffener.thickness)
            [AmatrixK,BmatrixK,DmatrixK,AshearK,Qreduced,Qbar]=...
                StiffenerCompositeParallel(Mat,Stiffener,layer,ebar);
            
            Amatrix=Amatrix+AmatrixK;
            Bmatrix=Bmatrix+BmatrixK;
            Dmatrix=Dmatrix+DmatrixK;
            Ashear=Ashear+AshearK;
        end
        bs=Stiffener.width;
        hs=Stiffener.height;
        A11=Amatrix(1,1);A16=Amatrix(1,2);A66=Amatrix(2,2);
        B11=Bmatrix(1,1);B16=Bmatrix(1,2);B66=Bmatrix(2,2);
        D11=Dmatrix(1,1);D16=Dmatrix(1,2);D66=Dmatrix(2,2);
        A55=Ashear;
        Ds=[bs*A11,     bs*A16,      0,      bs*B11,         bs*B16;
            bs*A16,     bs*A66,      0,      bs*B16,         bs*B66;
            0,           0,        bs*A55,      0,              0;
            bs*B11,     bs*B16,      0,      bs*D11,        bs*D16;
            bs*B16,     bs*B66,      0,      bs*D16,    bs*D66-1/12*bs^3*A55];
        % reduced GJ
        Ds(5,5)=1/6*(Mat.G12+Mat.G13)*Stiffener.width^3*hs;
        
        % ================================================================
        % ---------- Material Stiffness of the stiffener -----------------
        
    case 'perpendicular' %perpendicular'
        
        disp('------------ Stiffener PERPENDICULAR -------------');
        
        Stiffener.thickness=Stiffener.width/Stiffener.layer*ones(1,Stiffener.layer);%% lamina along width
        Amatrix=zeros(2,2);
        Bmatrix=zeros(2,2);
        Dmatrix=zeros(2,2);
        Ashear=zeros(1,1);
        for layer=1:length(Stiffener.thickness)
            
            [AmatrixK,BmatrixK,DmatrixK,AshearK,Qreduced,Qbar,Q_bend_new,Q_shear_new]=...
                StiffenerCompositePerpendicular(Mat,Stiffener,layer);
            
            Amatrix=Amatrix+AmatrixK;
            Bmatrix=Bmatrix+BmatrixK;
            Dmatrix=Dmatrix+DmatrixK;
            Ashear =Ashear+AshearK;
            
        end
        
        hs=Stiffener.height;
        A11=Amatrix(1,1);A16=Amatrix(1,2);A66=Amatrix(2,2);
        B11=Bmatrix(1,1);B16=Bmatrix(1,2);B66=Bmatrix(2,2);
        D11=Dmatrix(1,1);D16=Dmatrix(1,2);D66=Dmatrix(2,2);
        
        % A55 includes shear correction factors
        A55=Ashear;
        
        Ds=[hs*A11,     0,          hs*A16,         hs*ebar*A11,                hs*B16;
            0,         hs*A55,         0,                  0,                    hs*ebar*A55;
            hs*A16,      0,         hs*A66,       hs*ebar*A16,                    hs*B66;
            hs*ebar*A11, 0,         hs*ebar*A16,   (hs^3/12+hs*ebar^2)*A11,       hs*ebar*B16;
            hs*B16,     hs*ebar*A55,    hs*B66,     hs*ebar*B16,       ((hs^3/12+hs*ebar^2)*A55+hs*D66)];
        
end


if GJflag==1
    Ds(5,5)=1/6*(Mat.G12+Mat.G13)*Stiffener.width^3*Stiffener.height;
elseif GJflag==0
    Ds(5,5)=0;
elseif GJflag==3
    Ds(5,5)=Mat.E1*Ds(5,5);
end

Stiffene.Amatrix=Amatrix;
Stiffene.Bmatrix=Bmatrix;
Stiffene.Dmatrix=Dmatrix;
Stiffene.Ashear=Ashear;
Stiffene.Qreduced=Qreduced;
Stiffene.QbarBend=Q_bend_new;
Stiffene.QbarShear=Q_shear_new;

