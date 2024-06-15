function Jstiffener=StiffenerTorsionConstant1D(Stiffener)

a=Stiffener.height;
b=Stiffener.width;

if a/b>=1 && a/b<1.5
    
    Jstiffener=0.141*a*b^3;
elseif a/b>=1.5 && a/b<2
    
    Jstiffener=0.196*a*b^3;
    
elseif a/b>=2 && a/b<2.5
    
    Jstiffener=0.229*a*b^3;
elseif a/b>=2.5 && a/b<3
    
    Jstiffener=0.249*a*b^3;
elseif a/b>=3 && a/b<4
    
    Jstiffener=0.263*a*b^3;
    
elseif a/b>=4 && a/b<5
    
    Jstiffener=0.281*a*b^3;
    
elseif a/b>=5 && a/b<6
    
    Jstiffener=0.291*a*b^3;
    
elseif a/b>=6 && a/b<10
    
    Jstiffener=0.299*a*b^3;
else
    Jstiffener=1/3*a*b^3;
end