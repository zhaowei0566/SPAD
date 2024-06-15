%% Curvilinearly Stiffened Laminated Plate -- Buckling Code & Statics Code
% This code is only applicable to the blade-stiffener
% - static
% - buckling
% - vibration
% - PrestressModes
% - August, 7th, 2014
%% ===========Modification History ======================
% - May, 5 2015 Code for SciTech 2016, Vibration code added
% - May, 21 2015: Modify the alpha in the stiffener beam differential
%  stiffness for axial stress in the beam
% - May, 24 2015 : Modify alpha in stiffener geometric stiffness
%% ===============================================================
clear all;warning off;format long;
global Mat;
global FEM;
global Stru;
global Stiffener;
global Laminate;
% Solver='buckling';
%% NASTRUAL TO PHYSICAL

Natural2Physical;


%% beta=a/b  is aspect ratio;
beta1=6/8; % panel length-to-width ratio
ratioTW=0.01; %% plate thickness/plate length;

width_panel=0.8;

numberstiffenerelement=15;
%% =============== PANEL GEOMETRY================================
%%%==============================================================
eee=1;
beta=beta1(eee);
FEM.nodesCord(:,2)=platenodeX';
FEM.nodesCord(:,3)=platenodeY';
%---------Coordinates of node------------
Xcoord=FEM.nodesCord(:,2);
Ycoord=FEM.nodesCord(:,3);
Zcoord=FEM.nodesCord(:,4);
%% ============plot========================
% [xlength,ylength]=meshgrid(0:beta*width_panel/((length(find(Xcoord==0))-1)/2):beta*width_panel,...
%     0:width_panel/((length(find(Xcoord==0))-1)/2):width_panel);
% zlength=xlength.*0+ylength.*0;
% figure
% mesh(xlength,ylength,zlength);hold on;view(2)
% plot(Xcoord,Ycoord,'k.');hold on;
% title('Nodes of Stiffened Plate','FontSize',15);axis image;hold on;
% hold off; %% comment hold off if continue to plot stiffener
%-----------------------------------------
Stru.length=max(abs(Xcoord));
Stru.width=max(abs(Ycoord));
Stru.thickness=0.6/100;%ratioTW(1)*Stru.length; %factort*Stru.length;
%     Stru.thickness=1.04e-3; % Ali's example
%-----------D.O.F for each Node--------------
FEM.PlateNodeDof=5;%% for plate
Stiffener.nodedof=5; %% for stiffener
FEM.nodeCoordinates=FEM.nodesCord(:,2:3);

FEM.GDof=FEM.PlateNodeDof*size(FEM.nodeCoordinates,1);
FEM.typeplate='CQUAD8';
FEM.Dimension='2D';
switch FEM.typeplate
    case 'CQUAD4'
        FEM.GaussPointShear='1by1';
        FEM.GaussPointBend='2by2';
    case 'CQUAD8'
        FEM.GaussPointShear='2by2';
        FEM.GaussPointBend='3by3';
end
%---------------------------------------
FEM.NodeNumber=size(FEM.nodesCord,1);
FEM.elementNumber=size(FEM.elementNodes,1);

%% ----- Stiffener Elements and Nodes, Shape and Functions ----------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stiffenerelementnumber=numberstiffenerelement(eee);
stiffenernodesnumber=2*stiffenerelementnumber+1;
Stiffener.nodenum=stiffenernodesnumber;
Stiffener.elemnum=stiffenerelementnumber;
Stiffener.GDof=Stiffener.nodenum*Stiffener.nodedof;

% ===============================================================
% calculate Natural coordinates and save it in parametric studies
% ===============================================================
for stiffnum=[1 2]
    
    stiffener_Gauss_points;
end
%% ===========-Composite Material Properties ==========================

% MatProperty;
psipa=6894.75729;
% psipa=1;
Mat.kappa=5/6;
Mat.E1=19.2e6*psipa;
Mat.E2=1.56e6*psipa;
Mat.G12=0.82e6*psipa;
Mat.G13=Mat.G12;
Mat.G23=0.49e6*psipa;
Mat.v12=0.24;
Mat.v21=Mat.E2/Mat.E1*Mat.v12;
Mat.density=1800;
Mat.alpha1=-0.04e-6;
Mat.alpha2=16.7e-6;
%%  =================== Fiber ply orientation ==========================
Orientation=[0:1:4 5:5:90];
depthratio=[12]%]%1 5 12 20 30];
deltaT=300;

buckparameter=zeros(length(Orientation),length(depthratio));
TotalWeight=zeros(length(Orientation),length(depthratio));

for AngleNo=1:length(Orientation)
    
    theta=Orientation(AngleNo);
    
    % PlyAngle=[ 0 0 0 0 0 0 0 0];
    
%     PlyAngle=[ 0 90 0 90 90 0 90 0];
%         PlyAngle=[45 -45 45 -45 -45 45 -45 45];
        PlyAngle=[theta -theta theta -theta -theta theta -theta theta];
    
    PlyAnglePlate=PlyAngle;%
    PlyAngleStiffener=[ 0 0 0 0 0 0 0 0];%PlyAngle;%[theta theta theta theta theta theta theta theta];
    %% =====================================================================
    %    laminate thickness - layer thickness or tp/NL
    %%=====================================================================
    %       Laminates for PLATE, cross-ply laminates
    disp('----------Laminates for Plate----------');
    
    Laminate.theta=PlyAnglePlate;
    Laminate.layer=length(Laminate.theta);
    % %================== Laminate thickness ===========================
    Laminate.thickness=Stru.thickness/Laminate.layer*ones(1,Laminate.layer);%% Uniform thickness
    
    %
    Amatrix = zeros(3,3);
    Bmatrix = zeros(3,3);
    Dmatrix = zeros(3,3);
    Ashear =zeros(2,2);
    TempStress=zeros(3,1);
    for layer=1:length(Laminate.thickness)
        [AmatrixK,DmatrixK,AshearK,BmatrixK,QthetaK,ThermalExpCoeff]=LaminatedComposite(Mat,Stru,Laminate,layer);
        Amatrix=Amatrix+AmatrixK;
        Bmatrix=Bmatrix+BmatrixK;
        Dmatrix=Dmatrix+DmatrixK;
        Ashear=Ashear+AshearK;
        
        % % =========== thermal stress =============
        TempStress=TempStress+AmatrixK*ThermalExpCoeff'*deltaT;
    end
    ConstitutiveMatrix=zeros(6,6);
    ConstitutiveMatrix(1:3,1:3)=Amatrix;
    ConstitutiveMatrix(4:6,4:6)=Dmatrix;
    ConstitutiveMatrix(1:3,4:6)=Bmatrix;
    ConstitutiveMatrix(4:6,1:3)=Bmatrix;
    Laminate.A=Amatrix;
    Laminate.B=Bmatrix;
    Laminate.D=Dmatrix;
    Laminate.Ashear=Ashear;
    Laminate.ConstitutiveMatrix=ConstitutiveMatrix;
    
    %% ----LOADS-------
    FEM.p=1e4;% p is the uniform pressure.
    FEM.Mx=0;
    FEM.My=0;
    FORCE=LinearForceMatrix(FEM);
    %     pointload=1;
    %     FORCE=DiscreteForceMatrix(FEM,Stru,pointload);
    
    %% study the ratio of stress_ratio=Nx/Ny
    stress_ratio=1;
    sxx0=1e5*0;
    syy0=1E5*0; % positive for compression; negative for tensile
    txy0=1e5*stress_ratio;
    load=1;
    FEM.stress0=[sxx0,txy0;txy0,syy0]*load;
    %     FEM.stress0=[TempStress(1), TempStress(3);TempStress(3), TempStress(2)];
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% ------PLATE STIFFNESS Matrix --------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    FEM.typeplate='CQUAD8';
    [Kplate,Kelemp]=LinearStiffnessLaminatedPlate(FEM,Laminate);
    
    [KGplate]=GeometryStiffnessPlate(FEM,Stru); %% plate
    
    [MassPlate,MelemPlate]=LinearMassLaminatedPlate(FEM,Mat,Stru);
    
    %% Parameteric study for stiffener depth ratio
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % depth_ratio=Stiffener.height/Stiffener.width;
    
    for depthNO=1:length(depthratio)
        
        FEM.typestiff='CBAR3';
        
        KGStiffener=zeros(FEM.GDof,FEM.GDof); %% Timoshenko beam
        Kstiffener=zeros(FEM.GDof,FEM.GDof);
        MassStiffener=zeros(FEM.GDof,FEM.GDof);
        
        stiffenerWeight=0;
        
        disp(['----case #' num2str(depthNO) '-----'])
        
        %% ------------Stiffener-----------------------
        for stiffno=[ 1 2] % constistent with the number in the code: stiffener_Gauss_points;
            
            disp(['----Stiffener #' num2str(stiffno) '-----'])
            
            xstiffener=XXstiffener(stiffno,:);
            ystiffener=YYstiffener(stiffno,:);
            
            Stiffener.PointsCoord=zeros(stiffenernodesnumber,4);
            Stiffener.PointsCoord(:,1)=[5:4+stiffenernodesnumber]';
            Stiffener.PointsCoord(:,2)=[xstiffener]';
            Stiffener.PointsCoord(:,3)=[ystiffener]';
            
            Stiffener.element=zeros(stiffenerelementnumber,4);
            Stiffener.element(:,1)=[101:100+stiffenerelementnumber]';
            Stiffener.element(:,2)=[5:2:2*stiffenerelementnumber+3]';
            Stiffener.element(:,3)=[7:2:2*stiffenerelementnumber+5]';
            Stiffener.element(:,4)=[6:2:2*stiffenerelementnumber+4]';
            
            % specify gaussian points for each stiffener
            
            Stiffener.TargetPlate=TargetPlateID(stiffno,:);
            eval(['Stiffener.TargetGaussian=TargetGaussianNum' num2str(stiffno) ';'])
            
            %% ==========Laminates for stiffeners==========================
            %%%% ============Laminates for Stiffener
            %%%%========================================================
            display('----------Laminates for Stiffeners-----------');
            
            Stiffener.theta=PlyAngleStiffener;
            Stiffener.layer=length(Stiffener.theta);
            
            Stiffener.width=Stru.thickness;
            %             Stiffener.width=length(Stiffener.theta)*130e-6;
            Stiffener.height=depthratio(depthNO)*Stiffener.width;
            
            % % %  Stiffener Eccentricity % % %
            ebar=(Stru.thickness+Stiffener.height)/2;
            %                        ebar=0;
            
            Stiffener.CrossArea=Stiffener.width*Stiffener.height;
            Stiffener.Istiffener=1/12*Stiffener.width*Stiffener.height^3+ebar^2*Stiffener.CrossArea;
            %
            FEM.Stiffener=Stiffener;
            
            StiffenerElemNodeNum=size(Stiffener.element(1,:),2)-1;
            Stiffener.nodeID=unique(Stiffener.element(:,2:StiffenerElemNodeNum+1));
            Stiffener.nodecord= FEM.nodesCord(Stiffener.nodeID,:);
            
            %% Stiffener stiffness calculations
            
            %%%% Perpendicular or parallel laminate %%%%%%%%%%%%%%%
            Stiffener.orientation='perpendicular'; % 1-parallel, 2-perpendicular
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            flag=1; %% 1- effective GJ; 0 - GJ=0; 2- originial GJ
            % ------------Calculate Constitutive Matrix [Ds] for stiffener----
            [Ds,Qbar,Stiffene]=LaminatedBeamDs(Stiffener,FEM,ebar,Mat,flag);
            
            % ------------Stiffeness for each stiffener --------------
            StiffenerGaussPoint='2';%%
            
            [Kbeam,Kelem,StiffenerLength1]=...
                LaminatedStiffenerStiffness5Dof(FEM,Ds,Stiffener,StiffenerGaussPoint);
            
            %--------- Transform local displacement fields into global ones ----
            
            [Tmatrix,alpha]=TransformationMatrix(FEM,Stiffener);
            Kbeam=Tmatrix'*Kbeam*Tmatrix;
            
            % -------------- stiffener {d} to plate {d} -------------------
%             [ShapeMatrixSP,ShapeMatrixElem]=PlateShapeStiffenerv2(FEM,Stiffener);
                        eval(['ShapeMatrixSP=ShapeMatrixSP' num2str(stiffno) ';']);
            % ----------------- Sum stiffness -------------------------------
            Kstiffener=Kstiffener+ShapeMatrixSP'*Kbeam*ShapeMatrixSP;
            
            % ----------Geometric Stiffness of the stiffener -----------------
            KGbeam=GeometryStiffnessCurvedBeam(FEM,ebar,Stiffener);
            
            KGbeam=Tmatrix'*KGbeam*Tmatrix;
            
            KGStiffener=KGStiffener+ShapeMatrixSP'*KGbeam*ShapeMatrixSP;
            
            % ------------Mass matrix for each stiffener ------------
            [StiffenerMass,StiffenerMelem]=LinearMassLaminatedStiffener(FEM,Mat,ebar,Stiffener,StiffenerGaussPoint);
            Mbeam=Tmatrix'*StiffenerMass*Tmatrix;
            MassStiffener=MassStiffener+ShapeMatrixSP'*Mbeam*ShapeMatrixSP;
            
            % stiffeners total weight
            stiffenerWeight=stiffenerWeight+Stiffener.height*Stiffener.width*StiffenerLength1*Mat.density;
            
        end
        hold off;
        %% Total Material Stiffness and Gometric Stiffness, Mass matrix;
        Ktotal=Kplate+Kstiffener;
        KGtotal=KGplate+KGStiffener;
        Mtotal=MassPlate+MassStiffener;
        
        %%%======================================================
        %%               Buckling/Vibration ANALYSIS
        %%%=======================================================
        disp('---------EigenValue Computation (buckling eigenvalue or vibration eigenvalue)-------')
        %  KtotalBuck=Ktotal(1:FEM.GDof/5*3,1:FEM.GDof/5*3);
        %  KG=KG(1:FEM.GDof/5*3,1:FEM.GDof/5*3);
        
        SOL_flag='106';% 102 - pure plate; 105 - buckling; 103 - free vibration; 106 - presstressed vibration
        
        switch SOL_flag
            
            case '101'
                %% -----------STATIC ANALYSIS ------------- %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % FORCE
                % % distribute forces
                %%Boundary Conditions
                FEM.BCtype='SSSS-3'; %% Four simple supported sides
                [ActiveDof,SideDof]=EssentialBCPlate5Dof(FEM);
                FEM.displacement=zeros(1,FEM.GDof);
                %                 FEM.displacement(ActiveDof)=Ktotal(ActiveDof,ActiveDof)\FORCE(ActiveDof);
                FEM.displacement(ActiveDof)=Kplate(ActiveDof,ActiveDof)\FORCE(ActiveDof);
                %
                % % platetest; % static analysis of pure simply supported plate
                postprocessISOTROPIC;
                
            case '102'
                FEM.BCtype='SSSS-3'; %% Four simple supported sides
                [ActiveDof,SideDof]=EssentialBCPlate5Dof(FEM);
                Solver='vibration';
                [V,D]=eigs(Kplate(ActiveDof,ActiveDof),MassPlate(ActiveDof,ActiveDof),20,'sm');
                [DD,modeNo]=sort(diag(D));
                loadfactor=DD(1:10);
                VVsort=V(:,modeNo);
                cycleFreq=sqrt(loadfactor);
                frequency=sqrt(DD)/2/pi
                vibrparameter(AngleNo,depthNO)=cycleFreq(1)*Stru.length^2/Stru.thickness*sqrt(Mat.density/Mat.E2);
                vibrparameter2(AngleNo,depthNO)=cycleFreq(2)*Stru.length^2/Stru.thickness*sqrt(Mat.density/Mat.E2);
                vibrparameter3(AngleNo,depthNO)=cycleFreq(3)*Stru.length^2/Stru.thickness*sqrt(Mat.density/Mat.E2);
                vibrparameter4(AngleNo,depthNO)=cycleFreq(4)*Stru.length^2/Stru.thickness*sqrt(Mat.density/Mat.E2);
                
            case '105'
                FEM.BCtype='SSSS-3'; %% Four simple supported sides
                [ActiveDof,SideDof]=EssentialBCPlate5Dof(FEM);
                Solver='buckling';
                [V,D]=eigs(Ktotal(ActiveDof,ActiveDof),KGtotal(ActiveDof,ActiveDof),10,'sm');
                % % -----  buckling analysis of pure plate  ------
%                 [V,D]=eigs(Kplate(ActiveDof,ActiveDof),KGplate(ActiveDof,ActiveDof),10,'sm');
                [DD,modeNo]=sort((diag(D)));
                loadfactor=DD(1:10);
                VVsort=V(:,modeNo);
                frequency=loadfactor;
                %-------------minimal load factor---------------
                mineiglabel=find(min(abs(DD))==-DD);
                if isempty(mineiglabel)==1
                    mineiglabel=find(min(abs(DD))==DD);
                end
                plotmodenNo=mineiglabel;
%                 modeshapeplot;
                critical_loadfactor=DD(mineiglabel);
                
                buckparameter(AngleNo,depthNO)=critical_loadfactor*1e5*Stru.length^2/Mat.E2/Stru.thickness^2
                buckLF=critical_loadfactor
                buckLoad(AngleNo,depthNO)=critical_loadfactor*1e5*Stru.width*Stru.thickness;
                buckLoadNXX(AngleNo,depthNO)=critical_loadfactor*1e5*Stru.thickness/1e3 % 1e3 for mm
                
            case '103'
                %% -----------Boundary Conditions-----------
                FEM.BCtype='SSSS-3'; %% Four simple supported sides
                [ActiveDof,SideDof]=EssentialBCPlate5Dof(FEM);
                Solver='vibration';
                [V,D]=eigs(Ktotal(ActiveDof,ActiveDof),Mtotal(ActiveDof,ActiveDof),10,'sm');
                [DD,modeNo]=sort(diag(D));
                eigenvalue=DD(1:10);
                VVsort=V(:,modeNo);
                cycleFreq=sqrt(eigenvalue);
                frequency=cycleFreq/2/pi
                
                vibrparameter(AngleNo,depthNO)=cycleFreq(1)*Stru.length^2/Stru.thickness*sqrt(Mat.density/Mat.E2);
                vibrparameter2(AngleNo,depthNO)=cycleFreq(2)*Stru.length^2/Stru.thickness*sqrt(Mat.density/Mat.E2);
                vibrparameter3(AngleNo,depthNO)=cycleFreq(3)*Stru.length^2/Stru.thickness*sqrt(Mat.density/Mat.E2);
                vibrparameter4(AngleNo,depthNO)=cycleFreq(4)*Stru.length^2/Stru.thickness*sqrt(Mat.density/Mat.E2);
                %cycleFreq(1:4)*Stru.length^2/Stru.thickness*sqrt(Mat.density/Mat.E2)
                
                mineiglabel=find(min(abs(DD))==-DD);
                if isempty(mineiglabel)==1
                    mineiglabel=find(min(abs(DD))==DD);
                end
                
                plotmodenNo=mineiglabel;
                
                modeshapeplot;
                
            case '106'
                FEM.BCtype='SSSS-3'; %% Four simple supported sides
                [ActiveDof,SideDof]=EssentialBCPlate5Dof(FEM);
                Solver='buckling';
                [V,D]=eigs(Ktotal(ActiveDof,ActiveDof),KGtotal(ActiveDof,ActiveDof),10,'sm');
                % % -----  buckling analysis of pure plate  ------
                %                  [V,D]=eigs(Kplate(ActiveDof,ActiveDof),KGP(ActiveDof,ActiveDof),10,'sm');
                [DD,modeNo]=sort(diag(D));
                loadfactor=DD(1:10);
                VVsort=V(:,modeNo);
                frequency=DD;
                %-------------minimal load factor
                mineiglabel=find(min(abs(DD))==-DD);
                if isempty(mineiglabel)==1
                    mineiglabel=find(min(abs(DD))==DD);
                end
                plotmodenNo=mineiglabel;
                %                 modeshapeplot;
                critical_loadfactor=DD(mineiglabel);
                
                buckparameter(AngleNo,depthNO)=critical_loadfactor*1e5*Stru.length^2/Mat.E2/Stru.thickness^2
                buckLF=critical_loadfactor
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Free vibration
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [V,D]=eigs(Ktotal(ActiveDof,ActiveDof),Mtotal(ActiveDof,ActiveDof),10,'sm');
                [DD,modeNo]=sort(diag(D));
                eigenvalue=DD(1:10);
                VVsort=V(:,modeNo);
                cycleFreq=sqrt(eigenvalue);
                frequency=cycleFreq/2/pi
                
                vibrparameterFree(AngleNo,depthNO)=cycleFreq(1)*Stru.length^2/Stru.thickness*sqrt(Mat.density/Mat.E2);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                Solver='prestressed_vibr';
                
                jj=1
                for lambda_b_ratio=[-0.75]
                    %
                    Keff=Ktotal+lambda_b_ratio*buckLF*KGtotal; %% elastic stiffness + geometric stiffness
                    
                    [V,D]=eigs(Keff(ActiveDof,ActiveDof),Mtotal(ActiveDof,ActiveDof),10,'sm');
                    
                    [DD,modeNo]=sort(diag(D));
                    
                    loadfactor=DD(1:10);
                    cycleFreq=sqrt(loadfactor);
                    VVsort=V(:,modeNo);
                    frequency=cycleFreq/2/pi;
%                     vibrparameter(jj,1:4)=cycleFreq(1:4)*Stru.length^2/Stru.thickness*sqrt(Mat.density/Mat.E2);
                    jj=jj+1;
                    vibrparameter(AngleNo,depthNO)=cycleFreq(1)*Stru.length^2/Stru.thickness*sqrt(Mat.density/Mat.E2);
                    
                    if lambda_b_ratio==-0.75  % show mode shape at load factor ratio 0.9
                        mineiglabel=find(min(abs(DD))==-DD);
                        
                        if isempty(mineiglabel)==1
                            mineiglabel=find(min(abs(DD))==DD);
                        end
                        
                        %  plotmodenNo=mineiglabel;
                        plotmodenNo=1;
                        modeshapeplot;
                    end
                end
        end
    end
    %                 plotmodenNo=1;
    %         modeshapeplot;
    % Total weight of stiffened composite panel
    TotalWeight(AngleNo,depthNO)=stiffenerWeight+Stru.length*Stru.width*Stru.thickness*Mat.density;
    
end

%% plot buckparameter via stiffener depth ratio and ply orientation
%
% figure;
% [TTheta,DepthRatio]=meshgrid(depthratio,Orientation);
% surf(DepthRatio,TTheta,vibrparameter);
% xlabel(['Ply orientation (degree)'],'FontSize',14);
% ylabel(['Stiffener depth ratio (h_s/b_s)'],'FontSize',14);
% zlabel(['Buckling parameter'],'FontSize',14);

%% backup code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -----------STATIC ANALYSIS ------------- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORCE
% % distribute forces
% %%Boundary Conditions
% FEM.BCtype='SSSS-3'; %% Four simple supported sides
% [ActiveDof,SideDof]=EssentialBCPlate5Dof(FEM);
% FEM.displacement=zeros(1,FEM.GDof);
% FEM.displacement(ActiveDof)=Ktotal(ActiveDof,ActiveDof)\FORCE(ActiveDof);
% %
% % % platetest; % static analysis of pure simply supported plate
% postprocessISOTROPIC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%------Vibration analysis of PURE COMPOSITE PANEL-----
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     FEM.BCtype='SSSS-3'; %% Four simple supported sides
%     [ActiveDof,SideDof]=EssentialBCPlate5Dof(FEM);
%     disp('-----------Pure plate---Vibration Eigenvalue Computation-----------');
%     [V,D]=eig(Kplate(ActiveDof,ActiveDof),MassPlate(ActiveDof,ActiveDof));
%     [DD,modeNo]=sort(diag(D));
%     cycle_freq=DD(1:10);
%     VVsort=V(:,modeNo);
%     frequency=DD;
%     plotmodenNo=2;
%     Solver='vibration';
%     modeshapeplot;
%
%     disp('Natural Frequency of pure composite panel: ')
%     frequency=sqrt(cycle_freq)/2/pi
%% Plot stiffener height
% figure(8);hold on;
% plot(depthratio,abs(buckparameter));hold off;

%% Plot prestressed vibr mode
%
% vibrparameter=real(vibrparameter);
% Appliedloadfactor=[0:0.1:1];
% figure(2);hold on;
% plot(Appliedloadfactor,vibrparameter(:,1)/vibrparameter(1,1));
% plot(Appliedloadfactor,vibrparameter(:,2)/vibrparameter(1,2));hold off;
% % plot(Appliedloadfactor,vibrparameter(3,:));hold on;
% % plot(Appliedloadfactor,vibrparameter(4,:));hold off;
% xlabel('\lambda_b/\lambda_{b,cr}');
% ylabel('\lambda_n/\lambda_{n,i}');
%
%    vib_Nyy1=[1.000000000000000
%    0.959629884015866
%    0.916071876422738
%    0.868593113828158
%    0.816170470486641
%    0.757311204311431
%    0.689703948557392
%    0.609449122493775
%    0.508998009690508
%    0.369262431839685
%    0.000002148288225];
%
% vib_Nyy2=[1.000000000000000
%    0.962077863459251
%    0.921777354823483
%    0.878641469435184
%    0.832062631378748
%    0.781207632487016
%    0.724887559698307
%    0.661313289669127
%    0.587580527928434
%    0.498383236499411
%    0.381710535695889];
% hold on;
% plot(Appliedloadfactor,vib_Nyy1);hold on;
% plot(Appliedloadfactor,vib_Nyy2);hold off;

