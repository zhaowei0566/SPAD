
% postprocess for solution from stiffenedPlate

deformUZ=FEM.displacement(1:FEM.GDof/FEM.PlateNodeDof);
deformBx=FEM.displacement(1+FEM.GDof/FEM.PlateNodeDof:2*FEM.GDof/FEM.PlateNodeDof);
deformBy=FEM.displacement(1+2*FEM.GDof/FEM.PlateNodeDof:3*FEM.GDof/FEM.PlateNodeDof);


deformUX=FEM.displacement(1+3*FEM.GDof/FEM.PlateNodeDof:4*FEM.GDof/FEM.PlateNodeDof);
deformUY=FEM.displacement(1+4*FEM.GDof/FEM.PlateNodeDof:5*FEM.GDof/FEM.PlateNodeDof);




% s1 = find(FEM.nodeCoordinates_label(:,3) == Stru.width/2)
% 
% temp2 = FEM.nodeCoordinates_label(s1,:)
% 
% s2 = find(temp2(:,2) == Stru.length/2);
% 
% temp3 = temp2(s2,:)
% 
% deformUX(temp3(1))
% deformUY(temp3(1))





% Spline deflection and angle
dx=min(Xcoord):(max(Xcoord)-min(Xcoord))/100:max(Xcoord);
dy=min(Ycoord):(max(Ycoord)-min(Ycoord))/100:max(Ycoord);
[x3,y3]=meshgrid(dx,dy,0);

refined_UX =  griddata(Xcoord,Ycoord,deformUX,x3,y3,'v4');
refined_UY =  griddata(Xcoord,Ycoord,deformUY,x3,y3,'v4');
refined_UZ =  griddata(Xcoord,Ycoord,deformUZ,x3,y3,'v4');
refined_BetaX=griddata(Xcoord,Ycoord,deformBx,x3,y3,'v4');
refined_BetaY=griddata(Xcoord,Ycoord,deformBy,x3,y3,'v4');

figure
surf(x3,y3,refined_UX,'FaceColor','interp',...
    'EdgeColor','none',...
    'FaceLighting','phong');colormap(jet);view(2);axis image;
title('In-plane displacement along x-axis, u');colorbar
set(gcf,'color','w');
set(gca,'FontSize',16)


figure
surf(x3,y3,refined_UY,'FaceColor','interp',...
    'EdgeColor','none',...
    'FaceLighting','phong');colormap(jet);view(2);axis image;
title('In-plane displacement along y-axis, v');colorbar
set(gcf,'color','w');
set(gca,'FontSize',16)

figure
surf(x3,y3,refined_UZ,'FaceColor','interp',...
    'EdgeColor','none',...
    'FaceLighting','phong');colormap(jet);view(2);axis image;
title('Out-of-plane displacement, w');colorbar
set(gcf,'color','w');
set(gca,'FontSize',16)


figure
surf(x3,y3,refined_BetaX,'FaceColor','interp',...
    'EdgeColor','none',...
    'FaceLighting','phong');colormap(jet);view(2);axis image;
title('Rotational angle about y-axis, \beta_x');colorbar
set(gcf,'color','w');
set(gca,'FontSize',16)

figure
surf(x3,y3,refined_BetaY,'FaceColor','interp',...
    'EdgeColor','none',...
    'FaceLighting','phong');colormap(jet);view(2);axis image;
title('Rotational angle about x-axis, \beta_y');colorbar
set(gcf,'color','w');
set(gca,'FontSize',16)
% griddata(Xcoord,Ycoord,deformUX,Stru.length/2,Stru.width/2,'v4')

% griddata(Xcoord,Ycoord,deformUY,Stru.length/2,Stru.width/2,'v4')

% %% Plot in-plane deflections
% figure(111);hold on;
% % plot(Xcoord,Ycoord,'k.');hold on;
% title('Nodes of Stiffened Plate','FontSize',15);axis image;hold on;
% SF=0.1*Stru.length/max(abs([deformUX deformUY deformUZ]));
% 
% 
% FEM.nodeCoordinates_label_deformed = FEM.nodeCoordinates_label;
% FEM.nodeCoordinates_label_deformed(:,2) = Xcoord'+SF*deformUX;
% FEM.nodeCoordinates_label_deformed(:,3) = Ycoord'+SF*deformUY;
% 
% % deformed body
% patch_plot(FEM.elementNodes,FEM.nodeCoordinates_label_deformed,200,'mesh')
%%  Transverse Deflection Plot
% % % figure
% % % surf(x3,y3,z3,'FaceColor','interp',...
% % %     'EdgeColor','none',...
% % %     'FaceLighting','phong');
% % % if ebar==0
% % %     title(['Transverse Deflection of stiffened plate under uniform load, Concentric'],'FontSize',14);
% % % else
% % %     title(['Transverse Deflection of stiffened plate under uniform load, Eccentric'],'FontSize',14);
% % % end
% % % 
% % % xlabel('Length of the plate, a /m');
% % % ylabel('Width of the plate, b/m');
% % % zlabel('Transverse deflection /m');
% % % colorbar('location','SouthOutside');
% % % colormap('jet')
% % % view(2);
% % % axis image;

disp('The Largest Displacement UX is:');
[value,id] = max(abs(deformUX));
deformUX(id)

disp('The Largest Displacement UY is:');
[value,id] = max(abs(deformUY));
deformUY(id)



disp('The Largest Transverse Displacement is:');
[value,id] = max(abs(deformUZ));
deformUZ(id)



disp('The Largest Transverse Displacement is:');
[value,id] = max(abs(deformUZ));
deformUZ(id)


disp('The Largest Rotations BetaX is:');
[value,id] = max(abs(deformBx));
deformBx(id)

disp('The Largest Rotations BetaY is:');
[value,id] = max(abs(deformBy));
deformBy(id)


