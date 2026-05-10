function draw_groundwave_geometry_2d()
clc; close all;


Rx = 1.8;
Ry = 1.0;
A  = [-Rx-0.6, 0];
B  = [ Rx+0.6, 0.45];
P1 = [-0.6, -0.15];
P2 = [ 0.5, -0.05];
alpha_deg = 20;


set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');


t = linspace(0, 2*pi, 400);
x = Rx*cos(t); y = Ry*sin(t);

figure('Color','w'); hold on; axis equal;
fill(x, y, [0.75 0.75 0.75], ...
    'EdgeColor',[0.2 0.2 0.2], ...
    'LineStyle','--','LineWidth',1.5);

text(0, -0.05*Ry, 'S_p','HorizontalAlignment','center','FontSize',12,'FontWeight','bold');

Lx = 2.5; Ly = 2.0;
plot([-Lx Lx],[0 0],'k-','LineWidth',1.5);
plot([0 0],[-Ly Ly],'k-','LineWidth',1.5);

text(Lx*0.98, 0.03,'X','FontSize',12,'FontWeight','bold');
text(0.03, Ly*0.98,'Y','FontSize',12,'FontWeight','bold');


plot(A(1),A(2),'ro','MarkerFaceColor','r','MarkerSize',6);
text(A(1)-0.1,A(2)+0.08,'A (Tx)','Color','r','FontSize',12,'FontWeight','bold');

plot(B(1),B(2),'bo','MarkerFaceColor','b','MarkerSize',6);
text(B(1)+0.05,B(2)+0.05,'B (Rx)','Color','b','FontSize',12,'FontWeight','bold');

plot(P1(1),P1(2),'k.','MarkerSize',18);
text(P1(1)+0.05,P1(2)-0.08,'P_1','FontSize',11,'FontWeight','bold');

plot(P2(1),P2(2),'k.','MarkerSize',18);
text(P2(1)+0.05,P2(2)-0.08,'P_2','FontSize',11,'FontWeight','bold');


plot([A(1) B(1)], [A(2) B(2)], 'k--','LineWidth',1.5);


r = 0.35;
phi = linspace(0, deg2rad(alpha_deg), 60);

xarc = B(1) + r*cos(phi); 
yarc = B(2) + r*sin(phi);

plot(xarc, yarc, 'k','LineWidth',1.5);

text(B(1)+r*0.8*cosd(alpha_deg/2), ...
     B(2)+r*0.8*sind(alpha_deg/2)+0.05, ...
     '\alpha','FontSize',12,'FontWeight','bold');

plot([B(1) B(1)+r], [B(2) B(2)], 'k','LineWidth',1.5);
plot([B(1) B(1)+r*cosd(alpha_deg)], ...
     [B(2) B(2)+r*sind(alpha_deg)], 'k','LineWidth',1.5);

xlim([-Lx Lx]); ylim([-Ly Ly]);

xlabel('x','FontSize',12,'FontWeight','bold');
ylabel('y','FontSize',12,'FontWeight','bold');

title('2D schematic: vertical current at A over isolated surface S_p', ...
    'FontSize',13,'FontWeight','bold');

grid on;


ax = gca;
ax.LineWidth = 1.8;   
ax.Box = 'on';        


ax.TickDir = 'in';        
ax.TickLength = [0.015 0.015];
ax.FontSize = 11;
ax.FontWeight = 'bold';

end