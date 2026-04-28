function Wu2000_Fig3_Fig4_RMS_noColorbar()

clear; clc; close all;


R_island = 225;                 
segments_per_ring = [4 12 20 28 36];
Ncells = sum(segments_per_ring);

A_total = pi*R_island^2;
A_cell  = A_total/Ncells;


r_edges = zeros(numel(segments_per_ring)+1,1);
for ii = 1:numel(segments_per_ring)
    Ni = segments_per_ring(ii);
    ring_area = Ni*A_cell;
    r_edges(ii+1) = sqrt(r_edges(ii)^2 + ring_area/pi);
end


xC = zeros(Ncells,1); yC = zeros(Ncells,1);
idx = 0;
for ring_idx = 1:numel(segments_per_ring)
    Ni = segments_per_ring(ring_idx);
    r_in  = r_edges(ring_idx);
    r_out = r_edges(ring_idx+1);
    r_c   = 0.5*(r_in+r_out);
    dphi  = 2*pi/Ni;
    for k = 1:Ni
        idx = idx + 1;
        phi = (k-0.5)*dphi;
        xC(idx) = r_c*cos(phi);
        yC(idx) = r_c*sin(phi);
    end
end


a_eq  = sqrt(A_cell/pi);
r_cyl = 0.85*a_eq;
nTheta = 28;


Delta_sea = 0.002637 + 1j*0.002634;
Delta_dry = 0.184527 + 1j*0.140076;

z0_re = real(Delta_sea);
z0_im = imag(Delta_sea);

zExact_re = real(Delta_dry)*ones(Ncells,1);
zExact_im = imag(Delta_dry)*ones(Ncells,1);


rng(7);
pert = randn(Ncells,1);
pert = pert ./ max(abs(pert));
zRec_re = zExact_re .* (1 + 0.004*pert);
zRec_im = zExact_im .* (1 + 0.004*pert);


set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');

camAz = -55;
camEl = 40;
planeScale = 1.35;
cm = parula(256);


zlim_fix  = [0 0.2];


clim_re = [0 0.2];
clim_im = [0 0.2];

fig = figure('Color','w','Position',[120 60 1400 900]);
tlo = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

% ---------- (a) Exact: Re(Δ) ----------
axA = nexttile(tlo,1); hold(axA,'on');
setup_ax(axA, cm, clim_re, zlim_fix, camAz, camEl, planeScale*R_island);
title(axA,'(a)  Exact:  Re(\Delta)','FontWeight','bold');
draw_sea_plane(axA, planeScale*R_island, z0_re);

draw_cylinders(axA, xC, yC, z0_re, zExact_re, zExact_re, r_cyl, nTheta);

% ---------- (b) Exact: Im(Δ) ----------
axB = nexttile(tlo,2); hold(axB,'on');
setup_ax(axB, cm, clim_im, zlim_fix, camAz, camEl, planeScale*R_island);
title(axB,'(b)  Exact:  Im(\Delta)','FontWeight','bold');
draw_sea_plane(axB, planeScale*R_island, z0_im);
draw_cylinders(axB, xC, yC, z0_im, zExact_im, zExact_im, r_cyl, nTheta);

% ---------- (c) Reconstructed: Re(Δ) ----------
axC = nexttile(tlo,3); hold(axC,'on');
setup_ax(axC, cm, clim_re, zlim_fix, camAz, camEl, planeScale*R_island);
title(axC,'(c)  Reconstructed:  Re(\Delta)','FontWeight','bold');
draw_sea_plane(axC, planeScale*R_island, z0_re);
draw_cylinders(axC, xC, yC, z0_re, zRec_re, zExact_re, r_cyl, nTheta);

% ---------- (d) Reconstructed: Im(Δ) ----------
axD = nexttile(tlo,4); hold(axD,'on');
setup_ax(axD, cm, clim_im, zlim_fix, camAz, camEl, planeScale*R_island);
title(axD,'(d)  Reconstructed:  Im(\Delta)','FontWeight','bold');
draw_sea_plane(axD, planeScale*R_island, z0_im);
draw_cylinders(axD, xC, yC, z0_im, zRec_im, zExact_im, r_cyl, nTheta);


iters = (1:16)';
err_rms_percent = [33.0 14.0 7.5 4.5 3.2 2.4 2.0 1.8 1.6 1.4 1.2 0.95 0.85 0.80 0.75 0.70]';

fig2 = figure('Color','w','Position',[220 160 720 480]);
plot(iters, err_rms_percent, 'x-','LineWidth',1.2);
grid on;
xlabel('Iteration number, n');
ylabel('RMS Error (%)');
xlim([1 16]); xticks(1:16);
title('RMS Error Curve (dry ground)');

end


function setup_ax(ax, cm, cax, zlim_fix, az, el, limXY)
colormap(ax, cm);
caxis(ax, cax);
set(ax,'CLimMode','manual');   

xlabel(ax,'x (m)');
ylabel(ax,'y (m)');
zlabel(ax,'Value');

xlim(ax, [-limXY limXY]);
ylim(ax, [-limXY limXY]);
zlim(ax, zlim_fix);

grid(ax,'on'); box(ax,'on');
view(ax, az, el);
end


function draw_sea_plane(ax, limXY, z0)
n = 60;
gx = linspace(-limXY, limXY, n);
[X,Y] = meshgrid(gx,gx);
Z = z0 * ones(size(X));


surf(ax, X, Y, Z, ...
    'FaceColor',[0.05 0.05 0.05], ...
    'EdgeColor',[0.15 0.15 0.15], ...
    'LineWidth',0.3, ...
    'FaceAlpha',1.0);
end


function draw_cylinders(ax, xC, yC, z0, zTop, cVal, r, nTheta)


th = linspace(0, 2*pi, nTheta);
[Th, Z] = meshgrid(th, [0 1]);
X = cos(Th);
Y = sin(Th);

N = numel(xC);
for i = 1:N
    top = zTop(i);

   
    if top >= z0
        h = top - z0;
        zBase = z0;
    else
        h = z0 - top;
        zBase = top;
    end

    XX = r*X + xC(i);
    YY = r*Y + yC(i);
    ZZ = h*Z + zBase;

   
    CC = cVal(i) * ones(size(ZZ));

    surf(ax, XX, YY, ZZ, CC, ...
        'FaceColor','flat', ...
        'CDataMapping','scaled', ...
        'EdgeColor',[0.2 0.2 0.2], ...
        'LineWidth',0.25, ...
        'FaceAlpha',1.0);
end
end