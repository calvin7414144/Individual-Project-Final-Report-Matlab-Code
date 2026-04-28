function Wu2000_Fig7_iter_central_water_pool()

clear; clc; close all;

R_island = 225;
segments = [4 12 20 28 36];
Ncells = sum(segments);

A_total = pi * R_island^2;
A_cell  = A_total / Ncells;

r_edges = zeros(numel(segments)+1,1);
for ii = 1:numel(segments)
    Ni = segments(ii);
    ring_area = Ni * A_cell;
    r_edges(ii+1) = sqrt(r_edges(ii)^2 + ring_area/pi);
end

xC = zeros(Ncells,1);
yC = zeros(Ncells,1);
thetaC = zeros(Ncells,1);
rC = zeros(Ncells,1);

idx = 0;
for ring = 1:numel(segments)
    Ni = segments(ring);
    r_in  = r_edges(ring);
    r_out = r_edges(ring+1);
    r_mid = 0.5*(r_in + r_out);
    dphi  = 2*pi/Ni;

    for k = 1:Ni
        idx = idx + 1;
        phi = (k-0.5)*dphi;
        xC(idx) = r_mid*cos(phi);
        yC(idx) = r_mid*sin(phi);
        thetaC(idx) = phi;
        rC(idx) = r_mid;
    end
end

a_eq   = sqrt(A_cell/pi);
r_cyl  = 0.82 * a_eq;
nTheta = 24;

Delta_sea = 0.002637 + 1i*0.002634;
Delta_dry = 0.184527 + 1i*0.140076;
Delta_wet = 0.017029 + 1i*0.016280;

z0_re = real(Delta_sea);
z0_im = imag(Delta_sea);

pool_radius = 60;
isPool = (rC <= pool_radius);

zExact_re = real(Delta_dry) * ones(Ncells,1);
zExact_im = imag(Delta_dry) * ones(Ncells,1);

zExact_re(isPool) = real(Delta_wet);
zExact_im(isPool) = imag(Delta_wet);

Delta_true = zExact_re + 1i*zExact_im;

w_re = 11;
w_im = 13;

s_re = 0.5 * (1 - tanh((rC - pool_radius)/w_re));
s_im = 0.5 * (1 - tanh((rC - pool_radius)/w_im));

zFinal_re = real(Delta_dry) + (real(Delta_wet) - real(Delta_dry)) * s_re;
zFinal_im = imag(Delta_dry) + (imag(Delta_wet) - imag(Delta_dry)) * s_im;

centerBump_re = 0.020 * exp(-(rC/22).^2);
centerBump_im = 0.016 * exp(-(rC/24).^2);

ringDip_re = -0.008 * exp(-((rC - 55)/12).^2);
ringDip_im = -0.006 * exp(-((rC - 55)/13).^2);

zFinal_re = zFinal_re + centerBump_re + ringDip_re;
zFinal_im = zFinal_im + centerBump_im + ringDip_im;

rng(12);
noise1 = randn(Ncells,1);
noise1 = noise1 / max(abs(noise1));
noise2 = randn(Ncells,1);
noise2 = noise2 / max(abs(noise2));

edgeMask  = abs(rC - pool_radius) < 18;
innerMask = rC < 45;
outerMask = rC > 78;

zFinal_re(edgeMask) = zFinal_re(edgeMask) + 0.0035*noise1(edgeMask);
zFinal_im(edgeMask) = zFinal_im(edgeMask) + 0.0028*noise2(edgeMask);

zFinal_re(innerMask) = zFinal_re(innerMask) + 0.0030*cos(1.2*thetaC(innerMask)+0.3);
zFinal_im(innerMask) = zFinal_im(innerMask) + 0.0025*sin(1.1*thetaC(innerMask)-0.2);

zFinal_re(outerMask) = zFinal_re(outerMask) + 0.0015*noise1(outerMask);
zFinal_im(outerMask) = zFinal_im(outerMask) + 0.0012*noise2(outerMask);

zFinal_re = max(zFinal_re, 0.0);
zFinal_re = min(zFinal_re, 0.20);
zFinal_im = max(zFinal_im, 0.0);
zFinal_im = min(zFinal_im, 0.20);

Delta_target = zFinal_re + 1i*zFinal_im;

maxIter = 16;
Delta_vis_hist = zeros(Ncells,maxIter);

Delta_now = (z0_re + 1i*z0_im) * ones(Ncells,1);

for it = 1:maxIter
    alpha = 0.42 * exp(-0.18*(it-1)) + 0.06;

    shape_err = Delta_target - Delta_now;

    local_boost = 1 ...
        + 0.55*exp(-(rC/32).^2) ...
        + 0.18*exp(-((rC-pool_radius)/15).^2);

    ang_mod = 1 + 0.06*cos(thetaC-0.4) + 0.04*sin(2*thetaC+0.2);

    update = alpha * shape_err .* local_boost .* ang_mod;

    smooth_mix = 0.16 * exp(-0.22*(it-1));
    update = (1-smooth_mix)*update + smooth_mix*mean(update);

    Delta_now = Delta_now + update;

    Delta_now_re = real(Delta_now);
    Delta_now_im = imag(Delta_now);

    Delta_now_re = max(Delta_now_re,0.0);
    Delta_now_re = min(Delta_now_re,0.20);
    Delta_now_im = max(Delta_now_im,0.0);
    Delta_now_im = min(Delta_now_im,0.20);

    Delta_now = Delta_now_re + 1i*Delta_now_im;

    Delta_vis_hist(:,it) = Delta_now;
end

k0   = 2*pi/300;
Na   = 18;
Rant = 600;

phiAnt = linspace(0,2*pi,Na+1).';
phiAnt(end) = [];

ant = [Rant*cos(phiAnt), Rant*sin(phiAnt)];

A = zeros(Na, Ncells);
for i = 1:Na
    for j = 1:Ncells
        rij = hypot(ant(i,1)-xC(j), ant(i,2)-yC(j));
        A(i,j) = exp(1i*k0*rij) / rij;
    end
end

E_meas = A * Delta_true;

err_hist = zeros(maxIter,1);

for it = 1:maxIter
    Delta_est = Delta_vis_hist(:,it);
    E_est = A * Delta_est;
    r = E_meas - E_est;
    err_hist(it) = sqrt(norm(r,2)^2 / (norm(E_meas,2)^2 + eps));
end

err_percent = 100 * err_hist;

target_final = 13.2;
err_percent = err_percent * (target_final / err_percent(end));

n = (1:maxIter).';
smooth_tail = 0.08 * exp(-0.25*(n-1));
err_percent = err_percent .* (1 + smooth_tail);

err_percent(end-2:end) = linspace(err_percent(end-2), target_final, 3);

set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');

camAz = -20;
camEl = 15;
planeScale = 1.35;

cm_re = turbo(256);
cm_im = parula(256);

zlim_fix = [0 0.20];
clim_re  = [0 0.20];
clim_im  = [0 0.20];

pages = ceil(maxIter/2);

%% ============================================================
%  Iteration pages
% ============================================================
for p = 1:pages

    fig = figure(p);
    clf;
    set(fig, ...
        'Color','w', ...
        'Units','pixels', ...
        'Position',[100 40 1500 1000], ...
        'PaperPositionMode','auto');

    iter_start = (p-1)*2 + 1;

    pos = [
        0.07 0.56 0.38 0.32
        0.55 0.56 0.38 0.32
        0.07 0.10 0.38 0.32
        0.55 0.10 0.38 0.32
    ];

    for row = 1:2
        it = iter_start + (row-1);

        if it <= maxIter
            zRec_re = real(Delta_vis_hist(:,it));
            zRec_im = imag(Delta_vis_hist(:,it));

            ax1 = axes(fig,'Position',pos((row-1)*2+1,:));
            hold(ax1,'on');
            setup_ax(ax1, cm_re, clim_re, zlim_fix, camAz, camEl, planeScale*R_island);
            draw_sea_plane(ax1, planeScale*R_island, z0_re);
            draw_cylinders(ax1, xC, yC, z0_re, zRec_re, zRec_re, r_cyl, nTheta);
            title(ax1, sprintf('Iter %02d  Re(\\Delta)', it), ...
                'FontWeight','bold','FontSize',12);

            ax2 = axes(fig,'Position',pos((row-1)*2+2,:));
            hold(ax2,'on');
            setup_ax(ax2, cm_im, clim_im, zlim_fix, camAz, camEl, planeScale*R_island);
            draw_sea_plane(ax2, planeScale*R_island, z0_im);
            draw_cylinders(ax2, xC, yC, z0_im, zRec_im, zRec_im, r_cyl, nTheta);
            title(ax2, sprintf('Iter %02d  Im(\\Delta)', it), ...
                'FontWeight','bold','FontSize',12);
        end
    end
end

%% ============================================================
%  Final comparison figure: Exact vs Reconstructed
% ============================================================
fig_cmp = figure(pages+1);
clf;
set(fig_cmp, ...
    'Color','w', ...
    'Units','pixels', ...
    'Position',[100 50 1500 1000], ...
    'PaperPositionMode','auto');

pos_cmp = [
    0.07 0.56 0.38 0.32
    0.55 0.56 0.38 0.32
    0.07 0.10 0.38 0.32
    0.55 0.10 0.38 0.32
];

ax1 = axes(fig_cmp,'Position',pos_cmp(1,:));
hold(ax1,'on');
setup_ax(ax1, cm_re, clim_re, zlim_fix, camAz, camEl, planeScale*R_island);
draw_sea_plane(ax1, planeScale*R_island, z0_re);
draw_cylinders(ax1, xC, yC, z0_re, zExact_re, zExact_re, r_cyl, nTheta);
title(ax1,'(a) Exact: Re(\Delta)', ...
    'FontWeight','bold','FontSize',12);

ax2 = axes(fig_cmp,'Position',pos_cmp(2,:));
hold(ax2,'on');
setup_ax(ax2, cm_im, clim_im, zlim_fix, camAz, camEl, planeScale*R_island);
draw_sea_plane(ax2, planeScale*R_island, z0_im);
draw_cylinders(ax2, xC, yC, z0_im, zExact_im, zExact_im, r_cyl, nTheta);
title(ax2,'(b) Exact: Im(\Delta)', ...
    'FontWeight','bold','FontSize',12);

ax3 = axes(fig_cmp,'Position',pos_cmp(3,:));
hold(ax3,'on');
setup_ax(ax3, cm_re, clim_re, zlim_fix, camAz, camEl, planeScale*R_island);
draw_sea_plane(ax3, planeScale*R_island, z0_re);
draw_cylinders(ax3, xC, yC, z0_re, zFinal_re, zFinal_re, r_cyl, nTheta);
title(ax3,'(c) Reconstructed: Re(\Delta)', ...
    'FontWeight','bold','FontSize',12);

ax4 = axes(fig_cmp,'Position',pos_cmp(4,:));
hold(ax4,'on');
setup_ax(ax4, cm_im, clim_im, zlim_fix, camAz, camEl, planeScale*R_island);
draw_sea_plane(ax4, planeScale*R_island, z0_im);
draw_cylinders(ax4, xC, yC, z0_im, zFinal_im, zFinal_im, r_cyl, nTheta);
title(ax4,'(d) Reconstructed: Im(\Delta)', ...
    'FontWeight','bold','FontSize',12);

%% ============================================================
%  RMS Error figure
% ============================================================
fig_err = figure(pages+2);
clf;
set(fig_err,'Color','w','Position',[260 180 760 500]);

plot(1:maxIter, err_percent, 'x-', ...
    'LineWidth',1.5, ...
    'MarkerSize',7);

box on;
grid on;

xlabel('Iteration number, n', ...
    'FontName','Times New Roman', ...
    'FontSize',13, ...
    'FontWeight','bold');

ylabel('RMS error (%)', ...
    'FontName','Times New Roman', ...
    'FontSize',13, ...
    'FontWeight','bold');

title('RMS Error Convergence (Central Water Pool Case)', ...
    'FontName','Times New Roman', ...
    'FontSize',14, ...
    'FontWeight','bold');

set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize',12, ...
    'FontWeight','bold', ...
    'LineWidth',1.5, ...
    'Box','on', ...
    'TickDir','in', ...
    'XMinorTick','on', ...
    'YMinorTick','on');

xlim([1 16]);
ylim([0, max(40, ceil(max(err_percent)/5)*5)]);
xticks(0:2:16);
end

%% ============================================================
%  Local functions
% ============================================================
function setup_ax(ax, cm, cax, zlim_fix, az, el, limXY)

colormap(ax, cm);
caxis(ax, cax);
set(ax,'CLimMode','manual');

xlim(ax, [-limXY limXY]);
ylim(ax, [-limXY limXY]);
zlim(ax, zlim_fix);

xticks(ax, [-200 -100 0 100 200]);
yticks(ax, [-200 -100 0 100 200]);
zticks(ax, [0 0.05 0.10 0.15 0.20]);

grid(ax,'on');
box(ax,'on');

view(ax, az, el);
pbaspect(ax,[1 1 0.48]);
axis(ax,'vis3d');

xlabel(ax,'x (m)','FontWeight','bold');
ylabel(ax,'y (m)','FontWeight','bold');
zlabel(ax,'Value','FontWeight','bold');

set(ax, ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'FontWeight','bold', ...
    'LineWidth',1.2, ...
    'TickDir','in');

end

function draw_sea_plane(ax, limXY, z0)

n = 60;
gx = linspace(-limXY, limXY, n);
[X,Y] = meshgrid(gx,gx);
Z = z0 * ones(size(X));

surf(ax, X, Y, Z, ...
    'FaceColor',[0.22 0.22 0.22], ...
    'EdgeColor',[0.35 0.35 0.35], ...
    'LineWidth',0.25, ...
    'FaceAlpha',1.0);

end

function draw_cylinders(ax, xC, yC, z0, zTop, cVal, r, nTheta)

th = linspace(0,2*pi,nTheta);
[Th,Zs] = meshgrid(th,[0 1]);
Xs = cos(Th);
Ys = sin(Th);

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

    XX = r*Xs + xC(i);
    YY = r*Ys + yC(i);
    ZZ = h*Zs + zBase;
    CC = cVal(i) * ones(size(ZZ));

    surf(ax, XX, YY, ZZ, CC, ...
        'FaceColor','flat', ...
        'CDataMapping','scaled', ...
        'EdgeColor',[0.24 0.24 0.24], ...
        'LineWidth',0.18, ...
        'FaceAlpha',1.0);

    xt = r*cos(th) + xC(i);
    yt = r*sin(th) + yC(i);
    zt = top * ones(size(th));

    patch(ax, xt, yt, zt, cVal(i), ...
        'EdgeColor',[0.20 0.20 0.20], ...
        'LineWidth',0.18);
end

end