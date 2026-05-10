clear; clc; close all;


R0 = 210;
ring_radii = linspace(0, R0, 6);
segments = [4, 12, 20, 28, 36];
N = sum(segments);

Px = []; Py = [];

for k = 1:5
    r = (ring_radii(k) + ring_radii(k+1)) / 2;
    th = linspace(0, 2*pi, segments(k)+1);
    th(end) = [];
    Px = [Px, r*cos(th)];
    Py = [Py, r*sin(th)];
end


Tx_list = [
    300, 0;
    -200, 200;
    500, 0;
];

Tx_names = {"Tx (300,0)", "Tx (-200,200)", "Tx (500,0)"};


f = 1e6;
c = 3e8;
lambda = c/f;
k = 2*pi/lambda;


set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');


for case_id = 1:3

    Tx = Tx_list(case_id,:);


    E = zeros(1, N);

    for j = 1:N
        r = hypot(Px(j)-Tx(1), Py(j)-Tx(2));
        E(j) = exp(-1i * k * r) / r;
    end

    Eabs = abs(E);

    %% ----------  ----------
    [~, idx_max] = max(Eabs);

    Px2 = Px; Py2 = Py; Eabs2 = Eabs;
    Px2(idx_max) = [];
    Py2(idx_max) = [];
    Eabs2(idx_max) = [];

   
    fig = figure('Color','w','Position',[80 80 1400 520]);
    tlo = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

  
    ax1 = nexttile(tlo,1);

    scatter(ax1, Px, Py, 120, Eabs, 'filled'); hold(ax1,'on');

    plot(ax1, Tx(1), Tx(2), 'rp', ...
        'MarkerSize',16,'MarkerFaceColor','r');

    text(ax1, Tx(1)+10, Tx(2), Tx_names{case_id}, ...
        'Color','r','FontSize',12,'FontWeight','bold');

    colormap(ax1, turbo);
    cb = colorbar(ax1);
    cb.FontWeight = 'bold';
    cb.LineWidth = 1.5;

    axis(ax1,'equal'); grid(ax1,'on');

    xlabel(ax1,'x (m)','FontWeight','bold');
    ylabel(ax1,'y (m)','FontWeight','bold');

    title(ax1, ['(a)100 Cells distribution of |E|' Tx_names{case_id}], ...
        'FontWeight','bold');

    
    set(ax1,'FontName','Times New Roman',...
        'FontSize',11,...
        'FontWeight','bold',...
        'LineWidth',1.8,...
        'Box','on',...
        'TickDir','in',...
        'TickLength',[0.015 0.015]);

 
    ax2 = nexttile(tlo,2);
    hold(ax2,'on');

    scatter(ax2, 1:length(Eabs2), Eabs2, 80, Eabs2, 'filled');

    colormap(ax2, turbo);
    cb2 = colorbar(ax2);
    cb2.FontWeight = 'bold';
    cb2.LineWidth = 1.5;

    grid(ax2,'on');

    xlabel(ax2,'Index','FontWeight','bold');
    ylabel(ax2,'|E(P)|','FontWeight','bold');

    title(ax2, ['(b) |E(P)| vs Index (Tx Outside) ' ...
         Tx_names{case_id}], ...
        'FontWeight','bold');

     
    set(ax2,'FontName','Times New Roman',...
        'FontSize',11,...
        'FontWeight','bold',...
        'LineWidth',1.8,...
        'Box','on',...
        'TickDir','in',...
        'TickLength',[0.015 0.015]);

end