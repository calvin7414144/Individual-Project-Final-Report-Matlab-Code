clear; clc; close all;

N  = 18;
R0 = 210;

theta_full = linspace(0, 2*pi, N+1);
theta_full(end) = [];
theta = theta_full;

Px = R0*cos(theta);
Py = R0*sin(theta);
Pxy = [Px(:), Py(:)];

f = 1e6;
c = 3e8;

lambda = c/f;
k = 2*pi/lambda;

if lambda <= 0
    error('lambda must be positive.');
end

alpha_s = 0.25;

eps_w = 0.35;
q     = 0.7;

if eps_w < 0
    warning('eps_w is negative.');
end
if q <= 0
    warning('q is not positive.');
end

lw = 2.0;
ms = 3;

Tx_list = {
    [ 300,  300],  'Tx (300, 300)';
    [ 300,    0],  'Tx (300, 0)';
    [-300,  300],  'Tx (-300, 300)';
    [   0, -300],  'Tx (0, -300)';
    [-300, -300],  'Tx (-300, -300)';
    [   0,  300],  'Tx (0, 300)'
};

nCases = size(Tx_list,1);

for case_id = 1:nCases

    Tx        = Tx_list{case_id,1};
    case_name = Tx_list{case_id,2};

    if numel(Tx) ~= 2
        error('Tx must be a 1x2 vector.');
    end

    theta_Tx = atan2(Tx(2), Tx(1));

    E0 = zeros(1, N);

    for j = 1:N

        rx = Px(j);
        ry = Py(j);

        dx = rx - Tx(1);
        dy = ry - Tx(2);
        r  = hypot(dx, dy);

        if r == 0
            r = eps;
        end

        phaseTerm = exp(-1i*k*r);
        ampTerm   = 1 / r;

        E0(j) = phaseTerm * ampTerm;
    end

    E0_abs = abs(E0);

    Es = zeros(1, N);

    for j = 1:N

        dtheta = theta(j) - theta_Tx;

        s = abs(sin(dtheta));

        inside = s.^2 + eps_w^2;
        w_ang  = (sqrt(inside)).^q;

        Es(j) = alpha_s * w_ang * E0(j);
    end

    [~, idx_min] = min(abs(theta - theta_Tx));
    Es(idx_min) = 0;

    Es_abs = abs(Es);

    E_total     = E0 + Es;
    E_total_abs = abs(E_total);

    fig = figure('Color','w','Position',[80 80 1400 520]);
    tlo = tiledlayout(fig, 1, 2, 'TileSpacing','compact','Padding','compact');

    ax1 = nexttile(tlo, 1);

    sc = scatter(ax1, Px, Py, 120, E_total_abs, 'filled');
    hold(ax1,'on');

    plot(ax1, Tx(1), Tx(2), 'rp', 'MarkerSize', 16, 'MarkerFaceColor','r');

    for j = 1:N
        text(Px(j)+6, Py(j)+6, num2str(j), 'FontSize',10);
    end

    pad = 0.15;
    lim = 300 * (1 + pad);
    xlim(ax1, [-lim lim]);
    ylim(ax1, [-lim lim]);

    colormap(ax1, turbo);
    colorbar(ax1);
    axis(ax1, 'equal');
    grid(ax1, 'on');
    xlabel(ax1, 'x (m)');
    ylabel(ax1, 'y (m)');
    title(ax1, ['(a) Spatial |E_{total}| — ' case_name]);

    ax2 = nexttile(tlo, 2);
    hold(ax2, 'on');

    h1 = plot(ax2, 1:N, E_total_abs, '-o', 'LineWidth', lw, 'MarkerSize', ms);
    h2 = plot(ax2, 1:N, E0_abs,      '-s', 'LineWidth', lw, 'MarkerSize', ms);
    h3 = plot(ax2, 1:N, Es_abs,      '-d', 'LineWidth', lw, 'MarkerSize', ms);

    h1.MarkerFaceColor = h1.Color;
    h2.MarkerFaceColor = h2.Color;
    h3.MarkerFaceColor = h3.Color;

    xlim(ax2, [0.5 18.5]);
    ylim(ax2, [-0.001 0.018]);

    set(ax2, 'YScale', 'linear');

    grid(ax2, 'on');
    xlabel(ax2, 'Receiver Index');
    ylabel(ax2, 'Field Amplitude (linear)');
    legend(ax2, '|E_{total}|','|E_0|','|E_s|', 'Location','best');
    title(ax2, ['(b) Field comparison — ' case_name]);

end
