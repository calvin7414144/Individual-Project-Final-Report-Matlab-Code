
%  forward_model


clear; clc;


f  = 1e6;                
c0 = 3e8;
lambda0 = c0 / f;
k0 = 2*pi/lambda0;         
mu0 = 4*pi*1e-7;
eps0 = 1/(mu0*c0^2);
eta0 = sqrt(mu0/eps0);



eta_dry = 0.05;     
eta_sea = 0;       


ring_segs = [4, 12, 20, 28, 36];
R_surface = 150;         

cells = generate_cells_circular(R_surface, ring_segs);
Ncells = length(cells);


R_ant = 300;              
theta = 0;                
tx_pos = [R_ant*cos(theta), R_ant*sin(theta)];


eta_cells = eta_dry * ones(Ncells,1);


E1 = zeros(Ncells,1);
for i = 1:Ncells
    Pi = [cells(i).x, cells(i).y];
    E1(i) = compute_E1_simple(tx_pos, Pi, k0);
end


A = eye(Ncells);
b = E1;

for i = 1:Ncells
    Pi = [cells(i).x, cells(i).y];
    for j = 1:Ncells
        Pj = [cells(j).x, cells(j).y];
        aj = cells(j).a;


        I_ij = compute_Iij_11b_11c(Pi, Pj, aj, k0, i, j);

        delta_eta = eta_cells(j) - eta_sea;


        A(i,j) = A(i,j) - I_ij * delta_eta;
    end
end


E = A \ b;


figure('Color','w','Position',[100 100 1100 450]);


ax1 = subplot(1,2,1);

scatter(ax1, [cells.x], [cells.y], 50, real(E), 'filled');

axis(ax1,'equal');


xlim(ax1, [-200 200]);
ylim(ax1, [-200 200]);

colorbar(ax1);
title(ax1,'Real(E_z)');


set(ax1,'Box','on','LineWidth',1.8,...
    'FontWeight','bold',...         
    'FontName','Times New Roman'); 


ax2 = subplot(1,2,2);

scatter(ax2, [cells.x], [cells.y], 50, imag(E), 'filled');

axis(ax2,'equal');

 
xlim(ax2, [-200 200]);
ylim(ax2, [-200 200]);

colorbar(ax2);
title(ax2,'Imag(E_z)');

 
set(ax2,'Box','on','LineWidth',1.8,...
    'FontWeight','bold',...
    'FontName','Times New Roman');
disp('FWP (1)–(14)');

function cells = generate_cells_circular(R, segs)
N = sum(segs);
cells(N) = struct('x',0,'y',0,'a',0);

Atotal = pi*R^2;
Acell  = Atotal/N;
a_eq   = sqrt(Acell/pi);

r_in = 0;
idx  = 1;

for k = 1:length(segs)
    nseg = segs(k);
    A_ring = nseg * Acell;
    r_out = sqrt(r_in^2 + A_ring/pi);
    r_mid = (r_in + r_out)/2;

    dtheta = 2*pi/nseg;
    for m = 1:nseg
        th = (m-0.5)*dtheta;
        cells(idx).x = r_mid*cos(th);
        cells(idx).y = r_mid*sin(th);
        cells(idx).a = a_eq;
        idx = idx + 1;
    end
    r_in = r_out;
end
end

function E1 = compute_E1_simple(tx, P, k0)
R = sqrt((tx(1)-P(1))^2 + (tx(2)-P(2))^2);
E1 = exp(-1i*k0*R)/R;   % 
end

function I_ij = compute_Iij_11b_11c(Pi, Pj, a, k0, i, j)
Rij = sqrt( (Pi(1)-Pj(1))^2 + (Pi(2)-Pj(2))^2 );

if i == j

    I_ij = (2*pi/(1i*k0)) * (1 - exp(-1i*k0*a));
else

    if Rij < 1e-6, Rij = 1e-6; end
    I_ij = (exp(-1i*k0*Rij)/Rij) * (2*pi/k0^2) * (k0*a)*besselj(1, k0*a);
end
end


function EsB = compute_Es_15b(B, cells, E, eta_cells, eta_ref, k0)


N = length(cells);
EsB = 0;

for n = 1:N

    Pn = [cells(n).x, cells(n).y];
    En = E(n);
    an = cells(n).a;

    rhoBn = sqrt((B(1)-Pn(1))^2 + (B(2)-Pn(2))^2);

    delta_eta = eta_cells(n) - eta_ref;

    Fw2n = delta_eta;

    term = En * exp(-1i*k0*rhoBn) / (1i*k0*rhoBn) ...
           * delta_eta ...
           * Fw2n * (k0*an)*besselj(1, k0*an);

    EsB = EsB + term;
end
end

B = [300, 0];    

 
Es = compute_Es_15b(B, cells, E, eta_cells, eta_sea, k0);

E1_B = compute_E1_simple(tx_pos, B, k0);

E_total = E1_B + Es;

fprintf("E_total(B) = %.4e + j%.4e\n", real(E_total), imag(E_total));


x_range = linspace(-500, 500, 200);
E_total_array = zeros(size(x_range));

for idx = 1:length(x_range)
    B = [x_range(idx), 0];   
    
    
    Es_temp = compute_Es_15b(B, cells, E, eta_cells, eta_sea, k0);

   
    E1_temp = compute_E1_simple(tx_pos, B, k0);

    
    E_total_array(idx) = E1_temp + Es_temp;
end

figure;
plot(x_range, real(E_total_array), 'LineWidth', 2); hold on;
plot(x_range, imag(E_total_array), 'LineWidth', 2);
legend('Real(E total)', 'Imag(E total)');
xlabel('x (m)');
ylabel('Field amplitude');
title('Total Field E_{total}(B) along the x-axis');
grid on;

x = -500:20:500;
y = -500:20:500;

[E_X, E_Y] = meshgrid(x, y);
E_map = zeros(size(E_X));

for ix = 1:length(x)
    for iy = 1:length(y)
        B = [x(ix), y(iy)];
        Es_temp = compute_Es_15b(B, cells, E, eta_cells, eta_sea, k0);
        E1_temp = compute_E1_simple(tx_pos, B, k0);
        E_map(iy, ix) = abs(E1_temp + Es_temp);
    end
end

figure;
imagesc(x, y, E_map);
axis equal; colorbar;
title('|E_{total}(x,y)| 2D Distribution');



E_total_surface = zeros(Ncells,1);
for i = 1:Ncells
    Pi = [cells(i).x, cells(i).y];
    Es_i = compute_Es_15b(Pi, cells, E, eta_cells, eta_sea, k0);
    E1_i = compute_E1_simple(tx_pos, Pi, k0);
    E_total_surface(i) = E1_i + Es_i;
end

