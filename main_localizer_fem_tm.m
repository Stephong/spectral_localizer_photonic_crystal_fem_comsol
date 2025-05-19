clear 
close all
clc
%%
base_folder = "";
addpath(base_folder)
addpath(base_folder + "base_code")
%%
eps0_const = 8.8542*1e-12; % F/m
mu0_const = 1.2566*1e-6; % H/m
c_const = 299792458; % m/s
%%
% init_comsol_server();
%%

%% Get comsol data

disp('===== Get comsol data =====')
%%
a = 1e-6
%%
% need to import since not using commands in matlab terminal
import com.comsol.model.*
import com.comsol.model.util.*

ModelUtil.clear
[model, ~] = mphopen('main_design.mph', 'model')
%%
% launch Comsol GUI corresponding to model
% mphlaunch(model)
%%
[mphmatrix_str, mphxmeshinfo_str] = get_mphmatrix_str_model(model);
%%
import com.comsol.model.*
import com.comsol.model.util.*
ModelUtil.clear

clear model 
%%

%% Get lattice xtended mesh

disp('===== Get lattice xtended mesh =====')
%%
lattice_xnode = transpose(mphxmeshinfo_str.dofs.coords);
N_dof = size(lattice_xnode, 1)
%%

%% Get model matrix

disp('===== Get model matrix =====')
%%
[Ec, Dc, Kc, Ud, Null, Nullf, Uscale] = get_mphmatrix(mphmatrix_str);
%%
freq_norm = 2*pi*c_const/a;
get_H_eff = @(omega) (-1j*omega*freq_norm)^2 * Ec - (-1j*omega*freq_norm) * Dc + Kc;
%%

%% Localizer pre

disp('===== Localizer pre =====')
%%
[x0, y0, E0] = deal(0*a, 0*a, 0.37)
%%
X = transpose(Nullf) * spdiags(lattice_xnode(:,1), 0, N_dof, N_dof) * Null;
Y = transpose(Nullf) * spdiags(lattice_xnode(:,2), 0, N_dof, N_dof) * Null;
Iden = transpose(Nullf) * speye(size(N_dof,1)) * Null;
N_dofc = size(X, 1)
%%
% H_eff = H_eff_herm + 1j*H_eff_nonherm
H_eff = get_H_eff(E0);
H_eff_herm = 0.5 * (H_eff + H_eff');
H_eff_nonherm = -1j * 0.5 * (H_eff - H_eff');
ratio_herm = norm(H_eff_herm, 1) / norm(H_eff_nonherm, 1)
ishermitian(H_eff)
%%
X_unit = norm(X,1)
H_unit = norm(H_eff,1)/1e4
kappa_unit = H_unit/X_unit*1
%%

%% Localizer (x,y)

disp('===== Localizer (x,y) =====')
%%
N_x_lambd_0 = 10;
% N_y_lambd_0 = 100;
% tab_x_lambd_0 = linspace(min(lattice_xnode(:,1)), max(lattice_xnode(:,1)), N_x_lambd_0);
% tab_y_lambd_0 = linspace(min(lattice_xnode(:,2)), max(lattice_xnode(:,2)), N_y_lambd_0);
tab_E_lambd_0 = [E0];

tab_x_lambd_0 = linspace(min(lattice_xnode(:,1)), x0, N_x_lambd_0);
% tab_x_lambd_0 = tab_x_lambd_0(end);
tab_y_lambd_0 = [y0];
i_y0 = 1;

kappa_0 = 1.5 * kappa_unit;
N_spec = 10;

H_eff = sparse( get_H_eff(tab_E_lambd_0(1)) );
H = 0.5 * (H_eff + H_eff');

localizer_eval_0 = get_localizer_comsol_eval(tab_x_lambd_0, tab_y_lambd_0, tab_E_lambd_0*0, Iden, X, Y, H, kappa_0, N_spec=N_spec);

tab_L_spec_0 = localizer_eval_0.tab_L_spec;
tab_L_gap_0 = localizer_eval_0.tab_L_gap;
tab_C_L_0 = localizer_eval_0.tab_C_L;
%%
figure
colororder({'b','r'})
yyaxis left
for i_spec=1:size(tab_L_spec_0,4)
    plot(tab_x_lambd_0, tab_L_spec_0(1, :, i_y0, i_spec)/H_unit, 'b.')
    hold on
end
yline(0)
ylabel('\sigma(L_\lambda) (in H unit)')
% ylim([-1, 1])
ylim([-1, 1]*0.4)
yyaxis right
yline(0)
plot(tab_x_lambd_0, squeeze(tab_C_L_0(1,:,i_y0)), 'ro-')
yticks([0:1:1])
ylabel('C_L')
% ylim([0, 1])
title("y=" + tab_y_lambd_0(i_y0) + ", E=" + tab_E_lambd_0(1) + ", \kappa=" + kappa_0/kappa_unit + "(in kappa unit)")
xlabel('Position, x')
xlim([tab_x_lambd_0(1), tab_x_lambd_0(end)])
hold off

figure
plot(tab_x_lambd_0, tab_L_gap_0/H_unit, 'bx')
ylim([0, inf])
xlabel('Position, x')
ylabel('\mu_\lambda^C (in H unit)')
%%

%% Localizer (E)

disp('===== Localizer (E) =====')
%%
N_E_lambd_1 = 5;
tab_x_lambd_1 = [x0];
tab_y_lambd_1 = [y0];
tab_E_lambd_1 = linspace(0.05, 0.5, N_E_lambd_1);

kappa_1 = 1.5 * kappa_unit;
N_spec = 10;

tab_L_spec_1 = zeros(N_E_lambd_1, N_spec);
tab_L_gap_1 = zeros(N_E_lambd_1, 1);
tab_C_L_1 = zeros(N_E_lambd_1, 1);

tic
for ii=1:N_E_lambd_1
    fprintf([repmat('\b', 1, 200), '%f '], ii/N_E_lambd_1);
    H_eff = sparse( get_H_eff(tab_E_lambd_1(ii)) );
    H = 0.5 * (H_eff + H_eff');
    localizer_eval_1 = get_localizer_comsol_eval(tab_x_lambd_1, tab_y_lambd_1, tab_E_lambd_1(ii)*0, Iden, X, Y, H, kappa_0, N_spec=N_spec, ...
         Doclean=1);
    tab_L_spec_1(ii, :) = localizer_eval_1.tab_L_spec;
    tab_L_gap_1(ii) = localizer_eval_1.tab_L_gap;
    tab_C_L_1(ii) = localizer_eval_1.tab_C_L;
end
toc
%%
figure
colororder({'b','r'})
yyaxis left
for i_spec=1:size(tab_L_spec_1,2)
    plot(tab_E_lambd_1, tab_L_spec_1(:, i_spec)/H_unit, 'b.')
    hold on
end
yline(0)
ylabel('\sigma(L_\lambda) (in H unit)')
ylim([-1, 1]*0.1)
% ylim([-1, 1]*1*1e-3)
yyaxis right
yline(0)
plot(tab_E_lambd_1, tab_C_L_1, 'ro-')
% yticks([0:1:1])
ylabel('C_L')
% ylim([0, 1])
title("x=" + tab_x_lambd_1(1) + ", y=" + tab_y_lambd_1(1) + ", \kappa=" + kappa_1/kappa_unit + "(in kappa unit)")
xlabel('Energy, E')
xlim([tab_E_lambd_1(1), tab_E_lambd_1(end)])
hold off

figure
plot(tab_E_lambd_1, tab_L_gap_1/H_unit, 'bx')
ylim([0, inf])
title("x=" + tab_x_lambd_1(1) + ", y=" + tab_y_lambd_1(1) + ", \kappa=" + kappa_1/kappa_unit + "(in kappa unit)")
xlabel('Energy, E')
ylabel('\mu_\lambda^C (in H unit)')
%%

%%

%%

%%

%%

%%

%% Save

[file_pathstr, file_name, file_ext] = fileparts(matlab.desktop.editor.getActiveFilename);
save(file_name+".mat")
%%