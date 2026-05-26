%% =========================================================
%% Backward bifurcation diagram for the BPD model (Final Version)
%% =========================================================
clear; clc; close all;

%% ================= PARAMETERS =================
p.Lambda = 12;
p.alpha1 = 0.05;
p.alpha2 = 0.027;
p.pi1    = 1;
p.pi2    = 1;
p.eta    = 0.9;
p.gamma  = 0.6;
p.kappa  = 0.8;
p.theta  = 0.8;
p.psi    = 0.5;
p.mu     = 0.01;
p.delta  = 0.1;

% Fix beta2, vary beta1
p.beta2 = 0.0;

%% ================= THRESHOLD =================
K1 = p.alpha1 + p.mu;
K2 = p.alpha2 + p.mu;
K4 = p.gamma + p.mu;
K5 = p.kappa + p.theta + p.psi + p.mu + p.delta;
K7 = p.mu + p.delta;

Sc0 = p.Lambda / K1;
Sp0 = p.alpha1 * Sc0 / K2;

beta1_crit = (K4*K5*K7) / (p.gamma*p.theta*(Sc0+Sp0));
fprintf('Critical beta1 for R0 = 1: %.8f\n', beta1_crit);

beta1_vals = linspace(0, 3*beta1_crit, 1500);

%% ================= STORAGE =================
R0_vals = zeros(size(beta1_vals));

stable_endemic_R = [];
stable_endemic_I = [];

unstable_endemic_R = [];
unstable_endemic_I = [];

%% ================= LOOP =================
for k = 1:length(beta1_vals)
    beta1 = beta1_vals(k);

    R0 = reproduction_number(beta1, p);
    R0_vals(k) = R0;

    Iroots = find_roots(beta1, p);

    for j = 1:length(Iroots)
        Ieq = Iroots(j);
        Xeq = reconstruct_eq(Ieq, beta1, p);
        J = jacobian_bpd(Xeq, beta1, p);

        if max(real(eig(J))) < -1e-8
            stable_endemic_R(end+1) = R0;
            stable_endemic_I(end+1) = Ieq;
        else
            unstable_endemic_R(end+1) = R0;
            unstable_endemic_I(end+1) = Ieq;
        end
    end
end

%% ================= DFE =================
I_dfe = zeros(size(R0_vals));
stable_dfe = R0_vals < 1;
unstable_dfe = R0_vals >= 1;

%% ================= SORT =================
[stable_endemic_R, idx1] = sort(stable_endemic_R);
stable_endemic_I = stable_endemic_I(idx1);

[unstable_endemic_R, idx2] = sort(unstable_endemic_R);
unstable_endemic_I = unstable_endemic_I(idx2);

%% ================= Rc =================
Rc = min([stable_endemic_R, unstable_endemic_R]);
fprintf('Estimated Rc = %.6f\n', Rc);

%% ================= COLORS =================
col_stable_dfe   = [0.0 0.5 0.0];       % green
col_unstable_dfe = [0.85 0.325 0.098];  % red-orange
col_stable_end   = [0.0 0.447 0.741];   % blue
col_unstable_end = [0.0 0.0 0.0];       % black

%% ================= PLOT =================
figure('Color','w','Position',[100 100 900 650]);
hold on; box on;
plot(R0_vals(stable_dfe), I_dfe(stable_dfe), '-', 'Color', col_stable_dfe, 'LineWidth', 5);
plot(R0_vals(unstable_dfe), I_dfe(unstable_dfe), '-', 'Color', col_unstable_dfe, 'LineWidth', 5);
plot(stable_endemic_R, stable_endemic_I, '-', 'Color', col_stable_end, 'LineWidth', 5);
plot(unstable_endemic_R, unstable_endemic_I, '-', 'Color', col_unstable_end, 'LineWidth', 5);

% Threshold lines
xline(1, 'k:', 'LineWidth', 3);
xline(Rc, 'k--', 'LineWidth', 3);

% Labels
xlabel('$\mathcal{R}_0$', 'Interpreter','latex');
ylabel('Infected cocoa pods, $I^*$', 'Interpreter','latex');
%title('Backward bifurcation diagram for the BPD model','FontSize',18);

legend({'Stable DFE','Unstable DFE','Stable endemic equilibrium','Unstable endemic equilibrium'}, ...
        'Location','northwest');

xlim([0 1.15]);
ylim([0 1.05*max([stable_endemic_I,unstable_endemic_I,1])]);

set(gca,'FontSize',14,'LineWidth',2);
grid on;

%% =========================================================
%% ================= FUNCTIONS =============================
%% =========================================================

function R0 = reproduction_number(beta1, p)
    K1 = p.alpha1 + p.mu;
    K2 = p.alpha2 + p.mu;
    K4 = p.gamma + p.mu;
    K5 = p.kappa + p.theta + p.psi + p.mu + p.delta;
    K7 = p.mu + p.delta;
    K8 = p.mu + p.delta;

    Sc0 = p.Lambda / K1;
    Sp0 = p.alpha1 * Sc0 / K2;

    R0 = p.gamma * ( ...
        p.theta*K8*(Sc0+Sp0)*beta1 + ...
        p.psi*K7*(Sc0+Sp0)*p.beta2 ) / (K4*K5*K7*K8);
end

function Iroots = find_roots(beta1, p)
    f = @(I) scalar_eq(I, beta1, p);

    Igrid = linspace(1e-8, 50, 6000);
    F = arrayfun(f, Igrid);

    Iroots = [];
    for i = 1:length(Igrid)-1
        if F(i)*F(i+1) < 0
            r = fzero(f,[Igrid(i),Igrid(i+1)]);
            if r > 1e-8
                Iroots(end+1) = r;
            end
        end
    end

    Iroots = unique(round(Iroots,8));
end

function val = scalar_eq(I, beta1, p)
    K1 = p.alpha1 + p.mu;
    K2 = p.alpha2 + p.mu;
    K4 = p.gamma + p.mu;
    K5 = p.kappa + p.theta + p.psi + p.mu + p.delta;
    K7 = p.mu + p.delta;
    K8 = p.mu + p.delta;

    Is = (p.theta/K7)*I;
    Ip = (p.psi/K8)*I;

    lambda1 = beta1*Is*(1 + p.pi1*Is);
    lambda2 = p.beta2*Ip*(1 + p.pi2*Ip);
    L = lambda1 + lambda2;

    Sc = p.Lambda/(K1 + L);
    Sp = p.alpha1*Sc/(K2 + L);

    val = (Sc + Sp)*L - (K4*K5/p.gamma)*I;
end

function X = reconstruct_eq(I, beta1, p)
    K1 = p.alpha1 + p.mu;
    K2 = p.alpha2 + p.mu;
    K3 = p.eta + p.mu;
    K5 = p.kappa + p.theta + p.psi + p.mu + p.delta;
    K7 = p.mu + p.delta;
    K8 = p.mu + p.delta;

    Is = (p.theta/K7)*I;
    Ip = (p.psi/K8)*I;

    lambda1 = beta1*Is*(1 + p.pi1*Is);
    lambda2 = p.beta2*Ip*(1 + p.pi2*Ip);
    L = lambda1 + lambda2;

    Sc = p.Lambda/(K1 + L);
    Sp = p.alpha1*Sc/(K2 + L);
    Sr = p.alpha2*Sp/K3;
    E  = K5*I/p.gamma;
    R  = (p.eta*Sr + p.kappa*I)/p.mu;

    X = [Sc; Sp; Sr; E; I; R; Is; Ip];
end

function J = jacobian_bpd(X, beta1, p)
    Sc = X(1); Sp = X(2);
    Is = X(7); Ip = X(8);

    lam1 = beta1*Is*(1+p.pi1*Is);
    lam2 = p.beta2*Ip*(1+p.pi2*Ip);

    dlam1 = beta1*(1+2*p.pi1*Is);
    dlam2 = p.beta2*(1+2*p.pi2*Ip);

    K1 = p.alpha1 + p.mu;
    K2 = p.alpha2 + p.mu;
    K3 = p.mu + p.eta;
    K4 = p.gamma + p.mu;
    K5 = p.kappa + p.theta + p.psi + p.mu + p.delta;
    K7 = p.mu + p.delta;
    K8 = p.mu + p.delta;

    J = zeros(8);

    J(1,1)=-(lam1+lam2)-K1;
    J(1,7)=-Sc*dlam1; J(1,8)=-Sc*dlam2;

    J(2,1)=p.alpha1;
    J(2,2)=-(lam1+lam2)-K2;
    J(2,7)=-Sp*dlam1; J(2,8)=-Sp*dlam2;

    J(3,2)=p.alpha2; J(3,3)=-K3;

    J(4,1)=lam1+lam2; J(4,2)=lam1+lam2;
    J(4,4)=-K4;
    J(4,7)=(Sc+Sp)*dlam1;
    J(4,8)=(Sc+Sp)*dlam2;

    J(5,4)=p.gamma; J(5,5)=-K5;

    J(6,3)=p.eta; J(6,5)=p.kappa; J(6,6)=-p.mu;

    J(7,5)=p.theta; J(7,7)=-K7;
    J(8,5)=p.psi;   J(8,8)=-K8;
end