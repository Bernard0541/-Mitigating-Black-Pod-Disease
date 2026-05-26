%% =========================================================
%% Optimal control for the biologically consistent BPD model
%% Forward-Backward Sweep Method with switchable strategies
%% =========================================================
clear; clc; close all;

%% ===================== CHOOSE STRATEGY =====================
% 'A' = u1 and u2
% 'B' = u1 and u3
% 'C' = u2 and u3
% 'D' = u1, u2, and u3
strategy = 'A';

switch upper(strategy)
    case 'A'
        active_u1 = 1; active_u2 = 1; active_u3 = 0;
    case 'B'
        active_u1 = 1; active_u2 = 0; active_u3 = 1;
    case 'C'
        active_u1 = 0; active_u2 = 1; active_u3 = 1;
    case 'D'
        active_u1 = 1; active_u2 = 1; active_u3 = 1;
    otherwise
        error('Unknown strategy. Choose A, B, C, or D.');
end

fprintf('Running Strategy %s\n', upper(strategy));
fprintf('Active controls: u1=%d, u2=%d, u3=%d\n', active_u1, active_u2, active_u3);

%% ===================== PARAMETERS =====================
p.Lambda   = 12;
p.alpha1   = 0.05;
p.alpha2   = 0.027;
p.beta1    = 0.05;
p.beta2    = 0.05;
p.pi1      = 1.0;
p.pi2      = 1.0;
p.eta      = 0.9;
p.gamma    = 0.6;
p.kappa    = 0.8;
p.theta    = 0.8;
p.psi      = 0.5;
p.mu       = 0.01;
p.delta    = 0.1;

% control-related parameters
p.eps2     = 0.8;   % efficacy of treatment in reducing shedding
p.rho      = 0.7;   % efficacy of rouging on infected pods
p.xi       = 0.5;   % efficacy of rouging on inoculum classes

% objective weights
p.A1       = 300;   % weight on exposed pods
p.A2       = 500;   % weight on infected pods

% control costs
p.B1       = 30;
p.B2       = 30;
p.B3       = 30;

%% ===================== TIME GRID =====================
T  = 100;
N  = 1200;
t  = linspace(0,T,N+1);

%% ===================== INITIAL CONDITIONS =====================
% X = [Sc Sp Sr E I R Is Ip]
X0 = [80; 60; 20; 8; 5; 0; 3; 2];

%% ===================== INITIAL CONTROL GUESS =====================
u1 = 0.2*ones(size(t)) * active_u1;
u2 = 0.2*ones(size(t)) * active_u2;
u3 = 0.2*ones(size(t)) * active_u3;

%% ===================== NUMERICAL SETTINGS =====================
maxIter = 200;
tol     = 1e-5;
relax   = 0.4;

state_opts = odeset('RelTol',1e-7,'AbsTol',1e-9,'NonNegative',1:8);
adj_opts   = odeset('RelTol',1e-7,'AbsTol',1e-9);

Jhist = nan(maxIter,1);

fprintf('Starting forward-backward sweep...\n');

for iter = 1:maxIter

    u1_old = u1;
    u2_old = u2;
    u3_old = u3;

    %% -------------------------------------------------
    %% 1. FORWARD SOLVE FOR STATES
    %% -------------------------------------------------
    udata.t  = t;
    udata.u1 = u1;
    udata.u2 = u2;
    udata.u3 = u3;

    state_fun = @(tt,xx) state_rhs(tt,xx,p,udata);
    [ts, Xsol] = ode15s(state_fun, t, X0, state_opts);

    X = interp1(ts, Xsol, t, 'pchip')';
    X = max(X,0);

    if any(~isfinite(X(:))) || max(X(:)) > 1e10
        error('Forward state solve became unstable. Check parameters.');
    end

    %% -------------------------------------------------
    %% 2. BACKWARD SOLVE FOR ADJOINTS
    %% -------------------------------------------------
    L_T = zeros(8,1);

    tau = t;
    adj_tau_fun = @(tau,ll) -adjoint_rhs(T - tau, ll, p, t, X, udata);

    [tau_sol, Ltau] = ode15s(adj_tau_fun, tau, L_T, adj_opts);

    Ltau_interp = interp1(tau_sol, Ltau, tau, 'linear')';
    L = fliplr(Ltau_interp);

    %% -------------------------------------------------
    %% 3. CONTROL UPDATE
    %% -------------------------------------------------
    u1_new = zeros(size(t));
    u2_new = zeros(size(t));
    u3_new = zeros(size(t));

    for k = 1:length(t)
        Sc = X(1,k); Sp = X(2,k);
        I  = X(5,k); Is = X(7,k); Ip = X(8,k);

        lamSc = L(1,k); lamSp = L(2,k); lamE  = L(4,k);
        lamI  = L(5,k); lamR  = L(6,k);
        lamIs = L(7,k); lamIp = L(8,k);

        lambda1 = p.beta1 * Is * (1 + p.pi1*Is);
        lambda2 = p.beta2 * Ip * (1 + p.pi2*Ip);
        Phi = lambda1 + lambda2;

        % u1*
        u1_star = ( ...
            Phi*((Sc+Sp)*lamE - Sc*lamSc - Sp*lamSp) + ...
            p.theta*(1 - p.eps2*u2(k))*I*lamIs + ...
            p.psi  *(1 - p.eps2*u2(k))*I*lamIp ) / p.B1;

        % u2*
        u2_star = ( ...
            p.eps2*I*((1-u1(k))*p.theta*lamIs + (1-u1(k))*p.psi*lamIp ...
            - (p.theta + p.psi)*lamI) ) / p.B2;

        % u3*
        u3_star = ( ...
            p.rho*I*(lamI - lamR) + p.xi*Is*lamIs + p.xi*Ip*lamIp ) / p.B3;

        % turn off controls according to strategy
        u1_star = active_u1 * u1_star;
        u2_star = active_u2 * u2_star;
        u3_star = active_u3 * u3_star;

        % project to [0,1]
        u1_new(k) = min(1,max(0,u1_star));
        u2_new(k) = min(1,max(0,u2_star));
        u3_new(k) = min(1,max(0,u3_star));
    end

    % optional smoothing
    if active_u1, u1_new = smoothdata(u1_new,'movmean',11); end
    if active_u2, u2_new = smoothdata(u2_new,'movmean',11); end
    if active_u3, u3_new = smoothdata(u3_new,'movmean',11); end

    u1_new = min(1,max(0,u1_new)) * active_u1;
    u2_new = min(1,max(0,u2_new)) * active_u2;
    u3_new = min(1,max(0,u3_new)) * active_u3;

    % relaxed update
    u1 = relax*u1_old + (1-relax)*u1_new;
    u2 = relax*u2_old + (1-relax)*u2_new;
    u3 = relax*u3_old + (1-relax)*u3_new;

    % enforce off-switch exactly
    u1 = active_u1 * u1;
    u2 = active_u2 * u2;
    u3 = active_u3 * u3;

    %% -------------------------------------------------
    %% 4. OBJECTIVE FUNCTION
    %% -------------------------------------------------
    integrand = p.A1*X(4,:) + p.A2*X(5,:) + 0.5*(p.B1*u1.^2 + p.B2*u2.^2 + p.B3*u3.^2);
    Jhist(iter) = trapz(t, integrand);

    err = max([norm(u1-u1_old,inf), norm(u2-u2_old,inf), norm(u3-u3_old,inf)]);

    fprintf('Iter %3d: J = %.8e, error = %.3e\n', iter, Jhist(iter), err);

    if err < tol
        Jhist = Jhist(1:iter);
        fprintf('Converged in %d iterations.\n', iter);
        break;
    end

    if iter == maxIter
        Jhist = Jhist(1:iter);
        fprintf('Maximum iterations reached.\n');
    end
end

%% ===================== UNCONTROLLED COMPARISON =====================
udata0.t  = t;
udata0.u1 = zeros(size(t));
udata0.u2 = zeros(size(t));
udata0.u3 = zeros(size(t));

state_fun0 = @(tt,xx) state_rhs(tt,xx,p,udata0);
[ts0, Xsol0] = ode15s(state_fun0, t, X0, state_opts);

Xun = interp1(ts0, Xsol0, t, 'pchip')';
Xun = max(Xun,0);

%% ===================== PLOTS =====================

figure('Color','w');
plot(t, Xun(4,:), 'k--', 'LineWidth', 4); hold on;
plot(t, X(4,:),   'b-',  'LineWidth', 4);
xlabel('Time'); ylabel('Exposed cocoa pods');
legend('E uncontrolled','E controlled','Location','best');
%title(['Strategy ', upper(strategy), ': Controlled vs Uncontrolled Exposed Pods']);
xlim([0 50])
set(gca, 'LineWidth', 2)
grid on;

figure('Color','w');
plot(t, Xun(5,:), 'r--', 'LineWidth', 4); hold on;
plot(t, X(5,:),   'g-',  'LineWidth', 4);
xlabel('Time'); 
ylabel('Infected cocoa pods');
legend('I uncontrolled','I controlled','Location','best');
%title(['Strategy ', upper(strategy), ': Controlled vs Uncontrolled Infected Pods']);
xlim([0 50])
set(gca, 'LineWidth', 2)
grid on;

figure('Color','w');
plot(t, u1, 'r--', 'LineWidth', 3); 
hold on;
plot(t, u2, 'b--', 'LineWidth', 4);
%plot(t, u3, 'k--', 'LineWidth', 4);
xlabel('Time'); ylabel('Control profile');
%'$u_1$',
legend('$u_1$','$u_2$','Interpreter','latex','Location','best');
%title(['Strategy ', upper(strategy), ': Optimal Controls']);
ylim([0 1.05]);
set(gca, 'LineWidth', 2)
box on
grid on;

%% ===================== SAVE RESULTS =====================
results.strategy = upper(strategy);
results.t   = t;
results.X   = X;
results.Xun = Xun;
results.L   = L;
results.u1  = u1;
results.u2  = u2;
results.u3  = u3;
results.J   = Jhist(end);

disp('Final objective value:')
disp(results.J)

%% =========================================================
%% LOCAL FUNCTIONS
%% =========================================================

function dx = state_rhs(tt, x, p, udata)
    u1 = interp1(udata.t, udata.u1, tt, 'linear', 'extrap');
    u2 = interp1(udata.t, udata.u2, tt, 'linear', 'extrap');
    u3 = interp1(udata.t, udata.u3, tt, 'linear', 'extrap');

    Sc = x(1); Sp = x(2); Sr = x(3);
    E  = x(4); I  = x(5); R  = x(6);
    Is = x(7); Ip = x(8);

    lambda1 = p.beta1 * Is * (1 + p.pi1*Is);
    lambda2 = p.beta2 * Ip * (1 + p.pi2*Ip);
    Phi = lambda1 + lambda2;

    dx = zeros(8,1);

    dx(1) = p.Lambda - (1-u1)*Phi*Sc - (p.alpha1 + p.mu)*Sc;
    dx(2) = p.alpha1*Sc - (1-u1)*Phi*Sp - (p.alpha2 + p.mu)*Sp;
    dx(3) = p.alpha2*Sp - (p.mu + p.eta)*Sr;
    dx(4) = (1-u1)*Phi*(Sc + Sp) - (p.gamma + p.mu)*E;
    dx(5) = p.gamma*E - (p.kappa + p.mu + p.delta + ...
             p.theta*(1 - p.eps2*u2) + p.psi*(1 - p.eps2*u2) + p.rho*u3)*I;
    dx(6) = p.eta*Sr + p.kappa*I + p.rho*u3*I - p.mu*R;
    dx(7) = (1-u1)*p.theta*(1 - p.eps2*u2)*I - (p.mu + p.delta + p.xi*u3)*Is;
    dx(8) = (1-u1)*p.psi *(1 - p.eps2*u2)*I - (p.mu + p.delta + p.xi*u3)*Ip;
end

function dl = adjoint_rhs(tt, l, p, tgrid, X, udata)
    x = interp1(tgrid, X', tt, 'linear', 'extrap')';

    u1 = interp1(udata.t, udata.u1, tt, 'linear', 'extrap');
    u2 = interp1(udata.t, udata.u2, tt, 'linear', 'extrap');
    u3 = interp1(udata.t, udata.u3, tt, 'linear', 'extrap');

    Sc = x(1); Sp = x(2); Sr = x(3);
    Is = x(7); Ip = x(8);

    lSc = l(1); lSp = l(2); lSr = l(3); lE = l(4);
    lI  = l(5); lR  = l(6); lIs = l(7); lIp = l(8);

    lambda1 = p.beta1 * Is * (1 + p.pi1*Is);
    lambda2 = p.beta2 * Ip * (1 + p.pi2*Ip);
    Phi = lambda1 + lambda2;

    dlam1 = p.beta1 * (1 + 2*p.pi1*Is);
    dlam2 = p.beta2 * (1 + 2*p.pi2*Ip);

    dl = zeros(8,1);

    dl(1) = (1-u1)*Phi*(lSc - lE) + (p.alpha1 + p.mu)*lSc - p.alpha1*lSp;
    dl(2) = (1-u1)*Phi*(lSp - lE) + (p.alpha2 + p.mu)*lSp - p.alpha2*lSr;
    dl(3) = (p.mu + p.eta)*lSr - p.eta*lR;
    dl(4) = -p.A1 + (p.gamma + p.mu)*lE - p.gamma*lI;
    dl(5) = -p.A2 ...
          + (p.kappa + p.mu + p.delta + p.theta*(1-p.eps2*u2) + p.psi*(1-p.eps2*u2) + p.rho*u3)*lI ...
          - (p.kappa + p.rho*u3)*lR ...
          - (1-u1)*p.theta*(1-p.eps2*u2)*lIs ...
          - (1-u1)*p.psi *(1-p.eps2*u2)*lIp;
    dl(6) = p.mu*lR;
    dl(7) = (1-u1)*dlam1*(Sc*(lSc-lE) + Sp*(lSp-lE)) + (p.mu + p.delta + p.xi*u3)*lIs;
    dl(8) = (1-u1)*dlam2*(Sc*(lSc-lE) + Sp*(lSp-lE)) + (p.mu + p.delta + p.xi*u3)*lIp;
end