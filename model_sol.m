clearvars; clc
% Parameters
p1  = 12;     % Lambda
p2  = 0.05;   % alpha1
p3  = 0.027;  % alpha_2
p4  = 0.05;   % beta_1
p5  = 0.05;   % beta_2
p6  = 1;      % pi_1
p7  = 1;      % pi_2
p8  = 0.9;   % eta
p9  = 0.6;    % gamma
p10 = 0.8;    % kappa
p11 = 0.8;   % theta
p12 = 0.5;   % psi
p13 = 0.01;   % mu
p14 = 0.1; % delta

% Initial conditions
Sc0=2000; Sp0=100; Sr0=30; E0=10; I0=10; R0=0; Is0=5; Ip0=5;

Sc = p1/(p2 + p13);
Sp = (p1*p2)/((p3 + p13)*(p2 + p13));

R_nought = (p9*(p12*(p13+p14)*Sc*p5+p12*(p13+p14)*Sp*p5+p11*(p13+p14)*Sc*p4+p11*(p13+p14)*Sp*p5))...
    /((p9+p13)*(p10+p11+p12+p13+p14)*(p13+p14));
disp(['R_0 = ',num2str(R_nought)])


% ODE model
odemodel = @(t, y) [p1 - (p4*y(7)*(1 + p6*y(8)) + p5*y(8)*(1 + p7*y(8)))*y(1) - (p2 + p13)*y(1);
    p2*y(1) - (p4*y(7)*(1 + p6*y(8)) + p5*y(8)*(1 + p7*y(8)))*y(2) - (p3 + p13)*y(2);
    p3*y(2) - (p13 + p8)*y(3);
    (p4*y(7)*(1 + p6*y(8)) + p5*y(8)*(1 + p7*y(8)))*(y(1)+y(2)) - (p9 + p13)*y(4);
    p9*y(4) - (p10 + p11 + p12 + p13 + p14)*y(5);
    p8*y(3) + p10*y(5) - p13*y(6);
    p11*y(5) - (p13 + p14)*y(7);
    p12*y(5) - (p13 + p14)*y(8)];

% Initial condition vector
y0 = [Sc0 Sp0 Sr0 E0 I0 R0 Is0 Ip0];

% Time span
tspan = [0 365];

% Solve the differential equations using ode45
[tsol,xsol] = ode45(odemodel, tspan, y0);

% Plot the results
figure(1);
%plot(tsol,xsol(:,1),'b','linewidth',3);
hold on
plot(tsol,xsol(:,2),'k','linewidth',3);
%plot(tsol,xsol(:,3),'g','linewidth',3);
plot(tsol,xsol(:,4),'c','linewidth',3);
plot(tsol,xsol(:,5),'m','linewidth',3);
%plot(tsol,xsol(:,6),'y','linewidth',3);
%plot(tsol,xsol(:,7),'r--','linewidth',3);
%plot(tsol,xsol(:,8),'r-','linewidth',3);
%legend('S_c', 'S_p', 'S_r', 'E', 'I', 'R', 'I_s', 'I_p')
legend('S_p', 'E', 'I')
xlabel('Time (Days)');
ylabel('Number of pods');
xlim([0, 35]);
box on
set(gca, 'LineWidth', 2)
hold off;

% figure(2)
% subplot(2, 4, 1);
% hold on
% plot(tsol, xsol(:,1), 'b', 'LineWidth', 2);
% legend('S_c', 'Location', 'best');
% xlabel('Time (Days)');
% ylabel('Number of pods');
% xlim([0, 35]);
% box on
% set(gca, 'LineWidth', 2)
% subplot(2, 4, 2);
% plot(tsol, xsol(:,2), 'r', 'LineWidth', 2);
% legend('S_p', 'Location', 'best');
% xlabel('Time (Days)');
% ylabel('Number of pods');
% xlim([0, 35]);
% box on
% set(gca, 'LineWidth', 2)
% subplot(2, 4, 3);
% plot(tsol, xsol(:,3), 'r--', 'LineWidth', 2);
% legend('S_r', 'Location', 'best');
% xlabel('Time (Days)');
% ylabel('Number of pods');
% xlim([0, 35]);
% box on
% set(gca, 'LineWidth', 2)
% subplot(2, 4, 4);
% plot(tsol, xsol(:,4), 'c', 'LineWidth', 2);
% legend('E', 'Location', 'best');
% xlabel('Time (Days)');
% ylabel('Number of pods');
% xlim([0, 35]);
% box on
% set(gca, 'LineWidth', 2)
% subplot(2, 4, 5);
% plot(tsol, xsol(:,5), 'm', 'LineWidth', 2);
% legend('I', 'Location', 'best');
% xlabel('Time (Days)');
% ylabel('Number of pods');
% xlim([0, 35]);
% box on
% set(gca, 'LineWidth', 2)
% subplot(2, 4, 6);
% plot(tsol, xsol(:,6), 'y', 'LineWidth', 2);
% legend('R', 'Location', 'best');
% xlabel('Time (Days)');
% ylabel('Number of pods');
% xlim([0, 35]);
% box on
% set(gca, 'LineWidth', 2)
% subplot(2, 4, 7);
% plot(tsol, xsol(:,7), 'k', 'LineWidth', 2);
% legend('I_s', 'Location', 'best');
% xlabel('Time (Days)');
% ylabel('Number of pods');
% xlim([0, 35]);
% box on
% set(gca, 'LineWidth', 2)
% subplot(2, 4, 8);
% plot(tsol, xsol(:,8), 'g', 'LineWidth', 2);
% legend('I_p', 'Location', 'best');
% xlabel('Time (Days)');
% ylabel('Number of individuals');
% xlim([0, 35]);
% box on
% set(gca, 'LineWidth', 2)
% hold off;
