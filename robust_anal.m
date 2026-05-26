% Black Pod Disease SEIR Model Robustness Analysis in MATLAB
clc;
clear;

% Parameters (fixed values)
beta2 = 0.05;    % Primary infection rate
theta = 0.8;     % Spore release rate
psi = 0.5;       % Spore shedding rate
alpha1 = 0.05;   % Growth rate from Cherelles to young pods
alpha2 = 0.03;   % Growth rate from young pods to ripe pods
Lambda = 12;     % Recruitment rate of new pods
mu = 0.01;       % Natural death rate
eta = 0.9;       % Harvesting rate of healthy pods
%kappa = 0.8;     % Harvesting rate of infected pods
delta = 0.1;     % Death rate of infected pods
gamma = 0.6;     % Progression rate 

% Initial conditions
Sc0=2000; Sp0=50; Sr0=10; E0=0; I0=0; R0=0; Is0=30; Ip0=30;
%Sc0 = 100; Sp0 = 200; Sr0 = 200; E0 = 50; I0 = 50; R0 = 0; Is0 = 50; Ip0 = 50;
y0 = [Sc0, Sp0, Sr0, E0, I0, R0, Is0, Ip0];  % Initial compartments

% Time vector (160 days)
tspan = 0:1:160;

% Define ranges for beta1 and gamma
beta1_range = linspace(0.01, 1.0, 5);%100);   % Secondary infection rate range
kappa_range = linspace(0.01, 1.0, 5);%100);     % Progression rate range

% Initialize results matrices
peak_infections = zeros(length(beta1_range), length(kappa_range));
time_to_peak = zeros(length(beta1_range), length(kappa_range));
final_size = zeros(length(beta1_range), length(kappa_range));

% Define the SEIR model equations as an anonymous function
bpd_seir_model = @(t, y, beta1, kappa) [
    Lambda - (beta1 * y(7) * (1 + y(7)) + beta2 * y(8) * (1 + y(8))) * y(1) - (alpha1 + mu) * y(1);
    alpha1 * y(1) - (beta1 * y(7) * (1 + y(7)) + beta2 * y(8) * (1 + y(8))) * y(2) - (alpha2 + mu) * y(2);
    alpha2 * y(2) - (mu + eta) * y(3);
    (beta1 * y(7) * (1 + y(7)) + beta2 * y(8) * (1 + y(8))) * (y(1) + y(2)) - (gamma + mu) * y(4);
    gamma * y(4) - (kappa + theta + psi + mu + delta) * y(5);
    eta * y(3) + kappa * y(5) - mu * y(6);
    theta * y(5) - (mu + delta) * y(7);
    psi * y(5) - (mu + delta) * y(8);
];

% Loop over beta1 and gamma ranges
for i = 1:length(beta1_range)
    for j = 1:length(kappa_range)
        beta1 = beta1_range(i);
        kappa = kappa_range(j);
        
        % Solve the ODE system using ode45
        [t, Y] = ode45(@(t, y) bpd_seir_model(t, y, beta1, kappa), tspan, y0);
        
        % Extract compartments
        I = Y(:, 5)+Y(:, 7)+Y(:, 8);  % Infected compartment
        R = Y(:, 6);  % Removed compartment
        
        % Analyze outcomes
        peak_infections(i, j) = max(I);         % Peak infections
        [~, idx] = max(I);                      % Index of peak infection
        time_to_peak(i, j) = t(idx);            % Time to peak infection
        final_size(i, j) = R(end);              % Final size of epidemic
    end
end

% Plotting results
figure;

% Plot Peak Infections
subplot(1, 3, 1);
imagesc(kappa_range, beta1_range, peak_infections);
set(gca, 'YDir', 'normal');
colorbar;
xlabel('kappa');
ylabel('Beta1');
title('Peak Infections');

% Plot Time to Peak Infection
subplot(1, 3, 2);
imagesc(kappa_range, beta1_range, time_to_peak);
set(gca, 'YDir', 'normal');
colorbar;
xlabel('kappa');
ylabel('Beta1');
title('Time to Peak Infection');

% Plot Final Size of Epidemic
subplot(1, 3, 3);
imagesc(kappa_range, beta1_range, final_size);
set(gca, 'YDir', 'normal');
colorbar;
xlabel('kappa');
ylabel('Beta1');
title('Final Size of Epidemic');

% Adjust layout
sgtitle('Robustness Analysis of Black Pod Disease Model');
