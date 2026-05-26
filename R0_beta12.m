% Parameters
p1  = 12;     % Lambda
p2  = 0.05;   % alpha_1
p3  = 0.027;  % alpha_2
%p4  = 0.05;   % beta_1
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
Sc0=1500; Sp0=50; Sr0=10; E0=0; I0=0; R0=0; Is0=30; Ip0=30;

p4_range = linspace(0.0, 1.0, 100);
p5_range = linspace(0.0, 1.0, 100);

%R_nought = zeros(size(p4_range));
R_nought = zeros(length(p4_range),length(p5_range));

for i = 1:length(p4_range)
    for j = 1:length(p5_range)
        Sc = p1/(p2 + p13);
        Sp = (p1*p2)/((p3 + p13)*(p2 + p13));
        p4 = p4_range(i);
        p5 = p5_range(j);

        R_nought(i,j) = (p9*(p12*(p13+p14)*Sc*p5+p12*(p13+p14)*Sp*p5+p11*(p13+p14)*Sc*p4+p11*(p13+p14)*Sp*p5))...
            /((p9+p13)*(p10+p11+p12+p13+p14)*(p13+p14));
        %disp(['R_0 = ',num2str(R_nought)])
    end
end

%figure(1)
%plot(p4_range, R_nought, 'Color', 'k', 'LineWidth', 4)
%hold on
%xlabel('\beta_1');
%ylabel('$\mathcal{R}_0$');
%box on
%set(gca, 'LineWidth', 2)
%hold off

% Contour plot
cont=figure;
figure(2);
[C, h] = contourf(p4_range, p4_range, R_nought);
clabel(C, h);
xlabel('\beta_1');
ylabel('\beta_2');
box on
set(gca, 'LineWidth', 4)
colormap(jet);
colorbar;


% 3D Surface Plot
[Beta1, Beta2] = meshgrid(p4_range, p5_range);

figure(3);
surf(Beta1, Beta2, R_nought);
xlabel('\beta_1');
ylabel('\beta_2');
zlabel('$\mathcal{R}_0$');
set(gca, 'LineWidth', 4)
box on
colormap(jet); % Color map
colorbar; % Color bar
grid off
view(0,90)