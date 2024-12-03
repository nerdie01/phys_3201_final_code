function [time, theta, theta_dot, theta_ddot, thrust] = simulate(dry_mass, tx_shift, M, L, g, rho, c_D_motor, c_D_rod, r, b, R, eff, theta0, theta_dot0, filename)
    % data processing
    data = readtable(filename);
    time = data.t;        
    thrust = data.thrust * eff;    
    mass = data.mass + dry_mass;
    
    [time, thrust, mass] = shift.start(time, thrust, mass, tx_shift);
    
    % space discretization
    N = length(time);
    
    theta = zeros(N,1);
    theta_dot = zeros(N,1);
    theta_ddot = zeros(N,1);
    
    % initial conditions
    theta(1) = theta0;
    theta_dot(1) = theta_dot0;
    
    % precompute drag coefficient
    D_coef = ...
        (1/2) * rho * c_D_motor * pi * r^2 * L^3 + ...
        (1/4) * rho * c_D_rod * R * L^4;

    % forward euler method loop
    for n = 1:N-1
        dt = time(n+1) - time(n);
        dm_dt = (mass(n+1) - mass(n)) / dt;
        
        mg_term = ((1/2)*M + mass(n)) * (g/L) * sin(theta(n));
        I = (1/3)*M + mass(n);
        
        % drag force term
        F_drag = (D_coef * theta_dot(n)^2 + b * theta_dot(n)) / L^2;

        theta_ddot(n) = ( (thrust(n)/L) - dm_dt * theta_dot(n) - mg_term - F_drag ) / I;
        theta_dot(n+1) = theta_dot(n) + theta_ddot(n) * dt;
        theta(n+1) = theta(n) + theta_dot(n) * dt;
    end
        
    theta_ddot(N) = theta_ddot(N-1);
end

clc; clear all;

% constants
dry_mass = 25.0e-3; % dry mass of the rocket (kg)
tx_shift = 0; % shift in theoretical data start time (s)
M = 0.05387; % mass of rod (kg)
L = 0.59; % length of rod (m)
g = 9.81; % gravitational acceleration (m/s^2)
rho = 1.25744; % density of air (~ 1.2 kg/m^3)
c_D_motor = 0.83; % drag coefficient of motor (~0.83)
c_D_rod = 1.17; % drag coefficient of rod (~1.17)
r = 24e-3; % radius of motor (m)
R = 15.5e-3; % radius of rod (m)
b = 0.01; % frictional torque coefficient (k*r m)
eff = 0.7; % thrust efficiency multiplier

% initial conditions
theta0 = 0; % initial angle (rad)
theta_dot0 = 0; % initial angular velocity (rad/s)

% single burn rate coefficient data to focus on in graph
focus_csv = 'data/sixth_nakka_burn.csv';

% experimental data
exp_data = readtable('data/exp_data.csv');
time_exp = exp_data.t;        
theta_exp = exp_data.theta - exp_data.theta(1);

% simple visualization
figure;
hold on;

[time, theta, theta_dot, theta_ddot, thrust] = simulate(dry_mass, tx_shift, M, L, g, rho, c_D_motor, c_D_rod, r, b, R, eff, theta0, theta_dot0, focus_csv);

plot(time, theta, 'DisplayName', 'Theoretical (1/6 burn coefficient)');
plot(time_exp, theta_exp, '-r', 'DisplayName', 'Experimental');
plot(time, thrust, 'DisplayName', 'Thrust');

ylabel('\theta(t) (rad)');
xlabel('Time (s)');
legend('show');
legend('Location', 'best');

grid on;
hold off;

% comparative visualization
csv_file_labels = {'data/half_nakka_burn.csv', 'data/quarter_nakka_burn.csv', ...
             'data/sixth_nakka_burn.csv', 'data/eighth_nakka_burn.csv'};
burn_rate_labels = {'1/2 burn rate', '1/4 burn rate', '1/6 burn rate', '1/8 burn rate'};

figure;
hold on;

plot(time_exp, theta_exp, '-r', 'DisplayName', 'Experimental Data');

for i = 1:length(csv_file_labels)
    filename = csv_file_labels{i};

    [time, theta, ~, ~, thrust] = simulate(dry_mass, tx_shift, M, L, g, rho, ...
                                            c_D_motor, c_D_rod, r, b, R, eff, ...
                                            theta0, theta_dot0, filename);

    plot(time, theta, 'DisplayName', ['Theoretical for ', burn_rate_labels{i}]);
    plot(time, thrust, '--', 'DisplayName', ['Thrust for ', burn_rate_labels{i}]);
end

ylabel('\theta(t) (rad) / Thrust (N)');
xlabel('Time (s)');
legend('show');
legend('Location', 'best');
grid on;
title('Comparison of Thrust over Time for Different Burn Rate Coefficients');
hold off;