function [time, theta, theta_dot, theta_ddot, thrust] = simulate(dry_mass, tx_shift, M, L, g, theta0, theta_dot0, filename)
    % data processing
    data = readtable(filename);
    time = data.t;        
    thrust = data.thrust;    
    mass = data.mass + dry_mass;

    [time, thrust, mass] = shift.start(time, thrust, mass, tx_shift);

    % time discretization
    N = length(time);

    theta = zeros(N,1);
    theta_dot = zeros(N,1);
    theta_ddot = zeros(N,1);

    % apply initial conditions
    theta(1) = theta0;
    theta_dot(1) = theta_dot0;

    % forward euler method loop
    for n = 1:N-1
        dt = time(n+1) - time(n);
        dm_dt = (mass(n+1) - mass(n)) / dt;

        mg_term = ((1/2)*M + mass(n)) * (g/L) * sin(theta(n));
        I = (1/3)*M + mass(n);

        theta_ddot(n) = ((thrust(n)/L) - dm_dt*theta_dot(n) - mg_term ) / I;
        theta_dot(n+1) = theta_dot(n) + theta_ddot(n) * dt;
        theta(n+1) = theta(n) + theta_dot(n) * dt;
    end

    theta_ddot(N) = theta_ddot(N-1);
end

clc; clear all;

% constants
dry_mass = 15.5 / 1000; % dry mass of the rocket (kg)
tx_shift = 0; % shift in theoretical data start time (s)
M = 0.05387; % mass of rod (kg)
L = 0.59; % length of rod (m)
g = 9.81; % gravitational acceleration (m/s^2)
theta0 = 0; % initial angle
theta_dot0 = 0; % initial angular velocity

% import openmotor data
filename = 'data/sixth_nakka_burn.csv';
[time, theta, theta_dot, theta_ddot, thrust] = simulate(dry_mass, tx_shift, M, L, g, theta0, theta_dot0, filename);

% import experimental data
exp_data = readtable('data/exp_data.csv');
time_exp = exp_data.t;        
theta_exp = exp_data.theta - exp_data.theta(1); 

% visualization
figure;
hold on;

plot(time, theta, 'DisplayName', 'Theoretical');
plot(time_exp, theta_exp, '-r', 'DisplayName', 'Experimental', 'Color', 'Red');
plot(time, thrust, 'DisplayName', 'Thrust (N)');

ylabel('\theta (rad)');
xlabel('Time (s)');
legend('show');
legend('Location', 'best');
grid on;
title('Frictionless No-Drag Predicted Angle over Time');
% set(gca, 'YScale', 'log');
hold off;
