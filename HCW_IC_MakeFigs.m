clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 1.2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% FIGURE 1: x0 VARIATION
clear all
clc

h = 590e3; % [m] altitude of the target 
R = 6378.137e3; % [m] equatorial radius of Earth 
mu = 3.986004418e14; % [m^3/s^2] gravitational costant 
orbital_radius = R+h; % [m] orbital radius of the target 
n = sqrt(mu/(orbital_radius)^3); % [rad/s] mean motion(angular rate) of target 
tau =(2*pi)/n; % [s] orbital period of the chief 

A = [zeros(3) eye(3);3*n^2 0 0 0 2*n 0;0 0 0 -2*n 0 0;0 0 -n^2 0 0 0]; % state Matrix
phi=@(t) expm(A*t); % state transition matrix function

N = 3; % number of orbital period we want to visualize
time = linspace(0,N*tau,1000); % [s] time discretization to visualize the orbit of the deputy 

% initial relative state of deputy with respect to the target
initial_conditions.x = [-500 -100 -50 -10 0 10 50 100 500]; % [m] possible values of initial x-component position
initial_conditions.y = 0; % [m]
initial_conditions.z = 0; % [m]
initial_conditions.xdot = 0; % [m/s]
initial_conditions.ydot = 0; % [m/s]
initial_conditions.zdot = 0; % [m/s]

initial_conditions_tot = [0; initial_conditions.y; initial_conditions.z; initial_conditions.xdot; initial_conditions.ydot; initial_conditions.zdot]; % initial conditions as a 6-elements column matrix

actual_conditions = zeros(6,length(time)); % preallocation of the matrix with the condition of the deputy at i-th time in the i-th column

for j=1:length(initial_conditions.x)
    for i=1:length(time)
        initial_conditions_tot(1) = initial_conditions.x(j);
        actual_conditions(:,i) = phi(time(i))*initial_conditions_tot;
    end
    figure(1)
    grafico_1 = plot3(actual_conditions(2,:),actual_conditions(1,:),actual_conditions(3,:), 'LineWidth', 2, 'DisplayName', ['$x_{0}=', num2str(initial_conditions.x(j)), ' [m]$']);
    axis([-3000 3000 -2000 2000 -1000 1000]);
    grid on;
    hold on;
    view(2)
    xlabel('$Y [m]$', 'FontSize', 20, 'Interpreter','latex');
    ylabel('$X [m]$', 'FontSize', 20, 'Interpreter','latex');
    title('Spontaneous motion of deputy relative to target');
end
k = legend(gca, 'show'); 
set(k, 'Interpreter', 'latex');

%% FIGURE2: DOT_x0 VARIATION
clear all
clc

h = 590e3; % [m] altitude of the target 
R = 6378.137e3; % [m] equatorial radius of Earth 
mu = 3.986004418e14; % [m^3/s^2] gravitational costant 
orbital_radius = R+h; % [m] orbital radius of the target 
n = sqrt(mu/(orbital_radius)^3); % [rad/s] mean motion(angular rate) of target 
tau =(2*pi)/n; % [s] orbital period of the chief 

A = [zeros(3) eye(3);3*n^2 0 0 0 2*n 0;0 0 0 -2*n 0 0;0 0 -n^2 0 0 0]; % state Matrix
phi=@(t) expm(A*t); % state transition matrix function

N = 3; % number of orbital period we want to visualize
time = linspace(0,N*tau,1000); % time discretization to visualize the orbit of the deputy [s]

% initial relative state of deputy with respect to the target
initial_conditions.x = 0; % [m]
initial_conditions.y = 0; % [m]
initial_conditions.z = 0; % [m]
initial_conditions.xdot = [-50 -10 -5 -1 -0.5 0 0.5 1 5 10 50]; % [m/s] possible values of initial x-component of velocity
initial_conditions.ydot = 0; % [m/s]
initial_conditions.zdot = 0; % [m/s]

initial_conditions_tot = [initial_conditions.x; initial_conditions.y; initial_conditions.z; 0; initial_conditions.ydot; initial_conditions.zdot]; % initial conditions as a 6-elements column matrix

actual_conditions = zeros(6,length(time)); % preallocation of the matrix with the condition of the deputy at i-th time in the i-th column

for j=1:length(initial_conditions.xdot)
    for i=1:length(time)
        initial_conditions_tot(4) = initial_conditions.xdot(j);
        actual_conditions(:,i) = phi(time(i))*initial_conditions_tot;
    end
    figure(2)
    grafico_1 = plot3(actual_conditions(2,:),actual_conditions(1,:),actual_conditions(3,:), 'LineWidth', 2, 'DisplayName', ['$\dot x_{0}=', num2str(initial_conditions.xdot(j)), ' [m/s]$']);
    axis([-3000 3000 -2000 2000 -1000 1000]);
    grid on;
    hold on;
    view(2)
    xlabel('$Y [m]$', 'FontSize', 20, 'Interpreter','latex');
    ylabel('$X [m]$', 'FontSize', 20, 'Interpreter','latex');
    title('Spontaneous motion of deputy relative to target');
end
k = legend(gca, 'show'); 
set(k, 'Interpreter', 'latex');

%% FIGURE 3: DOT_y0 VARIATION
clear all
clc

h = 590e3; % [m] altitude of the target 
R = 6378.137e3; % [m] equatorial radius of Earth 
mu = 3.986004418e14; % [m^3/s^2] gravitational costant 
orbital_radius = R+h; % [m] orbital radius of the target 
n = sqrt(mu/(orbital_radius)^3); % [rad/s] mean motion(angular rate) of target 
tau =(2*pi)/n; % [s] orbital period of the chief 

A = [zeros(3) eye(3);3*n^2 0 0 0 2*n 0;0 0 0 -2*n 0 0;0 0 -n^2 0 0 0]; % state Matrix
phi=@(t) expm(A*t); % state transition matrix function

N = 3; % number of orbital period we want to visualize
time = linspace(0,N*tau,1000); % time discretization to visualize the orbit of the deputy [s]

% initial relative state of deputy with respect to the target
initial_conditions.x = 0; % [m]
initial_conditions.y = 0; % [m]
initial_conditions.z = 0; % [m]
initial_conditions.xdot = 0; % [m/s]
initial_conditions.ydot = [-50 -10 -5 -1 -0.5 -0.25 -0.113 0 0.113 0.25 0.5 1 5 10 50]; % [m/s] possible values of initial y-component of velocity
initial_conditions.zdot = 0; % [m/s]

initial_conditions_tot = [initial_conditions.x; initial_conditions.y; initial_conditions.z; initial_conditions.xdot; 0 ; initial_conditions.zdot]; % initial conditions as a 6-elements column matrix

actual_conditions = zeros(6,length(time)); % preallocation of the matrix with the condition of the deputy at i-th time in the i-th column

for j=1:length(initial_conditions.ydot)
    for i=1:length(time)
        initial_conditions_tot(5) = initial_conditions.ydot(j);
        actual_conditions(:,i) = phi(time(i))*initial_conditions_tot;
    end
    figure(3)
    grafico_1 = plot3(actual_conditions(2,:),actual_conditions(1,:),actual_conditions(3,:), 'LineWidth', 2, 'DisplayName', ['$\dot y_{0}=', num2str(initial_conditions.ydot(j)), ' [m/s]$']);
    axis([-3000 3000 -2000 2000 -1000 1000]);
    grid on;
    hold on;
    view(2)
    xlabel('$Y [m]$', 'FontSize', 20, 'Interpreter','latex');
    ylabel('$X [m]$', 'FontSize', 20, 'Interpreter','latex');
    title('Spontaneous motion of deputy relative to target');
end
k = legend(gca, 'show'); 
set(k, 'Interpreter', 'latex');

%% FIGURE 4: x0 AND DOT_x0 VARIATIONS
clear all
clc

h = 590e3; % [m] altitude of the target 
R = 6378.137e3; % [m] equatorial radius of Earth 
mu = 3.986004418e14; % [m^3/s^2] gravitational costant 
orbital_radius = R+h; % [m] orbital radius of the target 
n = sqrt(mu/(orbital_radius)^3); % [rad/s] mean motion(angular rate) of target 
tau =(2*pi)/n; % [s] orbital period of the chief 

A = [zeros(3) eye(3);3*n^2 0 0 0 2*n 0;0 0 0 -2*n 0 0;0 0 -n^2 0 0 0]; % state Matrix
phi=@(t) expm(A*t); % state transition matrix function

N = 1; % number of orbital period we want to visualize
time = linspace(0,N*tau,1000); % time discretization to visualize the orbit of the deputy [s]

% initial relative state of deputy with respect to the target
initial_conditions.x = input("Enter initial x-position of Deputy with respect to the target:"); % [m] initial x-position decided by the user as a keyboard input
initial_conditions.y = 0; % [m]
initial_conditions.z = 0; % [m]
initial_conditions.xdot = [-50 -10 -5 -1 -0.5 0 0.5 1 5 10 50]; % [m/s] possible values of initial x-component of velocity
initial_conditions.ydot = 0; % [m/s]
initial_conditions.zdot = 0; % [m/s]

initial_conditions_tot = [initial_conditions.x; initial_conditions.y; initial_conditions.z; 0; initial_conditions.ydot; initial_conditions.zdot]; % initial conditions as a 6-elements column matrix

actual_conditions = zeros(6,length(time)); % preallocation of the matrix with the condition of the deputy at i-th time in the i-th column

for j=1:length(initial_conditions.xdot)
    for i=1:length(time)
        initial_conditions_tot(4) = initial_conditions.xdot(j);
        actual_conditions(:,i) = phi(time(i))*initial_conditions_tot;
    end
    figure(4)
    grafico_1 = plot3(actual_conditions(2,:),actual_conditions(1,:),actual_conditions(3,:), 'LineWidth', 2, 'DisplayName', ['$\dot x_{0}=', num2str(initial_conditions.xdot(j)), ' [m/s]$']);
    axis([-3000 3000 -2000 2000 -1000 1000]);
    grid on;
    hold on;
    view(2)
    xlabel('$Y [m]$', 'FontSize', 20, 'Interpreter','latex');
    ylabel('$X [m]$', 'FontSize', 20, 'Interpreter','latex');
    title('Spontaneous motion of deputy relative to target');
end
k = legend(gca, 'show'); 
set(k, 'Interpreter', 'latex');

%% FIGURE 5 x0 AND DOT_y0 VARIATION
clear all
clc

h = 590e3; % [m] altitude of the target 
R = 6378.137e3; % [m] equatorial radius of Earth 
mu = 3.986004418e14; % [m^3/s^2] gravitational costant 
orbital_radius = R+h; % [m] orbital radius of the target 
n = sqrt(mu/(orbital_radius)^3); % [rad/s] mean motion(angular rate) of target 
tau =(2*pi)/n; % [s] orbital period of the chief 

A = [zeros(3) eye(3);3*n^2 0 0 0 2*n 0;0 0 0 -2*n 0 0;0 0 -n^2 0 0 0]; % state Matrix
phi=@(t) expm(A*t); % state transition matrix function

N = 3; % number of orbital period we want to visualize
time = linspace(0,N*tau,1000); % time discretization to visualize the orbit of the deputy [s]

% initial relative state of deputy with respect to the target
initial_conditions.x = input("Enter initial x-position of Deputy with respect to the target:"); % [m] initial x-position decided by the user as a keyboard input
initial_conditions.y = 0; % [m]
initial_conditions.z = 0; % [m]
initial_conditions.xdot = 0; % [m/s]
initial_conditions.ydot = [-24 -20 -2.5 -1.3 -0.109 -0.176 0.058 1.05 2.2 12 24]; % [m/s] possible values of initial y-component of velocity
initial_conditions.zdot = 0; % [m/s]

initial_conditions_tot = [initial_conditions.x; initial_conditions.y; initial_conditions.z; initial_conditions.xdot; 0 ; initial_conditions.zdot]; % initial conditions as a 6-elements column matrix

actual_conditions = zeros(6,length(time)); % preallocation of the matrix with the condition of the deputy at i-th time in the i-th column

for j=1:length(initial_conditions.ydot)
    for i=1:length(time)
        initial_conditions_tot(5) = initial_conditions.ydot(j);
        actual_conditions(:,i) = phi(time(i))*initial_conditions_tot;
    end
    figure(5)
    grafico_1 = plot3(actual_conditions(2,:),actual_conditions(1,:),actual_conditions(3,:), 'LineWidth', 2, 'DisplayName', ['$\dot y_{0}=', num2str(initial_conditions.ydot(j)), ' [m/s]$']);
    axis([-3000 3000 -2000 2000 -1000 1000]);
    grid on;
    hold on;
    view(2)
    xlabel('$Y [m]$', 'FontSize', 20, 'Interpreter','latex');
    ylabel('$X [m]$', 'FontSize', 20, 'Interpreter','latex');
    title('Spontaneous motion of deputy relative to target');
end
k = legend(gca, 'show'); 
set(k, 'Interpreter', 'latex');


%% FIGURE 6 trajectory on z
clear all
clc

h = 590e3; % [m] altitude of the target 
R = 6378.137e3; % [m] equatorial radius of Earth 
mu = 3.986004418e14; % [m^3/s^2] gravitational costant 
orbital_radius = R+h; % [m] orbital radius of the target 
n = sqrt(mu/(orbital_radius)^3); % [rad/s] mean motion(angular rate) of target 
tau =(2*pi)/n; % [s] orbital period of the chief 

A = [zeros(3) eye(3);3*n^2 0 0 0 2*n 0;0 0 0 -2*n 0 0;0 0 -n^2 0 0 0]; % state Matrix
phi=@(t) expm(A*t); % state transition matrix function

N = 3; % number of orbital period we want to visualize
time = linspace(0,N*tau,1000); % time discretization to visualize the orbit of the deputy [s]

% initial relative state of deputy with respect to the target
initial_conditions.x = input("Enter initial x-position of Deputy with respect to the target:"); % [m]
initial_conditions.y = input("Enter initial y-position of Deputy with respect to the target:"); % [m]
initial_conditions.z = input("Enter initial z-position of Deputy with respect to the target:"); % [m]
initial_conditions.xdot = input("Enter initial x-velocity of Deputy with respect to the target:"); % [m/s]
initial_conditions.ydot = input("Enter initial y-velocity of Deputy with respect to the target:"); % [m/s]
initial_conditions.zdot = input("Enter initial z-velocity of Deputy with respect to the target:"); % [m/s]

initial_conditions_tot = [initial_conditions.x; initial_conditions.y; initial_conditions.z; initial_conditions.xdot; initial_conditions.ydot ; initial_conditions.zdot]; % initial conditions as a 6-elements column matrix

final_conditions = zeros(6,length(time)); % preallocation of the matrix with the condition of the deputy at i-th time in the i-th column

for i=1:length(time)
    final_conditions(:,i) = phi(time(i))*initial_conditions_tot;
end
figure(6)
grafico_1 = plot3(final_conditions(2,:),final_conditions(1,:),final_conditions(3,:), 'LineWidth', 2);
grid on;
hold on;
view(3)
xlabel('$Y [m]$', 'FontSize', 20, 'Interpreter','latex');
ylabel('$X [m]$', 'FontSize', 20, 'Interpreter','latex');
% k = legend(gca, 'show'); 
% set(k, 'Interpreter', 'latex');


