function [Min_impulse, time_min] = MIN_2IMP (initial_position, initial_velocity)

    % This function gives out the minimum total impulse and the time of the 
    % randezvous in that case, for a Two-impulse linear rendezvous maneuver, among all
    % the total impulse needed to perform the rendezvous with different times 
    % between 1 second and half orbital period
    %
    % Syntax 
    % [Min_impulse, time_min] = MIN_2IMP (initial_position, initial_velocity)
    %
    % Input Arguments
    % initial_position --> (1,3) vector with initial relative position of
    % the deputy wrt target in [m]
    % initial_velocity --> (1,3) vector with initial relative velocity of
    % the deputy wrt target in [m/s]
    %
    % Output Arguments
    % Min_impulse --> A scalar value [m/s] representing the minimum total impulse
    % Time_min --> A scalar value [s] representing the time necessary to perform the Min_impulse

    initial_conditions.p = initial_position;
    initial_conditions.v = initial_velocity;
    
    h = 590e3; % [m] altitude of the target
    R = 6378.137e3; % [m] equatorial radius of Earth
    mu = 3.986004418e14; % [m^3/s^2] gravitational costant
    orbital_radius = R+h; % [m] orbital radius of the target
    n = sqrt(mu/(orbital_radius)^3); % [rad/s] mean motion(angular rate) of target
    tau =(2*pi)/n; % [s] orbital period of the chief
    
    A = [zeros(3) eye(3);3*n^2 0 0 0 2*n 0;0 0 0 -2*n 0 0;0 0 -n^2 0 0 0]; % state Matrix
    phi=@(t) expm(A*t); % state transition matrix function
    
    N = 0.5; % number of orbital period we want to visualize
    time = linspace(1,N*tau,10000); % time discretization to visualize the orbit of the deputy [s]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rendezvous velocity and Impulses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % preallocation of the matrix we will use
    rendezvous_conditions.v = zeros(3,length(time)); % relative velocity vector of deputy with respect to the target
    impulse_1.vector = zeros(3,length(time)); % first relative impulse(delta_v) vector we give to the deputy to reach the target
    impulse_2.vector = zeros(3,length(time)); % second relative impulse(delta_v) vector we give to the deputy to stop it once arrived to the target
    for i=1:length(time)
        M = phi(time(i)); %state transition matrix given the time
    
        rendezvous_conditions.v(:,i) = M(1:3,4:6)\([0 0 0]' - M(1:3,1:3)*initial_conditions.p');
    
        % First relative impulse as the different between the needed rendezvous
        % velocity and the velocity of the deputy at the time of we give the impulse
        impulse_1.vector(:,i) = rendezvous_conditions.v(:,i)-initial_conditions.v';
        impulse_1.magnitude(i) = norm(impulse_1.vector(:,i));
    
        % Second relative impulse as the opposite of the velocity vector once reached
        % the target to stop it and cancel its relative velocity wrt the target
        rndv_initial_conditions = [initial_conditions.p'; rendezvous_conditions.v(:,i)]; % conditons at the start of the randezvous
        impulse_2.vector(:,i) = -M(4:6,:)*rndv_initial_conditions;
        impulse_2.magnitude(i) = norm(impulse_2.vector(:,i));
    
        % Calcolo l'impulso totale
        impulse_tot(i) = impulse_1.magnitude(i) + impulse_2.magnitude(i);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % We find the minimun impulse among all the randezvous times we have computed
    [Min_impulse, time_min_ind] = min(impulse_tot)
    time_min = time(time_min_ind)

end
