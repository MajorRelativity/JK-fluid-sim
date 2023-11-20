%% Testing Section
main()

%% Main
function main()
    
    close all
    simulate_fluid(.0001,100,1000,.5,.001,.6,400)

end

%% Simulation Runner:
% Main Function:
function simulate_fluid(dt,sim_time,num_elements,size,f_f,n_f,g)
% The primary function used to run the simulation
% Takes:
%   dt: The time step between each iteration of the simulation
%   num_elements: The number of elements you want in the simulation
%   size: The size of each element in the simulation
%   sim_time: The amount of seconds you want to run the simulation for 

    % Run Setup:
    Data = spawn_elements([-20;40],size, num_elements, 20);

    plot_obj = create_plot(Data(:,:,1),size,'c', 200);

    % Simulation Loop:
    for t = 0:dt:sim_time
        % Forward Walk:
        Data = forward_walk(Data, num_elements, dt);

        % Molecule Collisions
        Data = run_element_collisions(Data,size);

        % Wall Collisions
        [u_collisions, l_collisions] = detect_wall_interaction(Data(:,:,1),size);
        Data = run_wall_collisions(Data,l_collisions,true,f_f,n_f,size,dt);
        Data = run_wall_collisions(Data,u_collisions,false,f_f,n_f,size,dt);

        % Add Gravity:
        Data(2,:,3) = Data(2,:,3) - g;

        % Update Plot:
        if mod(t, 25 * dt) == 0
            delete(plot_obj)
            plot_obj = create_plot(Data(:,:,1),size);
            drawnow
        end

    end

end

%% Data Manipulation:
% Setup
function Data = spawn_elements(center, size, num_elements, elements_wide)
% Creates the data matrix that will be used to store and keep track of the
% elements, spawning them at a default position. The velocities are random
% between a given speed interval.
% Takes:
%   center (vector): The spawning coordinates of the balls
%   num_elements (int): Gives the number of balls
%   elements_wide (int): How many elements wide the spawning box can be
% Returns:
%   data: a 3D matrix with dimension rows, element number as columns, and
%   position derivatives along the depth (r, v, a)

    % Positions:
    Data(:,:,1) = zeros(2,num_elements);

    length_x = elements_wide * 2 * size;
    c_point = center - length_x/2; % Set first point
    left_x = c_point(1);
    right_x = center(1) + length_x/2;
    
    c = 1; % Counter Variable
    while c <= num_elements
        if c_point(1) >= right_x
            c_point(2) = c_point(2) + 2 * size;
            c_point(1) = left_x;
        end

        Data(:,c,1) = c_point;
        c_point(1) = c_point(1) + 2 * size;
        c = c + 1;
    end

    % Random Velocities:
    Data(:,:,2) = zeros(2,num_elements);

    % Accelerations:
    Data(:,:,3) = zeros(2,num_elements);

end

% Forward Walk:
function Data = forward_walk(Data,num_elements,dt)
% Takes the data (as defined in spawn_elements) and iterates forward using
% the Euler method over dt:
% Takes:
%   data (3D Matrix): As defined in spawn_elements
%   dt: The interval to walk forward by
% Returns:
%   data (3D Matrix): Updated

% Velocity:
Data(:,:,2) = Data(:,:,2) + dt * Data(:,:,3);
%Data(:,:,2)
%Data(:,:,3)

% Position:
Data(:,:,1) = Data(:,:,1) + dt * Data(:,:,2); 

% Reset Acceleration:
Data(:,:,3) = zeros(2,num_elements);

end

% Plotting:
function plot_obj = create_plot(data, size, mode, resolution)
% Creates the plot that will hold the current state of the model
% Takes:
%   data: Only contains the position data
%   size: The radius of each circle
%   mode: 'c' = create. 'u' = Update
%   resolution: Resolution of the wall line data
% Returns:
%   plot_obj: Contains the handles for all of the elements on the plot
    
    plot_obj = viscircles(data',size,Color="blue");
    if nargin == 4 && mode == 'c'
        % Create Wall Data:
        x_wall = linspace(-40,20,resolution);
        y_u_wall = u_wall(x_wall);
        y_l_wall = l_wall(x_wall);
        
        % Create Plot Object
        hold on
        plot(x_wall,y_u_wall,x_wall,y_l_wall); 
        hold off
    
        % Plot Size:
        xlim([-30 30])
        ylim([-30 60])
    end

end

%% Wall Collision Detection and Correction:

% Manage Wall Collisions:
function Data = run_wall_collisions(Data,collisions,l,f_f,n_f,size,dt)
% Manages the wall collisions for each individual element
% Takes:
%   Data (3D Matrix): Contains all of the important data
%   collisions: Row vector indicating which elements collide with the
%       upper wall or lower wall
%   l: True if colliding with lower wall
%   size: Radius of the element
%   f_f: The friction factor
%   n_f: The proportion applied to the NORMAL acceleration
%   dt: Time step
% Returns:
%   Data: Updated
    
    % Load Data and Correct Position:
    r = correct_element_position_wall(Data(:,collisions,1),size,l);
    v = Data(:,collisions,2);
    
    num_collisions = width(r);
    aP = zeros(2,num_collisions);

    % Run Collisions:
    for i = 1:num_collisions
        
        aP(:,i) = wall_collision_force(r(:,i),v(:,i),l,f_f,n_f,dt);

    end
    
    Data(:,collisions,1) = r;
    Data(:,collisions,3) = Data(:,collisions,3) + aP;


end

% Wall Functions:
function y = u_wall(x)
% Defines the boundary you cannot find the element above!
% Takes:
%   x: The x values you desire
% Returns:
%   y: The y value associated with that x value
    
    % Initialize
    y = zeros(1,length(x));

    % Basic:
    %y = 50 .* ones(1,length(x));
    
    % Upper Tube:
    cond = x < -3;
    y(cond) = -5 * (x(cond) + 3) + 11;

    cond = x < 5 & x >= -3;
    y(cond) = -2 * x(cond) + 5;

    cond = x >= 5;
    y(cond) = -5;

end
function y = l_wall(x)
% Defines the boundary you cannot find the element below!
% Takes:
%   x: The x values you desire
% Returns:
%   y: The y value associated with that x value
    
    % Initialize:
    y = zeros(1,length(x));
    
    % Basic:
    %y = (x.^2) ./ 30 - 20;

    % Lower Tube:
    cond = x < -20;
    y(cond) = (1/8) * (x(cond) + 20).^2 + 29;

    cond = x < 4 & x >= -20;
    y(cond) = -(3/2) * x(cond) - 1;

    cond = x < 10 & x >= 4;
    y(cond) = -7;

    cond = x >= 10;
    y(cond) = -9;

end

% Detection, Collision, and Correction:
function [u_collisions, l_collisions] = detect_wall_interaction(data,size)
% Determines if the elements are colliding with the wall or not. If they
% are outside the wall this will register as a collision
% Takes:
%   data: Each column is a position of a different element. X positions
%   will be the row 1
%   size: The radius of each element
% Returns:
%   u_collisions: A Row logical vector depicting upper wall collisions
%   l_collisions: A row logical vector depicting lower wall collisions

    x = data(1,:);
    y = data(2,:);
    
    u_collisions = y + size >= u_wall(x);
    l_collisions = y - size <= l_wall(x);


end
function data = correct_element_position_wall(data,size,l)
% Corrects an element's position with respect to a wall. This is different
% from it colliding with an element and should be prioritized.
% Takes:
%   data: The position vectors of the elements
%   size: The radius of the element
%   l (bool): The wall interaction is lower
% Returns:
%   data: The corrected positions of the elements

    if l
        data(2,:) = size + l_wall(data(1,:)); 
    else
        data(2,:) = u_wall(data(1,:)) - size; 
    end

end
function aP = wall_collision_force(r, v, l, friction_factor, norm_factor, dt)
% Calculates the applied acceleration on an element interacting with a wall
% Takes:
%   r (vector): The position of the element
%   v (vector): The velocity of the element
%   l (bool): True if the interaction is with a lower wall
%   friction_factor (float): The proportion of vertical acceleration
%       applied to the collision as friction.
%   norm_factor: The scale applied to the normal acceleration (should be .5 to 1)
%   dt: The time that the collision happens over
% Returns:
%   aP (vector): Applied acceleration to the element

    % Wall Normal Vector:
    dx = .01; % Used to calculate the tangent
    x = r(1);
    
    if l
        T = [dx; l_wall(x + dx) - l_wall(x)];
        cross_vector = [0; 0; -1]; % Cross product used to determine direction of normal vector
    else
        T = [dx; u_wall(x + dx) - u_wall(x)];
        cross_vector = [0; 0; 1];
    end
    
    N = cross([T;0],cross_vector); % Normal Vector
    n = N / norm(N); % Normal Vector normalized

    n = n(1:2,1); % Return to 2D


    % Calculate Friction:
    norm_T = T / norm(T);
    projection_coeff = -1 * dot(v,norm_T); % We want project to face opposite the velocity
    
    f = projection_coeff * friction_factor * norm_T ./ dt;

    % Normal Acceleration Magnitude:
    a_magnitude = -2 * dot(v,n) * norm_factor / dt;

    % Full Acceleration:
    % This can be derived by combining the equation for reflection with the
    % equation for acceleration
    % Includes Friction
    aP = ( a_magnitude .* n) + f;
end

%% Molecule Collision Detection and Correction:
% Manage Collisions:
function Data = run_element_collisions(Data,size)
% Manages the collisions one by one for each pair of elements
% Takes:
%   Data (3D Matrix): Contains all of the position, velocity, and
%       acceleration data
% Returns:
%   Data: Updated

        
        % Assign Values:
        data = Data(:,:,1);
        max_R = 4;

        % Run Corrections and Collisions
        [R, dR] = calculate_distances(data, size, max_R);
        F = element_force(R,max_R);
        aP = build_accelerations(dR,F);

        % Update Acceleration:
        Data(:,:,3) = Data(:,:,3) + aP;

end

% Detection, Collision, and Correction
function F = element_force(R,max_R)
% Determines the force on an element based on the distance from the center
% of that element.
% Takes:
%   R: The magnitude of the distance between elements. In the form of an
%   upper traingular matrix
% Returns:
%   F: The magnitude of the force between elements with respect to the
%   element generating the force. Positive is away, negative is towards!
%   Shound mirror the format of R
    
    F = zeros(length(R));
    
    % For less than than 1:
    index = R <= 1 & R ~= 0; % Ignore anything that is 0
    M = 10000;
    F(index) = parabola(R(index), M, -1, 1);

    % For greater than than 1
    index = R > 1;
    M = -100;
    F(index) = parabola(R(index), M, 1, max_R);

    function f = parabola(r, extrema, root1, root2)
        % Defines a parabola based on extrema and roots
        r_avg = (root1 + root2) / 2;
        A = extrema / ((r_avg - root1).*(r_avg - root2));
        f = A.*(r - root1).*(r - root2);

    end


end
function [R, dR] = calculate_distances(data, size, max_R)
% Returns a upper triangular matrix containing the distances between
% particles in the row and column
% Takes:
%   data: A matrix of 2D position vectors where each column is a 
%       different point
%   size: The size of each molecule
%   max_radii_away: The maximum distance in radii that particles can be
%       from each other
% Returns:
%   R: An upper triangular matrix of distances
%   dR: A normalized 3D upper traingular matrix. The first layer in depth
%   is x and the second is y
    
    w = width(data);
    R = zeros(w);
    dR = zeros(w, w, 2);
    
    
    for j = 1:(width(data) - 1) % Do not check the last piece of data
    
        dr = data(:, (j + 1):end ) - data(:,j);
        r = sqrt(sum(dr .* dr));
        r(r == 0) = 1; % To prevent divide by 0 issue
        R(j,(j+1):end) = r;
        
        dR(j,(j + 1):end, 1) = dr(1,:) ./ r; % Save Normalized Versions
        dR(j,(j + 1):end, 2) = dr(2,:) ./ r; 

    end

    % Turn R into distance in radii. Change all radii more tha max_radii to
    % 0 as this will represent no interaction.

    R = R / size;
    R(R >= max_R) = 0;

    
end
function a = build_accelerations(dR,F)
% Takes the forces and the unit 3D matrix of unit vectors and sums them to
% create the accelerations for each particle!
% Takes:
%   F: A 2D matrix of all of the force magnitudes between element #row and
%   element #column
%   dR: The unit vectors corrosponding to the above where the depth is
%   differenet dimensions of position
% Returns:
%   a: A matrix where each column is a different vector acceleration for an
%   element.

    A = F .* dR; % Turn all of the unit vectors into full force vectors
    
    % To get the full force on a element n, you must take the sum of column
    % n and row n and add them together. This accounts for the forces on
    % all particles.
    
    a = sum(A,1) + permute(sum(-1 .* A,2),[2,1,3]);
    a = permute(a, [3 2 1]); % Switches the depth and the rows1

end


%% Defunct Functions
function [rP1, rP2] = correct_element_position(r1, v1, r2, v2, size)
% If the particles significantly overlap, then this function will move the
% particles backwards in time until they no longer overlap.
% Takes:
%   r1, r2: The position vectors
%   v1, v2: The velocity vectors
%   size: The radius of the particle
% Returns:
%   rP1, rP2: The corrected position vectors of the particles

    % Calculate Necessary Variables:
    dr = r1 - r2;
    dv = v1 - v2;
    R_sqr = size^2;
    
    % dt: Can be solved for using a quadratic equation
    dt = roots([dot(dv,dv),-2 * dot(dr,dv), dot(dr,dr) - (4 * R_sqr)]);
    
    % Choose the positive root:
    condition = dt > 0;
    if sum(condition) == 2
        dt = dt(1);
    elseif sum(condition) == 1
        dt = dt(condition);
    else
        dt = 0;
    end
    
    % Make Correction:
    rP1 = r1 - (dt * v1);
    rP2 = r2 - (dt * v2);

end
function [aP1, aP2] = element_collision_force(r1, v1, r2, v2, dt)
% Only should be used on verified collisions. This function takes in the
% vector velocities of the two colliding elements and returns the
% accelerations that should be applied by this collision to the element!
% Assumes that the masses of the two elements are the same
% Takes:
%   r1, r2: The positions of the two colliding elements
%   v1, v2: The velocities of the two colliding elements
%   dt: The time over which the collision should happen
% Returns:
%   aP1, aP2: The added accelerations to the elements.

    % Calculate Final velocity (only one is necessary):
    dv = v1 - v2; % Note direction (1 - 2)
    dr = r1 - r2;
    
    dr_mag = dot(dr,dr);

    if dr_mag ~= 0
        % This means no "NaN" issues
        vP1 = v1 - (dot(dv,dr) ./ dr_mag) .* dr;
    else
        vP1 = v1;
    end
    
    % Calculate Force on 1:
    aP1 = (vP1 - v1) ./ dt;

    % Force on 2:
    % This will be in the opposite direction of it is on 1 because of
    % Newton's 3rd law
    aP2 = -1 * aP1;



end
function collisions = detect_element_interaction(data, size)
% Returns a two column matrix where each value in the row
% denote the identity of the two particles that are colliding.
% Takes:
%   data: A matrix of 2D position vectors where each column is a 
%       different point
%   size: The size of each molecule
% Returns:
%   collisions: A two column matrix with indices in each of the columns

    collisions = boolean(zeros(width(data))); % Logical array
    size_sqr = size^2;
    
    for j = 1:(width(data) - 1) % Do not check the last piece of data
    
        dr = data(:, (j + 1):end ) - data(:,j);
        L_sqr = sum(dr .* dr);
        collisions(j,(j+1):end) = L_sqr <= size_sqr;
    
    end

    [row,col] = find(collisions);
    collisions = [row col];
    
end