%% Testing Section
%main()
close all
r1 = [0; 0]
r2 = [.5; 0]
v1 = [1; 0]
v2 = [1; 0]
size = .5

circles([r1(1) r2(1)],[r1(2) r2(2)],size)
xlim([-2 2])
ylim([-1.5 1.5])

figure
[r1, r2] = correct_element_position(r1, v1, r2, v2, size)
circles([r1(1) r2(1)],[r1(2) r2(2)],size)
xlim([-2 2])
ylim([-1.5 1.5])
dr = r1 - r2
sqrt(dot(dr,dr))



%% Main
function main()
    
    close
    simulate_fluid(.1,100,10,.4,.50,9.81)

end

%% Simulation Runner:
% Main Function:
function simulate_fluid(dt,sim_time,num_elements,size,v_loss,g)
% The primary function used to run the simulation
% Takes:
%   dt: The time step between each iteration of the simulation
%   num_elements: The number of elements you want in the simulation
%   size: The size of each element in the simulation
%   sim_time: The amount of seconds you want to run the simulation for 

    % Run Setup:
    Data = spawn_elements([0;0],num_elements,3);

    x = Data(1,:,1);
    y = Data(2,:,1);
    plot_obj = create_plot(x, y, 100);

    plot_obj(1).XDataSource = 'x';
    plot_obj(1).YDataSource = 'y';

    % Simulation Loop:
    for t = 0:dt:sim_time
        % Forward Walk:
        Data = forward_walk(Data, num_elements, dt);

        % Molecule Collisions
        collisions = detect_element_interaction(Data(:,:,1),size);
        Data = run_element_collisions(Data,collisions,size,dt);

        % Wall Collisions
        [u_collisions, l_collisions] = detect_wall_interaction(Data(:,:,1),size);
        Data = run_wall_collisions(Data,l_collisions,true,v_loss,size,dt);
        Data = run_wall_collisions(Data,u_collisions,false,v_loss,size,dt);

        % Add Gravity:
        Data(2,:,3) = Data(2,:,3) - g;

        % Update Plot:
        pause(.1)
        x(:) = Data(1,:,1);
        y(:) = Data(2,:,1);
        refreshdata(plot_obj,'caller')
        drawnow

    end

end

%% Data Manipulation:
% Setup
function Data = spawn_elements(center, num_elements, vmax)
% Creates the data matrix that will be used to store and keep track of the
% elements, spawning them at a default position. The velocities are random
% between a given speed interval.
% Takes:
%   center (vector): The spawning coordinates of the balls
%   num_elements (int): Gives the number of balls
%   vmax (int): The Maximum possible initial speed
% Returns:
%   data: a 3D matrix with dimension rows, element number as columns, and
%   position derivatives along the depth (r, v, a)

    % Positions:
    Data(:,:,1) = repmat(center,1,num_elements);

    % Random Velocities:
    Data(:,:,2) = rand([2 num_elements]) .* randi(vmax,[1 num_elements]);

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

% Position:
Data(:,:,1) = Data(:,:,1) + dt * Data(:,:,2); 

% Reset Acceleration:
Data(:,:,3) = zeros(2,num_elements);

end

% Plotting:
function plot_obj = create_plot(x, y , resolution)
% Creates the plot that will hold the current state of the model
% Takes:
%   x: Starting x data
%   y: Starting y data
%   resolution: Resolution of the wall line data

    % Create Wall Data:
    x_wall = linspace(-10,10,resolution);
    y_u_wall = u_wall(x_wall);
    y_l_wall = l_wall(x_wall);
    
    % Create Plot Object
    plot_obj = plot(x,y,'bo',x_wall,y_u_wall,x_wall,y_l_wall); 
    
    % Modify Element Plot Data
    plot_obj(1).MarkerSize = 20;

    % Plot Size:
    xlim([-5 5])
    ylim([-5 5])

end

%% Wall Collision Detection and Correction:

% Manage Wall Collisions:
function Data = run_wall_collisions(Data,collisions,l,v_loss,size,dt)
% Manages the wall collisions for each individual element
% Takes:
%   Data (3D Matrix): Contains all of the important data
%   u_collisions: Row vector indicating which elements collide with the
%       upper wall or lower wall
%   l: True if colliding with lower wall
%   size: Radius of the element
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
        
        aP(:,i) = wall_collision_force(r(:,i),v(:,i),l,v_loss,dt);

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

    y = -1 * x.^2 + 5;

end
function y = l_wall(x)
% Defines the boundary you cannot find the element below!
% Takes:
%   x: The x values you desire
% Returns:
%   y: The y value associated with that x value

    y = x.^2 - 5;

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
function aP = wall_collision_force(r, v, l, v_loss, dt)
% Calculates the applied acceleration on an element interacting with a wall
% Takes:
%   r (vector): The position of the element
%   v (vector): The velocity of the element
%   l (bool): True if the interaction is with a lower wall
%   v_loss (float): The proportion of velocity lost from the wall collision
%       (0 to 1)
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
    
    N = cross([T;0],cross_vector);
    n = N / norm(N);

    n = n(1:2,1); % Return to 2D

    % Acceleration:
    % This can be derived by combining the equation for reflection with the
    % equation for acceleration
    aP = (-2 .* dot(v,n) .* n .* v_loss) ./ dt;
end

%% Molecule Collision Detection and Correction:
% Manage Collisions:
function Data = run_element_collisions(Data,collisions,size,dt)
% Manages the collisions one by one for each pair of elements
% Takes:
%   Data (3D Matrix): Contains all of the position, velocity, and
%       acceleration data
%   collisions (Matrix): Contains two columns with column indices for the
%       Data matrix indicating where there is a collision. For example
%       [1 2] would mean element 1 is colliding with element 2
% Returns:
%   Data: Updated

    for i = 1:height(collisions)
        
        % Assign Values:
        E1 = collisions(i,1);
        E2 = collisions(i,2);

        r1 = Data(:,E1,1);
        r2 = Data(:,E2,1);

        v1 = Data(:,E1,2);
        v2 = Data(:,E2,2);

        % Run Corrections and Collisions
        [r1, r2] = correct_element_position(r1,v1,r2,v2,size);
        [aP1, aP2] = element_collision_force(r1, v1, r2, v2, dt);

        % Update Position and Acceleration:
        Data(:,E1,1) = r1;
        Data(:,E2,1) = r2;

        Data(:,E1,3) = Data(:,E1,3) + aP1;
        Data(:,E2,3) = Data(:,E2,3) + aP2;
    
    end

end

% Detection, Collision, and Correction
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