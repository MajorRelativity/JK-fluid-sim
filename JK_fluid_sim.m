%% Testing Section
main()


%% Main
function main()

end

%% Simulation Runner:
% Main Function:
function simulate_fluid(dt,num_elements,size)
% The primary function used to run the simulation



end


%% Data Manipulation:
% Setup
function data = spawn_elements(center, num_elements)
% Creates the data matrix that will be used to store and keep track of the
% elements, spawning them at a default position.
% Takes:
%   center (vector): The spawning coordinates of the balls
%   num_elements (int): Gives the number of balls
% Returns:
%   data: a 3D matrix with dimension rows, element number as columns, and
%   position derivatives along the depth (r, v, a)

    % Prealocate Matrix:
    data(:,:,1) = repmat(center,1,num_elements);
    data(:,:,2:3) = zeros(2,num_elements,2);

end

% Forward Walk:
function data = forward_walk(data,dt)
% Takes the data (as defined in spawn_elements) and iterates forward using
% the Euler method over dt:
% Takes:
%   data (3D Matrix): As defined in spawn_elements
%   dt: The interval to walk forward by
% Returns:
%   data (3D Matrix): Updated

% Velocity:
data(:,:,2) = data(:,:,2) + dt * data(:,:,3);

% Position:
data(:,:,1) = data(:,:,1) + dt * data(:,:,2); 

end

%% Wall Collision Detection and Correction:
% Wall Functions:
function y = u_wall(x)
% Defines the boundary you cannot find the element above!
% Takes:
%   x: The x values you desire
% Returns:
%   y: The y value associated with that x value

    y = x + 10;

end
function y = l_wall(x)
% Defines the boundary you cannot find the element below!
% Takes:
%   x: The x values you desire
% Returns:
%   y: The y value associated with that x value

    y = x - 10;

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
function collisions = detect_element_interaction(data, size)
% Returns an upper triangular logical matrix where the row and column
% denote the identity of the two particles, and the value is true if these
% two particles have collides.
% Takes:
%   data: A matrix of 2D position vectors where each column is a 
%       different point
%   size: The size of each molecule
% Returns:
%   collisions: A logical upper triangular matrix of identified collisions

    collisions = boolean(zeros(width(data))); % Logical array
    size_sqr = size^2;
    
    for j = 1:(width(data) - 1) % Do not check the last piece of data
    
        dr = data(:, (j + 1):end ) - data(:,j);
        L_sqr = sum(dr .* dr);
        collisions(j,(j+1):end) = L_sqr <= size_sqr;
    
    end
    
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
    if dt(1) > 0
        dt = dt(1);
    else
        dt = dt(2);
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

    vP1 = v1 - (dot(dv,dr) ./ dot(dr,dr)) .* dr;
    
    % Calculate Force on 1:
    aP1 = (vP1 - v1) ./ dt;

    % Force on 2:
    % This will be in the opposite direction of it is on 1 because of
    % Newton's 3rd law
    aP2 = -1 * aP1;



end