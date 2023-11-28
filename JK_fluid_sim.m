%% Testing Section
result = main();

%% Main
function f_obj = main()
    
    close all
    f_obj = simulate_fluid("test");

end

%% Simulation Runner:
% Main Function:
function f_obj = simulate_fluid(preset)
% The primary function used to run the simulation
% Takes:
%   preset: Which preset the fluid object shoulds use!

    % Run Setup:
    f_obj = fluid_obj(preset);
    f_obj = f_obj.spawn_elements();
    [plot_obj, arrow_obj] = create_plot(f_obj,'c', 200);

    % Initialize Frame Recording
    f_obj.frames = getframe();
    frame_count = 2;
    iterations = 0;

    dt_per_frame = floor((1/f_obj.dt) / f_obj.fps);

    % Simulation Loop:
    for t = 0:f_obj.dt:f_obj.sim_time
        % Forward Walk:
        iterations = iterations + 1;
        f_obj = f_obj.forward_walk();

        % Molecule Collisions
        f_obj.Data = run_element_collisions(f_obj,f_obj.Data,f_obj.e_radius);

        % Wall Collisions
        [u_collisions, l_collisions] = detect_wall_interaction(f_obj,f_obj.Data(:,:,1));
        f_obj.Data = run_wall_collisions(f_obj,f_obj.Data,l_collisions,true);
        f_obj.Data = run_wall_collisions(f_obj,f_obj.Data,u_collisions,false);

        % Add Gravity:
        f_obj.Data(2,:,3) = f_obj.Data(2,:,3) - f_obj.g;

        % Record:
        f_obj = f_obj.record_velocities(iterations);

        % Update Plot:
        if mod(t, dt_per_frame * f_obj.dt) == 0
            delete(plot_obj)
            delete(arrow_obj)
            [plot_obj, arrow_obj] = create_plot(f_obj);
            disp("t = " + string(t))
            drawnow
            f_obj.frames(frame_count) = getframe();

            frame_count = frame_count + 1;
        end

    end

end

%% Data Manipulation:
% Setup

% Plotting:
function [plot_obj, arrow_obj] = create_plot(f_obj, mode, resolution)
% Creates the plot that will hold the current state of the model
% Takes:
%   data: Only contains the position data
%   size: The radius of each circle
%   mode: 'c' = create. 'u' = Update
%   resolution: Resolution of the wall line data
% Returns:
%   plot_obj: Contains the handles for all of the elements on the plot
%   arrow_obj: Returns the handles for all of the arrows on the plot
    
    % Initialize:
    data = f_obj.Data(:,:,1);
    hold on

    % Main Code:
    plot_obj = viscircles(data',f_obj.e_radius,Color="blue");

    % Make Arrows for Regions:
    x = permute(f_obj.rec_v_region(1,3,:), [1 3 2]);
    y = permute(f_obj.rec_v_region(2,3,:), [1 3 2]);
    u = permute(f_obj.rec_v_current(1,1,:), [1 3 2]);
    v = permute(f_obj.rec_v_current(2,1,:), [1 3 2]);

    u(isnan(u)) = 0;
    v(isnan(v)) = 0;
    
    arrow_obj = quiver(x,y,u,v,'r');

    if nargin == 3 && mode == 'c'
        % Create Wall Data:
        x_wall = linspace(f_obj.x_axis(1),f_obj.x_axis(2),resolution);
        y_u_wall = f_obj.u_wall(x_wall);
        y_l_wall = f_obj.l_wall(x_wall);
        
        % Create Plot Object
        
        plot(x_wall,y_u_wall,x_wall,y_l_wall); 
    
        % Create Rectangle Regions:
        for depth = 1:size(f_obj.rec_v_region,3)
            
            pos = [f_obj.rec_v_region(:,1,depth)',...
                f_obj.rec_v_region(:,2,depth)' - f_obj.rec_v_region(:,1,depth)'];
            rectangle('Position',pos,'EdgeColor', f_obj.rec_v_color(:,:,depth))

        end
        

        % Plot Size:
        xlim(f_obj.x_axis)
        ylim(f_obj.y_axis)
    end

    hold off

end

%% Wall Collision Detection and Correction:

% Manage Wall Collisions:
function Data = run_wall_collisions(f_obj,Data,collisions,l)
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
    
    % Initialize:
    f_f = f_obj.friction_factor;
    n_f = f_obj.normal_factor;
    size = f_obj.e_radius;
    dt = f_obj.dt;
    
    % Load Data and Correct Position:
    r = correct_element_position_wall(f_obj,Data(:,collisions,1),size,l);
    v = Data(:,collisions,2);
    
    num_collisions = width(r);
    aP = zeros(2,num_collisions);

    % Run Collisions:
    for i = 1:num_collisions
        
        aP(:,i) = wall_collision_force(f_obj,r(:,i),v(:,i),l,f_f,n_f,dt);

    end
    
    Data(:,collisions,1) = r;
    Data(:,collisions,3) = Data(:,collisions,3) + aP;


end

% Detection, Collision, and Correction:
function [u_collisions, l_collisions] = detect_wall_interaction(f_obj,data)
% Determines if the elements are colliding with the wall or not. If they
% are outside the wall this will register as a collision
% Takes:
%   data: Each column is a position of a different element. X positions
%   will be the row 1
%   size: The radius of each element
% Returns:
%   u_collisions: A Row logical vector depicting upper wall collisions
%   l_collisions: A row logical vector depicting lower wall collisions
    
    % Initialize:
    size = f_obj.e_radius;
    
    % Code:
    x = data(1,:);
    y = data(2,:);
    
    u_collisions = y + size >= f_obj.u_wall(x);
    l_collisions = y - size <= f_obj.l_wall(x);


end
function data = correct_element_position_wall(f_obj,data,size,l)
% Corrects an element's position with respect to a wall. This is different
% from it colliding with an element and should be prioritized.
% Takes:
%   data: The position vectors of the elements
%   size: The radius of the element
%   l (bool): The wall interaction is lower
% Returns:
%   data: The corrected positions of the elements

    if l
        data(2,:) = size + f_obj.l_wall(data(1,:)); 
    else
        data(2,:) = f_obj.u_wall(data(1,:)) - size; 
    end

end
function aP = wall_collision_force(f_obj,r, v, l, friction_factor, norm_factor, dt)
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
        T = [dx; f_obj.l_wall(x + dx) - f_obj.l_wall(x)];
        cross_vector = [0; 0; -1]; % Cross product used to determine direction of normal vector
    else
        T = [dx; f_obj.u_wall(x + dx) - f_obj.u_wall(x)];
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
function Data = run_element_collisions(f_obj,Data,size)
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
        F = element_force(f_obj,R,max_R);
        aP = build_accelerations(dR,F);

        % Update Acceleration:
        Data(:,:,3) = Data(:,:,3) + aP;

end

% Detection, Collision, and Correction
function F = element_force(f_obj,R,max_R)
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
    index = R <=2 & R ~= 0; % Ignore anything that is 0
    M = f_obj.e_repulse;
    F(index) = parabola(R(index), M, -2, 2);

    % For greater than than 1
    index = R > 2;
    M = f_obj.e_attract;
    F(index) = parabola(R(index), M, 2, max_R);

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
