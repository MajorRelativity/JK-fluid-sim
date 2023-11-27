classdef fluid_obj
    %FLUID_OBJ Contains all of the properties and data related to the fluid
    %simulation.
    %
    % Properties:
    %   preset(str): Keeps track of what preset is used for the model
    %   sim_time (int): The amount of seconds you want to run the simulation for
    %   dt (float): The time step between each iteration of the simulation
    %   Data (3D Matrix): a 3D matrix with dimension rows, element number as columns, and
    %       position derivatives along the depth (r, v, a)
    %
    %   e_num (int): The number of elements you want in the simulation
    %   e_radius (float): The radius of each element
    %   e_repulse (float): The maximum repulsion at the center of each
    %       element
    %   e_attract (float): The minimum attraction (negative) for each particle
    %   spawn_region ([x1 x2;y1 y2]): A matrix of two vectors designating
    %       the bounds of the spawn region
    %
    %   friction_factor (float): Controls the friction between elements and the
    %       wall
    %   normal_factor (float): Controls the elasticity of the collision
    %
    %   
    %   g (float): Gravity constant
    %   x_axis (Row Vector): Beginning and end of plotted x axis
    %   y_axis (Row Vector): same as above but for y
    %
    %   frames (Vector): will contain all of the frames from the simulation
    %   fps (int): number of frames per second in final movie (approx)
    %
    % Units:
    %   t = seconds
    %   d = mm
    properties
        preset
        sim_time
        dt
        Data

        e_num
        e_radius
        e_repulse
        e_attract
        spawn_region
        
        friction_factor
        normal_factor

        g
        x_axis
        y_axis

        frames
        fps
    end
    
    methods
        % Setup:
        function obj = fluid_obj(preset)
            %FLUID_OBJ Construct an instance of this class based on which
            %preset is chosen:
            %   General: Ran in all instances
            %   test: Used to test new parameter combinations
            %   tube: Used to analyze how changing pipe size effects flow
            %       rate
            
            % Presets:
            switch preset
                case "test"
                
                    obj.dt = .0001;
                    obj.sim_time = 3;
    
                    obj.e_num = 100;
                    obj.e_radius = .5;
                    obj.e_repulse = 100000;
                    obj.e_attract = -1000;
                    obj.spawn_region = [-10 10; 15 20];
    
                    obj.friction_factor = .001;
                    obj.normal_factor = .6;
    
                    obj.g = 980;
                    obj.x_axis = [-30 30];
                    obj.y_axis = [-20 40];

                case "tube"
                    obj.dt = .0001;
                    obj.sim_time = 1;
    
                    obj.e_num = 1000;
                    obj.e_radius = .5;
                    obj.e_repulse = 10^5;
                    obj.e_attract = -1000;
                    obj.spawn_region = [-200 0;0 200];
    
                    obj.friction_factor = 0;
                    obj.normal_factor = .6;
    
                    obj.g = 981;
                    obj.x_axis = [-40 40];
                    obj.y_axis = [-20 60];

            end
            
            % General:
            obj.preset = preset;
            obj.fps = 30;

            obj.Data(:,:,1) = zeros(2,obj.e_num); % Position
            obj.Data(:,:,2) = zeros(2,obj.e_num); % Velocity
            obj.Data(:,:,3) = zeros(2,obj.e_num); % Acceleration

        end
        function obj = spawn_elements(obj)
            % Modifies the Data matrix to apply spawn positions based on
            % the parameters defined in fluid_obj
            % Takes:
            %   obj.spawn_center (vector): The spawning coordinates of the balls
            %   obj.e_num (int): Gives the number of balls
            %   obj.elements_wide (int): How many elements wide the
            %   spawning box can be
            
            c = 1;
            c_point = obj.spawn_region(:,1);
            while c <= obj.e_num && c_point(2) <= obj.spawn_region(2,2)
                if c_point(1) > obj.spawn_region(1,2) 
                    % If x is over max spot, reset and increase y
                    c_point(1) = obj.spawn_region(1,1);
                    c_point(2) = c_point(2) + 2 * obj.e_radius; % add to y
                end
                
                cond = obj.u_wall(c_point(1)) > c_point(2) && obj.l_wall(c_point(1)) < c_point(2);
                
                if cond
                    % Place point
                    obj.Data(:,c,1) = c_point;
                    c = c + 1;
                end
                
                c_point(1) = c_point(1) + 2 * obj.e_radius; % add to x

            end
            disp("Spawned " + string(c - 1) + " elements")

        end

        % Modification:
        function obj = forward_walk(obj)
            % Takes the data (as defined in spawn_elements) and iterates forward using
            % the Euler method over dt:
            % Takes:
            %   data (3D Matrix): As defined in spawn_elements
            %   dt: The interval to walk forward by
            % Returns:
            %   data (3D Matrix): Updated
            
            % Velocity:
            obj.Data(:,:,2) = obj.Data(:,:,2) + obj.dt * obj.Data(:,:,3);
            
            % Position:
            obj.Data(:,:,1) = obj.Data(:,:,1) + obj.dt * obj.Data(:,:,2); 
            
            % Reset Acceleration:
            obj.Data(:,:,3) = zeros(2,obj.e_num);
        
        end


        % Wall Functions:
        function y = u_wall(obj,x)
            % Defines the boundary you cannot find the element above!
            % Takes:
            %   x: The x values you desire
            % Returns:
            %   y: The y value associated with that x value
            
            % Initialize
            y = zeros(1,length(x));
            
            % Upper Tube:
            switch obj.preset
                case "test"
                    y = 35 .* ones(1,length(x));
                case "tube"
                    cond = x < 0;
                    y(cond) = (1/5) * x(cond).^2 + 4;
                
                    cond = x >= 0;
                    y(cond) = 4;
            end
        
        end
        function y = l_wall(obj,x)
            % Defines the boundary you cannot find the element below!
            % Takes:
            %   x: The x values you desire
            % Returns:
            %   y: The y value associated with that x value
            
            % Initialize:
            y = zeros(1,length(x));
            
            % Basic:
            switch obj.preset
                case "test"
                    y = (x.^2) ./ 30;

                case "tube"
                    cond = x < 0;
                    y(cond) = (1/15) * x(cond).^2 + 2;
                
                    cond = x < 15 & x >= 0;
                    y(cond) = 2;
                
                    cond = x < 25 & x >= 15;
                    y(cond) = -(1/5) * (x(cond) - 15) + 2;

                    cond = x >= 25;
                    y(cond) = 0;
            end
        
        end

    end
end

