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
    %   spawn_center ([x;y]): A vector showing where the center of the
    %       elements spawn should be
    %   elements_wide (int): How many elements wide the spawn rectangle
    %       should be
    %
    %   friction_factor (float): Controls the friction between elements and the
    %       wall
    %   normal_factor (float): Controls the elasticity of the collision
    %
    %   
    %   g (float): Gravity constant
    %   
    %   
    properties
        preset
        sim_time
        dt
        Data

        e_num
        e_radius
        spawn_center
        elements_wide
        
        friction_factor
        normal_factor

        g
    end
    
    methods
        % Setup:
        function obj = fluid_obj(preset)
            %FLUID_OBJ Construct an instance of this class based on which
            %preset is chosen:
            %   General: Ran in all instances
            %   test: Used to test new parameter combinations
            
            % Presets:
            if preset == "test"
                
                obj.dt = .0001;
                obj.sim_time = 100;

                obj.e_num = 100;
                obj.e_radius = .5;
                obj.spawn_center = [-20;40];
                obj.elements_wide = 20;

                obj.friction_factor = .001;
                obj.normal_factor = .6;

                obj.g = 400;

            end
            
            % General:
            obj.preset = preset;

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
        
            % Positions:
        
            length_x = obj.elements_wide * 2 * obj.e_radius;
            c_point = obj.spawn_center - length_x/2; % Set first point
            left_x = c_point(1);
            right_x = obj.spawn_center(1) + length_x/2;
            
            c = 1; % Counter Variable
            while c <= obj.e_num
                if c_point(1) >= right_x
                    c_point(2) = c_point(2) + 2 * obj.e_radius;
                    c_point(1) = left_x;
                end
        
                obj.Data(:,c,1) = c_point;
                c_point(1) = c_point(1) + 2 * obj.e_radius;
                c = c + 1;
            end

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
        function y = u_wall(~,x)
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
        function y = l_wall(~,x)
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

    end
end

