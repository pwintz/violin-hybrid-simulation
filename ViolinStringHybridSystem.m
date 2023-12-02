classdef ViolinStringHybridSystem < HybridSystem
    % A bouncing ball modeled as a HybridSystem subclass.

    % Define parameters.
    properties
        linear_density = 0.1;
        tension = 100;
        string_length = 10;

        vel_damping_coef = 1;
        
        bow_normal_force = 1000
        bow_speed = 200;
        bow_grid_ndx = 6;  % Index of the bow.
        bow_static_friction_coeff = 1;
        bow_kinetic_friction_coeff = 0.8;
    end

    % Define constant properties that cannot be modified (i.e., "immutable").
    properties(SetAccess = immutable) 
        n_grid;

        % string_pos is the height of the string at each grid point, excluding
        % the end points, which are fixed.
        string_pos_indices;
        string_vel_indices;

        % The state q indicates whether the bow is in a stick mode (q=0) where
        % the interaction the bow and string is governed by static friction or
        % slip mode (q=1) where the interaction is governed by kinetic friction.
        q_index;
    end

    properties(Constant)
        STICK_MODE = 0;
        SLIP_MODE = 1;
    end

    methods 
        function this = ViolinStringHybridSystem(n_grid)
            % Constructor

            % One state component to store the mode, and two for each spatial
            % discretization (position+velocity).            % 
            state_dim = 2*n_grid + 1;
            this = this@HybridSystem(state_dim);

            this.n_grid = n_grid;

            this.string_pos_indices = 1:n_grid;
            this.string_vel_indices = n_grid + this.string_pos_indices;
            this.q_index = state_dim;
        end

        % To define the data of the system, we implement 
        % the abstract functions from HybridSystem.m

        function xdot = flowMap(this, x, t, j)
            % Extract the state components.
            string_pos = x(this.string_pos_indices);
            string_vel = x(this.string_vel_indices);
            q = x(this.q_index);

            % Divide by one less then the number of grid points because the
            % demoninator should be the number of gaps between grid points.
            dx = this.string_length / (this.n_grid - 1);
            
            % string_pos_padded = [0; string_pos; 0];
            string_pos_left = [0; string_pos(1:end-1)];
            string_pos_right = [string_pos(2:end); 0];
            pos_second_derivative ...
                = (string_pos_right - 2*string_pos + string_pos_left) / dx^2;

            % Define the value of f(x). 
            xdot = zeros([this.state_dimension, 1]);
            xdot(this.string_pos_indices) = string_vel;

            % Wave velocity
            % wave_vel = sqrt(this.tension / this.linear_density);
            wave_vel_sq = this.tension / this.linear_density;
            
            xdot(this.string_vel_indices) = wave_vel_sq * pos_second_derivative ...
                                            - this.vel_damping_coef * string_vel;

            if q == this.STICK_MODE
                xdot(this.string_pos_indices(this.bow_grid_ndx)) = this.bow_speed;
                xdot(this.string_vel_indices(this.bow_grid_ndx)) = 0;
            elseif q == this.SLIP_MODE
                % Apply the force of the bow
                bow_ndx_speed = string_vel(this.bow_grid_ndx);
                vel_diff = this.bow_speed - bow_ndx_speed;

                bow_force = sign(vel_diff) * this.bow_normal_force * this.bow_kinetic_friction_coeff;
                string_portion_mass = 1;
                xdot(this.string_vel_indices(this.bow_grid_ndx)) ...
                    = xdot(this.string_vel_indices(this.bow_grid_ndx)) + bow_force / string_portion_mass;

                % xdot(this.string_vel_indices(this.bow_grid_ndx)) = this.bow_speed;
            end

        end

        function xplus = jumpMap(this, x)
            % Define the value of g(x). 
            q = x(this.q_index);
            xplus = x; 
            qplus = 1 - q;
            if qplus == this.STICK_MODE
                % If we are transitioning to stick mode, then set the string
                % speed at the bow index to the bow speed. 
                xplus(this.string_vel_indices(this.bow_grid_ndx)) = this.bow_speed;
            end
            xplus(this.q_index) = qplus;
        end
        
        function inC = flowSetIndicator(this, x) %#ok<INUSD>
            inC = 1;
        end

        function inD = jumpSetIndicator(this, x)
            q = x(this.q_index);
            left_index_pos = x(this.string_pos_indices(this.bow_grid_ndx-1));
            bow_index_pos = x(this.string_pos_indices(this.bow_grid_ndx));
            right_index_pos = x(this.string_pos_indices(this.bow_grid_ndx+1));
            bow_index_vel = x(this.string_vel_indices(this.bow_grid_ndx));
            max_static_friction_force = this.bow_normal_force * this.bow_kinetic_friction_coeff;

            % Vertical component of the tensile force at the bow.
            tensile_force = this.tension *(left_index_pos - bow_index_pos) ...
                + this.tension *(right_index_pos - bow_index_pos);
            if q == this.STICK_MODE
                % Jump to the slip mode if tensile force is stronger than the
                % tensile force. THIS IS NOT QUITE RIGHT BECUASE IT NEGLECTS
                % MOMENTUM.
                inD = abs(tensile_force) >= max_static_friction_force;
            elseif q == this.SLIP_MODE
                % Jump to stick mode if the bow is matching the speed of the
                % string AND the static friction would not be immediately
                % overcome.
                % inD = abs(bow_index_vel - this.bow_speed) < 5e-2*this.bow_speed ...
                inD = bow_index_vel >= this.bow_speed ...
                    && abs(tensile_force) <= 0.9*max_static_friction_force;
            end
        end
    end
end