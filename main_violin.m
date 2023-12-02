
n_grid = 40;
violin = ViolinStringHybridSystem(n_grid);

grid = linspace(0, violin.string_length, n_grid+2)';
grid_interior = grid(2:end-1);
h0 = sin(pi*grid_interior / violin.string_length);
v0 = 0*h0;
q0 = ViolinStringHybridSystem.STICK_MODE;
x0 = [h0; v0; q0];

tspan = [0, 10];
jspan = [0, 1000];
sol = violin.solve(x0, tspan, jspan)
t_grid = linspace(0, sol.t(end), 400);
sol = sol.interpolateToHybridArc(t_grid);

figure(1);
for i = 1:size(sol.x, 1)
    clf
    ymax_in_sol = max(max(abs(sol.x(:, violin.string_pos_indices))));
    ylim(1.5*ymax_in_sol*[-1, 1])
    hold on
    plot(grid, [0; sol.x(i, violin.string_pos_indices)'; 0])
    q = sol.x(i, violin.q_index);

    switch q
        case violin.STICK_MODE 
            color = 'red';
        case violin.SLIP_MODE
            color = 'blue';
    end

    plot(grid(1+violin.bow_grid_ndx), sol.x(i, violin.string_pos_indices(violin.bow_grid_ndx))', '*', 'Color', color)
    plot(grid, 0.01*[0; sol.x(i, violin.string_vel_indices)'; 0])
    title(sprintf('Time: %.2f seconds, q=%d', sol.t(i), q))
    drawnow
    pause(0.1)
end
disp('End')
% plot(grid, [0; sol.xf(violin.string_pos_indices); 0], 'r')
