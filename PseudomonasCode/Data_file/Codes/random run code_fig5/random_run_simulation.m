clear; close all; clc;

total_time     = 8;  dt = 0.001;      n_steps = total_time/dt;
inflow_rate    = 0;    speed = 18;    runtime_lambda = 0.9;
x_max = 800;   y_max = 800;

C_min =  0; C_max = 100;
D_r = 0.03;           % rad²/s
D_t = 0.4;            % µm²/s
L_cell = 1.5;         % µm
L_tail=5               %µm

particles_x = {}; particles_y = {};
particles_theta = {}; particles_runtime = {};
particles_active = {};

current_run_x = {}; current_run_y = {};
runs_x = {}; runs_y = {};

for k = 1:50
    px = rand()*x_max;  py = rand()*y_max;  th = 2*pi*rand();
    particles_x{end+1} = px;
    particles_y{end+1} = py;
    particles_theta{end+1} = th;
    particles_runtime{end+1} = exprnd(1/runtime_lambda);
    particles_active{end+1} = true;
    current_run_x{end+1} = px;
    current_run_y{end+1} = py;
end

figure; hold on; grid on; title('Bacterial motion');
xlabel('X (µm)'); ylabel('Y (µm)');
xlim([-50 x_max+50]); ylim([-10 y_max+10]);
plot([0 x_max],[0 0],'k', [0 x_max],[y_max y_max],'k', [x_max x_max],[0 y_max],'k','LineWidth',2)
plot(0,35,'g>','MarkerFaceColor','g','MarkerSize',10);
hPts = scatter([],[],'o','filled'); drawnow;

for step = 1:n_steps
    for j = 1:numel(particles_x)
        if ~particles_active{j}
            continue;
        end

        x = particles_x{j}(end);
        y = particles_y{j}(end);
        th = particles_theta{j};
        rt = particles_runtime{j} - dt;

        if rt <= 0
            runs_x{end+1} = current_run_x{j};
            runs_y{end+1} = current_run_y{j};
            current_run_x{j} = [];
            current_run_y{j} = [];
            particles_active{j} = false;
            continue;
        end

        C = C_min + (C_max-C_min)*(x/x_max);
        v_dp = 1; % diff velocity
        omega = (-v_dp/L_cell)*sin(th);
        th = th + omega*dt + sqrt(2 * D_r * dt) * randn();
        x = x + speed*cos(th)*dt + sqrt(2 * D_t * dt) * randn() + v_dp*dt;
        y = y + speed*sin(th)*dt + sqrt(2 * D_t * dt) * randn();

        particles_x{j} = [particles_x{j} x];
        particles_y{j} = [particles_y{j} y];
        particles_theta{j} = mod(th, 2*pi);
        particles_runtime{j} = rt;

        current_run_x{j}(end+1) = x;
        current_run_y{j}(end+1) = y;
    end

    set(hPts,'XData',cellfun(@(z)z(end),particles_x), ...
             'YData',cellfun(@(z)z(end),particles_y));
    drawnow;
end

for j = 1:numel(current_run_x)
    if ~isempty(current_run_x{j})
        runs_x{end+1} = current_run_x{j};
        runs_y{end+1} = current_run_y{j};
    end
end

nRuns = numel(runs_x);
straightness = zeros(nRuns,1);
for k = 1:nRuns
    x = runs_x{k}; y = runs_y{k};
    dx = x(end) - x(1); dy = y(end) - y(1);
    displacement = sqrt(dx^2 + dy^2);
    path_length = sum(sqrt(diff(x).^2 + diff(y).^2));
    straightness(k) = displacement / max(path_length, eps);
end

straightness_norm = (straightness);

cmap = cool(32);%color map
colors = cmap(round(1 + 31 * straightness_norm), :);

figure; hold on; grid on; axis equal;
xlabel('\Deltax (µm)'); ylabel('\Deltay (µm)');

for k = 1:nRuns
    dx = runs_x{k} - runs_x{k}(1);
    dy = runs_y{k} - runs_y{k}(1);
    z = zeros(size(dx));
    c = straightness(k) * ones(size(dx));  % Use raw straightness values directly

    surface([dx; dx], [dy; dy], [z; z], [c; c], ...
        'FaceColor', 'none', ...
        'EdgeColor', 'interp', ...
        'LineWidth', 2);
end

colormap(cmap);
caxis([min(straightness) max(straightness)]);  % Color limits set by raw data range

cb = colorbar;
cb.Label.String = 'Straightness';

% Ticks dynamically spaced in your data range
cb.Ticks = linspace(0.4, 0.5, 5);
cb.TickLabels = arrayfun(@(v) sprintf('%.2f', v), cb.Ticks, 'UniformOutput', false);