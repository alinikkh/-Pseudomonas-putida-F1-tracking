close all
forced_index_test = [28];  %  % test_code
total_runs = length(all_runs);
num_forced = length(forced_index_test);
num_to_plot = min(50, total_runs);

% Exclude the forced indices from the random pool
available_idx = setdiff(1:total_runs, forced_index_test);

% Number of remaining runs to randomly choose
num_random = num_to_plot - num_forced;

% Randomly select the remaining indices
random_other_idx = randsample(available_idx, num_random);

% Combine forced + random
random_idx = [forced_index_test, random_other_idx];

num_to_plot = min(50, length(all_runs));
% rng(1);  % for reproducibility
%random_idx = randperm(length(all_runs), num_to_plot);
%122
% Prepare colormap
cmap = cool(32);  % You can use 'jet', 'parula', 'hot', etc.

figure; hold on; grid on; axis equal;
xlabel('\Deltax (\mum)', 'FontSize', 18);
ylabel('\Deltay (\mum)', 'FontSize', 18);

set(gca, 'FontSize', 18, 'LineWidth', 1.5);  % Axis font and line styling
xlim([-40 40]);
ylim([-40 40]);

for i = 1:num_to_plot
    run = all_runs{random_idx(i)};
    dx = run(:,1) - run(1,1);  % Shift to origin
    dy = run(:,2) - run(1,2);
    z = zeros(size(dx));
    c = all_straightness(random_idx(i)) * ones(size(dx));  % constant color per run

    % Plot using surface trick
    surface([dx'; dx'], [dy'; dy'], [z'; z'], [c'; c'], ...
        'FaceColor', 'none', ...
        'EdgeColor', 'interp', ...
        'LineWidth', 2);
end

% Color settings
colormap(cmap);
caxis([min(all_straightness) max(all_straightness)]);

cb = colorbar;
cb.Label.String = 'Straightness';
cb.Label.FontSize = 18;
cb.Ticks = linspace(min(all_straightness), max(all_straightness), 6);
cb.TickLabels = arrayfun(@(v) sprintf('%.2f', v), cb.Ticks, 'UniformOutput', false);

box on;
%%


% Randomly sample 50 runs (or fewer if less available)
num_to_plot = min(50, length(all_runs));
% rng(1);  % for reproducibility
random_idx = randperm(length(all_runs), num_to_plot);

% Prepare colormap
cmap = cool(32);  % You can use 'jet', 'parula', 'hot', etc.

figure; hold on; grid on; axis equal;
xlabel('\Deltax (\mum)', 'FontSize', 18);
ylabel('\Deltay (\mum)', 'FontSize', 18);

set(gca, 'FontSize', 18, 'LineWidth', 1.5);  % Axis font and line styling
xlim([-40 40]);
ylim([-40 40]);

for i = 1:num_to_plot
    run = all_runs{random_idx(i)};
    dx = run(:,1) - run(1,1);  % Shift to origin
    dy = run(:,2) - run(1,2);
    z = zeros(size(dx));
    c = all_straightness(random_idx(i)) * ones(size(dx));  % constant color per run

    % Plot using surface trick
    surface([dx'; dx'], [dy'; dy'], [z'; z'], [c'; c'], ...
        'FaceColor', 'none', ...
        'EdgeColor', 'interp', ...
        'LineWidth', 2);
end

% Color settings
colormap(cmap);
caxis([min(all_straightness) max(all_straightness)]);

cb = colorbar;
cb.Label.String = 'Straightness';
cb.Label.FontSize = 18;
cb.Ticks = linspace(min(all_straightness), max(all_straightness), 6);
cb.TickLabels = arrayfun(@(v) sprintf('%.2f', v), cb.Ticks, 'UniformOutput', false);

box on;
