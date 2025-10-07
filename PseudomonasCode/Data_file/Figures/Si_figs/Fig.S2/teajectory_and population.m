clear all; close all; clc;
tic;

%% ---------------- PARAMETERS ----------------
total_time = 175;      % seconds
dt = 1e-2;             % timestep
n_steps = round(total_time/dt);

speed = 20;            % µm/s
runtime_lambda = 0.82; % s
x_max = 800;           % main channel length
y_max = 60;            % main channel height
D_r = 0.03;            % rotational diffusion
D_t = 0.4;             % translational diffusion
L_cell = 1.5;          % µm (effective length for dp drift)
D=1500                 %um^/s %salt diffusion
%% ---------------- RESERVOIR PARAMETERS ----------------
x_res_min = -800; x_res_max = 0;
y_res_min = 0;   y_res_max = 60;
N_res = 200;

%% ---------------- Run PARAMETERS ----------------
N_MC = 3;                 % independent runs
nBac_init = 200;          % initial main channel bacteria

%% ---------------- DIFFUSIOPHORESIS PARAMETERS ----------------
M = 100;                  % mobility
c_min = 1; c_max = 100;  % salt concentration at channel ends

%% ---------------- TRAJECTORY PARAMETERS ----------------
nTrack = 2000;               % max number of bacteria to track
time_windows = [1 25; 75 100; 150 175]; % [start end] in seconds
traj_time_idx = round(time_windows/dt);

%% ---------------- STORAGE ----------------
populations_all = zeros(N_MC, n_steps);
traj_results = cell(N_MC,1);

for mc = 1:N_MC
    fprintf(' Randome run %d/%d\n', mc, N_MC);

    %% --- Initialize main channel bacteria
    px = rand(nBac_init,1) * x_max;
    py = rand(nBac_init,1) * y_max;
    pth = 2*pi*rand(nBac_init,1);
    prt = -runtime_lambda*log(rand(nBac_init,1));  % exponential

    % Assign persistent IDs
    track_ids = (1:nBac_init)';

    % Initially tracked bacteria
    nBac_current = numel(track_ids);
    tracked_ids = [];
    if nBac_current > 0
        nAdd = min(nTrack, nBac_current);
        tracked_ids = track_ids(randperm(nBac_current,nAdd))';
    end

    %% --- Initialize reservoir bacteria
    res_x = x_res_min + (x_res_max-x_res_min)*rand(N_res,1);
    res_y = y_res_min + (y_res_max-y_res_min)*rand(N_res,1);
    res_th = 2*pi*rand(N_res,1);
    res_rt = -runtime_lambda*log(rand(N_res,1));

    %% --- Trajectory storage
    traj_x = NaN(numel(tracked_ids), n_steps);
    traj_y = NaN(numel(tracked_ids), n_steps);
    if ~isempty(tracked_ids)
        [~, ia, ib] = intersect(track_ids, tracked_ids);
        traj_x(ib,1) = px(ia);
        traj_y(ib,1) = py(ia);
    end

    %% --- Simulation loop
    for step = 1:n_steps
        %% --- Update reservoir bacteria
        res_rt = res_rt - dt;
        reset_idx = (res_rt <= 0);
        res_th(reset_idx) = 2*pi*rand(sum(reset_idx),1);
        res_rt(reset_idx) = -runtime_lambda*log(rand(sum(reset_idx),1));
        %res_th = res_th + sqrt(2*D_r*dt).*randn(size(res_th));
t_current = step*dt;

% Diffusiophoresis for reservoir bacteria near channel entrance 
 dc_dx_max = (c_max - c_min) / (2 * sqrt(pi * D * t_current)); 
% max gradient
 v_dp_max = M * dc_dx_max / (c_min+9.5); 
 
% v_dp at entrance
%v_dp_max = M * slope0 ./ (cmin + 9.0);  % same +90 as you used for normalization
near_channel_idx = (res_x > -100);   % reservoir region to apply entrance slope
v_dp_res = zeros(size(res_x));
v_dp_res(near_channel_idx) = v_dp_max;

        omega = (-v_dp_res / L_cell) .* sin(res_th);
        res_th = res_th + omega*dt + sqrt(2*D_r*dt).*randn(numel(res_th),1);

        % Update reservoir positions with v_dp
        res_x_new = res_x + speed*cos(res_th)*dt + sqrt(2*D_t*dt).*randn(size(res_x)) + v_dp_res*dt;
        res_y_new = res_y + speed*sin(res_th)*dt + sqrt(2*D_t*dt).*randn(size(res_y));


        % Reflective x
        left_idx = res_x_new < x_res_min;
        res_x_new(left_idx) = 2*x_res_min - res_x_new(left_idx);
        res_th(left_idx) = pi - res_th(left_idx);

        % Reflective y
        low_idx = res_y_new < y_res_min;
        res_y_new(low_idx) = 2*y_res_min - res_y_new(low_idx);
        res_th(low_idx) = -res_th(low_idx);

        high_idx = res_y_new > y_max;
        res_y_new(high_idx) = 2*y_max - res_y_new(high_idx);
        res_th(high_idx) = -res_th(high_idx);

        res_x = res_x_new;
        res_y = res_y_new;

        %% --- Bacteria entering channel
        enter_idx = (res_x >=0);
        n_enter = sum(enter_idx);
        if n_enter > 0
            new_ids = max(track_ids) + (1:n_enter)';  % unique IDs
            px = [px; zeros(n_enter,1)];
            py = [py; res_y(enter_idx)];
            pth = [pth; 2*pi*rand(n_enter,1)];
            prt = [prt; -runtime_lambda*log(rand(n_enter,1))];
            track_ids = [track_ids; new_ids];

            n_tracked = numel(tracked_ids);
            if n_tracked < nTrack
                n_add = min(nTrack - n_tracked, n_enter);
                if n_add > 0
                    % randomly pick n_add from entering bacteria and make column vector
                    new_tracked = new_ids(randperm(n_enter, n_add));
                    tracked_ids = [tracked_ids(:); new_tracked(:)];  % ensure both are column
                end
            end
        end

        % Remove entered bacteria from reservoir
        res_x(enter_idx) = [];
        res_y(enter_idx) = [];
        res_th(enter_idx) = [];
        res_rt(enter_idx) = [];

        %% --- Update main channel bacteria
        prt = prt - dt;
        tumble_idx = (prt <= 0);
        pth(tumble_idx) = 2*pi*rand(sum(tumble_idx),1);
        prt(tumble_idx) = -runtime_lambda*log(rand(sum(tumble_idx),1));

        % Position-dependent diffusiophoresis
        % --- Use semi-infinite erf profile for diffusiophoresis
        if ~isempty(px)
            t_current = step*dt;  % current simulation time
            % erf-based concentration profile
            c_x = c_max + (c_min-c_max) * erfc(px / (2*sqrt(D*t_current)));

            % gradient dc/dx = -(c_min-c_max)/(2*sqrt(pi*D*t)) * exp(-x^2/(4*D*t))
            dc_dx = (c_max - c_min) ./ (1*sqrt(pi*D*t_current)) .* exp(-px.^2 ./ (4*D*t_current));

            % v_dp = M * (dc/dx)/c
            v_dp = M * dc_dx ./ (c_x+9.5);
        else
            v_dp = [];
        end

        omega = (-v_dp / L_cell) .* sin(pth);
        pth = pth + omega*dt + sqrt(2*D_r*dt).*randn(numel(pth),1);

        x_new = px + speed.*cos(pth)*dt + sqrt(2*D_t*dt).*randn(numel(px),1) + v_dp*dt;
        y_new = py + speed.*sin(pth)*dt + sqrt(2*D_t*dt).*randn(numel(py),1);

        % Reflective y
        below_idx = (y_new < 0);
        y_new(below_idx) = -y_new(below_idx);
        pth(below_idx & sin(pth)<0) = -pth(below_idx & sin(pth)<0);

        above_idx = (y_new > y_max);
        y_new(above_idx) = 2*y_max - y_new(above_idx);
        pth(above_idx & sin(pth)>0) = -pth(above_idx & sin(pth)>0);

        % Exit back to reservoir
        back_idx = (x_new < 0);
        if any(back_idx)
            res_x = [res_x; x_res_max*ones(sum(back_idx),1)];
            res_y = [res_y; y_new(back_idx)];
            res_th = [res_th; 2*pi*rand(sum(back_idx),1)];
            res_rt = [res_rt; -runtime_lambda*log(rand(sum(back_idx),1))];
        end

        % Reflective closed channel end
        right_idx = (x_new > x_max);
        x_new(right_idx) = 2*x_max - x_new(right_idx);
        pth(right_idx & cos(pth)>0) = pi - pth(right_idx & cos(pth)>0);

        % Update survivors
        keep_idx = ~back_idx;
        px = x_new(keep_idx);
        py = y_new(keep_idx);
        pth = pth(keep_idx);
        prt = prt(keep_idx);
        track_ids = track_ids(keep_idx);

        %% --- Save trajectories
        if ~isempty(tracked_ids) && ~isempty(px)
            [~, ia, ib] = intersect(track_ids, tracked_ids);
            traj_x(:,step) = NaN;
            traj_y(:,step) = NaN;
            traj_x(ib,step) = px(ia);
            traj_y(ib,step) = py(ia);
        end

        %% --- Save population
        populations_all(mc,step) = numel(px);
    end

    traj_results{mc} = struct('x',traj_x,'y',traj_y);
end

%% ---------------- PLOT TRAJECTORIES SEPARATELY ----------------
colors_max = 1000; % max distinct colors
for mc = 1:N_MC
    traj = traj_results{mc};
    nBac_plot = size(traj.x,1);
    cmap = lines(min(nBac_plot,colors_max));

    for w = 1:size(time_windows,1)
        t_start = round(time_windows(w,1)/dt);
        t_end   = round(time_windows(w,2)/dt);

        figure; hold on; grid on;
        figure; hold on; grid on;

        % plot your trajectories
        for b = 1:nBac_plot
            plot(traj.x(b,t_start+1:t_end), traj.y(b,t_start+1:t_end), ...
                'Color', cmap(mod(b-1,size(cmap,1))+1,:), 'LineWidth', 1.5);
        end

        xlabel('x (\mum)','FontSize',14);
        ylabel('y (\mum)','FontSize',14);
        xlim([0 x_max]); ylim([0 y_max]);

        axis equal  % makes units in x and y the same
        title(sprintf('Trajectories MC=%d, Time %d-%d s', mc, time_windows(w,1), time_windows(w,2)),'FontSize',16);
        ;
    end
end

%% ---------------- PLOT POPULATION ----------------
time_axis = (0:n_steps-1)*dt;
pop_mean = mean(populations_all,1)/nBac_init;
pop_std  = std(populations_all,0,1)/nBac_init;

figure; hold on; grid on;
fill([time_axis fliplr(time_axis)], [pop_mean+pop_std fliplr(pop_mean-pop_std)], ...
    'b','FaceAlpha',0.2,'EdgeColor','none');
plot(time_axis, pop_mean,'b','LineWidth',2);
xlabel('Time (s)','FontSize',18);
ylabel('Normalized cell population','FontSize',18);
title('Population Dynamics','FontSize',16);


toc
%%
%% ---------------- PLOT ANOTHER 50 RANDOM TRAJECTORIES ----------------
for mc = 1:N_MC
    traj = traj_results{mc};
    nBac_total = size(traj.x,1);

    %rand_ids = randperm(nBac_total, min(50,nBac_total));

    for w = 1:size(time_windows,1)

        rand_ids = randperm(nBac_total, min(50,nBac_total));

        t_start = round(time_windows(w,1)/dt);
        t_end   = round(time_windows(w,2)/dt);

        figure; hold on; grid on;
        cmap = lines(numel(rand_ids));

        for k = 1:numel(rand_ids)
            b = rand_ids(k);
            plot(traj.x(b,t_start+1:t_end), traj.y(b,t_start+1:t_end), ...
                'Color', cmap(mod(k-1,size(cmap,1))+1,:), 'LineWidth', 1.5);
        end

% draw walls
plot([0 x_max], [0 0], 'k-', 'LineWidth', 2);        % bottom
plot([0 x_max], [y_max y_max], 'k-', 'LineWidth', 2);% top
plot([x_max x_max], [0 y_max], 'k-', 'LineWidth', 2);% right end
%plot([0 0], [0 y_max], 'k--', 'LineWidth', 1.5);     % reservoir entrance

xlabel('x (\mum)','FontSize',14);
ylabel('y (\mum)','FontSize',14);
xlim([-10 x_max+10]);   % extend to show reservoir region if you want
ylim([-10 60]);
axis equal

        title(sprintf(' (MC=%d, %d-%d s)', ...
            mc, time_windows(w,1), time_windows(w,2)), 'FontSize',16);
    end
end
