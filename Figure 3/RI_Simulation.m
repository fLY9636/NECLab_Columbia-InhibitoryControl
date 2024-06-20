function SIM_subj_events = RI_Simulation(loopIDX, num_iters)
tic
SIM_subj_events = cell(1, length(loopIDX));
cnt = 0;
for i = loopIDX(:).'
    cnt = cnt+1;
    SIM_session_events = cell(num_iters, 1);
    for j = 1:num_iters
        rng('shuffle')
        % generate random tone durations
        CDF_values = [0.1 0.08];
        points = [5 7.5];
        lambda = 1/(points(2)-points(1));
        scale = 1/CDF_values(1);
        random_numbers1 = scale*exprnd(lambda, 80, 1);
        random_numbers1 = random_numbers1+5;
        random_numbers1(random_numbers1>12|random_numbers1<5) = [];
        num_trials = length(random_numbers1);
        % generate random free-period durations
        lower_bound = 12;
        upper_bound = 17;
        random_numbers2 = unifrnd(lower_bound, upper_bound, num_trials, 1);
        segments = zeros(1, 2*num_trials);
        segments(1:2:end) = random_numbers2;
        segments(2:2:end) = random_numbers1;
        endpoints = zeros(1, 2*num_trials+1);
        for n = 2:2*num_trials+1
            endpoints(n) = endpoints(n-1)+segments(n-1);
        end
        SIM_session_events{j} = endpoints;
    end
    SIM_subj_events{cnt} = SIM_session_events;
end
return
toc