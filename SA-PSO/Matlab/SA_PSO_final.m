clc; clear; close all;

% Load external TSP data
cC = load('TSP_2D_734.txt');
numCities = size(cC,1);

% holds all solutions per cooling loop for n instances of SA for all 10
% swarm particles
total_x1 = [];
total_x2 = [];
total_x3 = [];
total_x4 = [];
total_x5 = [];
total_x6 = [];
total_x7 = [];
total_x8 = [];
total_x9 = [];
total_x10 = [];

% holds number of equilibrium loops for each particle for every cooling loop
% for all n instances (1001 loops x (10 particles x 100 instances)
% sorry for this mess
all_xx = [];

total_endTime = [];

for record=1:1
    
    % Graphing
    total_solutions = [];
    total_xx = [];
    total_time = [];
    
    %PSO parameters
    ww = 0.3;
    c1 = 0.1;
    c2 = 0.3;
    numSwarm = 10;
    
    p_best = zeros([numSwarm, 1]); % value of personal best
    p_params = []; % parameters of personal best
    g_best = 0; % value of global best
    g_params = []; % parameters of global best
    
    xx = []; % xx = numEquilibriumLoops
    vv = [];
    
    %SA parameters
    numCoolingLoops = 1000;
    initMin = 0;
    initMax = 1000;
    pEnd = 0.00001;
    pStart = 0.3;        % Probability of accepting worse solution at the start
    pEnd = 0.00001;        % Probability of accepting worse solution at the end
    tStart = -1.0/log(pStart); % Initial temperature
    tEnd = -1.0/log(pEnd);     % Final temperature
    frac = (tEnd/tStart)^(1.0/(numCoolingLoops-1.0));% Fractional reduction per cycle
    
    % initialize PSO
    
    % randomize swarm
    for i = 1:numSwarm
        xx = [xx; randi(initMax-initMin)+initMin];
    end
    xx(1) = initMax; % set particle 1 to max
    
    vv = -0.1*xx;
    p_params = xx; % initialize
    display(xx)
    
    % Initialize SA
    cityRoutes_i = [];
    for z = 1:numSwarm
        cityRoutes_i = [cityRoutes_i; randperm(numCities)];
    end
    
    cityRoutes_i = repmat(randperm(numCities), numSwarm, 1) % same initial route instead of all random (above)
    
    cityRoutes_j = cityRoutes_i; % current route
    cityRoutes_k = cityRoutes_i; % new route
    cityRoutes_o = cityRoutes_i; % optimal route
    
    D_j = calcDist_vect(numCities, numSwarm, cityRoutes_i, cC, ones(numSwarm)); % current distance
    D_k = D_j; % new distance
    D_o = D_j; % optimal distance
    D = [D_j]; % route distance for each temperature setting
    
    numAcceptedSolutions = ones([numSwarm 1]); % number of changed solutions
    tCurrent = tStart;         % Current temperature = initial temperature
    DeltaE_avg = zeros([numSwarm, 1]);   % DeltaE Average
    tic;
    
    total_solutions = [total_solutions; D_j'];
    total_xx = [total_xx; xx'];
    
    for ii=1:numCoolingLoops
        disp(['Cycle: ',num2str(ii),' starting temperature: ',num2str(tCurrent), ' maxEq: ', num2str(max(xx))])
        startIter = tic;
        startInner = tic;
        numSolVect = ones([numSwarm 1]);
        
        innerTime = zeros(numSwarm, 1);
        stdMat = [];
        for jj=1:max(xx)
            %display(cityRoutes_j)
            maxIters = jj<=xx;
            cityRoutes_k = perturbRoute_vect(numCities, cityRoutes_j); % new perturbed route from old route
            D_k = calcDist_vect(numCities, numSwarm, cityRoutes_k, cC, maxIters) + D_k.*~maxIters;
            DeltaE = abs(D_k-D_j);
            if (ii==1 && jj==1) DeltaE_avg = DeltaE; end
            worseArr = D_k > D_j;
            p = exp(-DeltaE./(DeltaE_avg .* tCurrent));
            pacceptArr = p > rand();
            acceptArr = ~worseArr | pacceptArr & worseArr;
            numSolVect = numSolVect + acceptArr.*maxIters;
            cityRoutes_j = (cityRoutes_k.*acceptArr + cityRoutes_j.*~acceptArr).*maxIters + cityRoutes_j.*~maxIters;
            D_j = (D_k.*acceptArr + D_j.*~acceptArr).*maxIters + D_j.*~maxIters;
            
            numAcceptedSolutions = numAcceptedSolutions + acceptArr;
            DeltaE_avg = ((DeltaE_avg .* (numAcceptedSolutions-1.0) + DeltaE) ./ numAcceptedSolutions) .* acceptArr + DeltaE_avg.*~acceptArr;
            stdMat = [stdMat, D_j];
            innerTime(:,1) = innerTime.*~maxIters + toc(startInner).*maxIters;
        end
        
        % calculate stdev while considering when each particle stopped
        % iterating
        fitness = [];
        check = stdMat';
        %clf;
        %hold all
        for j = 1:numSwarm
            fitness(j) = std(check(1:xx(j),j));
            %plot(check(1:xx(j),j));
            
        end
        %{ display graph of partlices while iterating
        pause(0.0);
        display(fitness)
        display(xx')
        %hold off
        %}
        
        % Update swarm
        fitness = fitness';
        p_newArr = fitness > p_best;
        p_best = fitness.*p_newArr + p_best.*~p_newArr;
        p_params = xx.*p_newArr + p_params.*~p_newArr;
        [fMax, iMax] = max(fitness);
        
        g_best = fMax;
        g_params = xx(iMax);
        
        vv = ww.*vv + c1.*rand().*(p_params-xx) + c2.*rand().*(g_params - xx);
        xx = xx + vv;
        
        UBArr = xx > initMax;
        LBArr = xx < initMin;
        xx = round(UBArr.*initMax + xx.*~UBArr.*~LBArr + LBArr.*initMin);
        display(xx')
        
        tCurrent = frac * tCurrent; % Lower the temperature for next cycle
        D = [D D_k]; %record the route distance for each temperature setting
        optimalArr = D_k < D_o;
        D_o = D_k.*optimalArr + D_o.*~optimalArr; % Update optimal distance
        cityRoutes_o = cityRoutes_k.*optimalArr + cityRoutes_k.*~optimalArr;  % Update optimal route at each cycle
        endIter = toc(startIter);
        fprintf ('time: %f ', endIter)
        
        total_solutions = [total_solutions; D_j'];
        total_xx = [total_xx; xx'];
        total_time = [total_time; innerTime'];
        
    end
    
    all_xx = [all_xx, total_xx];
    display(num2str(D_k))
    display(num2str(D_o))
    endTime = toc    
    
    total_x1 = [total_x1, total_solutions(:,1)];
    total_x2 = [total_x2, total_solutions(:,2)];
    total_x3 = [total_x3, total_solutions(:,3)];
    total_x4 = [total_x4, total_solutions(:,4)];
    total_x5 = [total_x5, total_solutions(:,5)];
    total_x6 = [total_x6, total_solutions(:,6)];
    total_x7 = [total_x7, total_solutions(:,7)];
    total_x8 = [total_x8, total_solutions(:,8)];
    total_x9 = [total_x9, total_solutions(:,9)];
    total_x10 = [total_x10, total_solutions(:,10)];
    
    total_endTime = [total_endTime; sum(total_time)];
end
