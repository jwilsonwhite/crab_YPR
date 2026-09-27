function dungy_YPR

% Yield-per-recruit analysis for Dungeness crab under climate stress
% Based on original code by Ben Walker
% Authored by Will White. September 2026

%%%%%%%%%%%%%%%%%%%%%%%
% Model parameters
legal_cw = 159; % legal size limit (CA Commercial)
legal_cw_values = [159*0.9, 159]; % range of possible alternatives

ages = 3:5;
widths = 1:200;

% Allometry
wgt_b = 3; 
wgt_legal = 900; % g
wgt_a = wgt_legal/(legal_cw^wgt_b); % solve for a coefficient 
wgt_legal_values = [700, 900, 1100]; % sensitivity analysis;

% Molting probability
molt_alpha = -21.5;
molt_beta = 0.13;

% Molt increment
molt_increment_baseline = 30; % mm
molt_inc = 24:35;
molt_inc_dist = normcdf(molt_inc(2:end),molt_increment_baseline,2.5)-normcdf(molt_inc(1:(end-1)),molt_increment_baseline,2.5); % weighting function
molt_inc_dist = molt_inc_dist./sum(molt_inc_dist); % integrate to one
molt_inc = molt_inc(2:end); % shrink by one dimension bc of the integration step

% Natural mortality
mort_gam = 1.0; % coefficient for exponential decline in mortality with size
mort0 = 35; % coefficient for baseline mortality function

% Molting mortality
mort_molt_base = 0.05;

% Harvest mortality
harvest_prop = 0.95; %repmat(0.95,[3,1]); % proportion of legal crabs harvested in one season


% starting CW distribution
mean_cw_start = 140;
sd_cw_start = 5; 
cw_dist_tmp = normcdf(widths(2:end),mean_cw_start,sd_cw_start)-normcdf(widths(1:(end-1)),mean_cw_start,sd_cw_start); % weighting function
cw_dist = [0, cw_dist_tmp]; % add a 0 so it's the same length as 'widths'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run scenarios

% How this is organized:
% Baseline
% Mm
% M0
% pM
% Mm + M0
% pM +  Mm
% pM + M0
% All

% Whether or not each term is adjusted
molt_prob_multiplier =      [0, 1, 0, 0, 0, 1,  1, 1]; % reductions in molt probability
molt_mortality_multiplier = [0, 0, 1, 0, 1, 1,  0, 1]; % increases in molting mortality
baseline_mort_multiplier =  [0, 0, 0, 1, 1, 0,  1, 1]; % increases in natural mortality

% The level of adjustment
range = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3]; % 

YPR = nan(length(molt_prob_multiplier),1);
Yield = nan(length(1:max(ages)),length(molt_prob_multiplier));
CW = nan(length(widths),length(1:max(ages)),length(molt_prob_multiplier));

% Loop over scenarios
for s = 1:length(molt_prob_multiplier)
    for r =  1:length(range)
        for c = 1:length(legal_cw_values)

    % Now call the specific adjustment for each one
    molt_prob_multiplier_tmp = (1 - molt_prob_multiplier(s) * range(r) );
    molt_mort_multiplier_tmp = (1 + molt_mortality_multiplier(s) * range(r));
    baseline_mort_multiplier_tmp = (1 + baseline_mort_multiplier(s) * range(r));


    [YPR_tmp, Yield_tmp, CW_tmp] = do_YPR(cw_dist,ages,widths,wgt_a,wgt_b,molt_alpha,molt_beta,mort0,mort_gam,mort_molt_base,legal_cw_values(c),...
        molt_prob_multiplier_tmp,baseline_mort_multiplier_tmp,molt_mort_multiplier_tmp,harvest_prop,...
        molt_inc,molt_inc_dist);


    %Store the results
    YPR(s,r,c) = YPR_tmp;
    Yield(:,s,r,c) = Yield_tmp(:);
    CW(:,:,s,r,c) = CW_tmp;

        end % end loop over legal limit
    end % end loop over range
end % end loop over scenarios

YPR = YPR./YPR(1,1,2); % rescale YPR to be proportional to the baseline


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the plots
% Figure 2: YPR vs parameter change
figure(1)
clf

set(gcf,'units','cent','position',[10,10,7,7])


hold on
ScenLab = {'0','p(molt)','M_m_o_l_t','M_b_a_s_e','M_m_o_l_t & M_b_a_s_e','p(molt) & M_m_o_l_t','p(molt) & M_b_a_s_e','All'};
Col = flipud(parula(size(YPR,1)-1));
Style = {':','-'};


for s = 2:size(YPR,1) % loop over scenarios
    for c = 1:2 % legal size limit scenarios

    plot(range,YPR(s,:,c),'-','linewidth',1.5,'color',Col(s-1,:),'linestyle',Style{c})
    
    if c == 2
    text(range(end)*1.01,YPR(s,r,c),ScenLab{s});
    end

    end % end loop over legal limit
end % end  loop over s

plot([0,range(end)],[1,1-range(end)],'k--')
xlim([0 0.35])
set(gca,'tickdir','out','ticklength',[0.015 0.015])
ylabel('Proportion of maximum YPR')
xlabel('Proportional change in parameter')

print('-depsc2','Fig2_YPR.eps')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3: Size distributions over successive ages in 3 representative
% scenarios
% For the manuscript we ran this twice for each of the legal size limits
figure(2)
clf

set(gcf,'units','cent','position',[10,10,7,21])

subplot(3,1,1) % plot the baseline
hold on

LW = [1, 2, 3];
CM = jet(3);

for a = 1:3 % plot each age class
plot(widths,CW(:,a,1,end,2),'linewidth',LW(a),'color',CM(a,:))
end % end loop over ages

xlim([120,200])
ylim([0, 0.02])
set(gca,'tickdir','out','ticklength',[0.015, 0.015],'ygrid','on')


subplot(3,1,2) % plot the impacted case
hold on

for a = 1:3 % plot each age class
plot(widths,CW(:,a,2,end,2),'linewidth',LW(a),'color',CM(a,:))
end % end loop over ages

xlim([120,200])
ylim([0, 0.02])
set(gca,'tickdir','out','ticklength',[0.015, 0.015],'ygrid','on')


subplot(3,1,3) % plot the impacted case
hold on

for a = 1:3 % plot each age class
    plot(widths,CW(:,a,4,end,2),'linewidth',LW(a),'color',CM(a,:))
end % end loop over ages

xlim([120,200])
ylim([0, 0.02])
set(gca,'tickdir','out','ticklength',[0.015, 0.015],'ygrid','on')

print('-depsc2','Fig3_size_cohorts.eps')

end % end dungy_YPR


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Helper Functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main model function: do_YPR
function [YPR, Yield, CW] = do_YPR(start_cw,ages,widths,wgt_a,wgt_b,molt_alpha,molt_beta,mort0,mortgam,mort_molt_base,...
                    legal_cw,molt_prob_multiplier,baseline_mort_multiplier,molt_mortality_multiplier,...
                    harvest_prop,molt_increment,molt_inc_dist)

CW = zeros(length(widths),length(ages)); % state variable: carapace width
CW(:,1) = start_cw(:); % Initial CW distribution of the first age class

Yield = zeros(1:length(max(ages)),1);


for a = 2:max(ages) % loop over ages

    % Order of operations: molt, growth, mortality, harvest

    % Molting
    p_molt = p_molt_baseline(CW(:,a-1),molt_alpha,molt_beta); % probability of molting
    p_molt = p_molt .* molt_prob_multiplier;
    p_molt = min(max(p_molt,0),1); % constrain probabilities

    % Split the cohort into molters and non-molters
    CW_grow = CW(:,a-1).*p_molt(:); % the part of the cohort that molts & grows
    CW_no_grow = CW(:,a-1).*(1-p_molt(:)); % the part of the cohort that does not molt & remains the same

    % Growth
    if isscalar(molt_increment) % constant molt_increment
        molt_increment = round(molt_increment);
        CW_grow = [zeros(molt_increment,1); CW_grow(1:(end-molt_increment))]; 

    else % variable molt increment

        CW_tmp = zeros(length(CW_grow),length(molt_increment)); % one column for each possible growth increment


        for m = 1:length(molt_increment)
        CW_tmp(:,m) = [zeros(molt_increment(m),1); CW_grow(1:(end-molt_increment(m)))].*molt_inc_dist(m);
        end % end loop over molt_increment

        CW_grow = sum(CW_tmp,2); % combine across growth increments

    end % end if isscalar(molt_increment)

    % Mortality
    M_base = baseline_mort(widths,mort0,mortgam).*baseline_mort_multiplier;
    M_mort = mort_molt_base .* molt_mortality_multiplier;

    CW_no_grow = CW_no_grow .* exp(-M_base(:)); % mortality of non-molters
    CW_grow = CW_grow .* exp(-(M_base(:) + M_mort(:)));

    CW(:,a) = CW_no_grow + CW_grow; % re-combine the two groups before applying harvest

    % Harvest
    Harvest = CW(:,a) .* (widths(:) >= legal_cw) .* harvest_prop; % legal-size crabs are harvested

    CW(:,a) = CW(:,a) - Harvest; % remove harvested crabs

    Yield(a) = sum( Harvest(:) .* cw_to_weight(widths(:),wgt_a,wgt_b));


end % end loop over age

YPR = sum(Yield);

end % end do_YPR

% Helper function: carapace width to weight
function W = cw_to_weight(cw,a,b)
W = a.*cw.^b;
end

% Helper function: molt probability
function P_molt = p_molt_baseline(cw,a,b)
P_molt = 1./(1 + exp(a + b.*cw));
end

% Helper function: baseline natural mortality
function M = baseline_mort(cw,mort0,mort_gam)
M = mort0 .* cw .^ (-mort_gam);
end