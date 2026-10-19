close all;
% local truncation error experiment for all four methods
% starts ON the analytical solution at t_ref, takes one step of size h,
% and compares that to the analytical solution at t_ref + h
% 100 step sizes from 1e-5 to 1e1 (even in log space), then fits
% error = k*h^p to get p. runs both test IVPs, only plots the first one

% figures for the report, saved one folder up so they don't end up in the repo
fig_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Assignment 3 Report Figures');
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% not a multiple of pi so the derivative isn't zero here
t_ref = 0.492;
h_list = logspace(-5,1,100);

% same order as the rows in the table
step_funcs = {@forward_euler_step,@backward_euler_step,@explicit_midpoint_step,@implicit_midpoint_step};
method_names = {'Forward Euler','Backward Euler','Explicit Midpoint','Implicit Midpoint'};
rate_funcs = {@rate_func01,@rate_func02};
solutions = {@solution01,@solution02};

% only fit the part where the error is really scaling like h^p. too big of
% a step and the higher order terms take over, too small and it's just
% round off noise
h_fit_max = 1e-2;
error_fit_min = 1e-12;

% row 1 is |X(t+h)-X(t)|, rows 2 to 5 are the methods, one column per test
p_list = zeros(5,2);

for test_num = 1:2
    rate_func = rate_funcs{test_num};
    solution = solutions{test_num};
    X_ref = solution(t_ref);

    diff_list = zeros(1,length(h_list));
    error_list = zeros(4,length(h_list));
    for i = 1:length(h_list)
        h = h_list(i);
        X_true = solution(t_ref+h);

        % how much X itself changes over the step
        diff_list(i) = norm(X_true - X_ref);

        % one step of each method from the same starting point
        for m = 1:4
            XB = step_funcs{m}(rate_func,t_ref,X_ref,h);
            error_list(m,i) = norm(XB - X_true);
        end
    end

    k_list = zeros(5,1);
    fit_spots = h_list <= h_fit_max & diff_list > error_fit_min;
    [p_list(1,test_num),k_list(1)] = generate_error_fit(h_list(fit_spots),diff_list(fit_spots));
    for m = 1:4
        fit_spots = h_list <= h_fit_max & error_list(m,:) > error_fit_min;
        [p_list(m+1,test_num),k_list(m+1)] = generate_error_fit(h_list(fit_spots),error_list(m,fit_spots));
    end

    % hang on to test 1 for the plots
    if test_num == 1
        diff_list01 = diff_list;
        error_list01 = error_list;
        k_list01 = k_list;
    end
end

% blue/red for the explicit ones, green/purple for the implicit ones
colors = [0.00, 0.25, 0.85; 0.2235, 1.0000, 0.0784; 0.90, 0.15, 0.05; 1, 1, 0];
diff_color = [0,0.6,0.6];

% plot 1: forward euler and explicit midpoint with fit lines, plus the
% difference |X(t+h)-X(t)| to compare against
fig1 = figure('Color','w','Position',[100,100,1000,700]);
loglog(h_list,diff_list01,'o','Color',diff_color,'MarkerFaceColor',diff_color,'MarkerSize',6, ...
    'DisplayName','$|X(t+h)-X(t)|$ data');
hold on;
loglog(h_list,k_list01(1)*h_list.^p_list(1,1),'-','Color','k','LineWidth',2, ...
    'DisplayName',sprintf('$|X(t+h)-X(t)|$ fit, $p$ = %.2f',p_list(1,1)));
for m = [1,3]
    loglog(h_list,error_list01(m,:),'o','Color',colors(m,:),'MarkerFaceColor',colors(m,:),'MarkerSize',6, ...
        'DisplayName',[method_names{m},' data']);
    loglog(h_list,k_list01(m+1)*h_list.^p_list(m+1,1),'-','Color','k','LineWidth',2, ...
        'DisplayName',sprintf('%s fit, $p$ = %.2f',method_names{m},p_list(m+1,1)));
end
grid on; box on;
xlabel('Step size, $h$ (-)','Interpreter','latex');
ylabel('Local truncation error, $\epsilon$ (-)','Interpreter','latex');
title('Local Truncation Error vs. Step Size for the Explicit Methods','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','FontSize',16);
% data goes from bottom left to top right, top left is always empty
legend('Interpreter','latex','Location','northwest','FontSize',13);
exportgraphics(fig1,fullfile(fig_dir,'local_truncation_explicit_fit.png'),'Resolution',300);

% plot 2: all four methods, data only. forward and backward euler land on
% top of each other (same h^2 term, opposite sign) so the implicit ones
% get smaller markers drawn on top so both show up
fig2 = figure('Color','w','Position',[100,100,1000,700]);
marker_list = {'o','s','o','s'};
size_list = [9,5,9,5];
for m = 1:4
    loglog(h_list,error_list01(m,:),marker_list{m},'Color',colors(m,:),'MarkerFaceColor',colors(m,:), ...
        'MarkerSize',size_list(m),'MarkerEdgeColor','k','DisplayName',method_names{m});
    hold on;
end
grid on; box on;
xlabel('Step size, $h$ (-)','Interpreter','latex');
ylabel('Local truncation error, $\epsilon$ (-)','Interpreter','latex');
title('Local Truncation Error vs. Step Size for All Four Methods','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','FontSize',16);
legend('Interpreter','latex','Location','northwest','FontSize',14);
exportgraphics(fig2,fullfile(fig_dir,'local_truncation_all_methods.png'),'Resolution',300);

% table for the report (p is dimensionless)
row_names = {'|X(t+h)-X(t)|','Forward Euler','Backward Euler','Explicit Midpoint','Implicit Midpoint'};
fprintf('\nTable: local truncation error scaling with step size, p from error = k*h^p\n');
fprintf('(t_ref = %.3f, fit over %.0e <= h <= %.0e, all values dimensionless (-))\n',t_ref,min(h_list),h_fit_max);
fprintf('%-20s %16s %16s\n','Method','Test 1: p (-)','Test 2: p (-)');
for r = 1:5
    fprintf('%-20s %16.4f %16.4f\n',row_names{r},p_list(r,1),p_list(r,2));
end

% also saved as a csv so it's easy to drop into the word table
local_table = table(row_names',p_list(:,1),p_list(:,2),'VariableNames',{'Method','Test1_p','Test2_p'});
writetable(local_table,fullfile(fig_dir,'local_truncation_p_table.csv'));
