close all;
% global truncation error experiment for all four methods
% integrates each test IVP over t = [0,5] for 40 step sizes (even in log
% space) and compares the last point to the analytical solution at t = 5
% fits error = k*h^p and error = k*(rate function calls)^p for each method
% runs both test IVPs, only plots the first one

% figures for the report, saved one folder up so they don't end up in the repo
fig_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Assignment 3 Report Figures');
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

tspan = [0,5];
h_ref_list = logspace(-4,0,40);

% same order as the rows in the table
step_funcs = {@forward_euler_step,@backward_euler_step,@explicit_midpoint_step,@implicit_midpoint_step};
method_names = {'Forward Euler','Backward Euler','Explicit Midpoint','Implicit Midpoint'};
rate_funcs = {@rate_func01,@rate_func02};
solutions = {@solution01,@solution02};

% same idea as the local experiment, only fit where it's scaling like h^p
h_fit_max = 1e-2;
error_fit_min = 1e-12;

% one row per method, one column per test
p_h = zeros(4,2);
p_calls = zeros(4,2);

for test_num = 1:2
    rate_func = rate_funcs{test_num};
    solution = solutions{test_num};
    X0 = solution(tspan(1));
    X_true_end = solution(tspan(2));

    h_avg_list = zeros(4,length(h_ref_list));
    calls_list = zeros(4,length(h_ref_list));
    error_list = zeros(4,length(h_ref_list));
    k_h = zeros(4,1);
    k_calls = zeros(4,1);

    for m = 1:4
        for i = 1:length(h_ref_list)
            [~,X_list,h_avg,num_evals] = fixed_step_integration(rate_func,step_funcs{m},tspan,X0,h_ref_list(i));

            % last row of X_list is the final value, turned back into a column
            error_list(m,i) = norm(X_list(end,:)' - X_true_end);
            h_avg_list(m,i) = h_avg;
            calls_list(m,i) = num_evals;
        end

        fit_spots = h_avg_list(m,:) <= h_fit_max & error_list(m,:) > error_fit_min;
        [p_h(m,test_num),k_h(m)] = generate_error_fit(h_avg_list(m,fit_spots),error_list(m,fit_spots));
        [p_calls(m,test_num),k_calls(m)] = generate_error_fit(calls_list(m,fit_spots),error_list(m,fit_spots));
    end

    % hang on to test 1 for the plots
    if test_num == 1
        h_avg_list01 = h_avg_list;
        calls_list01 = calls_list;
        error_list01 = error_list;
        k_h01 = k_h;
        k_calls01 = k_calls;
    end
end

% blue/red for the explicit ones, green/purple for the implicit ones
colors = [0.00, 0.25, 0.85; 0.2235, 1.0000, 0.0784; 0.90, 0.15, 0.05; 1, 1, 0];
marker_list = {'o','s','o','s'};
size_list = [9,5,9,5];

% plot 1: error vs step size, explicit methods with fit lines
fig1 = figure('Color','w','Position',[100,100,1000,700]);
for m = [1,3]
    loglog(h_avg_list01(m,:),error_list01(m,:),'o','Color',colors(m,:),'MarkerFaceColor',colors(m,:), ...
        'MarkerSize',7,'DisplayName',[method_names{m},' data']);
    hold on;
    loglog(h_avg_list01(m,:),k_h01(m)*h_avg_list01(m,:).^p_h(m,1),'-','Color','k','LineWidth',2, ...
        'DisplayName',sprintf('%s fit, $p$ = %.2f',method_names{m},p_h(m,1)));
end
grid on; box on;
xlabel('Step size, $h$ (-)','Interpreter','latex');
ylabel('Global truncation error, $\epsilon$ (-)','Interpreter','latex');
title('Global Truncation Error vs. Step Size for the Explicit Methods','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','FontSize',16);
legend('Interpreter','latex','Location','northwest','FontSize',14);
exportgraphics(fig1,fullfile(fig_dir,'global_truncation_step_size_explicit_fit.png'),'Resolution',300);

% plot 2: error vs step size, all four methods, no fit lines
fig2 = figure('Color','w','Position',[100,100,1000,700]);
for m = 1:4
    loglog(h_avg_list01(m,:),error_list01(m,:),marker_list{m},'Color',colors(m,:),'MarkerFaceColor',colors(m,:), ...
        'MarkerSize',size_list(m),'MarkerEdgeColor', 'k', 'DisplayName',method_names{m});
    hold on;
end
grid on; box on;
xlabel('Step size, $h$ (-)','Interpreter','latex');
ylabel('Global truncation error, $\epsilon$ (-)','Interpreter','latex');
title('Global Truncation Error vs. Step Size for All Four Methods','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','FontSize',16);
legend('Interpreter','latex','Location','northwest','FontSize',14);
exportgraphics(fig2,fullfile(fig_dir,'global_truncation_step_size_all_methods.png'),'Resolution',300);

% plot 3: error vs number of rate function calls, explicit methods with fit lines
% more calls means a smaller error, so the data goes from top left to
% bottom right and the top right corner is empty
fig3 = figure('Color','w','Position',[100,100,1000,700]);
for m = [1,3]
    loglog(calls_list01(m,:),error_list01(m,:),'o','Color',colors(m,:),'MarkerFaceColor',colors(m,:), ...
        'MarkerSize',7,'DisplayName',[method_names{m},' data']);
    hold on;
    loglog(calls_list01(m,:),k_calls01(m)*calls_list01(m,:).^p_calls(m,1),'-','Color','k','LineWidth',2, ...
        'DisplayName',sprintf('%s fit, $p$ = %.2f',method_names{m},p_calls(m,1)));
end
grid on; box on;
xlabel('Number of rate function calls, $N$ (-)','Interpreter','latex');
ylabel('Global truncation error, $\epsilon$ (-)','Interpreter','latex');
title('Global Truncation Error vs. Rate Function Calls for the Explicit Methods','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','FontSize',16);
legend('Interpreter','latex','Location','northeast','FontSize',14);
exportgraphics(fig3,fullfile(fig_dir,'global_truncation_rate_calls_explicit_fit.png'),'Resolution',300);

% plot 4: error vs number of rate function calls, all four methods, no fit lines
fig4 = figure('Color','w','Position',[100,100,1000,700]);
for m = 1:4
    loglog(calls_list01(m,:),error_list01(m,:),marker_list{m},'Color',colors(m,:),'MarkerFaceColor',colors(m,:), ...
        'MarkerSize',size_list(m),'DisplayName',method_names{m});
    hold on;
end
grid on; box on;
xlabel('Number of rate function calls, $N$ (-)','Interpreter','latex');
ylabel('Global truncation error, $\epsilon$ (-)','Interpreter','latex');
title('Global Truncation Error vs. Rate Function Calls for All Four Methods','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','FontSize',16);
legend('Interpreter','latex','Location','northeast','FontSize',14);
exportgraphics(fig4,fullfile(fig_dir,'global_truncation_rate_calls_all_methods.png'),'Resolution',300);

% tables for the report (p is dimensionless)
fprintf('\nTable: global truncation error scaling with step size, p from error = k*h^p\n');
fprintf('(t = [%g,%g], fit over h <= %.0e, all values dimensionless (-))\n',tspan(1),tspan(2),h_fit_max);
fprintf('%-20s %16s %16s\n','Method','Test 1: p (-)','Test 2: p (-)');
for m = 1:4
    fprintf('%-20s %16.4f %16.4f\n',method_names{m},p_h(m,1),p_h(m,2));
end

fprintf('\nTable: global truncation error scaling with rate function calls, p from error = k*N^p\n');
fprintf('%-20s %16s %16s\n','Method','Test 1: p (-)','Test 2: p (-)');
for m = 1:4
    fprintf('%-20s %16.4f %16.4f\n',method_names{m},p_calls(m,1),p_calls(m,2));
end

% average rate function calls per step for test 1, shows how much extra
% the newton solve costs for the implicit methods
fprintf('\naverage rate function calls per step (test 1):\n');
for m = 1:4
    fprintf('%-20s %.2f\n',method_names{m},mean(calls_list01(m,:)./((tspan(2)-tspan(1))./h_avg_list01(m,:))));
end

% also saved as csvs so they're easy to drop into the word tables
global_h_table = table(method_names',p_h(:,1),p_h(:,2),'VariableNames',{'Method','Test1_p','Test2_p'});
writetable(global_h_table,fullfile(fig_dir,'global_truncation_step_size_p_table.csv'));
global_calls_table = table(method_names',p_calls(:,1),p_calls(:,2),'VariableNames',{'Method','Test1_p','Test2_p'});
writetable(global_calls_table,fullfile(fig_dir,'global_truncation_rate_calls_p_table.csv'));
