close all;
% runs each integration method on the first test IVP (eqn 5) for a few
% step sizes and plots them against the analytical solution x(t) = cos(t)
% one figure per method, the numerical solution should land closer and
% closer to the dashed black line as h gets smaller

% figures for the report, saved one folder up so they don't end up in the repo
fig_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Assignment 3 Report Figures');
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

tspan = [0,10];
X0 = solution01(tspan(1));

% forward euler and explicit midpoint use their own integration functions
% from day 10, the implicit ones go through fixed_step_integration
fe_integration = @(h_ref) forward_euler_fixed_step_integration(@rate_func01,tspan,X0,h_ref);
em_integration = @(h_ref) explicit_midpoint_fixed_step_integration(@rate_func01,tspan,X0,h_ref);
be_integration = @(h_ref) fixed_step_integration(@rate_func01,@backward_euler_step,tspan,X0,h_ref);
im_integration = @(h_ref) fixed_step_integration(@rate_func01,@implicit_midpoint_step,tspan,X0,h_ref);

% forward euler goes unstable on this IVP once h > 0.4 (see the stability
% part) so its biggest step stays under that
plot_example(fe_integration,[0.435,0.42,0.1],tspan,'Forward Euler', ...
    fullfile(fig_dir,'forward_euler_example_integration.png'));
plot_example(em_integration,[0.42,0.41,0.1],tspan,'Explicit Midpoint', ...
    fullfile(fig_dir,'explicit_midpoint_example_integration.png'));

% implicit methods stay stable with big steps so they can start bigger
plot_example(be_integration,[2,0.7,0.05],tspan,'Backward Euler', ...
    fullfile(fig_dir,'backward_euler_example_integration.png'));
plot_example(im_integration,[1.5,0.8,0.05],tspan,'Implicit Midpoint', ...
    fullfile(fig_dir,'implicit_midpoint_example_integration.png'));


%runs one method for every h_ref in the list and plots it with the
%analytical solution, then saves the figure
function plot_example(integration_func,h_ref_list,tspan,method_name,file_name)
    % biggest step first, one color each
    colors = [0.93,0.69,0.13; 0.85,0.33,0.10; 0.47,0.67,0.19; 0.00,0.45,0.74];

    fig = figure('Color','w','Position',[100,100,1000,650]);
    hold on; box on; grid on;

    % smallest step gets drawn first so the bigger (worse) ones end up on
    % top and you can actually see them miss
    for i = length(h_ref_list):-1:1
        [t_list,X_list,h_avg] = integration_func(h_ref_list(i));
        num_plots(i) = plot(t_list,X_list,'o-','Color',colors(i,:),'MarkerFaceColor',colors(i,:), ...
            'MarkerSize',5,'LineWidth',1.5,'DisplayName',sprintf('%s, $h$ = %.3f',method_name,h_avg));

        % how far off the last point is, just to check it's shrinking
        fprintf('%s h = %.4f: error at t = %g is %.3e\n',method_name,h_avg,tspan(2), ...
            abs(X_list(end)-solution01(tspan(2))));
    end

    % analytical solution goes on last so it doesn't get covered up
    t_fine = linspace(tspan(1),tspan(2),1000);
    true_plot = plot(t_fine,solution01(t_fine),'k--','LineWidth',2.5, ...
        'DisplayName','Analytical solution, $x(t) = \cos(t)$');

    % cos(t) runs through every corner, so leave room above it for the legend
    xlim(tspan);
    ylim([-1.35,1.5]);
    xlabel('Time, $t$ (-)','Interpreter','latex');
    ylabel('Solution, $x(t)$ (-)','Interpreter','latex');
    title(sprintf('%s Approximations of $\\dot{x} = -5x + 5\\cos(t) - \\sin(t)$',method_name), ...
        'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex','FontSize',16);
    legend([true_plot,num_plots],'Interpreter','latex','Location','north','NumColumns',2,'FontSize',13);

    exportgraphics(fig,file_name,'Resolution',300);
end
