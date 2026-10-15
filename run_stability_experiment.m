% stability experiment, all four methods on the first test IVP (eqn 5)
% over t = [0,20] with a target step of 0.38 and then 0.45
% one figure per step size, four subplots stacked (one per method)
% for x' = -5x + ..., forward euler multiplies the error by (1 - 5h) every
% step, so it's only stable when |1 - 5h| < 1 (h < 0.4). 0.38 should be
% fine and 0.45 should blow up for the explicit methods

% figures for the report, saved one folder up so they don't end up in the repo
fig_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Assignment 3 Report Figures');
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

tspan = [0,20];
X0 = solution01(tspan(1));
h_ref_list = [0.38,0.45];
file_names = {'stability_h_ref_0_38.png','stability_h_ref_0_45.png'};

step_funcs = {@forward_euler_step,@backward_euler_step,@explicit_midpoint_step,@implicit_midpoint_step};
method_names = {'Forward Euler','Backward Euler','Explicit Midpoint','Implicit Midpoint'};
colors = [0.00,0.45,0.74; 0.47,0.67,0.19; 0.85,0.33,0.10; 0.49,0.18,0.56];

t_fine = linspace(tspan(1),tspan(2),1000);

for k = 1:2
    fig = figure('Color','w','Position',[50,50,1000,1300]);

    for m = 1:4
        [t_list,X_list,h_avg] = fixed_step_integration(@rate_func01,step_funcs{m},tspan,X0,h_ref_list(k));
        fprintf('h_ref = %.2f, %s: h = %.4f, max |x| = %.3e\n',h_ref_list(k),method_names{m},h_avg,max(abs(X_list)));

        subplot(4,1,m);
        hold on; box on; grid on;
        num_plot = plot(t_list,X_list,'o-','Color',colors(m,:),'MarkerFaceColor',colors(m,:), ...
            'MarkerSize',5,'LineWidth',1.5);
        true_plot = plot(t_fine,solution01(t_fine),'k--','LineWidth',2);

        xlim(tspan);
        title(sprintf('%s, $h$ = %.4f',method_names{m},h_avg),'Interpreter','latex');
        xlabel('Time, $t$ (-)','Interpreter','latex');
        ylabel('Solution, $x(t)$ (-)','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex','FontSize',14);

        if max(abs(X_list)) < 2
            % stable, cos(t) fills every corner so give it headroom and
            % put the legend in one row across the top
            ylim([-1.3,2.2]);
            legend([true_plot,num_plot],{'Analytical solution',method_names{m}}, ...
                'Interpreter','latex','Location','north','Orientation','horizontal','FontSize',12);
        else
            % blew up, it starts out fine so the bottom left is empty
            y_pad = 0.1*(max(X_list)-min(X_list));
            ylim([min(X_list)-y_pad,max(X_list)+y_pad]);
            legend([true_plot,num_plot],{'Analytical solution',method_names{m}}, ...
                'Interpreter','latex','Location','southwest','FontSize',12);
        end
    end

    sgtitle(sprintf('Explicit vs. Implicit Stability with a Target Step Size of $h_{ref}$ = %.2f',h_ref_list(k)), ...
        'Interpreter','latex','FontSize',18);
    exportgraphics(fig,fullfile(fig_dir,file_names{k}),'Resolution',300);
end
