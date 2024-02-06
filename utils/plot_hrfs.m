
% plot_hrfs.m plots the estimated and true hemodynamic response functions
%
% INPUT:
%   final_sol: matrix containing the estimated HRFs
%   H        : matrix containing the true HRFs
%   u        : time axis of the HRFs
%   cmap     : colormap used for visualization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_hrfs(final_sol,H,u,cmap)

M = size(final_sol,2);

fig = figure;
for m = 1:M
    subplot(M,1,m);
    final_sol(1,m) = 0;
    plot(u,final_sol(:,m)/max(final_sol(:,m)),'LineWidth',2,...
        'Color',cmap(m,:)); hold on;
    plot(u,H(:,m,1)/max(H(:,m,1)),'--','LineWidth',2,...
        'Color',cmap(m,:)); hold on;
    set(gca,'YTick',[]);
    legend(['Estimated HRF, m = ', num2str(m)],...
        ['True HRF, m = ', num2str(m)],'FontSize',15);
    legend boxoff;
    ax = gca; set(gca,'FontSize',15);
    if m == 1
        set(gca,'XTick',[]);
        ax.Position = [0.01 0.705 0.98 0.285];
    elseif m == 2
        set(gca,'XTick',[]);
        ax.Position = [0.01 0.405 0.98 0.285];
    else
        set(gca,'XTick',0:8);
        ax.Position = [0.01 0.105 0.98 0.285];
    end
end
xlabel('Time (s)');

% set(fig,'Units','Inches');
% pos = get(fig,'Position');
% set(fig,'PaperPositionMode','Auto','PaperUnits','Inches',...
%     'PaperSize',[pos(3), pos(4)])
% print(fig,'../Results/example_hrf_est','-dpdf','-r0')
