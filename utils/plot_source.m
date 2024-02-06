
% plot_source.m plots the estimated source signal across the given stimulus
% times
%
% INPUT:
%   ep_rec  : estimated source vector 
%   stim_on : on-times of stimuli (in seconds)
%   stim_off: off-times of stimuli (in seconds)
%   t_axis  : time axis (i.e., x-axis of the plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_source(ep_rec,stim_on,stim_off,t_axis)

fig = figure; 
ymin = min(ep_rec)-.4; ymax = max(ep_rec)+1.5; 
numReps = length(stim_on); % N = length(t_axis);

for i = 1:numReps
    p1 = fill([stim_on(i) stim_on(i) stim_off(i) stim_off(i)],...
        [ymin ymax ymax ymin],'r');
    set(p1,'facealpha',.1);
    set(p1,'EdgeColor','none');
    hold on;
end
p2 = plot(t_axis,ep_rec,'LineWidth',2,'Color','k');
xlim([0 t_axis(end)]); ylim([ymin ymax]); xlim([0 t_axis(end)]);
set(gca,'YTick',[]); set(gca,'FontSize',15); xlabel('Time (s)');
legend([p1 p2],{'Experimental Paradigm','Estimated Source Signal'},...
    'FontSize',15,'Location','northwest');
% thres = 1.4; hold on;
% plot(t_axis,thres*ones(1,N),'--k','DisplayName','Global Threshold');

% set(fig,'Units','Inches');
% pos = get(fig,'Position');
% set(fig,'PaperPositionMode','Auto','PaperUnits','Inches',...
%     'PaperSize',[pos(3), pos(4)])
% print(fig,'../Results/example_source_est','-dpdf','-r0')
