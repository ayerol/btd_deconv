function plot_source(ep_rec,stim_on,stim_off,t_axis)

figure; 
ymin = min(ep_rec)-.2; ymax = max(ep_rec)+.5; 
numReps = length(stim_on); N = length(t_axis);

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
thres = 1.4; hold on;
plot(t_axis,thres*ones(1,N),'--k','DisplayName','Global Threshold');