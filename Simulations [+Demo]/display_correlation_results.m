clear; close all; hold on;

cmap = lines(7);

mean_corr = []; std_corr = [];
mean_corr_fixedHrf = []; std_corr_fixedHrf = []; 
mean_corr_lowest_cost_Hrf = []; std_corr_lowest_cost_Hrf = []; 

for snr_val = [-20 -10 0 10 20]

    load(['../Results/corrs_' num2str(snr_val) 'dB.mat']);

    load(['../Results/corrs_fixedHrf_' num2str(snr_val) 'dB.mat']);

    load(['../Results/corrs_lowest_cost_Hrf_' num2str(snr_val) 'dB.mat']);

    std_corr = cat(2,std_corr,std(corrs));

    mean_corr = cat(2,mean_corr,mean(corrs));

    std_corr_fixedHrf = cat(2,std_corr_fixedHrf,std(corrs_fixedHrf));

    mean_corr_fixedHrf = cat(2,mean_corr_fixedHrf,mean(corrs_fixedHrf));

    std_corr_lowest_cost_Hrf = cat(2,std_corr_lowest_cost_Hrf,...
        std(corrs_lowest_cost_Hrf));

    mean_corr_lowest_cost_Hrf = cat(2,mean_corr_lowest_cost_Hrf,...
        mean(corrs_lowest_cost_Hrf));

end

x_axis = 1:5; fig = figure;
p2 = errorbar(x_axis,mean_corr_fixedHrf,std_corr_fixedHrf,'d--',...
    'MarkerSize',15,'Color',cmap(4,:),'LineWidth',2,'DisplayName',...
    'SNR = 0 dB'); hold on;
p3 = errorbar(x_axis,mean_corr_lowest_cost_Hrf,std_corr_lowest_cost_Hrf,...
    '*-.','MarkerSize',15,'Color',cmap(7,:),'LineWidth',2,'DisplayName',...
    'SNR = 0 dB'); hold on;
p1 = errorbar(x_axis,mean_corr,std_corr,'-s','MarkerSize',15,'Color',...
    cmap(1,:),'LineWidth',2,'DisplayName','SNR = 0 dB'); hold on;
% p3 = plot(0:5,4*ones(1,6),'--','Color',cmap(5,:),'LineWidth',2);
p4 = plot(0,0,'Marker','s','LineWidth',2,'MarkerSize',10,...
    'Color',cmap(1,:),'LineStyle','none'); hold on;
p5 = plot(0,0,'Marker','d','LineWidth',2,'MarkerSize',10,...
    'Color',cmap(4,:),'LineStyle','none'); hold on;
p6 = plot(0,0,'Marker','*','LineWidth',2,'MarkerSize',10,...
    'Color',cmap(7,:),'LineStyle','none'); hold on;

set(gca,'FontSize',15); ylim([-.2 1.2]); set(gca,'yTick',0:.2:1);
xlim([0.5 5.5]);  xticks(1:5); xticklabels({'-20','-10','0','10','20'}); 
xlabel('SNR (dB)'); label_h = ylabel('Correlation Score');
label_h.Position(1) = 0.12;

legend([p4,p6,p5],{'Our Method','Lowest-Cost Solution','Fixed HRF'},...
    'Orientation','horizontal','FontSize',15);

grid on; box on;

% set(fig,'Units','Inches');
% pos = get(fig,'Position');
% set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig,'Results/simulation_ep','-dpdf','-r0');