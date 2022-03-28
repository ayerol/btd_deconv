clear; close all; hold on;

cmap = lines(7);

med_iou = []; std_iou = [];
med_iou_fixedHrf = []; std_iou_fixedHrf = []; 

for snr_val = [0 5 10 20]

    load(['Results/ious_' num2str(snr_val) 'dB.mat']);

    load(['Results/ious_fixedHrf' num2str(snr_val) 'dB.mat']);

    std_iou = cat(2,std_iou,std(ious));

    med_iou = cat(2,med_iou,median(ious));

    std_iou_fixedHrf = cat(2,std_iou_fixedHrf,std(ious_fixedHrf));

    med_iou_fixedHrf = cat(2,med_iou_fixedHrf,median(ious_fixedHrf));

end

x_axis = 1:4; figure;
p1 = errorbar(x_axis,med_iou,std_iou,'-s','MarkerSize',15,'Color',...
    cmap(6,:),'LineWidth',2,'DisplayName','SNR = 0 dB'); hold on;
p2 = errorbar(x_axis,med_iou_fixedHrf,std_iou_fixedHrf,'*-.',...
    'MarkerSize',15,'Color',cmap(4,:),'LineWidth',2,'DisplayName',...
    'SNR = 0 dB'); hold on;
p3 = plot(0:5,4*ones(1,6),'--','Color',cmap(5,:),'LineWidth',2);
p4 = plot(0,0,'Marker','s','LineWidth',2,'MarkerSize',10,...
    'Color',cmap(6,:),'LineStyle','none'); hold on;
p5 = plot(0,0,'Marker','*','LineWidth',2,'MarkerSize',10,...
    'Color',cmap(4,:),'LineStyle','none'); hold on;

set(gca,'FontSize',15); ylim([0 5.5]); set(gca,'yTick',1:4);
xlim([0.5 4.5]);  xticks(1:4); xticklabels({'0','5','10','20'}); 
xlabel('SNR (dB)'); label_h = ylabel('IoU (s)');
label_h.Position(1) = 0.25;

legend([p4,p5,p3],{'Our method','Fixed HRF','Ideal Value'},'FontSize',15);

grid on; box on