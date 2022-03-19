function iou = compute_iou(ep_rec,t_axis,stim_on,stim_off)

pks = findpeaks(ep_rec);

thres = pks(1)-0.15;
binarized_ep = ep_rec;
binarized_ep(binarized_ep <= thres) = 0;
binarized_ep(binarized_ep > thres) = 1;

iou = zeros(1,length(stim_on));
est_start = t_axis(find(diff(binarized_ep)==1)+1);
est_end = t_axis(find(diff(binarized_ep)==-1)+1);

while length(est_start) ~= length(stim_on) % Adjust the threshold if 
                                                         % necessary
    if length(est_start) < length(stim_on)

        thres = thres - .1;

    else

        thres = thres + .01;

    end

    binarized_ep = ep_rec;
    binarized_ep(binarized_ep <= thres) = 0;
    binarized_ep(binarized_ep > thres) = 1;

    est_start = t_axis(find(diff(binarized_ep)==1)+1);
    est_end = t_axis(find(diff(binarized_ep)==-1)+1);

    if length(est_start) < 8 % If the threshold is too low

        binarized_ep = imbinarize(ep_rec);

        est_start = t_axis(find(diff(binarized_ep)==1)+1);
        est_end = t_axis(find(diff(binarized_ep)==-1)+1);
    
    end

end

for i = 1:length(stim_on)

    true_ep = [stim_on(i),stim_off(i)];
    est_ep = [est_start(i),est_end(i)];

    iou(i) = compute_iou_single_stimulus(true_ep,est_ep);

end

end