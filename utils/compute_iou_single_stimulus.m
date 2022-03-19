function iou = compute_iou_single_stimulus(rng1,rng2)

    min_val = min(rng1(1),rng2(1));

    if min_val < 0

        rng1 = rng1 + abs(min_val);
        rng2 = rng2 + abs(min_val);

    end
    
    union = max(rng1(2),rng2(2)) - min(rng1(1),rng2(1));
    intersection = max(0,min(rng1(2),rng2(2)) - max(rng1(1),rng2(1)));
    iou = intersection/union;

end