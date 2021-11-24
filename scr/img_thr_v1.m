function v = img_thr_v1(vol, min_val, max_val, norm)
    v = vol;
    v(v > max_val) = max_val;
    v(v < min_val) = min_val;
    v = v - min_val+1;
    v = v * (norm / (max_val-min_val));
end