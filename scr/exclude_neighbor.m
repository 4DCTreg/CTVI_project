function index_all = exclude_neighbor(pts_fix,pts_mov,spc)
 n_all = size(pts_fix,1);
 koef = repmat(spc, [n_all, 1]);
 index_fix = [];
for n_num = 1:n_all
    if ismember(n_num,index_fix)

        cur_fix = repmat(pts_fix(n_num,:), [n_all, 1]);

        pt_errs_phys = sqrt( sum((  (pts_fix - cur_fix).*koef  ).^2, 2) );
        index_cur = find(pt_errs_phys<2);
        index_cur(find(index_cur==n_num))=[];
        index_fix = [index_fix;index_cur];
    end

end

index_mov = [];
for n_num = 1:n_all
    if ismember(n_num,index_mov)
        cur_mov = repmat(pts_mov(n_num,:), [n_all, 1]);
        pt_errs_phys = sqrt( sum((  (pts_mov - cur_mov).*koef  ).^2, 2) );
        index_cur = find(pt_errs_phys<2);
        index_cur(find(index_cur==n_num))=[];
        index_mov = [index_mov;index_cur];
    end

end

index_all = [index_fix;index_mov];
index_all = unique(index_all);