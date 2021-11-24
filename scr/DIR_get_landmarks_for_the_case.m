function pts = DIR_get_landmarks_for_the_case(caseN, basepth) 
    cn = num2str(caseN);
    if caseN < 10
        path = [basepth, ...
            '/Subject_0', cn, '/'];
    else
        path = [basepth, ...
            '/Subject_', cn, '/'];
    end
        
    pts = struct();
    pts.extreme.in = read_DIR_points_file([path, 'T01.txt']);
    pts.extreme.ex = read_DIR_points_file([path, 'T05.txt']);
    
%     path2 = [basepth, ...
%             '4dCT/Case', cn, 'Pack/Case', cn, 'Pack/Sampled4D/'];
%     for i = 1 : 6
%         ni = num2str(i - 1);
%         pts.smp{i}.pts = read_DIR_points_file([path2, 'case', cn, ...
%             '_4D-75_T', ni, '0.txt']);
%     end
end