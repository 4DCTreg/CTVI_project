clear
points = [0 100  150  500  550   600   700   1000  2000 4000 ;
          0 0.27 0.31 0.32 0.319 0.318 0.315 0.31  0.29 0.285 ];
points(2,:) = points(2,:)-0.02;

points_homrf = [0 100  150  500   550   600   700   1000   2000 4000 ;
                0 0.28 0.31 0.315 0.314 0.312 0.310 0.305  0.29 0.285 ];

points_homrf(2,:) = points_homrf(2,:)-0.015;

a_homrf = points_homrf(1,:);
b_homrf = points_homrf(2,:);
values_homrf = spcrv(points_homrf,3);
p_points_homrf = [100, 500,1000,2000,];
p_value_homrf = zeros(1,4);
for id = 1:length(p_points_homrf)
    diff = abs(values_homrf(1,:)-p_points_homrf(id));
    n_id = find(diff == min(diff));
    p_value_homrf(id) = values_homrf(2,n_id(1));
    
end

a = points(1,:);
b = points(2,:);
values = spcrv(points,3);
p_points = [100, 500,1000,2000];
p_value = zeros(1,4);
for id = 1:length(p_points)
    diff = abs(values(1,:)-p_points(id));
    n_id = find(diff == min(diff));
    p_value(id) = values(2,n_id(1));
    
end
semilogx(values(1,:),values(2,:),'Color',[0.85,0.33,0.10],'LineWidth',2.5);

hold on
semilogx(values_homrf(1,:),values_homrf(2,:),'Color',[0.00,0.45,0.74],'LineWidth',2.5);
hold on
semilogx(p_points_homrf,p_value_homrf,'^','Color',[0.00,0.45,0.74],'LineWidth',2.5);
hold on
semilogx(p_points,p_value,'o','Color',[0.85,0.33,0.10],'LineWidth',2.5);
set(gca,'XLim',[500 5000])
grid on
legend('isoPTV','HOMRF');

%% legend
figure
plot([1,10],[1,1],'Color',[0.00,0.45,0.74],'LineWidth',2.5)
hold on 
plot(4.5,1,'^','Color',[0.00,0.45,0.74],'LineWidth',2.5)
hold on 
plot([1,10],[2,2],'Color',[0.85,0.33,0.10],'LineWidth',2.5)
hold on 
plot(4.5,2,'o','Color',[0.85,0.33,0.10],'LineWidth',2.5)
set(gca,'YLim',[0 3])



% hold on, plot(values(1,:),values(2,:)), hold off