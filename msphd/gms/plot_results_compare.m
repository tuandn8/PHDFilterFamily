function handles= plot_results_compare(model,truth,meas,est_phd, est_msphd)

[X_track,k_birth,k_death]= extract_tracks(truth.X,truth.track_list,truth.total_tracks);


%plot ground truths
limit= [ model.range_c(1,1) model.range_c(1,2) model.range_c(2,1) model.range_c(2,2) ];
figure; truths= gcf; hold on;
for i=1:truth.total_tracks
    Pt= X_track(:,k_birth(i):1:k_death(i),i); Pt=Pt([1 3],:);
    plot( Pt(1,:),Pt(2,:),'r-'); 
    plot( Pt(1,1), Pt(2,1), 'ro','MarkerSize',6);
    plot( Pt(1,(k_death(i)-k_birth(i)+1)), Pt(2,(k_death(i)-k_birth(i)+1)), 'r^','MarkerSize',6);
end

plot(model.S(1,1),model.S(2,1), 'v','Color','r','MarkerSize',10);
plot(model.S(1,2),model.S(2,2), 'v', 'Color', 'r','MarkerSize',10);

for i = 1 : truth.K
    %est_xy{i} = est_cphd.X{i}([1 3],:);
    est_phd_xy = est_msphd.X{i}([1 3], :);
    plot( est_phd_xy(1,:), est_phd_xy(2,:),'bo', 'MarkerSize', 3);
end

grid on; axis equal; axis(limit); title('Ground Truths'); 

%plot tracks and measurments in x/y
figure; tracking= gcf; hold on;

%plot x measurement
subplot(211); box on; 

for k=1:meas.K
    if ~isempty(meas.Z{k})
        hlined= line(k*ones(size(meas.Z{k,1},2),1),meas.Z{k,1}(1,:),'LineStyle','none','Marker','x','Markersize',5,'Color',0.7*ones(1,3));
        %hlined2= line(k*ones(size(meas.Z{k,2},2),1),meas.Z{k,2}(1,:),'LineStyle','none','Marker','x','Markersize',5,'Color','m');
    end   
end

%plot x track
for i=1:truth.total_tracks
    Px= X_track(:,k_birth(i):1:k_death(i),i); Px=Px([1 3],:);
    hline1= line(k_birth(i):1:k_death(i),Px(1,:),'LineStyle','-','Marker','none','LineWidth',1,'Color',0*ones(1,3));
end

%plot x estimate phd, msphd
for k=1:meas.K
    if ~isempty(est_phd.X{k})
        P_phd= est_phd.X{k}([1 3],:);
        P_msphd= est_msphd.X{k}([1 3],:);
        hline2= line(k*ones(size(est_phd.X{k},2),1),P_phd(1,:),'LineStyle','none','Marker','.','Markersize',8,'Color','r');
        hline3= line(k*ones(size(est_msphd.X{k},2),1),P_msphd(1,:),'LineStyle','none','Marker','o','Markersize',5,'Color','b');
    end
end

%plot y measurement
subplot(212); box on;
    
for k=1:meas.K
    if ~isempty(meas.Z{k})
        yhlined= line(k*ones(size(meas.Z{k},2),1),meas.Z{k}(2,:),'LineStyle','none','Marker','x','Markersize',5,'Color',0.7*ones(1,3));
    end
end

%plot y track
for i=1:truth.total_tracks
        Py= X_track(:,k_birth(i):1:k_death(i),i); Py=Py([1 3],:);
        yhline1= line(k_birth(i):1:k_death(i),Py(2,:),'LineStyle','-','Marker','none','LineWidth',1,'Color',0*ones(1,3));
end

%plot y estimate phd, msphd
for k=1:meas.K
    if ~isempty(est_phd.X{k}),
        P_phd= est_phd.X{k}([1 3],:);
        P_msphd= est_msphd.X{k}([1 3],:);
        yhline2= line(k*ones(size(est_phd.X{k},2),1),P_phd(2,:),'LineStyle','none','Marker','.','Markersize',8,'Color','r');
        yhline3= line(k*ones(size(est_msphd.X{k},2),1),P_msphd(2,:),'LineStyle','none','Marker','o','Markersize',5,'Color','b');
    end
end

subplot(211); xlabel('Time', 'FontSize',14); ylabel('x-coordinate (m)', 'FontSize',14);
set(gca, 'XLim',[1 truth.K]); set(gca, 'YLim',model.range_c(1,:));
x_legend=legend([hline3 hline2 hline1 hlined],'MS-PHD    ','PHD      ','True tracks','Measurements');
set(x_legend, 'FontSize',14);

subplot(212); xlabel('Time', 'FontSize',14); ylabel('y-coordinate (m)', 'FontSize',14);
set(gca, 'XLim',[1 truth.K]); set(gca, 'YLim',model.range_c(2,:));
y_legend=legend([yhline3   yhline2 yhline1 yhlined],'MS-PHD    ','PHD      ','True tracks','Measurements');
set(y_legend, 'FontSize',14);

%plot error
ospa_phd= zeros(truth.K,3);
ospa_msphd= zeros(truth.K,3);
ospa_c= 100;
ospa_p= 1;
for k=1:meas.K
    [ospa_phd(k,1),   ospa_phd(k,2),   ospa_phd(k,3)]  = ospa_dist(get_comps(truth.X{k},[1 3]),get_comps(est_phd.X{k},[1 3]),ospa_c,ospa_p);
    [ospa_msphd(k,1), ospa_msphd(k,2), ospa_msphd(k,3)]= ospa_dist(get_comps(truth.X{k},[1 3]),get_comps(est_msphd.X{k},[1 3]),ospa_c,ospa_p);
end

figure; ospa= gcf; hold on;
subplot(3,1,1); 
plot(1:meas.K,ospa_phd(:,1),'r', 1:meas.K,ospa_msphd(:,1), 'b'); grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Dist','FontSize',14);
dist_legend=legend('PHD', 'MSPHD', 'Location', 'northeast');
set(dist_legend, 'FontSize',14);

subplot(3,1,2); 
plot(1:meas.K,ospa_phd(:,2),'r',1:meas.K,ospa_msphd(:,2),'b'); grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Loc','FontSize',14);
loc_legend=legend('PHD', 'MSPHD', 'Location', 'northeast');
set(loc_legend, 'FontSize',14);

subplot(3,1,3); 
plot(1:meas.K,ospa_phd(:,3),'r', 1:meas.K,ospa_msphd(:,3),'b'); grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Card','FontSize',14);
xlabel('Time', 'FontSize',14);
card_legend=legend('PHD', 'MSPHD', 'Location', 'northeast');
set(card_legend, 'FontSize',14);

dist_error = mean(ospa_phd(:,1)), 

%plot cardinality
figure; cardinality= gcf; 
subplot(2,1,1); box on; hold on;
stairs(1:meas.K,truth.N,'k'); 
plot(1:meas.K,est_phd.N,'ko','Color', 'r');
grid on;
phd_card_legend=legend(gca,'True','PHD Estimated');
set(phd_card_legend, 'FontSize',14);
set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 max(truth.N)+1]);
xlabel('Time', 'FontSize',14); ylabel('Cardinality phd', 'FontSize',14);

subplot(2,1,2); box on; hold on;
stairs(1:meas.K,truth.N,'k'); 
plot(1:meas.K,est_msphd.N,'ko','Color','b');

grid on;
msphd_card_legend=legend(gca,'True','MSPHD Estimated');
set(msphd_card_legend, 'FontSize',14);
set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 max(truth.N)+1]);
xlabel('Time', 'FontSize',14); ylabel('Cardinality msphd', 'FontSize',14);

%return
handles=[ truths tracking ospa cardinality ];


function [X_track,k_birth,k_death]= extract_tracks(X,track_list,total_tracks)

K= size(X,1); 
x_dim= size(X{K},1); 
k=K-1; while x_dim==0, x_dim= size(X{k},1); k= k-1; end;
X_track= zeros(x_dim,K,total_tracks);
k_birth= zeros(total_tracks,1);
k_death= zeros(total_tracks,1);

max_idx= 0;
for k=1:K
    if ~isempty(X{k}),
        X_track(:,k,track_list{k})= X{k};
    end;
    if max(track_list{k})> max_idx, %new target born?
        idx= find(track_list{k}> max_idx);
        k_birth(track_list{k}(idx))= k;
    end;
    if ~isempty(track_list{k}), max_idx= max(track_list{k}); end;
    k_death(track_list{k})= k;
end;

function Xc= get_comps(X,c)

if isempty(X)
    Xc= [];
else
    Xc= X(c,:);
end
