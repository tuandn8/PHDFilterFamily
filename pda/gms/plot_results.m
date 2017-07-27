function handles= plot_results(model,truth,meas,est_pda,est_cphd)

[X_track,k_birth,k_death]= extract_tracks(truth.X,truth.track_list,truth.total_tracks);


%plot ground truths
limit= [ model.range_c(1,1) model.range_c(1,2) model.range_c(2,1) model.range_c(2,2) ];
figure; truths= gcf; hold on;

for i=1:truth.total_tracks
    Pt= X_track(:,k_birth(i):1:k_death(i),i); Pt=Pt([1 3],:);
    plot( Pt(1,:),Pt(2,:),'k-'); 
    plot( Pt(1,1), Pt(2,1), 'kx','MarkerSize',6);
    plot( Pt(1,(k_death(i)-k_birth(i)+1)), Pt(2,(k_death(i)-k_birth(i)+1)), 'k^','MarkerSize',6);
end
%est_xy = cell(truth.K);
for i = 1 : truth.K
    %est_xy{i} = est_cphd.X{i}([1 3],:);
    est_cphd_xy = est_cphd.X{i}([1 3], :);
    est_pda_xy  = est_pda.X{i}([1 3], :);
    plot( est_cphd_xy(1,:), est_cphd_xy(2,:),'ro', 'MarkerSize', 3);
    plot( est_pda_xy(1,:), est_pda_xy(2,:),'bx', 'MarkerSize', 8);
end


axis equal; axis(limit); title('Ground Truths');

%plot tracks and measurments in x/y
figure; tracking= gcf; hold on;

%plot x measurement
subplot(211); box on; 

for k=1:meas.K
    if ~isempty(meas.Z{k})
        hlined= line(k*ones(size(meas.Z{k},2),1),meas.Z{k}(1,:),'LineStyle','none','Marker','x','Markersize',5,'Color',0.7*ones(1,3));
    end   
end

%plot x track
for i=1:truth.total_tracks
    Px= X_track(:,k_birth(i):1:k_death(i),i); Px=Px([1 3],:);
    hline1= line(k_birth(i):1:k_death(i),Px(1,:),'LineStyle','-','Marker','none','LineWidth',1,'Color',0*ones(1,3));
end

%plot x estimate
for k=1:meas.K
    if ~isempty(est_cphd.X{k})
        P_cphd= est_cphd.X{k}([1 3],:);
        P_pda= est_pda.X{k}([1 3],:);
        hline2= line(k*ones(size(est_cphd.X{k},2),1),P_cphd(1,:),'LineStyle','none','Marker','o','Markersize',5,'Color','r');
        %hline3= line(k*ones(size(est_pda.X{k},2),1),P_pda(1,:),'LineStyle','none','Marker','.','Markersize',8,'Color','b');
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

%plot y estimate
for k=1:meas.K
    if ~isempty(est_cphd.X{k}),
        P_cphd= est_cphd.X{k}([1 3],:);
        P_pda= est_pda.X{k}([1 3],:);
        yhline2= line(k*ones(size(est_cphd.X{k},2),1),P_cphd(2,:),'LineStyle','none','Marker','o','Markersize',5,'Color','r');
        %yhline3= line(k*ones(size(est_pda.X{k},2),1),P_pda(2,:),'LineStyle','none','Marker','.','Markersize',8,'Color','b');
    end
end

subplot(211); xlabel('Time'); ylabel('x-coordinate (m)');
set(gca, 'XLim',[1 truth.K]); set(gca, 'YLim',model.range_c(1,:));
legend([ hline2 hline1 hlined],'Estimates cphd   ','True tracks','Measurements');

subplot(212); xlabel('Time'); ylabel('y-coordinate (m)');
set(gca, 'XLim',[1 truth.K]); set(gca, 'YLim',model.range_c(2,:));
legend([ yhline2 yhline1 yhlined],'Estimates cphd   ','True tracks','Measurements');

%plot error
ospa_cphd= zeros(truth.K,3);
ospa_pda= zeros(truth.K,3);
ospa_c= 50;
ospa_p= 1;
for k=1:meas.K
    [ospa_cphd(k,1), ospa_cphd(k,2), ospa_cphd(k,3)]= ospa_dist(get_comps(truth.X{k},[1 3]),get_comps(est_cphd.X{k},[1 3]),ospa_c,ospa_p);
    [ospa_pda(k,1),  ospa_pda(k,2),  ospa_pda(k,3)]= ospa_dist(get_comps(truth.X{k},[1 3]),get_comps(est_pda.X{k},[1 3]),ospa_c,ospa_p);
end

figure; ospa= gcf; hold on;
subplot(3,1,1); 
plot(1:meas.K,ospa_pda(:,1),'b', 1:meas.K,ospa_cphd(:,1), 'r'); grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Dist');
legend('PDA', 'CPHD', 'Location', 'northeast');

subplot(3,1,2); 
plot(1:meas.K,ospa_pda(:,2),'b',1:meas.K,ospa_cphd(:,2),'r'); grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Loc');
legend('PDA', 'CPHD', 'Location', 'northeast');
pda_loc = mean(ospa_pda(:,2));
cphd_loc = mean(ospa_cphd(:,2));

subplot(3,1,3); 
plot(1:meas.K,ospa_pda(:,3),'b', 1:meas.K,ospa_cphd(:,3),'r'); grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Card');
xlabel('Time');
legend('PDA', 'CPHD', 'Location', 'northeast');
xlabel('Time');

%plot cardinality
figure; cardinality= gcf; hold on;
subplot(2,1,1); 
stairs(1:meas.K,truth.N,'k'); hold on;
plot(1:meas.K,est_pda.N,'bo','MarkerSize', 5); 
grid on;
legend('PDA','Location', 'northeast');
subplot(2,1,2); 
stairs(1:meas.K,truth.N,'k'); hold on;
plot(1:meas.K, est_cphd.N,'ro','MarkerSize', 5);
grid on;
legend('CPHD', 'Location', 'northeast');

grid on;
legend(gca,'True','Estimated');
set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 max(truth.N)+1]);
xlabel('Time'); ylabel('Cardinality');

%return
%handles=[ truths tracking ospa cardinality cphd_loc pda_loc];
handles=[ ospa cardinality pda_loc cphd_loc];


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
