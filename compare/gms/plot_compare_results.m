function handles = plot_compare_results( model,truth,meas,est_phd,est_cphd )
%PLOT_COMPARE_RESULTS Summary of this function goes here
%   Detailed explanation goes here

[X_track,k_birth,k_death] = extract_tracks(truth.X, truth.track_list,truth.total_tracks);

%plot ground truths
limit= [ model.range_c(1,1) model.range_c(1,2) model.range_c(2,1) model.range_c(2,2) ];
figure; truths= gcf; hold on;
for i=1:truth.total_tracks
    Pt= X_track(:,k_birth(i):1:k_death(i),i); Pt=Pt([1 3],:);
    plot( Pt(1,:),Pt(2,:),'r-'); 
    plot( Pt(1,1), Pt(2,1), 'ro','MarkerSize',6);
    plot( Pt(1,(k_death(i)-k_birth(i)+1)), Pt(2,(k_death(i)-k_birth(i)+1)), 'r^','MarkerSize',6);
end
for i = 1 : truth.K
    %est_xy{i} = est_cphd.X{i}([1 3],:);
    est_phd_xy = est_cphd.X{i}([1 3], :);
    plot( est_phd_xy(1,:), est_phd_xy(2,:),'bo', 'MarkerSize', 3);
end
axis equal; axis(limit); title('Ground Truths');

% plot track and clutter in Oxy
figure; datatracking= gcf; hold on;
plot_clear_est = 0;
for k = 1:meas.K
    if plot_clear_est
        clf;
        plot(meas.Z{k}(1,:), meas.Z{k}(2,:), 'kx');   % measurement at k index
        hold on;
        %     P_track = squeeze(X_track([1 3],k,truth.track_list{k}));
        %     plot(P_track(1,:), P_track(2,:), 'ro', 'MarkerSize', 5);
        
        for i = 1:truth.total_tracks
            if k < k_death(i) && k > k_birth(i)
                Pt= X_track(:,k_birth(i):1:k,i); Pt=Pt([1 3],:);
                hold on;
                plot( Pt(1,:),Pt(2,:),'r-');  % truth track from birth to k
                plot( Pt(1,1), Pt(2,1), 'ko','MarkerSize',6); % track begin position
                plot( Pt(1,(k-k_birth(i)+1)), Pt(2,(k-k_birth(i)+1)), 'r^','MarkerSize',6); % track end position
            end
        end
        for i = 1 : est_cphd.N
            Et=est_cphd.X{k}([1 3], :);
            plot(Et(1,:), Et(2, :), 'go');
        end
        %pause(0.5);
    else
        measure= plot(meas.Z{k}(1,:), meas.Z{k}(2,:), 'kx');   % measurement at k index
        hold on;
        %     P_track = squeeze(X_track([1 3],k,truth.track_list{k}));
        %     plot(P_track(1,:), P_track(2,:), 'ro', 'MarkerSize', 5);
        
        for i = 1:truth.total_tracks
            if k < k_death(i) && k > k_birth(i)
                %Pt= X_track(:,k_birth(i):1:k,i); Pt=Pt([1 3],:);
                Pt= X_track(:,k,i); Pt=Pt([1 3],:);             
                plot( Pt(1,:),Pt(2,:),'r-', 'MarkerSize', 2);  % truth track from birth to k
                hold on;
                %plot( Pt(1,1), Pt(2,1), 'ko','MarkerSize',6); % track begin position
                %plot( Pt(1,(k-k_birth(i)+1)), Pt(2,(k-k_birth(i)+1)), 'r^','MarkerSize',6); % track end position
            end
        end
        for i = 1 : est_cphd.N
            Et=est_cphd.X{k}([1 3], :);
            plot(Et(1,:), Et(2, :), 'bo','MarkerSize', 2);
            hold on
        end
        %pause(0.5);
        delete(measure);
    end
end
axis equal; axis(limit); title('Target tracking');

% plot card distribution

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
    hline1= line(k_birth(i):1:k_death(i),Px(1,:),'LineStyle','-','Marker','none','LineWidth',1,'Color','r');
end

%plot x estimate
% for k=1:meas.K
%     if ~isempty(est_phd.X{k})
%         P_phd= est_phd.X{k}([1 3],:);
%         hline2= line(k*ones(size(est_phd.X{k},2),1),P_phd(1,:),'LineStyle','none','Marker','+','Markersize',5,'Color','b');
%     end
% end

for k=1:meas.K
    if ~isempty(est_cphd.X{k})
        P_cphd= est_cphd.X{k}([1 3],:);
        hline3= line(k*ones(size(est_cphd.X{k},2),1),P_cphd(1,:),'LineStyle','none','Marker','o','Markersize',5,'Color','b');
    end
end


%plot y measurement
subplot(212); box on;
    
for k=1:meas.K
    if ~isempty(meas.Z{k})
        yhlined= line(k*ones(size(meas.Z{k},2),1),meas.Z{k}(2,:),'LineStyle','none','Marker','x','Markersize',3,'Color',0.7*ones(1,3));
    end
end

%plot y track
for i=1:truth.total_tracks
        Py= X_track(:,k_birth(i):1:k_death(i),i); Py=Py([1 3],:);
        yhline1= line(k_birth(i):1:k_death(i),Py(2,:),'LineStyle','-','Marker','none','LineWidth',1,'Color','r');
end

%plot y estimate
% for k=1:meas.K
%     if ~isempty(est_phd.X{k}),
%         P_phd= est_phd.X{k}([1 3],:);
%         yhline2= line(k*ones(size(est_phd.X{k},2),1),P_phd(2,:),'LineStyle','none','Marker','+','Markersize',5,'Color','b');
%     end
% end

for k=1:meas.K
    if ~isempty(est_cphd.X{k}),
        P_cphd= est_cphd.X{k}([1 3],:);
        yhline3= line(k*ones(size(est_cphd.X{k},2),1),P_cphd(2,:),'LineStyle','none','Marker','o','Markersize',3,'Color','b');
    end
end

subplot(211); xlabel('Time','FontSize',14); ylabel('x-coordinate (m)','FontSize',14);
set(gca, 'XLim',[1 truth.K]); set(gca, 'YLim',model.range_c(1,:));
%x_legend=legend([hline3 hline2 hline1 hlined],'Estimates cphd         ','Estimates phd         ','True tracks','Measurements', 'Location', 'northeast');
x_legend=legend([hline3 hline1 hlined],'Estimates cphd         ','True tracks','Measurements', 'Location', 'northeast');
set(x_legend,'FontSize',14);

subplot(212); xlabel('Time','FontSize',14); ylabel('y-coordinate (m)','FontSize',14);
set(gca, 'XLim',[1 truth.K]); set(gca, 'YLim',model.range_c(2,:));
%y_legend=legend([yhline3 yhline2 yhline1 yhlined],'Estimates cphd        ', 'Estimates phd         ','True tracks','Measurements', 'Location', 'northeast');
y_legend=legend([yhline3 yhline1 yhlined],'Estimates cphd        ', 'True tracks','Measurements', 'Location', 'northeast');
set(y_legend,'FontSize',14);
%plot error
ospa_phd_vals= zeros(truth.K,3);
ospa_cphd_vals= zeros(truth.K,3);
ospa_c= 50;
ospa_p= 1;
for k=1:meas.K
    [ospa_phd_vals(k,1), ospa_phd_vals(k,2), ospa_phd_vals(k,3)]= ospa_dist(get_comps(truth.X{k},[1 3]),get_comps(est_phd.X{k},[1 3]),ospa_c,ospa_p);
    [ospa_cphd_vals(k,1), ospa_cphd_vals(k,2), ospa_cphd_vals(k,3)]= ospa_dist(get_comps(truth.X{k},[1 3]),get_comps(est_cphd.X{k},[1 3]),ospa_c,ospa_p);
end

figure; ospa= gcf; hold on;
subplot(3,1,1); 
plot(1:meas.K,ospa_phd_vals(:,1),'b', 1:meas.K,ospa_cphd_vals(:,1), 'r'); grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Dist','FontSize',14);
xlabel('Time','FontSize',14);
dist_legend=legend('PHD', 'CPHD', 'Location', 'northeast');
set(dist_legend,'FontSize',14);

subplot(3,1,2); 
plot(1:meas.K,ospa_phd_vals(:,2),'b',1:meas.K,ospa_cphd_vals(:,2),'r'); grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Loc','FontSize',14);
xlabel('Time','FontSize',14);
loc_legend=legend('PHD', 'CPHD', 'Location', 'northeast');
set(loc_legend,'FontSize',14);

subplot(3,1,3); 
plot(1:meas.K,ospa_phd_vals(:,3),'b', 1:meas.K,ospa_cphd_vals(:,3),'r'); grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Card','FontSize',14);
xlabel('Time','FontSize',14);
card_legend=legend('PHD', 'CPHD', 'Location', 'northeast');
set(card_legend,'FontSize',14);

dist_error = [mean(ospa_phd_vals(:,1)), mean(ospa_cphd_vals(:,1))];
loc_error = [mean(ospa_phd_vals(:,2)), mean(ospa_cphd_vals(:,2))];
card_error = [mean(ospa_phd_vals(:,3)), mean(ospa_cphd_vals(:,3))];

%plot cardinality
figure; cardinality= gcf; hold on;
subplot(2,1,1); 
stairs(1:meas.K,truth.N,'r'); hold on;
plot(1:meas.K,est_phd.N,'bo','MarkerSize', 3); 
grid on;
%phd_card_legend=legend('PHD','Location', 'northeast');
phd_card_legend=legend(gca,'True','Estimated PHD');
set(phd_card_legend,'FontSize',14);

subplot(2,1,2); 
stairs(1:meas.K,truth.N,'r'); hold on;
plot(1:meas.K, est_cphd.N,'bo', 'MarkerSize', 3);
grid on;
%cphd_card_legend=legend('CPHD', 'Location', 'northeast');
cphd_card_legend=legend(gca,'True','Estimated CPHD');
set(cphd_card_legend,'FontSize',14);
% 
% grid on;
% legend(gca,'True','Estimated');
% set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 max(truth.N)+1]);
% xlabel('Time'); ylabel('Cardinality');

%return
handles=[ truths datatracking tracking ospa cardinality dist_error loc_error card_error];
end



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
        idx= track_list{k}> max_idx;
        k_birth(track_list{k}(idx))= k;
    end;
    if ~isempty(track_list{k}), max_idx= max(track_list{k}); end;
    k_death(track_list{k})= k;
end;
end


function Xc= get_comps(X,c)

if isempty(X)
    Xc= [];
else
    Xc= X(c,:);
end
end


