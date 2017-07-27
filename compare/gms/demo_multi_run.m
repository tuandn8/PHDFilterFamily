% This is a script for the GM-CPHD/GM-PHD filter proposed in thesis
% Run simultaneously two algorithm to compare the results

% compare phd and cphd with mutiple run

nrun =20;
phd_card = zeros(100, nrun); % 100 = model.K
cphd_card = zeros(100, nrun); 

plot_flag= 0;

for i = 1:nrun
    model= gen_model;
    truth= gen_truth(model);
    meas=  gen_meas(model,truth);
    est_cphd =   run_cphd_filter(model,meas);
    est_phd  =   run_phd_filter(model, meas);
    if plot_flag
        handles_cphd_phd = plot_compare_results(model,truth,meas,est_phd,est_cphd);
    end
    
    phd_card(:, i) = est_phd.N;
    cphd_card(:, i) = est_cphd.N;
end

phd_cov_card = zeros(size(phd_card,1), 1);
cphd_cov_card = zeros(size(cphd_card,1), 1);
phd_mean_card = zeros(size(cphd_card,1), 1);
cphd_mean_card = zeros(size(cphd_card,1), 1);

for i = 1:size(phd_card,1)
    phd_cov_card(i) = std(phd_card(i,:));
    cphd_cov_card(i) = std(cphd_card(i,:));
    
    phd_mean_card(i) = mean(phd_card(i,:));
    cphd_mean_card(i) = mean(cphd_card(i,:));
end


figure; cardinality= gcf; hold on;
subplot(2,1,1); 
stairs(1:meas.K,truth.N,'k'); hold on;
plot(1:meas.K,phd_mean_card,'b.', 1:meas.K, phd_mean_card + phd_cov_card, 'b-', 1:meas.K, phd_mean_card - phd_cov_card,  'b-','MarkerSize', 5); 
grid on;
phd_legend=legend('PHD','Location', 'northeast');
set(phd_legend, 'FontSize',14);
subplot(2,1,2); 
stairs(1:meas.K,truth.N,'k'); hold on;
plot(1:meas.K,cphd_mean_card,'r.', 1:meas.K,cphd_mean_card + cphd_cov_card, 'r-', 1:meas.K, cphd_mean_card - cphd_cov_card,  'r-','MarkerSize', 5); 
grid on;
cphd_legend=legend('CPHD', 'Location', 'northeast');
set(cphd_legend, 'FontSize',14);
% 
% if nrun > 1
%     figure; hold on;
%     subplot(3,1,1);
%     plot(pD,squeeze( error(1,1,:)), 'bo-', pD, squeeze(error(2,1,:)),'ro'); grid on; set(gca, 'XLim',[pD(1) pD(size(pD, 2))]); set(gca, 'YLim',[0 50]); ylabel('OSPA Dist');
%     legend('PHD', 'CPHD');
%     
%     subplot(3,1,2);
%     plot(pD, squeeze(error(1,2,:)), 'bo-', pD, squeeze(error(2,2,:)),'ro'); grid on; set(gca, 'XLim',[pD(1) pD(size(pD, 2))]); set(gca, 'YLim',[0 50]); ylabel('OSPA Loc');
%     legend('PHD', 'CPHD');
%     
%     subplot(3,1,3);
%     plot(pD, squeeze(error(1,3,:)), 'bo-', pD, squeeze(error(2,3,:)),'ro'); grid on; set(gca, 'XLim',[pD(1) pD(size(pD, 2))]); set(gca, 'YLim',[0 20]); ylabel('OSPA Card');
%     legend('PHD', 'CPHD');
% end