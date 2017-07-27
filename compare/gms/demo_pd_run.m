% This is a script for the GM-CPHD/GM-PHD filter proposed in thesis
% Run simultaneously two algorithm to compare the results

% compare phd and cphd with variant detect prob

pD=[0.92 0.94 0.96 0.98 1];
error = zeros(2, 3, size(pD,2));

nrun = size(pD, 2);
nrun = 1;

for i = 1:nrun
    model= gen_model;
    model.P_D = pD(i);
    truth= gen_truth(model);
    meas=  gen_meas(model,truth);
    est_cphd =   run_cphd_filter(model,meas);
    est_phd  =   run_phd_filter(model, meas);
    handles_cphd_phd = plot_compare_results(model,truth,meas,est_phd,est_cphd);
    
    error(:,:,i) = reshape(handles_cphd_phd(6:11), 2, 3);
end


if nrun > 1
    figure; hold on;
    subplot(3,1,1);
    plot(pD,squeeze( error(1,1,:)), 'bo-', pD, squeeze(error(2,1,:)),'ro'); grid on; set(gca, 'XLim',[pD(1) pD(size(pD, 2))]); set(gca, 'YLim',[0 50]); ylabel('OSPA Dist');
    legend('PHD', 'CPHD');
    
    subplot(3,1,2);
    plot(pD, squeeze(error(1,2,:)), 'bo-', pD, squeeze(error(2,2,:)),'ro'); grid on; set(gca, 'XLim',[pD(1) pD(size(pD, 2))]); set(gca, 'YLim',[0 50]); ylabel('OSPA Loc');
    legend('PHD', 'CPHD');
    
    subplot(3,1,3);
    plot(pD, squeeze(error(1,3,:)), 'bo-', pD, squeeze(error(2,3,:)),'ro'); grid on; set(gca, 'XLim',[pD(1) pD(size(pD, 2))]); set(gca, 'YLim',[0 20]); ylabel('OSPA Card');
    legend('PHD', 'CPHD');
end