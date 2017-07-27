
pD=[0.7 0.75 0.80 0.85 0.90 0.95 1];
loc_error = zeros(2, size(pD,2));

nrun = size(pD, 2);
nrun = 1;

for i = 1:nrun
    model= gen_model;
    %model.P_D = pD(i);
    truth= gen_truth(model);
    meas=  gen_meas(model,truth);
    est_cphd= run_cphd_filter(model,truth, meas);
    est_pda = run_pda_filter(model, truth, meas);
    handles=  plot_results(model,truth,meas,est_pda, est_cphd);
    
    loc_error(1, i) = handles(3);
    loc_error(2, i) = handles(4);
end