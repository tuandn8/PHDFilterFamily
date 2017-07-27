% This is a script for the GM-CPHD/GM-PHD filter proposed in thesis
% Run simultaneously two algorithm to compare the results


model= gen_model;
truth= gen_truth(model);
meas=  gen_meas(model,truth);
est_cphd =   run_cphd_filter(model,meas);
handles_cphd = plot_results(model,truth,meas,est_cphd);


