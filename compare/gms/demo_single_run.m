% compare phd and cphd


model= gen_model;
truth= gen_truth(model);
meas=  gen_meas(model,truth);
est_cphd =   run_cphd_filter(model,meas);
est_phd  =   run_phd_filter(model, meas);
handles_cphd_phd = plot_compare_results(model,truth,meas,est_phd,est_cphd);
