
    model= gen_model;
    %model.P_D = pD(i);
    truth= gen_truth(model);
    meas=  gen_meas(model,truth);
    est_phd = run_phd_filter(model, meas);
    est_msphd=   run_msphd_filter(model,meas);
    handles= plot_results_compare(model, truth, meas, est_phd, est_msphd);
    %handles= plot_results(model,truth,meas,est);
