function meas= gen_meas(model,truth)

%variables
meas.K= truth.K;
meas.Z= cell(truth.K, model.sensors);

%generate measurements
for k=1:truth.K
    for s=1:model.sensors
        if truth.N(k) > 0
            idx= find( rand(truth.N(k),1) <= calculate_pd(model, truth.X{k}, model.S(:,s)));                                            %detected target indices
            %idx= find( rand(truth.N(k),1) <= model.P_D); 
            meas.Z{k,s}= gen_observation_fn(model,truth.X{k}(:,idx),'noise');                          %single target observations if detected 
        end
        N_c= poissrnd(model.lambda_c);                                                               %number of clutter points
        C= repmat(model.range_c(:,1),[1 N_c])+ diag(model.range_c*[ -1; 1 ])*rand(model.z_dim,N_c);  %clutter generation
        meas.Z{k,s}= [ meas.Z{k,s} C ];                                                                  %measurement is union of detections and clutter
    end
end
end

function pd= calculate_pd(model, X, S)
    %pds =  exp(-0.5 * dot((X([1 3], :)- repmat(S, 1, size(X,2))), model.inv_RD*(X([1 3], :)-repmat(S, 1, size(X,2)))));
    pd = -2*10^(-4) * sqrt(power(X(1,:) - repmat(S(1),1, size(X,2)), 2) + power(X(2,:) - repmat(S(2),1, size(X,2)), 2)) + 1;
    pd = pd';
end




























    