function est = run_cphd_filter(model,meas)

% This is the MATLAB code for the GMCPHD filter 

%=== Setup

%output variables
est.X= cell(meas.K,1);
est.N= zeros(meas.K,1);
est.L= cell(meas.K,1);

%filter parameters
filter.L_max= 100;                  %limit on number of Gaussians
filter.elim_threshold= 1e-5;        %pruning threshold
filter.merge_threshold= 4;          %merging threshold

filter.N_max= 20;                   %maximum cardinality number (for cardinality distribution)

filter.P_G= 0.999;                               %gate size in percentage
filter.gamma= chi2inv(filter.P_G,model.z_dim);   %inv chi^2 dn gamma value
filter.gate_flag= 1;                             %gating on or off 1/0

filter.run_flag= 'disp';            %'disp' or 'silence' for on the fly output

est.filter= filter;

%=== Filtering 

%initial prior
% w_update(1)= eps;
% m_update(:,1)= [0.1;0;0.1;0];
% P_update(:,:,1)= diag([1 1 1 1]).^2;
% L_update = 1;
 cdn_update= [1; zeros(filter.N_max,1)];

tt_cphd_update=cell(0,1);                       % track table for cphd (cell array of structs for invidual tracks)
% tt_cphd_update{1}.w=eps;
% tt_cphd_update{1}.m=[0.1;0;0.1;0];
% tt_cphd_update{1}.P=diag([1 1 1 1]).^2;
% tt_cphd_update{1}.l=[1,1];


%figure; 
%recursive filtering
for k=1:meas.K
    % ------------------------------------------------------------------------
    % B1. Calculate prediction step which include intensity prediction v(k|k-1)
    %     and cardainality prediction rho(k|k-1)
    % ------------------------------------------------------------------------
    
    
    % B1.1 Calculate intensity prediction: include birth component from model and
    % surviving component from existed component in previous step
    
    % for birth component
    tt_cphd_birth=cell(model.L_birth,1);                                                    %number of birth component
    for tabidx=1:model.L_birth
        tt_cphd_birth{tabidx}.w= model.w_birth(tabidx);
        tt_cphd_birth{tabidx}.m= model.m_birth(:,tabidx);
        tt_cphd_birth{tabidx}.P= model.P_birth(:,:,tabidx);
        tt_cphd_birth{tabidx}.l= [k;tabidx];
    end
    
    % for surviving component
    tt_cphd_survive= cell(length(tt_cphd_update), 1);
    for tabidx=1:length(tt_cphd_update)
        [mtemp_predict, Ptemp_predict]= kalman_predict_multiple(model, tt_cphd_update{tabidx}.m, tt_cphd_update{tabidx}.P);
        tt_cphd_survive{tabidx}.w = model.P_S * tt_cphd_update{tabidx}.w;
        tt_cphd_survive{tabidx}.m = mtemp_predict;
        tt_cphd_survive{tabidx}.P = Ptemp_predict;
        tt_cphd_survive{tabidx}.l = tt_cphd_update{tabidx}.l;
    end
    
    % merge surviving and birth components
    tt_cphd_predict=cat(1,tt_cphd_birth, tt_cphd_survive);                  L_predict = length(tt_cphd_predict);
    
    
    % B1.2 Calculate cardinality prediction: equal to sum from zero to N_max of product
    % birth's cardinality component with binominal of j existed component
    % surviving cardinality distribution
    
    % terms = sum_l_j_inf(binominal_C_j_l * surv_prob^j * (1 - surv_prob)^(l-j) * cdn_update_k_minus_1)
    % cdn_pre = sum_j_zero_n(card_birth(n-j) * terms)
    survive_cdn_predict = zeros(filter.N_max+1,1);
    % First calculate sum_l_j_inf(binominal_C_j_l * surv_prob^j * (1 - surv_prob)^(l-j))
    % remember that log(binominal_c_j_l) = sum(1..l) - sum(1..j) - sum(l-j..l)   
    for j=0:filter.N_max
        idxj=j+1;
        terms= zeros(filter.N_max+1,1);
        for ell=j:filter.N_max
            idxl= ell+1;
            terms(idxl) = exp(sum(log(1:ell))-sum(log(1:j))-sum(log(1:ell-j))+j*log(model.P_S)+(ell-j)*log(model.Q_S))*cdn_update(idxl);
        end
        survive_cdn_predict(idxj) = sum(terms);
    end

    % predicted cardinality= convolution of birth and surviving cardinality distribution
    % remember that cdn_birth is iid cluster RFS
    % terms = 1/[(n-j)! * sum(w_birth)] * sum(w_brith)^(n-j) * survive_cdn_predict
    cdn_predict = zeros(filter.N_max+1,1);
    for n=0:filter.N_max
        idxn=n+1;
        terms= zeros(filter.N_max+1,1);
        for j=0:n
            idxj= j+1;
            terms(idxj)= exp(-sum(model.w_birth)+(n-j)*log(sum(model.w_birth))-sum(log(1:n-j)))*survive_cdn_predict(idxj);
        end
        cdn_predict(idxn) = sum(terms);
    end

    %normalize predicted cardinality distribution
    cdn_predict = cdn_predict/sum(cdn_predict);
    
    % -----------------------------------------------------------------------------
    % B2. Perform calculate gate to reduce the complex. It's normal step in
    %     tracking problem
    % -----------------------------------------------------------------------------
    if filter.gate_flag
        m_tracks = [];
        P_tracks = [];
        for tabidx = 1:length(tt_cphd_predict)
            m_tracks = cat(2, m_tracks, tt_cphd_predict{tabidx}.m);
            P_tracks = cat(3, P_tracks, tt_cphd_predict{tabidx}.P);
        end       
        meas.Z{k}= gate_meas_gms(meas.Z{k},filter.gamma,model,m_tracks,P_tracks);
        %meas.Z{k}= gate_meas_gms(meas.Z{k},filter.gamma,model,m_predict,P_predict);        
    end
        
    % ----------------------------------------------------------------------------
    % B3. Calculate intensity update v(k|k) and cardinality update rho(k|k)
    % ----------------------------------------------------------------------------
    %number of measurements
    m= size(meas.Z{k},2);
    % number of update components
    tt_update = cell((m+1)*length(tt_cphd_predict),1);                          L_posterior = length(tt_update);
    
    qz_update = zeros(length(tt_cphd_predict), m);
    w_predict = zeros(length(tt_cphd_predict), 1);
    for tabidx = 1:length(tt_cphd_predict)
        w_predict(tabidx, 1) = tt_cphd_predict{tabidx}.w;
        tt_update{tabidx} = tt_cphd_predict{tabidx};
    end
    
  
    % using Kalman filter to calculate state, cov, likelihood of each
    % Gaussian component
    for emm = 1:m
        for tabidx = 1:length(tt_cphd_predict)
            stoidx = emm * length(tt_cphd_predict) + tabidx;
            [qz_temp, m_temp, P_temp] = kalman_update_multiple(meas.Z{k}(:,emm), model, tt_cphd_predict{tabidx}.m, tt_cphd_predict{tabidx}.P);
            qz_update(tabidx, emm) = qz_temp;
            tt_update{stoidx}.m = m_temp; 
            tt_update{stoidx}.P = P_temp;
            tt_update{stoidx}.l = tt_cphd_predict{tabidx}.l;
        end
    end

    % pre calculation for elementary symmetric functions   
    XI_vals = zeros(m,1); 
    for ell=1:m
       XI_vals(ell) = model.P_D * w_predict' * qz_update(:, ell)/model.pdf_c;
    end
    
    esfvals_E = esf(XI_vals);                   %calculate esf for entire observation set
    esfvals_D = zeros(m,m);                     %calculate esf with each observation index removed one-by-one
    for ell=1:m
        esfvals_D(:,ell) = esf([XI_vals(1:ell-1); XI_vals(ell+1:m)]);
    end
    
    %pre calculation for upsilons
    upsilon0_E = zeros(filter.N_max+1,1);
    upsilon1_E = zeros(filter.N_max+1,1);
    upsilon1_D = zeros(filter.N_max+1,m);
    
    for n=0:filter.N_max
        idxn= n+1;
        
        % calculate upsilon0_E(idxn)
        % Remember that formular
        % ternsu_E(idxj) = (|Z|-j)! * p_K_k(|Z|-j) * permu_n_(j+u) * (1-pD)^(n-j-u) / sum(w) * sym_ej  
        % upsilon_k_u[w,Z](n) = sum_0_min(|Z|,n) (termu_E) 
        terms0_E= zeros(min(m,n)+1,1);  
        for j=0:min(m,n)
            idxj= j+1;
            terms0_E(idxj) = exp(-model.lambda_c+(m-j)*log(model.lambda_c)+sum(log(1:n))-sum(log(1:n-j))+(n-j)*log(model.Q_D)-j*log(sum(w_predict)))*esfvals_E(idxj);
        end
        upsilon0_E(idxn)= sum(terms0_E);
        
        terms1_E= zeros(min(m,n)+1,1);  %calculate upsilon1_E(idxn)
        for j=0:min(m,n)
            idxj= j+1;
            if n>=j+1
                terms1_E(idxj) = exp(-model.lambda_c+(m-j)*log(model.lambda_c)+sum(log(1:n))-sum(log(1:n-(j+1)))+(n-(j+1))*log(model.Q_D)-(j+1)*log(sum(w_predict)))*esfvals_E(idxj);
            end
        end
        upsilon1_E(idxn)= sum(terms1_E);
        
        if m~= 0                        %calculate upsilon1_D(idxn,:) if m>0
            terms1_D= zeros(min((m-1),n)+1,m);
            for ell=1:m
                for j=0:min((m-1),n)
                    idxj= j+1;
                    if n>=j+1
                        terms1_D(idxj,ell) = exp(-model.lambda_c+((m-1)-j)*log(model.lambda_c)+sum(log(1:n))-sum(log(1:n-(j+1)))+(n-(j+1))*log(model.Q_D)-(j+1)*log(sum(w_predict)))*esfvals_D(idxj,ell);
                    end
                end
            end
            upsilon1_D(idxn,:)= sum(terms1_D,1);
        end
    end
    
    %---cardinality update
    cdn_update= upsilon0_E.*cdn_predict;
    cdn_update= cdn_update/sum(cdn_update);
            
    % missed detection term
    for tabidx = 1:length(tt_cphd_predict)
        tt_update{tabidx}.w = model.Q_D * (upsilon1_E'*cdn_predict)/(upsilon0_E'*cdn_predict) * tt_cphd_predict{tabidx}.w;
        %tt_update{tabidx}
    end
  
    % detected term
    if m~=0
        %m detection terms 
        %Kalman update precalculated and stored
        for ell=1:m
            for tabidx = 1 : length(tt_cphd_predict)
                stoidx = ell * length(tt_cphd_predict) + tabidx;
                %disp('qz_update '); qz_update(:,ell)
                %disp('w_predict '); w_predict
%                 disp('1 '); (upsilon1_D(:,ell)'*cdn_predict)/(upsilon0_E'*cdn_predict)
%                 disp('2 '); (upsilon1_D(:,ell)'*cdn_predict)/(upsilon0_E'*cdn_predict)*model.P_D.*qz_update(tabidx,ell)/model.pdf_c
%                 disp('3 '); w_predict(tabidx,1)
                tt_update{stoidx}.w = (upsilon1_D(:,ell)'*cdn_predict)/(upsilon0_E'*cdn_predict)*model.P_D.*qz_update(tabidx,ell)/model.pdf_c.*w_predict(tabidx,1);
            end
        end
    end 

    %---mixture management
   
    %clf;
    %stem3(m_update(1,:), m_update(2,:), w_update(:), 'ro');
    %hold on;
    %pause(0.5);
    
    
    %pruning, merging, capping   
    [tt_update] = gaus_prune_l(tt_update, filter.elim_threshold);
    [tt_update] = gaus_merge_l(tt_update, model, filter.merge_threshold);
    [tt_update] = gaus_cap_l(tt_update,   model, filter.L_max);
    
    tt_cphd_update = tt_update;
    L_clean = length(tt_cphd_update);
    %--- state extraction
    [est.X{k}, est.N(k),est.L{k}] = extract_estimate(tt_cphd_update, cdn_update, model);
    
   
    %---display diagnostics
    if ~strcmp(filter.run_flag,'silence')
        disp([' time= ',num2str(k),...
         ' #eap cdn=' num2str([0:filter.N_max]*cdn_update,4),...
         ' #var cdn=' num2str([0:filter.N_max].^2*cdn_update-([0:filter.N_max]*cdn_update)^2,4),...
         ' #est card=' num2str(est.N(k),4),...
         ' #gaus pred=',num2str(L_predict),...
         ' #gaus post=',num2str(L_posterior), ...
         ' #gaus updt=',num2str(L_clean)   ]);
    end
end
end


function [X,N,L]=extract_estimate(tt_cphd_update, cdn_update, model)
    %extract estimates via MAP cardinality and corresponding tracks
    [~,idx_max_cdn] = max(cdn_update);
    map_cdn = idx_max_cdn-1;
    N = min(length(tt_cphd_update),map_cdn);
    
    X = zeros(model.x_dim, N);
    L = zeros(2, N);
    
    w_update = zeros(length(tt_cphd_update), 1);
    for tabidx = 1 : length(tt_cphd_update)
        w_update(tabidx, 1)= tt_cphd_update{tabidx}.w;
    end
    [~, idxcmp] = sort(w_update, 'descend');
    for n = 1: N
        [~, idxtrk] = max(tt_cphd_update{idxcmp(n)}.w);
        X(:,n)=tt_cphd_update{idxcmp(n)}.m;
        L(:,n)=tt_cphd_update{idxcmp(n)}.l;
    end
    
end

            