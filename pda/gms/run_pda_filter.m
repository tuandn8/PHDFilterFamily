function est = run_pda_filter(model,truth, meas)
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

filter.P_G= 0.90;                               %gate size in percentage
filter.gamma= chi2inv(filter.P_G,model.z_dim);   %inv chi^2 dn gamma value
filter.gate_flag= 1;                             %gating on or off 1/0

filter.run_flag= 'disp';            %'disp' or 'silence' for on the fly output

est.filter= filter;

%=== Filtering 

%initial prior
m_update = zeros(model.x_dim, 0);
P_update = zeros(model.x_dim, model.x_dim, 0);

%recursive filtering
max_idx=0;
last_track_list= zeros(1, 0);
for k=1:meas.K
    % remove lost track in predicted state vector and covariance matrix
    if (size(last_track_list,2)) 
        [~, lostIdx] = setdiff(last_track_list, truth.track_list{k});
        [~, lostIdx] = setdiff(last_track_list, last_track_list(lostIdx));
        m_update = m_update(:, lostIdx);
        P_update = P_update(:,:,lostIdx);
    end
    last_track_list = truth.track_list{k};
    
    %---intensity prediction 
    if ~isempty(truth.X{k})
        if max(truth.track_list{k}) > max_idx % new born target 
            idx = find(truth.track_list{k} > max_idx);
            m_update= cat(2, m_update, truth.X{k}(:,idx));
            P_update= cat(3, P_update, reshape(repmat(model.P_init * model.P_init', 1, size(idx,2)), model.x_dim, model.x_dim, size(idx,2)));
        end
        if ~isempty(truth.track_list{k})
            max_idx= max(truth.track_list{k}); 
        end;
    end
    

    
    [m_predict,P_predict] = kalman_predict_multiple(model,m_update,P_update);                       %surviving components
       
    if ~isempty(m_predict)
        for i=1:size(m_predict,2)
            %---gatingZ
            Z = [];
            if filter.gate_flag
                Z= gate_meas_gms(meas.Z{k},filter.gamma,model,m_predict(:,i),P_predict(:,:,i));
            end
            
            m= size(Z,2);
            qz= zeros(1, m + 1);
            qz(1)=1-model.P_D*filter.P_G; 
            
            if m~=0
                [qz(2:m+1),~,~] = kalman_update_multiple(Z,model,m_predict(:,i),P_predict(:,:,i));
            end

            qz= qz/sum(qz);
            %-- update
            z0 = m_predict(:,i);
            z0 = z0([1 3]);
            Z = [z0 Z];
            [m_temp, P_temp] = pda_update_single(Z, model.H, model.R, qz, m_predict(:,i), P_predict(:,:,i));
            m_update(:,i) = m_temp;
            P_update(:,:,i) = P_temp;
        end
    end
    
     
    %--- state extraction
    est.N(k) = length(m_update);
    est.X{k} = m_update;
    
    %---display diagnostics
    if ~strcmp(filter.run_flag,'silence')
        disp([' time= ',num2str(k),...
         ' #est card=' num2str(length(m_update))]);
    end
end
end

function [m_update, P_update] = pda_update_single(Z, H, R, qz, m, P) 
    mu = H*m;
    S  = R+H*P*H'; Vs= chol(S); det_S= prod(diag(Vs))^2; inv_sqrt_S= inv(Vs); iS= inv_sqrt_S*inv_sqrt_S';
    K  = P*H'*iS;
    
    v_temp = Z - repmat(mu, [1 size(Z,2)]);
    v= bsxfun(@times, qz, v_temp);
    v= sum(v,2);
    
    VV = zeros(size(Z,1), size(Z,1));
    for i = 1:size(v_temp,2)
        VV_temp = v_temp(:,i) * v_temp(:,i)';
        VV = VV + VV_temp;
    end
    
    P_inov= K*(VV - v*v')*K';
          
    m_update = m + K*v;
    P_c = (eye(size(P))-K*H)*P;
    P_update = qz(1) * P + (1-qz(1))*P_c + P_inov;
end
            