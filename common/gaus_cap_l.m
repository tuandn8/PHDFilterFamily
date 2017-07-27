function [tt_update_new]= gaus_cap_l(tt_update,model, max_number)

w = zeros(length(tt_update), 1);
m = zeros(model.x_dim, length(tt_update));
P = zeros(model.x_dim, model.x_dim, length(tt_update));
for tabidx = 1 : length(tt_update)
    w(tabidx, 1)= tt_update{tabidx}.w;
    m(:, tabidx) = tt_update{tabidx}.m;
    P(:, :, tabidx) = tt_update{tabidx}.P;
end

if length(w) > max_number
    [~,idx]= sort(w,1,'descend');
    tt_update = tt_update(idx(1:max_number));
    w_new= w(idx(1:max_number)); w_new = w_new * (sum(w)/sum(w_new));
    x_new= m(:,idx(1:max_number));
    P_new= P(:,:,idx(1:max_number));
    
    for tabidx = 1 : length(tt_update)
        tt_update{tabidx}.w = w_new(tabidx);
        tt_update{tabidx}.m = x_new(:, tabidx);
        tt_update{tabidx}.P = P_new(:,:,tabidx);
    end
end
 tt_update_new = tt_update;


