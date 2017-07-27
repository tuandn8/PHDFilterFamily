function [tt_update_new]= gaus_prune_l(tt_update,elim_threshold)

Ij = [];
for tabidx = 1 : length(tt_update)
    if tt_update{tabidx}.w > elim_threshold
        Ij = [Ij tabidx];
    end
end

tt_update_new = tt_update(Ij');
return;

