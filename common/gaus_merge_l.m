function [tt_update_new]= gaus_merge_l(tt_update, model, threshold)

L= length(tt_update); x_dim= model.x_dim;
I= 1:L;
el= 1;

if L == 0
    tt_update_new = tt_update;
    return;
end



while ~isempty(I),
    w = zeros(length(tt_update), 1);
    m = zeros(model.x_dim, length(tt_update));
    P = zeros(model.x_dim, model.x_dim, length(tt_update));
    for tabidx = 1 : length(tt_update)
        w(tabidx, 1)= tt_update{tabidx}.w;
        m(:, tabidx) = tt_update{tabidx}.m;
        P(:, :, tabidx) = tt_update{tabidx}.P;
    end

    [~,j]= max(w); j= j(1);
    Ij= []; iPt= inv(P(:,:,j));
    tt_update_new{el}.w = 0; 
    tt_update_new{el}.m = zeros(x_dim,1); 
    tt_update_new{el}.P = zeros(x_dim,x_dim);
    
    for i= I
        val= (m(:,i)-m(:,j))'*iPt*(m(:,i)-m(:,j));
        if val <= threshold,
            Ij= [ Ij i ];
        end;
    end;
    
    
   tt_update_new{el}.w = sum(w(Ij));
   tt_update_new{el}.m = wsumvec(w(Ij),m(:,Ij),x_dim);
   tt_update_new{el}.P = wsummat(w(Ij),P(:,:,Ij),x_dim);

   tt_update_new{el}.m= tt_update_new{el}.m/tt_update_new{el}.w;
   tt_update_new{el}.P= tt_update_new{el}.P/tt_update_new{el}.w;
   tt_update_new{el}.l= tt_update{j}.l;
   
   I= setdiff(I,Ij);
   tt_update = tt_update(I);
   L= length(tt_update); 
   I= 1:L;
   
   el= el+1;
end

function out = wsumvec(w,vecstack,xdim)
    wmat = repmat(w',[xdim,1]);
    out  = sum(wmat.*vecstack,2);

function out = wsummat(w,matstack,xdim)
    w = reshape(w,[1,1,size(w)]);
    wmat = repmat(w,[xdim,xdim,1]);
    out = sum(wmat.*matstack,3);
    
