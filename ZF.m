function pet = ZF(P,T,Y,D,X_mod,V)
    % form Ys
    Yp = Y(:,1:P,:) + V(:,1:P,:);
    Yp1 = tenmat(Yp,1).data;
    
    Y = Y(:,P+1:T,:) + V(:,P+1:T,:);
    Y2 = tenmat(Y,2).data;
    
    % estimate H
    Gamma_h = khatrirao(D,X_mod(1:P,:)).';
    H_est = Yp1*pinv(Gamma_h);
    
    % estimate S
    Gamma_s = khatrirao(D,H_est).';
    pet = Y2*pinv(Gamma_s);
end