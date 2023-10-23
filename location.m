function H=location(L,N,K)

    % Area parameters
    ll = -400;
    ul = 400;
    
    % AP and UE location
    ap_lx = ll + ul*rand(L,1);
    ap_ly = ll + ul*rand(L,1);
    
    ue_lx = ll + ul*rand(K,1);
    ue_ly = ll + ul*rand(K,1);
    
    % calculate channel gain
    H = [];
    for i=1:L
        h_l = [];
        for j=1:K
            beta_db = -30.5 - 36.7*log10(((ap_lx(i) - ue_lx(j))^2 + (ap_ly(i) - ue_ly(j))^2 + 100)^(1/3)) + 4*randn(1,1);
            beta = 10.^(beta_db/20);
            h_l = [h_l, sqrt(beta)*randn(N,1)+1i*sqrt(beta)*randn(N,1)];
        end
        H = [H ; h_l];
    end
end

