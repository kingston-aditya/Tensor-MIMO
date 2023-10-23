function ans = ALS(P,T,K,Y1,Y2,D,xmod)

    % initialization
    Ht = randn(P,K);
    St = randn(T,K);
    St_prev = zeros(T,K);

    J = 200;
    while norm(St - St_prev) > 10^-5 || J>1
        % update H
        Ht = Y1*pinv(khatrirao(D,St).');

        % update S
        St_prev = St;
        St = Y2*pinv(khatrirao(D,Ht).');

        % update J
        J = J-1;
    end

    % calculate lambda
    lambda_s = diag(real(St(1,:))./real(xmod));
    % H_hat = Ht/lambda_s;
    S_hat = St*pinv(lambda_s);

    ans = S_hat;
end