function [S_hat , H_hat] = R_ALS(P,T,K,Y1,Y2,D,xmod)

        % initialization
        Ht = randn(P,K);
        St = randn(T,K);
        St_prev = zeros(T,K);
                
        lambd = 0;
        J = 200;
        while norm(St - St_prev) > 10^-5 || J>1
             % update H
             Ht = Y1*khatrirao(D,St)*pinv(khatrirao(D,St).'*khatrirao(D,St) + lambd*eye(K));
        
             % update S
             St_prev = St;
             St = Y2*khatrirao(D,Ht)*pinv(khatrirao(D,Ht).'*khatrirao(D,Ht) + lambd*eye(K));
             
             % update iterator
             J=J-1;
        end
            
            
        % calculate lambda
        lambda_s = diag(real(St(1,:))./real(xmod));
        H_hat = Ht*lambda_s;
        S_hat = St*pinv(lambda_s);
end