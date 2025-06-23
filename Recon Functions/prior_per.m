function pi = prior_per(w_sg, pi_sg,pi_out)
    % Your original density function
   % Scaling and shifting to achieve desired range [0.05, 0.95]
    pi_rest = (1-pi_sg-pi_out)./(size(w_sg,4)-2);
    pi = zeros(size(w_sg));
    pi(w_sg==0) = pi_rest;
    pi(w_sg==1) = pi_sg;
    pi(:,:,end) = pi_out;
end
