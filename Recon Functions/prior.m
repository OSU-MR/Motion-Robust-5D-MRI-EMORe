function pi = prior(w_sg, ratio_sg,ratio_out)
    bin_dim = ndims(w_sg);
    pi_rest = 1./(size(w_sg,bin_dim)-2+ratio_sg+ratio_out);
    pi_sg = ratio_sg.*pi_rest;
    pi_out = ratio_out.*pi_rest;
    pi = zeros(size(w_sg));
    pi(w_sg==0) = pi_rest;
    pi(w_sg==1) = pi_sg;
    pi(:,:,:,81) = pi_out;
end
