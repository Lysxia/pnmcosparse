function M = rand_sampling_matrix( m,d )

    M = randn(m,d);
    M_sum = sqrt(sum(M.^2,1));
    M = M./repmat(M_sum,m,1);

end

