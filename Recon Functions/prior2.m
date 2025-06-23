function wo = prior2(w_init, pf1,pf2)
    w_init = squeeze(w_init);
    w_init = reshape(w_init,[size(w_init,1) 20 4]);
    wo_11 = 1./(1+(size(w_init,3)-1)./pf1+(size(w_init(:,:),2)-4)/pf2);
    wo_10 = wo_11./pf1;
    wo_00 = wo_11./pf2;
    wo = zeros(size(w_init));
    % Iterate through each slice along the third dimension
    for k = 1:size(w_init,1)
        % Iterate through each row of the current slice
        for i = 1:size(w_init,2)
            % Check if there is a single '1' in the current row
            if sum(w_init(k,i, :) > 0)
                % Set the entire row to 1
                wo(k,i, :) = wo_10;
            end
        end
    end
    wo(w_init==1) = wo_11;
    wo(wo==0) = wo_00;
    wo = reshape(wo,[1 size(w_init,1) 1 size(w_init(:,:),2)]);
    % wo = wo./sum(wo,bin_dim);
end
