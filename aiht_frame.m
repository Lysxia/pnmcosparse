function [ X ] = aiht_frame( Xinit,gradient_step_handle,projection_handle,max_iteration )
%aiht_frame AIHT and AHTP frame function

    X = Xinit;
    
    for t = 1:max_iteration
        X = projection_handle(gradient_step_handle(X));
    end
    
end

