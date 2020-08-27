function [M_branch_PD,x_A_PD,M_branch_torus,x_A_torus,M_branch_unstable,x_A_unstable,M_branch_stable,x_A_stable] = branch_stability(N,x_branch,M_branch,floquet_mtp_branch,s_branch)

% Initialize
M_branch_PD = []; x_A_PD = [];
M_branch_torus = []; x_A_torus = [];
M_branch_unstable = []; x_A_unstable = [];
M_branch_stable = []; x_A_stable = [];
L = length(M_branch);
m = zeros(1,N);
%%%%
abs_tol = 1;
real_tol = 1e-2;
imag_tol = 1e-3;
% Search
for j=1:L
    % Check for unstable parts of branch
    % If any multiplier is outside the unit circle (Give a little margin
    % for error on multiplier at +1) and it is far enough in the branch
    if any(abs(floquet_mtp_branch{j}(:))>abs_tol)  
        jU = abs(floquet_mtp_branch{j}(:))>abs_tol;
        if any(abs(imag(floquet_mtp_branch{j}(jU)))>imag_tol)         % If there is any multiplier with nonzero imaginary part            
            % Torus part of branch
            M_branch_PD = [M_branch_PD; nan]; x_A_PD = [x_A_PD; nan(1,N)];
            M_branch_torus = [M_branch_torus; M_branch{j}]; for k=1:N, m(1,k) = max(x_branch{j}(:,k)); end; x_A_torus = [x_A_torus; m];
            M_branch_unstable = [M_branch_unstable; nan]; x_A_unstable = [x_A_unstable; nan(1,N)];
            M_branch_stable = [M_branch_stable; nan]; x_A_stable = [x_A_stable; nan(1,N)];
        elseif any(real(floquet_mtp_branch{j}(:))<-1)             % If there is at least one multiplier with Re(lambda) < -1
            % Period-doubled part of branch
            M_branch_PD = [M_branch_PD; M_branch{j}]; for k=1:N, m(1,k) = max(x_branch{j}(:,k)); end; x_A_PD = [x_A_PD; m];
            M_branch_torus = [M_branch_torus; nan]; x_A_torus = [x_A_torus; nan(1,N)];
            M_branch_unstable = [M_branch_unstable; nan]; x_A_unstable = [x_A_unstable; nan(1,N)];
            M_branch_stable = [M_branch_stable; nan]; x_A_stable = [x_A_stable; nan(1,N)];
        elseif any(real(floquet_mtp_branch{j}(:))>1+real_tol)
            % Unstable through pitchfork, fold or transcritical bifurcation
            % of cycles
            M_branch_PD = [M_branch_PD; nan]; x_A_PD = [x_A_PD; nan(1,N)];
            M_branch_torus = [M_branch_torus; nan]; x_A_torus = [x_A_torus; nan(1,N)];
            M_branch_unstable = [M_branch_unstable; M_branch{j}]; for k=1:N, m(1,k) = max(x_branch{j}(:,k)); end; x_A_unstable = [x_A_unstable; m];
            M_branch_stable = [M_branch_stable; nan]; x_A_stable = [x_A_stable; nan(1,N)];
        else % Stable part of the branch
            M_branch_PD = [M_branch_PD; nan]; x_A_PD = [x_A_PD; nan(1,N)];
            M_branch_torus = [M_branch_torus; nan]; x_A_torus = [x_A_torus; nan(1,N)];
            M_branch_unstable = [M_branch_unstable; nan]; x_A_unstable = [x_A_unstable; nan(1,N)];
            M_branch_stable = [M_branch_stable; M_branch{j}]; for k=1:N, m(1,k) = max(x_branch{j}(:,k)); end; x_A_stable = [x_A_stable; m];
        end
    else % Stable part of the branch
        M_branch_PD = [M_branch_PD; nan]; x_A_PD = [x_A_PD; nan(1,N)];
        M_branch_torus = [M_branch_torus; nan]; x_A_torus = [x_A_torus; nan(1,N)];
        M_branch_unstable = [M_branch_unstable; nan]; x_A_unstable = [x_A_unstable; nan(1,N)];
        M_branch_stable = [M_branch_stable; M_branch{j}]; for k=1:N, m(1,k) = max(x_branch{j}(:,k)); end; x_A_stable = [x_A_stable; m];
    end
end

end