function [ifpm_o, F] = MWT_PSO(A_RGB, B_RGB)
% Perform Single-Level 2D Multiwavelet decomposition on RGB channels of MRI image %  
X_Red = GHM(A_RGB(:,:,1));
X_Green = GHM(A_RGB(:,:,2));
X_Blue = GHM(A_RGB(:,:,3));

% Perform Single-Level 2D Multiwavelet decomposition on RGB channels of SPECT image % 
Y_Red = GHM(B_RGB(:,:,1));
Y_Green = GHM(B_RGB(:,:,2));
Y_Blue = GHM(B_RGB(:,:,3));

% Concatenate RGB channels
X = cat(3,X_Red,X_Green,X_Blue);
Y = cat(3,Y_Red,Y_Green,Y_Blue);

% % Plot LL, LH, HL, and HH frequency sub-bands %
% figure(1);
% imshow(X);
% figure(2);
% imshow(Y);

% Extract LL, LH, HL, and HH frequency sub-bands %
[N_A,M_A]=size(A_RGB(:,:,1));
[N_B,M_B]=size(B_RGB(:,:,1));

LL_X = X(1:N_A, 1:N_A, :);
LH_X = X(1:N_A, N_A+1:2*N_A, :);
HL_X = X(N_A+1:2*N_A, 1:N_A, :);
HH_X = X(N_A+1:2*N_A, N_A+1:2*N_A, :);

LL_Y = Y(1:N_B, 1:N_B, :);
LH_Y = Y(1:N_B, N_B+1:2*N_B, :);
HL_Y = Y(N_B+1:2*N_B, 1:N_B, :);
HH_Y = Y(N_B+1:2*N_B, N_B+1:2*N_B, :);

n_particles = 25; % number of particles to use in PSO
w = zeros(n_particles,8); % weights/particle positions
p = zeros(n_particles,9); % particles' best known positions + corresponding IFPMs
g = zeros(8); % swarm's best known position
ifpm_g = 0; % swarm's best IFPM value
v = zeros(n_particles,8); % particles' velocities
for i=1:n_particles
    % Generate a candidate solution (particle position) drawn from a uniform distribution %
    w(i,:) = rand(8,1); 
    
    % Initialize particle's best known position as this position %
    p(i,1:8) = w(i,:);
    
    % Fusion of sub-bands based on weights (from particle i) %
    LL_F = p(i,1)*LL_X + p(i,2)*LL_Y;
    LH_F = p(i,3)*LH_X + p(i,4)*LH_Y;
    HL_F = p(i,5)*HL_X + p(i,6)*HL_Y;
    HH_F = p(i,7)*HH_X + p(i,8)*HH_Y;
    
    % Concatenate resulting sub-bands %
    F = [LL_F,LH_F;HL_F,HH_F];    

    % Perform Inverse Multiwavelet transform on RGB channels %
    F_Red = IGHM(F(:,:,1));
    F_Green = IGHM(F(:,:,2));
    F_Blue = IGHM(F(:,:,3));
    
    % Concatenate RGB channels %
    F = cat(3,F_Red,F_Green,F_Blue);
    
    % Compute IFPM of candidate solution %
    ifpm_p = IFPM(im2uint8(A_RGB), im2uint8(B_RGB), im2uint8(F));
    p(i,9) = ifpm_p;
    
    % Update the swarm's best known position %
    if(ifpm_p > ifpm_g)
        ifpm_g = ifpm_p;
        g = p(i,1:8);
    end
    
    % Generate a particle velocity drawn from a uniform distribution %
    v(i,:) = -1 + (1+1)*rand(1,8);
end

%omegas = [-0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.1 0.2 0.3 0.4 0.5 0.6];
omegas = [-0.6 -0.4 -0.2 0.1 0.3 0.5];
phi_ps = [-1 -0.6 -0.15 0.5 2.1 2.5];
phi_gs = [0.6 1.33 2.2 2.6 3.4 4 4.9];
best_config = zeros(9);
max_iter = 10; % termination criterion

ifpm_g_last = 0;
count = 0;
term = 0;

for omega_i=1:size(omegas,2)
    for phi_p_i=1:size(phi_ps,2)
        for phi_g_i=1:size(phi_gs,2)
            
            for iter=1:max_iter
                for i=1:n_particles
                    for d=1:8
                        rp = rand(1);
                        rg = rand(1);
            
                        % Update particle's velocity %
                        v(i,d) = omegas(omega_i)*v(i,d) + phi_ps(phi_p_i)*rp*(p(i,d) - w(i,d)) + phi_gs(phi_g_i)*rg*(g(d) - w(i,d));
                    end
        
                    % Update particle's position %
                    w(i,:) = w(i,:) + v(i,:);
        
                    % Fusion of sub-bands based on new weights %
                    LL_F = w(i,1)*LL_X + w(i,2)*LL_Y;
                    LH_F = w(i,3)*LH_X + w(i,4)*LH_Y;
                    HL_F = w(i,5)*HL_X + w(i,6)*HL_Y;
                    HH_F = w(i,7)*HH_X + w(i,8)*HH_Y;
    
                    % Concatenate resulting sub-bands %
                    F = [LL_F,LH_F;HL_F,HH_F];    

                    % Perform Inverse Multiwavelet transform on RGB channels %
                    F_Red = IGHM(F(:,:,1));
                    F_Green = IGHM(F(:,:,2));
                    F_Blue = IGHM(F(:,:,3));
    
                    % Concatenate RGB channels %
                    F = cat(3,F_Red,F_Green,F_Blue);
    
                    % Compute IFPM of new solution %
                    ifpm_w = IFPM(im2uint8(A_RGB), im2uint8(B_RGB), im2uint8(F));
        
                    % Update the particle's best known position %
                    if(ifpm_w > p(i,9))
                        p(i,1:8) = w(i,:);
                        p(i,9) = ifpm_w;
                        % Update the swarm's best known position %
                        if(p(i,9) > ifpm_g)
                            ifpm_g = p(i,9);
                            g = p(i,1:8);
                        end
                    end
                end
                verbose = sprintf('Iteration: %d',iter);
                disp(verbose);
            end
            
            if(ifpm_g > best_config(9))
                best_config(9) = ifpm_g;
                best_config(1:8) = g;
            end
            
            if(ifpm_g - ifpm_g_last < 0.000001)
                count = count + 1;
            end
            ifpm_g_last = ifpm_g;
            
            verbose = sprintf('Configuration: %d %d %d ====> IFPM = %f',omega_i, phi_p_i, phi_g_i, best_config(9));
            disp(verbose);
            
            if(count == 5)
                term = 1;
                break
            end
        end
        if(term == 1)
            break
        end
    end
    if(term == 1)
        break
    end
end

% Fusion of sub-bands based on optimal weights %
LL_F = best_config(1)*LL_X + best_config(2)*LL_Y;
LH_F = best_config(3)*LH_X + best_config(4)*LH_Y;
HL_F = best_config(5)*HL_X + best_config(6)*HL_Y;
HH_F = best_config(7)*HH_X + best_config(8)*HH_Y;

% Concatenate resulting sub-bands %
F = [LL_F,LH_F;HL_F,HH_F];    

% Perform Inverse Multiwavelet transform on RGB channels %
F_Red = IGHM(F(:,:,1));
F_Green = IGHM(F(:,:,2));
F_Blue = IGHM(F(:,:,3));
    
% Concatenate RGB channels %
F = cat(3,F_Red,F_Green,F_Blue);
    
% Compute IFPM of optimal solution %
ifpm_o = best_config(9);

%%save('best_config.mat', 'best_config');
