% Read MRI image and convert it to RGB format %
[A,map] = imread("MRI.gif",1); 
A_RGB = ind2rgb(A,map);

% Read SPECT image and convert it to RGB format %
[B,map] = imread("SPECT.gif",1);
B_RGB = ind2rgb(B,map);

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

n_particles = 500; % number of particles to use in PSO
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


max_iter = 10; % termination criterion
omega = 0.3; % PSO behavioral parameter
phi_p = 2; % PSO behavioral parameter
phi_g = 1; % PSO behavioral parameter
for iter=1:max_iter
    for i=1:n_particles
        for d=1:8
            rp = rand(1);
            rg = rand(1);
            
            % Update particle's velocity %
            v(i,d) = omega*v(i,d) + phi_p*rp*(p(i,d) - w(i,d)) + phi_g*rg*(g(d) - w(i,d));
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

% Fusion of sub-bands based on optimal weights %
LL_F = g(1)*LL_X + g(2)*LL_Y;
LH_F = g(3)*LH_X + g(4)*LH_Y;
HL_F = g(5)*HL_X + g(6)*HL_Y;
HH_F = g(7)*HH_X + g(8)*HH_Y;

% Concatenate resulting sub-bands %
F = [LL_F,LH_F;HL_F,HH_F];    

% Perform Inverse Multiwavelet transform on RGB channels %
F_Red = IGHM(F(:,:,1));
F_Green = IGHM(F(:,:,2));
F_Blue = IGHM(F(:,:,3));
    
% Concatenate RGB channels %
F = cat(3,F_Red,F_Green,F_Blue);
    
% Compute IFPM of optimal solution %
ifpm_o = IFPM(im2uint8(A_RGB), im2uint8(B_RGB), im2uint8(F));
