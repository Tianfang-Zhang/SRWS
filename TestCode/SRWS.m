%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is an implement of the Self-Regualrized Weighted Sparse (SRWS)
% model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you use this code in your publications, please cite:
% Zhang T, Peng Z, Wu H, et al. Infrared small target detection via self-regularized weighted sparse model[J]. Neurocomputing, 420: 124-148.
% @article{zhang420infrared,
%   title={Infrared small target detection via self-regularized weighted sparse model},
%   author={Zhang, Tianfang and Peng, Zhenming and Wu, Hao and He, Yanmin and Li, Chaohai and Yang, Chunping},
%   journal={Neurocomputing},
%   volume={420},
%   pages={124--148},
%   publisher={Elsevier}
% }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you have any questions, please contact:
% Author: Tianfang Zhang
% Email: sparkcarleton@gmail.com
% Copyright:  University of Electronic Science and Technology of China
% Date: 2019/7/3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%* License: Our code is only available for non-commercial research use.

classdef SRWS

    properties
        len = 50;
        step = 10;
        
        p = 1;
        alpha = 2;
        ogsK = 4;
        
        betaL = 1;
        lambdaL = 1;
        gammaL = 0.09;   
        
        loss;
        Z;
        rstB;
        rstE;
        result;
        
    end

    methods
        
        function obj = process(obj, inImg)
            [m, n] = size(inImg);
            
            % % Construct image-patch
            weight = obj.getOGS(inImg, obj.ogsK, obj.p, obj.alpha);
            patchWeight = obj.image2patch(weight, obj.len, obj.step);
            patchImg = obj.image2patch(inImg, obj.len, obj.step);
            
            % % Iterate solution and the implement of FRR model
            beta = obj.betaL / sqrt(max(m, n));
            lambda = obj.lambdaL / sqrt(max(m, n));
            gamma = obj.gammaL / sqrt(max(m, n));
            [obj.Z, E, T, obj.loss] = obj.optimization(patchImg, patchWeight, beta, lambda, gamma);
            B = patchImg * obj.Z;
            
            % % Reconstruct target image and background image
            rstT = obj.patch2image(T, obj.len, obj.step, size(inImg));
            obj.rstB = obj.patch2image(B, obj.len, obj.step, size(inImg));
            obj.rstE = obj.patch2image(E, obj.len, obj.step, size(inImg));
            
            obj.result = rstT .* (rstT>0);
            obj.rstB = obj.rstB .* (obj.rstB>0);
            obj.rstE = obj.rstE .* (obj.rstE>0);
        end
        
        function ogsMat = getOGS(obj, inImg, K, p, alpha)

            h = ones(K, K);

            [fx, fy] = gradient(inImg);

            ogs_fx = imfilter(fx, h);
            ogs_fy = imfilter(fy, h);

            ogsMat = (ogs_fx.^p + ogs_fy.^p).^(1/p);
            ogsMat = ogsMat ./ max(ogsMat(:));

            minValue = min(ogsMat(:));
            maxValue = max(ogsMat(:));

            ogsMat = exp(alpha*((ogsMat-minValue)/(maxValue-minValue)));
        end
        
        function patchImg = image2patch(obj, inImg, len, step)

            [m, n] = size(inImg);
            patchImg = zeros(len*len, (length(1:step:m-len)*length(1:step:n-len)));
            counter = 1;

            for i = 1:step:m-len
                for j = 1:step:n-len
                    tmp_patch = inImg(i:i+len-1, j:j+len-1);
                    patchImg(:, counter) = reshape(tmp_patch, len*len, 1);
                    counter = counter + 1;
                end
            end

        end
        
        function [Z, E, T, loss] = optimization(obj, inImg, weight, beta, lambda, gamma)

            % % Parameters setting
            [m, n] = size(inImg);

            W_k = zeros(n, n);  W_kp1 = zeros(n, n);
            V_k = zeros(n, n);  V_kp1 = zeros(n, n);
            J_k = zeros(n, n);  J_kp1 = zeros(n, n);
            Z_k = zeros(n, n);  Z_kp1 = zeros(n, n);
            E_k = zeros(m, n);  E_kp1 = zeros(m, n);
            T_k = zeros(m, n);  T_kp1 = zeros(m, n);

            Y1_k = zeros(m, n); Y1_kp1 = zeros(m, n);
            Y2_k = zeros(1, n); Y2_kp1 = zeros(1, n);
            Y3_k = zeros(n, n); Y3_kp1 = zeros(n, n);
            mu_k = 1 / (5*std(inImg(:)));
            
            M = zeros(n, n);

            iterNum = 0;
            Converged = false;
            normIn = norm(inImg, 'fro'); 
            loss = [];

            % % Iteraton steps
            while ~Converged
                % Update V
                G = (Z_k + Z_k') / 2;
                L = diag(sum(G, 2)) - G;
                try
                    V_kp1 = lyap(2*(W_k')*W_k, beta*(L+L'), -2*W_k'*Z_k);
                catch
                    V_kp1 = V_k;
                end
                % Update W
                W_kp1 = Z_k * pinv(V_kp1);

                % Update J
                    for i = 1:n
                        for j = 1:n
                            M(i, j) = norm(V_kp1(:,i)-V_kp1(:,j), 2) ^ 2;
                        end
                    end
                J_kp1 = obj.softThreshold(Z_k+Y3_k/mu_k, 0.5*beta*M/mu_k);

                % Update Z
                one_v = ones(1, n);
                A = 2*eye(n) + mu_k*(inImg'*inImg+one_v'*one_v+eye(n));
                B = 2*W_kp1 + inImg'*(mu_k*inImg-mu_k*E_k-mu_k*T_k+Y1_k) + one_v'*(mu_k*one_v-Y2_k)...
                    + mu_k*J_kp1 - Y3_k;
                Z_kp1 = A^(-1) * B;

                % Update E
                E_kp1 = obj.solveL21(inImg-inImg*Z_kp1-T_k+Y1_k/mu_k, lambda/mu_k);

                % Update T
                T_kp1 = obj.softThreshold(inImg-inImg*Z_kp1-E_kp1+Y1_k/mu_k, gamma*weight/mu_k);

                % Update Y1, Y2, Y3 and mu
                Y1_kp1 = Y1_k + mu_k*(inImg-inImg*Z_kp1-E_kp1-T_kp1);
                Y2_kp1 = Y2_k + mu_k*(one_v*Z_kp1-one_v);
                Y3_kp1 = Y3_k + mu_k*(Z_kp1-J_kp1);
                mu_kp1 = 1.5 * mu_k;

                % Judge converged
                iterNum = iterNum + 1;
                stopCriterion = norm(inImg-inImg*Z_kp1-E_kp1-T_kp1, 'fro') / normIn;
                loss(iterNum) = stopCriterion;

                nonzeroNum_k = sum(T_k(:) ~= 0);
                nonzeroNum_kp1 = sum(T_kp1(:) ~= 0);

                if stopCriterion < 1e-7 || abs(nonzeroNum_k-nonzeroNum_kp1) <= 1
                    Converged = true;
                end

                % Display the status
%                 disp(['Iteration Number: ', num2str(iterNum), '    loss: ', num2str(stopCriterion)]);
%                 disp(['Non zero number: ', num2str(nonzeroNum_k)]);
%                 disp(' ');  

                % Assignment to continue
                V_k = V_kp1;    W_k = W_kp1;    J_k = J_kp1;
                Z_k = Z_kp1;    E_k = E_kp1;    T_k = T_kp1;    
                Y1_k = Y1_kp1;  Y2_k = Y2_kp1;  Y3_k = Y3_kp1;  
                mu_k = mu_kp1;

            end
            % % Outputs
            Z = Z_k;
            E = E_k;
            T = T_k;


        end
        
        function image = patch2image(obj, patchImg, len, step, imageSize)

            m = imageSize(1); n = imageSize(2);

            recons = zeros(m, n, 100);
            countMatrix = zeros(m, n);
            index = 1;

            for i = 1:step:m-len
                for j = 1:step:n-len
                    % Count the time each pixel used and record the value
                    repatch = reshape(patchImg(:,index), len, len);
                    countMatrix(i:i+len-1, j:j+len-1) = countMatrix(i:i+len-1, j:j+len-1) + 1;

                    % Record the value of each pixel
                    for ii = i:i+len-1
                        for jj = j:j+len-1
                            recons(ii, jj, countMatrix(ii,jj)) = repatch(ii-i+1, jj-j+1);
                        end
                    end

                    index = index + 1;
                end
            end

            image = zeros(m, n);

            for i = 1:m
                for j = 1:n
                    % Median
                    if countMatrix(i ,j) > 0
                        vector = recons(i, j, 1:countMatrix(i,j));

                        image(i, j) = mean(vector);
                    end
                end
            end

        end
        
        function output = softThreshold(obj, inImg, epsilon)
            output = zeros(size(inImg));
            posImg = inImg .* (inImg>epsilon);
            negImg = inImg .* (inImg<-epsilon);
            output = output + (posImg-epsilon).*(inImg>epsilon) + (negImg+epsilon) .* (inImg<-epsilon);
        end
        
        function output = solveL21(obj, inImg, alpha)

            [m, n] = size(inImg);
            output = zeros(m, n);

            for i = 1:n
                tmp = norm(inImg(:, i), 'fro');

                if tmp > alpha
                    coef = (tmp - alpha) / tmp;
                    output(:, i) = coef * inImg(:, i);
                end
            end


        end
        
    end
    
    
end