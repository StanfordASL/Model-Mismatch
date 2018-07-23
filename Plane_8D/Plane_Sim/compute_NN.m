function [N_nn,T_nn,C_nn] = compute_NN(nodes,r)

%ds: step-size
%r: connection radius
%nodes: N x 4

N = size(nodes,1);

N_nn = zeros(N,N);
T_nn = zeros(N,N);
C_nn = zeros(N,N);

%% First compute matrix of euclidean norm differences to prune

fprintf('Pruning...');
prune_mat = eye(N); %already ignore (i,i) pairs
for i = 1:N
   dists = norms(nodes(i+1:N,1:3)-repmat(nodes(i,1:3),N-i,1),2,2);
   prune_mat(i,i+1:N) = (dists > r)';
   prune_mat(i+1:N,i) = prune_mat(i,i+1:N)';
end
fprintf('Done\n');

%%

for i = 1:N
    for j = 1:N
        if prune_mat(i,j)
            continue;
        else
            if (mod(i,10)==0 && mod(j,100)==0)
                fprintf('(%d,%d)\n',i,j);
            end
            [~,T_ij,C_ij] = steer(nodes(i,:)',nodes(j,:)',0);
            if (C_ij <= r)
                N_nn(i,j) = 1;
                T_nn(i,j) = T_ij; C_nn(i,j) = C_ij;
            end
        end
    end
end

N_nn = sparse(N_nn);
T_nn = sparse(T_nn);
C_nn = sparse(C_nn);

end

