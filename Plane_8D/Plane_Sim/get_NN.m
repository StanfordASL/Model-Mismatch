function [nn_i, Cnn_i, Tnn_i, Pnn_i] = get_NN(i,nodes,N,r)

nn_i = zeros(1,N);
Cnn_i = zeros(1,N);
Tnn_i = zeros(1,N);
Pnn_i = cell(1,N);

%euclidean prune first
dists = norms(nodes(1:N,1:3)-repmat(nodes(i,1:3),N,1),2,2);
prune_mat = (dists > r);
prune_mat(i) = 1;

%now go through list
for j = 1:N
    if prune_mat(j)
        continue;
    else
        [~,T_ij,C_ij] = steer(nodes(i,:)',nodes(j,:)',0);
        if (C_ij <= r)
            nn_i(j) = 1; Cnn_i(j) = C_ij; Tnn_i(j) = T_ij;
            [Pnn_i{1,j},~,~] = steer(nodes(i,:)',nodes(j,:)',1);
        end
    end
end


end