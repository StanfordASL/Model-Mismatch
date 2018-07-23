function P_nn = get_steer_paths(nodes,N_nn)

N = size(nodes,1);

P_nn = cell(N,N);

for i = 1:N
    if mod(i,10)==0
        disp(i);
    end
    NN_i = find(N_nn(i,:));
    for j = 1:length(NN_i)
        [P_nn{i,NN_i(j)},~,~] = steer(nodes(i,:)',nodes(NN_i(j),:)',1);
    end
end

end