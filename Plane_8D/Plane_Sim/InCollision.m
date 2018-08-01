function coll = InCollision(path,teb,obs)

% global v_hat dT;

% ds = v_hat*dT; %path step size
% dS = 10; %step size
% skip = round(dS/ds);

n_obs = size(obs,3);

obs_cen = reshape(mean(obs),3,n_obs);
N = size(path,1);

%% simple check

%for each obstacle, find closest point along path
dists_obs = zeros(n_obs,2);
for i = 1:n_obs
   [dists_obs(i,1),dists_obs(i,2)] = min( norms(path(:,1:3) - repmat(obs_cen(:,i)',N,1),2,2) ); 
end

%ignore some
ignore = dists_obs(:,1) >= 15; %could be more precise but it'll do

coll = 0;
%check non-ignore
if (sum(ignore)~=n_obs)
    obs_check = find(~ignore);
    found_coll = 0; j=1;
    while(~found_coll && j<= length(obs_check))
        %get obstacle index 
        o_idx = obs_check(j);
        %get path index of closest approach
        p_idx = dists_obs(o_idx,2);
        %get position at closest approach
        pos = path(p_idx,1:3)';
        %now rotate & translate teb to this position
        teb_V = update_teb(teb,pos,path(p_idx,4),0);
        
        %check
        found_coll = GJK(teb_V,obs(:,:,o_idx),4);
        
        %increment
        j = j+1;
    end
    if (found_coll)
        coll = 1;
        return;
    end
end

%% path wise check
% for i = 1:skip:size(path,1)
%    %let's first do quick pruning
%    pos = path(i,1:3)';
%    mean_dist = norms(obs_cen-repmat(pos,1,n_obs),2,1);
%    
%    ignore = mean_dist >= 15; %could be more precise but it'll do
%    
%    if (sum(ignore)~=0)
%         %now rotate & translate teb
%         teb_V = update_teb(teb,pos,path(i,4),0);
%         %find obstacles needing to be checked
%         obs_check = find(~ignore);
%         found_coll = 0; j = 1;
%         while (~found_coll && j <= length(obs_check))
%             o_idx = obs_check(j);
%             found_coll = GJK(teb_V,obs(:,:,o_idx),4);
%             j = j+1;
%         end
%         if (found_coll)
%             coll = 1;
%             return;
%         end           
%    end  
%    
% end




end