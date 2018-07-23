function box_V = update_teb(teb,r,psi,do_plot,varargin)

R =   [cos(psi), -sin(psi), 0;
       sin(psi), cos(psi), 0;
       0,0,1];
   
box_V = ( R*teb.V' + repmat(r,1,teb.nV) )';

if (do_plot)
    figure(varargin{1})
    trimesh(teb.hull,box_V(:,1),box_V(:,2),box_V(:,3),'FaceColor','k','EdgeColor','k',...
                                                      'FaceAlpha',0.2);
end

end