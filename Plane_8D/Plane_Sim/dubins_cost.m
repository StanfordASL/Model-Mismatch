function [param, cost] = dubins_cost(p1,p2,r)

param = dubins_core(p1, p2, r);
if param.STATUS < 0
    cost = inf;
else
    cost = dubins_length(param);
end

end

function length = dubins_length(param)
length = param.SEG_param(1);
length = length + param.SEG_param(2);
length = length + param.SEG_param(3);
length = length * param.r;
end