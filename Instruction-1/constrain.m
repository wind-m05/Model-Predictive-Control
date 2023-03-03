function y = constrain(x,lb,ub)
if x <= lb
    y = lb;
elseif x >= ub
    y = ub;
else
    y=x;
end

