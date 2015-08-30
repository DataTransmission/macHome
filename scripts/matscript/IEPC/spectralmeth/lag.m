function lagpoly = lag(points, m, XI)
    lagpoly = ones(length(XI), 1);
    for x=1:length(XI)
        for j=1:size(points)
            if j ~= m
               lagpoly(x) = lagpoly(x) * (XI(x) - points(j)) / (points(m) - points(j));     
            end
        end
    end
end