function dlagpoly = dlag(points, m, XI)
    dlagpoly = zeros(length(XI), 1);
    for x=1:length(XI)
        for i=1:size(points)
            if i ~= m
                prod = 1 / (points(m) - points(i));
                for j=1:size(points)
                    if j ~= m && j ~= i
                       prod = prod * (XI(x) - points(j)) / (points(m) - points(j));     
                    end
                end
                dlagpoly(x) = dlagpoly(x) + prod;
            end
        end
    end
end