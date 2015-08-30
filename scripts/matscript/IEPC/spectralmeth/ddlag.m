function ddlagpoly = ddlag(points, m, XI)
    ddlagpoly = zeros(length(XI), 1);
    for x=1:length(XI)
        for i=1:size(points)
            for n=1:size(points)
                if i ~= m && n ~= m && i ~= n
                    prod = 1 / ((points(m) - points(i)) * (points(m) - points(n)));
                    for j=1:size(points)
                        if j ~= m && j ~= i && j ~= n
                           prod = prod * (XI(x) - points(j)) / (points(m) - points(j));     
                        end
                    end
                    ddlagpoly(x) = ddlagpoly(x) + prod;
                end
            end
        end
    end
end