function D = calcDistMat(numCities, numSwarm, routes, coords, loops)
% Return distance of route from (x, y) coords of numCities
    D = zeros([numSwarm, 1]);
    count = 0;
    for row = routes.'
        count = count + 1;
        if loops(count) == 1
            for ii=1:numCities-1
                D(count) = D(count) + sqrt((coords(row(ii),2) - coords(row(ii+1),2))^2 + (coords(row(ii),3)-coords(row(ii+1),3))^2);
            end
            D(count) = D(count) + sqrt((coords(row(numCities),2)-coords(row(1),2))^2 + (coords(row(numCities),3)-coords(row(1),3))^2);
        end
    end
    
    
end

