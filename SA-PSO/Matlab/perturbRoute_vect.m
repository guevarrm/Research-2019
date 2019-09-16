function [theCityRoute] = genRoute(numCities, theCityRoute)
    randIndex1 = randi(numCities);
    alreadyChosen = true;
    while alreadyChosen == true
        randIndex2 = randi(numCities);
        if randIndex2 ~= randIndex1
            alreadyChosen = false;
        end
    end
    theCityRoute(:,[randIndex1 randIndex2]) = theCityRoute(:,[randIndex2 randIndex1]);
end
