classdef UlResults
    properties
        nqIndex;
        trueSE;
        theoreticalSE;
        numericalSE;
    end
    methods
        function ulResults = UlResults(n_vec)
            if nargin ~= 0
                m = length(n_vec);
                ulResults = UlResults.empty(m,0);
                for i = 1:m
                    ulResults(i) = UlResults;
                    ulResults(i).nqIndex = n_vec(i);
                    ulResults(i).trueSE = zeros(7,1);
                    ulResults(i).theoreticalSE = zeros(7,1);
                    ulResults(i).numericalSE = zeros(7,1);
                end
            end
        end
    end
end