classdef DlResults
    properties
        nqIndex;
        trueSE;
        theoreticalSE;
        numericalSE;
    end
    methods
        function dlResults = DlResults(n_vec)
            if nargin ~= 0
                m = length(n_vec);
                dlResults = DlResults.empty(m,0);
                for i = 1:m
                    dlResults(i) = DlResults;
                    dlResults(i).nqIndex = n_vec(i);
                    dlResults(i).trueSE = zeros(7,1);
                    dlResults(i).theoreticalSE = zeros(7,1);
                    dlResults(i).numericalSE = zeros(7,1);
                end
            end
        end
    end
end