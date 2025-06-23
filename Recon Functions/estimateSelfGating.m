function [SG_pair, sgI, SGfirstOccurrence, SGlastOccurrence] = estimateSelfGating(PEInd_SG)
  % Returns the most frequent [phase,partition] pair (SG channel),
  % its repeat interval (sgI), and first/last occurrence indices.

  % Convert pairs to string for easy comparison
    pairStrings = string(PEInd_SG(:,1)) + "," + string(PEInd_SG(:,2));
    
    % Find unique pairs and their counts
    [uniquePairs, ~, indices] = unique(pairStrings, 'stable');
    counts = accumarray(indices, 1);
    
    % Find the pair with the highest count
    [~, maxIdx] = max(counts);
    SG_pair = uniquePairs(maxIdx);
    
    % Find all occurrences of this pair
    occurrences = find(pairStrings == SG_pair);
    % Compute sgI and first occurrence location
    if length(occurrences) > 1
        sgI = diff(occurrences(1:2)); % Distance between first two occurrences
    else
        error("The most frequent pair appears only once. No repeating pattern found.");
    end
    SGfirstOccurrence = occurrences(1);
    SGlastOccurrence = occurrences(end);
    
    % Display results
    disp("SG Location in k-space (ky, kz): " + SG_pair);
    % disp("Frequency: " + num2str(sgI));
    % disp("First Occurrence: " + num2str(SGfirstOccurrence));
    % disp("Last Occurrence: " + num2str(SGlastOccurrence));
end
