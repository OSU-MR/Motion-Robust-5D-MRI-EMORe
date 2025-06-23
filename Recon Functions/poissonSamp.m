function points = poissonSamp(N, startVal, stopVal, r)
    % Default values if not provided
    if nargin < 2, startVal = 0; end
    if nargin < 3, stopVal = 1000; end
    if nargin < 4
        % A simple default: roughly evenly spaced minimum distance.
        r = max(1, floor((stopVal - startVal) / (N)));
    end

    points = [];       % Array to store accepted points
    maxAttempts = 10000; % Limit to avoid infinite loops
    attempts = 0;

    while numel(points) < N && attempts < maxAttempts
        candidate = randi([startVal, stopVal]);
        % Accept candidate if it's at least distance r from all already chosen points
        if isempty(points) || all(abs(points - candidate) >= r)
            points = [points, candidate];
        end
        attempts = attempts + 1;
    end

    if numel(points) < N
        warning('Could not sample enough points with the given minimum distance r.');
    end
    points = sort(points); % Optional: sort the points in increasing order
end