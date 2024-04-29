function spins = SpinVerts(vertices, permno)
    % Compute designated # of permutations/spins of the input surface data
    % in FreeSurfer fsaverage5.

    % Initialize waitbar
    h = waitbar(0, 'Spinning vertices...');

    % Set up paths

    % Added 07/31/2020
    % Exclude the medial wall by labeling those vertices with NaN
    % Note: in the lh(rh).aparc.a2009s.annot for fsaverage5 data, 1644825 is
    % the label of vertices in the medial wall. Assign those vertices to NaN

    % Extract the corresponding sphere surface coordinates for rotation
    % gL = gifti(surface);
    % vertices = gL.vertices;

    % Initialize random number generator for reproducibility
    rng(0);

    % Initialize variables to save rotation
    spins = zeros(permno, size(vertices, 1));
    bl = vertices;
    t = zeros(permno,1);
    % Permutation starts
    for j = 1:permno
        tic
        % The updated uniform sampling procedure
        A = normrnd(0, 1, 3, 3);
        [TL, temp] = qr(A);
        TL = TL * diag(sign(diag(temp)));
        if (det(TL) < 0)
            TL(:, 1) = -TL(:, 1);
        end

        % Reflect across the Y-Z plane for right hemisphere
        I1 = eye(3, 3);
        I1(1, 1) = -1;
        TR = I1 * TL * I1;
        bl = bl * TL;

        % Find the pair of matched vertices with the min distance and reassign
        % values to the rotated surface.
        spins(j, :) = nearestneighbour(vertices', bl');

        t(j) = toc;

        tmean = mean(t(1:j));

        sec_remain = tmean*(permno-j);

        % Update waitbar and estimate time remaining
        waitbar(j / permno, h, sprintf('Spinning vertices... %d%% (%.2f seconds remaining)', round(j / permno * 100),sec_remain));

        % Uncomment the following lines if you want to save rotated data
        % bigrotl=[bigrotl; datal(spins)'];
        % bigrotr=[bigrotr; datar(Ir)'];
    end

    % Close waitbar when done
    close(h);

    % Save the workspace if needed
    % save(wsname,'bigrotl','bigrotr')
end