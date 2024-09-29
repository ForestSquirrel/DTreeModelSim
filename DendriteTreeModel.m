classdef DendriteTreeModel
    properties
        dendrites       % Array of dendrite objects
        numDendrites    % Number of dendrites
        Stimuli         % Cell array of stimuli
        Connectivity    % Struct array representing connections
        h               % Integration step size
        coupling_matrix % Matrix for sparse computations
        id_to_idx       % Mapping from ID to index
    end

    methods
        function obj = DendriteTreeModel()
            obj.dendrites = Dendrite.empty;
            obj.numDendrites = 0;
            obj.Stimuli = {};
            obj.Connectivity = struct('ID', {}, 'List_Proximal', {}, 'List_Distal', {});
            obj.h = [];
            obj.id_to_idx = [];
        end

        function obj = addDendrite(obj, dendrite)
            obj.numDendrites = obj.numDendrites + 1;
            obj.dendrites(obj.numDendrites) = dendrite;
            obj.Stimuli{obj.numDendrites} = [];
            obj.id_to_idx(dendrite.ID + 1) = obj.numDendrites; % Assuming IDs start from 0
            % Initialize connectivity for this dendrite
            obj.Connectivity(obj.numDendrites).ID = dendrite.ID;
            obj.Connectivity(obj.numDendrites).List_Proximal = [];
            obj.Connectivity(obj.numDendrites).List_Distal = [];
        end

        function idx = getDendriteIndex(obj, ID)
            idx = obj.id_to_idx(ID + 1); % Adjust for MATLAB indexing
        end

        function valid = isValidID(obj, ID)
            valid = any([obj.dendrites.ID] == ID);
        end

        function obj = addConnection(obj, ID, List_Proximal, List_Distal)
            idx = obj.getDendriteIndex(ID);
            if ~isempty(List_Proximal)
                obj.Connectivity(idx).List_Proximal = unique([obj.Connectivity(idx).List_Proximal, List_Proximal]);
                % Add reciprocal connections
                for pid = List_Proximal
                    pidx = obj.getDendriteIndex(pid);
                    obj.Connectivity(pidx).List_Distal = unique([obj.Connectivity(pidx).List_Distal, ID]);
                end
            end
            if ~isempty(List_Distal)
                obj.Connectivity(idx).List_Distal = unique([obj.Connectivity(idx).List_Distal, List_Distal]);
                % Add reciprocal connections
                for did = List_Distal
                    didx = obj.getDendriteIndex(did);
                    obj.Connectivity(didx).List_Proximal = unique([obj.Connectivity(didx).List_Proximal, ID]);
                end
            end
        end

        function obj = addConnectionStr(obj, connStr)
            % Parse the input connection string and add connections to the model
            % according to the specified rules, allowing groups anywhere in the string.

            % Remove any spaces from the input string
            connStr = strrep(connStr, ' ', '');

            % Check if the string starts with '-'
            if startsWith(connStr, '-')
                error('Invalid connection string: Soma cannot have "-" to its left.');
            end

            % Split the string into tokens, respecting groups
            tokens = splitConnectionString(connStr);

            % Parse tokens into lists of IDs
            tokenIDs = cell(1, length(tokens));
            for i = 1:length(tokens)
                token = tokens{i};
                tokenIDs{i} = parseToken(token);
            end

            % Validate IDs
            allIDs = unique([tokenIDs{:}]);
            for i = 1:length(allIDs)
                if ~obj.isValidID(allIDs(i))
                    error('Invalid ID %d in connection string.', allIDs(i));
                end
            end

            % Now, for each pair of adjacent tokens, create connections
            for i = 1:length(tokenIDs)-1
                idsA = tokenIDs{i};
                idsB = tokenIDs{i+1};
                % Create connections from every ID in idsA to every ID in idsB
                for a = idsA
                    for b = idsB
                        obj = obj.addConnection(a, [], b);
                    end
                end
            end

            % Nested function to split the connection string
            function tokens = splitConnectionString(str)
                str = char(str);
                tokens = {};
                currentToken = '';
                insideGroup = false;
                idx = 1;
                while idx <= length(str)
                    ch = str(idx);
                    if ch == '['
                        insideGroup = true;
                        currentToken = [currentToken, ch];
                    elseif ch == ']'
                        insideGroup = false;
                        currentToken = [currentToken, ch];
                    elseif ch == '-' && ~insideGroup
                        % End of token
                        tokens{end+1} = currentToken;
                        currentToken = '';
                    else
                        currentToken = [currentToken, ch];
                    end
                    idx = idx + 1;
                end
                if ~isempty(currentToken)
                    tokens{end+1} = currentToken;
                end
            end

            % Nested function to parse a token into IDs
            function ids = parseToken(token)
                if any(startsWith(token, '[')) && any(endsWith(token, ']'))
                    % It's a group
                    content = token(2:end-1);
                    if isempty(content)
                        error('Empty group in connection string.');
                    end
                    idStrs = strsplit(content, ',');
                    try
                        ids = cellfun(@str2double, idStrs);
                    catch
                        error('Invalid IDs in connection string.');
                    end
                else
                    % It's a single ID
                    try
                        ids = str2double(token);
                    catch
                        error('Invalid IDs in connection string.');
                    end
                end
            end
        end

        function obj = addStimuli(obj, ID, Signal)
            idx = obj.getDendriteIndex(ID);
            obj.Stimuli{idx} = Signal;
        end

        function initCon = getInitCon(obj)
            initCon = [];
            for i = 1:obj.numDendrites
                initCon = [initCon; obj.dendrites(i).InitialConditions];
            end
        end

        function obj = buildModel(obj, h)
            obj.h = h;
            obj = obj.preprocessStimuli();
            % The model is ready for simulation after this step
        end

        function dxdt = computeRHS(obj, t, x)
            dxdt = zeros(size(x));
            Len = length(x);
            for i = 1:2:Len
                idx = (i + 1) / 2; % Dendrite index
                dendrite = obj.dendrites(idx);
                u = x(i);
                v = x(i+1);
                % Compute coupling term
                coupling_term = 0;
                conn = obj.Connectivity(idx);
                % Proximal connections
                for pid = conn.List_Proximal
                    pidx = obj.getDendriteIndex(pid);
                    pu = x((pidx-1)*2 + 1);
                    coupling_term = coupling_term + dendrite.params.gc * (pu - u);
                end
                % Distal connections
                for did = conn.List_Distal
                    didx = obj.getDendriteIndex(did);
                    du = x((didx-1)*2 + 1);
                    coupling_term = coupling_term + dendrite.params.gc * (du - u);
                end
                % Get external stimuli if any
                Stim = 0;
                StimSignal = obj.Stimuli{idx};
                if ~isempty(StimSignal)
                    Stim = StimSignal(t);
                end
                % Compute RHS
                if dendrite.ID == 0 % Soma
                    [dUdt, dVdt] = soma_rhs(u, v, dendrite.params, coupling_term);
                else % Dendrite
                    [dUdt, dVdt] = dendrite_rhs(u, v, dendrite.params, coupling_term, Stim);
                end
                dxdt(i) = dUdt;
                dxdt(i+1) = dVdt;
            end
        end

        %%%%%% Sparse computations functions %%%%%%
        function dxdt = sparseRHS(obj, t, x)
            N = obj.numDendrites;
            U = x(1:2:end);
            V = x(2:2:end);

            % Pre-extract parameters
            params = [obj.dendrites.params];
            alpha = [params.alpha]';
            b = [params.b]';
            Tau = [params.Tau]';
            TauR = [params.TauR]';
            NaX = [params.NaX]';
            gc = [params.gc]';

            % Precompute coupling matrix if not already done
            if isempty(obj.coupling_matrix)
                obj.coupling_matrix = buildCouplingMatrix(obj);
            end
            C = obj.coupling_matrix;

            % Compute coupling term
            coupling_term = C * U - sum(C, 2) .* U;

            Stim = zeros(N, 1);
            for idx = 1:N
                StimSignal = obj.Stimuli{idx};
                if ~isempty(StimSignal)
                    Stim(idx) = StimSignal(t);
                end
            end

            % Compute derivatives
            dUdt = Tau .* (NaX .* (U .* (U - 1) .* (1 - alpha .* U) - V) + coupling_term + Stim);
            dVdt = TauR .* b .* U;

            % Combine derivatives
            dxdt = zeros(size(x));
            dxdt(1:2:end) = dUdt;
            dxdt(2:2:end) = dVdt;
        end

        function C = buildCouplingMatrix(obj)
            N = obj.numDendrites;
            % Estimate the total number of non-zero entries
            nonZeroEntries = 0;
            for i = 1:N
                conn = obj.Connectivity(i);
                nonZeroEntries = nonZeroEntries + length(conn.List_Proximal) + length(conn.List_Distal);
            end

            % Preallocate arrays for the non-zero entries
            rows = zeros(nonZeroEntries, 1);
            cols = zeros(nonZeroEntries, 1);
            vals = zeros(nonZeroEntries, 1);

            % Index for inserting into preallocated arrays
            idx = 1;
            for i = 1:N
                conn = obj.Connectivity(i);
                gc_i = obj.dendrites(i).params.gc;
                % Proximal connections
                for pid = conn.List_Proximal
                    pidx = obj.getDendriteIndex(pid);
                    rows(idx) = i;
                    cols(idx) = pidx;
                    vals(idx) = gc_i;
                    idx = idx + 1;
                end
                % Distal connections
                for did = conn.List_Distal
                    didx = obj.getDendriteIndex(did);
                    rows(idx) = i;
                    cols(idx) = didx;
                    vals(idx) = gc_i;
                    idx = idx + 1;
                end
            end

            % Create the sparse coupling matrix
            C = sparse(rows, cols, vals, N, N);
        end
        function obj = preprocessStimuli(obj)
            % Preprocess stimuli to ensure they can be evaluated efficiently during RHS computation
            for idx = 1:obj.numDendrites
                StimSignal = obj.Stimuli{idx};
                if ~isempty(StimSignal)
                    if isnumeric(StimSignal)
                        % Create a griddedInterpolant for efficient interpolation
                        time_vector = 0:1:(length(StimSignal)-1);
                        F = griddedInterpolant(time_vector.*obj.h, StimSignal, 'pchip', 'pchip');
                        obj.Stimuli{idx} = F;
                    end
                end
            end
        end

        %%%%%%% MISC %%%%%%%%%%%%%%%%
        function verify(obj)
            % Verify that all dendrites are connected to the soma (ID = 0)
            % If there are floating dendrites, output an error with their IDs

            % Get the list of all IDs
            allIDs = [obj.dendrites.ID];

            % Initialize visited set
            visitedIDs = [];

            % Initialize a queue for BFS
            queue = [];

            % Start from soma (ID = 0)
            if ~obj.isValidID(0)
                error('Soma (ID = 0) not found in the model.');
            end

            queue(end+1) = 0;
            visitedIDs(end+1) = 0;

            while ~isempty(queue)
                currentID = queue(1);
                queue(1) = []; % Remove the first element

                idx = obj.getDendriteIndex(currentID);
                conn = obj.Connectivity(idx);

                % Get connected IDs (both proximal and distal)
                connectedIDs = [conn.List_Proximal, conn.List_Distal];

                for i = 1:length(connectedIDs)
                    connectedID = connectedIDs(i);
                    if ~ismember(connectedID, visitedIDs)
                        visitedIDs(end+1) = connectedID;
                        queue(end+1) = connectedID;
                    end
                end
            end

            % Now compare visitedIDs with allIDs
            unvisitedIDs = setdiff(allIDs, visitedIDs);
            if ~isempty(unvisitedIDs)
                error('Floating dendrites detected with IDs: %s', mat2str(unvisitedIDs));
            end
            % If everything is connected, do nothing
        end

        %%%%%% Force-Directed 3D %%%%%%
        function visualizeFDGL(obj, maxIterations, eps)
            % visualizeFDGL (force-directed graph layout) performs a force-directed layout with
            % spherical initialization and convergence criteria.
            %
            % Inputs:
            % - maxIterations: Maximum number of iterations for the algorithm
            % - eps: Convergence threshold for the total energy

            % Constants for force calculation
            K_attraction = .2;    % Spring constant for attraction (tune as necessary)
            K_repulsion = .2;   % Constant for repulsion (tune as necessary)
            timestep = .2;     % Time step for updating positions (tune as necessary)

            % Assign initial positions using spherical initialization
            [nodePositions, nodeSizes] = obj.sphericalInitialization();

            % Get edges
            [sourceNodes, targetNodes] = obj.getEdges();

            % Iterative force calculation with convergence check
            for iter = 1:maxIterations
                % Initialize force accumulators
                forces = zeros(obj.numDendrites, 3);
                totalEnergy = 0; % Track total energy for convergence check

                % Calculate repulsive forces between all pairs of nodes
                for i = 1:obj.numDendrites
                    for j = i+1:obj.numDendrites
                        diff = nodePositions(i, :) - nodePositions(j, :);
                        dist = norm(diff);
                        if dist > 0
                            repulsiveForce = K_repulsion * (diff / dist^2);
                            forces(i, :) = forces(i, :) + repulsiveForce;
                            forces(j, :) = forces(j, :) - repulsiveForce;
                            totalEnergy = totalEnergy + K_repulsion / dist; % Repulsive energy
                        end
                    end
                end

                % Calculate attractive forces along the edges
                for e = 1:length(sourceNodes)
                    sourceIdx = obj.getDendriteIndex(sourceNodes(e));
                    targetIdx = obj.getDendriteIndex(targetNodes(e)); % Corrected line
                    diff = nodePositions(sourceIdx, :) - nodePositions(targetIdx, :);
                    dist = norm(diff);
                    if dist > 0 % Prevent division by zero
                        attractiveForce = -K_attraction * (diff / dist); % Hooke's law
                        forces(sourceIdx, :) = forces(sourceIdx, :) + attractiveForce;
                        forces(targetIdx, :) = forces(targetIdx, :) - attractiveForce;
                        totalEnergy = totalEnergy + 0.5 * K_attraction * dist^2; % Attractive energy
                    end
                end

                % Update node positions based on forces, but don't update the soma (ID = 0)
                for i = 2:obj.numDendrites  % Start from node 2 to skip the soma
                    nodePositions(i, :) = nodePositions(i, :) + timestep * forces(i, :);
                end

                % Keep the soma fixed at the origin
                nodePositions(1, :) = [0, 0, 0];

                % Check for convergence: if total energy is below epsilon, stop iterating
                if totalEnergy <= eps
                    fprintf('Converged after %d iterations with energy: %f\n', iter, totalEnergy);
                    break;
                end
            end

            % Prepare colors for nodes
            nodeColors = zeros(obj.numDendrites, 3); % RGB colors

            % Define colors
            somaColor = [1, 0, 0]; % Red for soma
            stimulatedColor = [0, 1, 0]; % Green for stimulated dendrites
            defaultColor = [0, 0, 1]; % Blue for other dendrites
            for i = 1:obj.numDendrites
                if obj.dendrites(i).ID == 0
                    % Soma
                    nodeColors(i, :) = somaColor;
                    nodeSizes(i) = nodeSizes(i) * 2; % Make soma larger
                elseif ~isempty(obj.Stimuli{i})
                    % Dendrites receiving stimuli
                    nodeColors(i, :) = stimulatedColor;
                else
                    % Default dendrites
                    nodeColors(i, :) = defaultColor;
                end
            end

            % Plot the resulting layout
            figure;
            scatter3(nodePositions(:, 1), nodePositions(:, 2), nodePositions(:, 3), nodeSizes, nodeColors, 'filled');
            hold on;

            % Add labels
            % Prepare labels
            labels = cell(obj.numDendrites, 1);
            for i = 1:obj.numDendrites
                ID = obj.dendrites(i).ID;
                if ID == 0
                    labels{i} = 'Soma';
                else
                    labels{i} = sprintf('Dendrite %d', ID);
                end
            end
            % Add labels
            for i = 1:obj.numDendrites
                text(nodePositions(i, 1), nodePositions(i, 2), nodePositions(i, 3), [' ' labels{i}], 'FontSize', 10);
            end

            % Plot edges
            for e = 1:length(sourceNodes)
                sourceIdx = obj.getDendriteIndex(sourceNodes(e));
                targetIdx = obj.getDendriteIndex(targetNodes(e));
                plot3([nodePositions(sourceIdx, 1), nodePositions(targetIdx, 1)], ...
                    [nodePositions(sourceIdx, 2), nodePositions(targetIdx, 2)], ...
                    [nodePositions(sourceIdx, 3), nodePositions(targetIdx, 3)], 'k-');
            end

            title('Force-Directed 3D Visualization of Dendritic Tree');
            grid on;
            axis equal;
            rotate3d on;
        end

        function [nodePositions, nodeSizes] = sphericalInitialization(obj)
            % Spherical initialization: evenly distribute nodes in 3D space around a sphere
            %
            % Outputs:
            % - nodePositions: Nx3 array of initial 3D positions for each node
            % - nodeSizes: Node sizes for visualization

            numNodes = obj.numDendrites;
            nodePositions = zeros(numNodes, 3);
            nodeSizes = zeros(numNodes, 1);

            % Parameters
            radius = 10; % Arbitrary radius for the sphere

            % Initialize positions for each node on a sphere (except soma)
            for i = 2:numNodes
                theta = pi * rand();   % Elevation angle
                phi = 2 * pi * rand(); % Azimuthal angle
                r = radius * (0.8 + 0.2 * rand()); % Random distance from center (slightly varied)

                nodePositions(i, 1) = r * sin(theta) * cos(phi); % x-coordinate
                nodePositions(i, 2) = r * sin(theta) * sin(phi); % y-coordinate
                nodePositions(i, 3) = r * cos(theta);            % z-coordinate
            end

            % Set the soma (ID = 0) at the center
            nodePositions(1, :) = [0, 0, 0];

            % Assign node sizes (optional, can be based on other properties)
            for i = 1:numNodes
                if obj.dendrites(i).ID == 0
                    nodeSizes(i) = 100; % Larger size for the soma
                else
                    nodeSizes(i) = 100*obj.dendrites(i).params.NaX; % Default size for other dendrites
                end
            end
        end

        %%%%%%% 3D tree %%%%%%%

        function [nodePositions, nodeSizes] = assignSpatialCoordinates(obj)
            % Initialize positions and sizes arrays
            numNodes = obj.numDendrites;
            nodePositions = zeros(numNodes, 3);
            nodeSizes = zeros(numNodes, 1); % For storing node sizes

            % Start recursive assignment from the soma (ID = 0)
            somaIdx = obj.getDendriteIndex(0);
            nodePositions(somaIdx, :) = [0, 0, 0]; % Soma at origin
            nodeSizes(somaIdx) = 100; % Base size for the soma

            % Define initial angle
            initialAngle = [pi/2; 0]; % Elevation; Azimuth in radians

            % Recursive function to assign positions
            function assignCoordinates(nodeID, parentPos, depth, angle)
                idx = obj.getDendriteIndex(nodeID);
                % Get child nodes (Distal connections)
                conn = obj.Connectivity(idx);
                childIDs = conn.List_Distal;
                numChildren = length(childIDs);
                if numChildren > 0
                    % Distribute angles for child branches
                    angleStep = pi / 4; % Adjust as needed
                    childAngles = linspace(-angleStep, angleStep, numChildren);
                    for i = 1:numChildren
                        childID = childIDs(i);
                        childIdx = obj.getDendriteIndex(childID);
                        % Retrieve dendrite parameters
                        dendrite = obj.dendrites(childIdx);
                        gc = dendrite.params.gc;
                        NaX = dendrite.params.NaX;

                        % Calculate node size based on NaX
                        nodeSizes(childIdx) = NaX * 100; % Scale factor for visualization

                        % Calculate branch length based on gc
                        baseLength = 1; % Base length
                        branchLength = baseLength / gc; % Inverse relation for visualization

                        % Calculate child angle
                        childAngle = angle + [pi/6; childAngles(i)]; % Elevation; Azimuth

                        % Calculate child position
                        x = parentPos(1) + branchLength * sin(childAngle(1)) * cos(childAngle(2));
                        y = parentPos(2) + branchLength * sin(childAngle(1)) * sin(childAngle(2));
                        z = parentPos(3) + branchLength * cos(childAngle(1));
                        childPos = [x, y, z];

                        % Store position
                        nodePositions(childIdx, :) = childPos;

                        % Recurse for child dendrites
                        assignCoordinates(childID, childPos, depth + 1, childAngle);
                    end
                end
            end

            % Start the recursion from the soma
            assignCoordinates(0, nodePositions(somaIdx, :), 0, initialAngle);
        end

        function [sourceNodes, targetNodes] = getEdges(obj)
            sourceNodes = [];
            targetNodes = [];
            for i = 1:obj.numDendrites
                currentID = obj.dendrites(i).ID;
                conn = obj.Connectivity(i);
                % Proximal connections (from current node to proximal nodes)
                for pid = conn.List_Proximal
                    sourceNodes(end+1) = currentID;
                    targetNodes(end+1) = pid;
                end
                % Distal connections (from current node to distal nodes)
                for did = conn.List_Distal
                    sourceNodes(end+1) = currentID;
                    targetNodes(end+1) = did;
                end
            end
        end
        function visualize3D(obj)
            % Assign spatial coordinates and node sizes
            [nodePositions, nodeSizes] = obj.assignSpatialCoordinates();

            % Get edges
            [sourceNodes, targetNodes] = obj.getEdges();

            % Create a mapping from IDs to indices
            idList = [obj.dendrites.ID];
            idToIndex = containers.Map(idList, 1:length(idList));

            % Prepare labels
            labels = cell(obj.numDendrites, 1);
            for i = 1:obj.numDendrites
                ID = obj.dendrites(i).ID;
                if ID == 0
                    labels{i} = 'Soma';
                else
                    labels{i} = sprintf('Dendrite %d', ID);
                end
            end

            % Prepare colors for nodes
            nodeColors = zeros(obj.numDendrites, 3); % RGB colors

            % Define colors
            somaColor = [1, 0, 0]; % Red for soma
            stimulatedColor = [0, 1, 0]; % Green for stimulated dendrites
            defaultColor = [0, 0, 1]; % Blue for other dendrites

            % Assign colors
            for i = 1:obj.numDendrites
                if idList(i) == 0
                    % Soma
                    nodeColors(i, :) = somaColor;
                    nodeSizes(i) = nodeSizes(i) * 2; % Make soma larger
                elseif ~isempty(obj.Stimuli{i})
                    % Dendrites receiving stimuli
                    nodeColors(i, :) = stimulatedColor;
                else
                    % Default dendrites
                    nodeColors(i, :) = defaultColor;
                end
            end

            % Plot nodes
            figure;
            scatter3(nodePositions(:,1), nodePositions(:,2), nodePositions(:,3), nodeSizes, nodeColors, 'filled');
            hold on;

            % Label nodes
            for i = 1:length(idList)
                text(nodePositions(i,1), nodePositions(i,2), nodePositions(i,3), sprintf('  %s', labels{i}), 'FontSize', 10);
            end

            % Plot edges
            numEdges = length(sourceNodes);
            for i = 1:numEdges
                sourceID = sourceNodes(i);
                targetID = targetNodes(i);
                sourceIdx = idToIndex(sourceID);
                targetIdx = idToIndex(targetID);
                sourcePos = nodePositions(sourceIdx, :);
                targetPos = nodePositions(targetIdx, :);
                % Adjust line width based on conductance (gc)
                gc = obj.dendrites(sourceIdx).params.gc;
                lineWidth = gc * 0.5; % Scale factor for visualization
                plot3([sourcePos(1), targetPos(1)], [sourcePos(2), targetPos(2)], [sourcePos(3), targetPos(3)], 'k-', 'LineWidth', lineWidth);
            end

            % Customize plot
            title('3D Visualization of Dendritic Tree');
            grid on;
            axis equal;
            rotate3d on;
        end

        %%%%% CPP stuff %%%%%%
        function [t, solution] = RK4CPP(obj, tmax)
            % RK4CPP prepares the data for the MEX function and runs the simulation
            %
            % Inputs:
            % - tmax: Maximum simulation time
            %
            % Outputs:
            % - t: Time vector
            % - solution: Simulated state variables over time

            % Prepare time vector
            h = obj.h;
            t = 0:h:tmax;
            numTimeSteps = length(t);

            % Get initial conditions
            X0 = obj.getInitCon(); % Should be a column vector

            % Prepare dendrite parameters
            numDendrites = obj.numDendrites;
            % Columns: [alpha, b, Tau, TauR, NaX, gc, ID]
            paramsArray = zeros(numDendrites, 7);

            for idx = 1:numDendrites
                dendrite = obj.dendrites(idx);
                paramsArray(idx, :) = [
                    dendrite.params.alpha, dendrite.params.b, dendrite.params.Tau, ...
                    dendrite.params.TauR, dendrite.params.NaX, dendrite.params.gc, dendrite.ID
                    ];
            end

            % Prepare coupling matrix data
            % Ensure the coupling matrix is built
            if isempty(obj.coupling_matrix)
                obj.coupling_matrix = obj.buildCouplingMatrix();
            end
            C = obj.coupling_matrix; % Sparse matrix

            % Convert the sparse matrix to COO format (row, col, value)
            [rowIndices, colIndices, values] = find(C);
            % Adjust indices for zero-based indexing in C++
            couplingData = [rowIndices - 1, colIndices - 1, values]; % Nx3 matrix

            % Prepare stimuli data
            % Identify dendrites with stimuli
            stimDendriteIDs = [];
            StimuliMatrix = [];
            stimDendriteIDtoIndex = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
            stimIdx = 1;

            for idx = 1:numDendrites
                StimSignal = obj.Stimuli{idx};
                if ~isempty(StimSignal)
                    dendriteID = obj.dendrites(idx).ID;
                    stimDendriteIDs(end+1) = dendriteID;
                    % Evaluate stimuli at all time steps
                    StimValues = zeros(1, numTimeSteps);
                    for k = 1:numTimeSteps
                        tt = t(k);
                        StimValue = StimSignal(tt);
                        if isnan(StimValue)
                            StimValue = 0;
                        end
                        StimValues(k) = StimValue;
                    end
                    % Append to StimuliMatrix
                    StimuliMatrix = [StimuliMatrix; StimValues];
                    % Map dendrite ID to row index in StimuliMatrix
                    stimDendriteIDtoIndex(dendriteID) = stimIdx;
                    stimIdx = stimIdx + 1;
                end
            end

            % Prepare mapping from dendrite IDs to indices (zero-based for C++)
            dendriteIDtoIndex = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
            for idx = 1:numDendrites
                dendriteID = obj.dendrites(idx).ID;
                % Subtract 1 for zero-based indexing in C++
                dendriteIDtoIndex(dendriteID) = idx - 1;
            end

            % Convert maps to arrays for passing to C++
            stimDendriteIDArray = stimDendriteIDs;
            dendriteIDArray = [obj.dendrites.ID];

            % Run the MEX function
            % Assuming the MEX function is named 'mexRK4Solver' and has the appropriate signature
            [t, solution] = mexRK4Solver(paramsArray, couplingData, StimuliMatrix, stimDendriteIDArray, X0, h, numTimeSteps, dendriteIDArray);

            % The solution returned should be an array of size [numStateVars x numTimeSteps]
            % where numStateVars = numDendrites * 2 (for u and v variables)

        end
    end
end
