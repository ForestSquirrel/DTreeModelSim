classdef DendriteTreeModel
    properties
        dendrites       % Array of dendrite objects
        numDendrites    % Number of dendrites
        Stimuli         % Cell array of stimuli
        Connectivity    % Struct array representing connections
        h               % Integration step size
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
                    if isa(StimSignal, 'function_handle')
                        Stim = StimSignal(t);
                    elseif isnumeric(StimSignal) % Stimulus is an array
                        Stim = interp_point_by_point(StimSignal, t, obj.h);
                    end
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
        %%%%%%% VISUALIZATION %%%%%%%
        function G = makeGraph(obj)
            % Method to create a digraph object representing the dendritic tree
            sourceNodes = [];
            targetNodes = [];

            % Loop through dendrites to extract connections
            for i = 1:obj.numDendrites
                currentID = obj.dendrites(i).ID;
                conn = obj.Connectivity(i);

                % Proximal connections (edges from connected dendrites to current dendrite)
                for pid = conn.List_Proximal
                    sourceNodes(end+1) = pid;
                    targetNodes(end+1) = currentID;
                end

                % Distal connections (edges from current dendrite to connected dendrites)
                for did = conn.List_Distal
                    sourceNodes(end+1) = currentID;
                    targetNodes(end+1) = did;
                end
            end

            % Create a digraph object
            G = digraph(sourceNodes, targetNodes);

            % Assign node labels
            nodeIDs = [obj.dendrites.ID];
            nodeLabels = arrayfun(@(id) sprintf('Dendrite %d', id), nodeIDs, 'UniformOutput', false);
            somaIdx = find(nodeIDs == 0);
            if ~isempty(somaIdx)
                nodeLabels{somaIdx} = 'Soma'; % Label the soma
            end
            G.Nodes.Name = nodeLabels';
        end

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
                        nodeSizes(childIdx) = NaX * 1000; % Scale factor for visualization

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

    end
end
