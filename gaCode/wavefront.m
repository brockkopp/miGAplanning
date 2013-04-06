% Wavefront example
close all; clc

m=1;
% Vehicle start and end position
startPos = startPt;
endPos = endPt;

%% Set up environment
map = obsGrid(50:80, 50:90);

% Region Bounds
posMinBound = [0 0];
posMaxBound = size(map);

% Plot environment
figure(1); clf;
hold on;
axis([posMinBound(1) posMaxBound(1) posMinBound(2) posMaxBound(2)]);
imagesc(1-map');
colormap('bone')

%% Graph setup
N = posMaxBound(1);
M = posMaxBound(2);

rx = posMaxBound(1)-posMinBound(1);
ry = posMaxBound(2)-posMinBound(2);

for i = 1:N
    for j = 1:M
        x(i,j) = posMinBound(1) + i/N*rx;
        y(i,j) = posMinBound(2) + j/M*ry;
    end
end

% Start node
[startx, starti] = min(abs(x(:,1)-startPos(1)));
[starty, startj] = min(abs(y(1,:)-startPos(2)));
start = (M-1)*starti + startj;

% End node
[finishx, finishi] = min(abs(x(:,1)-endPos(1)));
[finishy, finishj] = min(abs(y(1,:)-endPos(2)));
finish = M*(finishi-1) + finishj;

% Node numbering and graph connectivity
e = sparse(N*M,N*M);
% For each cell in map
for i=1:N
    for j = 1:M
        % If cell is empty
        if (map(i,j) == 0)
            cur = M*(i-1)+j;
            % Link up if empty
            if (i>1)
                if (map(i-1,j) == 0)
                    e(cur, cur-M) = 1;
                    e(cur-M,cur) = 1;
                end
            end
            % Link left
            if (j>1)
                if (map(i,j-1) == 0)
                    e(cur, cur-1) = 1;
                    e(cur-1,cur) = 1;
                end
            end    
            % Link down
            if (i<N)
                if (map(i+1,j) == 0)
                    e(cur, cur+M) = 1;
                    e(cur+M,cur) = 1;
                end
            end      
            % Link right
            if (j<M)
                if (map(i,j+1) == 0)
                    e(cur, cur+1) = 1;
                    e(cur+1,cur) = 1;
                end
            end   
            % Link up-left if empty
            if (i>1 && j>1)
                if (map(i-1,j-1) == 0)
                    e(cur, cur-M-1) = 1;
                    e(cur-M-1,cur) = 1;
                end
            end      
            % Link down-left if empty
            if (i<N && j>1)
                if (map(i+1,j-1) == 0)
                    e(cur, cur+M-1) = 1;
                    e(cur+M-1,cur) = 1;
                end
            end
            % Link down-right if empty
            if (i<N && j<M)
                if (map(i+1,j+1) == 0)
                    e(cur, cur+M+1) = 1;
                    e(cur+M+1,cur) = 1;
                end
            end
            % Link up-right if empty
            if (i>1 && j<M)
                if (map(i-1,j+1) == 0)
                    e(cur, cur-M+1) = 1;
                    e(cur-M+1,cur) = 1;
                end
            end
        end % if empty
    end % j
end % i

%% Wavefront
% Essentially a breadth-first search for the cells that can be reached

% Initialize open set (node cost)
O = [finish 0];
% Initialize closed set (same form as open set)
C = [];
done = 0;
t = 1;

dist = 0;
wave = 2*(N+M)*map;

while (~done)
    % Check end condition
    if (length(O)==0)
        done = 1;
        continue;
    end

    % Grab next node in open set
    curnode = O(1,:);
    
    % Move to closed set and save distance for plotting
    C = [C; curnode];
    O = O([2:end],:);
    curi = floor((curnode(1)-1)/M)+1;
    curj = mod(curnode(1),M);
    if (curj==0) curj=M; end
    wave(curi,curj) = curnode(2);

    % Get all neighbours of current node
    neigh = find(e(curnode(1),:)==1);
    
    % Process each neighbour
    for i=1:length(neigh)
        % If neighbour is already in closed list, skip it
        found = find(C(:,1)==neigh(i),1);
        if (length(found)==1)
            continue;
        end
        % If neighbour is already in open list, skip it
        found = find(O(:,1)==neigh(i),1);
        % Otherwise, add to open list at the bottom
        if (length(found)==0)
            O = [O; neigh(i) curnode(2)+1]; 
        end
    end
    if (curnode(2) > dist)
        figure(4); clf; hold on;
        axis([0 N 0 M]);
        colormap('default')
        imagesc(wave', [0 1.5*(M+N)])
        plot(finishi,finishj,'r*');
        plot(starti,startj,'b*');
        dist = curnode(2);
%         F(t) = getframe(gcf);
%         t = t+1;
    end
end

%% Shortest path

len = wave(starti,startj);
path = zeros(len,2);
path(1,:) = [starti startj];
for i=1:len
    options = [];
    if (path(i,1)>1) %left
        options = [options; [path(i,1)-1 path(i,2)]];
    end
    if (path(i,1)<N) %right
        options = [options; [path(i,1)+1 path(i,2)]];
    end
    if (path(i,2)>1) %down
        options = [options; [path(i,1) path(i,2)-1]];
    end
    if (path(i,2)<M) %up
        options = [options; [path(i,1) path(i,2)+1]];
    end
    if (path(i,1)>1 && path(i,2)>1) %%left-down
        options = [options; [path(i,1)-1 path(i,2)-1]];
    end
    if (path(i,1)>1 && path(i,2)<M) %%left-up
        options = [options; [path(i,1)-1 path(i,2)+1]];
    end
    if (path(i,1)<N && path(i,2)>1)  %%right-down
        options = [options; [path(i,1)+1 path(i,2)-1]];
    end
    if (path(i,1)<N && path(i,2)<M)  %%right-up
        options = [options; [path(i,1)+1 path(i,2)+1]];
    end
    oplen = length(options(:,1));
    for j = 1:oplen
        costs(j) = wave(options(j,1),options(j,2));
    end
    [dist, best] = min(costs);
    path(i+1,:) = options(best,:);
end

% figure(4); hold on;
figure(4); clf; hold on;
axis([0 N 0 M]);
colormap('default')
imagesc(wave', [0 1.5*(M+N)])
plot(finishi,finishj,'r*');
plot(starti,startj,'b*');
dist = curnode(2);
plot(path(:,1),path(:,2), 'rx-')

% for i=1:100
%     F(t) = getframe(gcf);
%     t=t+1;
% end

