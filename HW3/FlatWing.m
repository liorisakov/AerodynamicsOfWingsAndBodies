classdef FlatWing < handle
    properties
        span
        area
        aspectRatio
        
        mac
        macY
        macLeadingEdgeX
        
        chordLeadingEdges % (xr, yr; x1, y1; ...; xt, yt)
        chordLengths
        spanSectionCount
        spanSectionLengths
        
        N
        tileDivision
        tiles
        controlPointsX
        controlPointsY
    end
    
    methods
        function obj = FlatWing(chordLeadingEdges, chordLengths, N)
            if any(N ~= floor(N)) || any(length(N) > 2)
                error("Argument 'N' must be a 1- or 2-length integer vector")
            elseif size(chordLeadingEdges, 2) ~= 2
                msg1 = "Argument 'chordLeadingEdges' must be an Mx2 matrix,";
                msg2 = ' where M is the number of specified leading edges';
                msg = strcat(msg1, msg2);
            elseif size(chordLeadingEdges, 1) ~= length(chordLengths)
                msg1 = "Arguments 'chordLeadingEdges' and 'chordLengths' ";
                msg2 = 'sizes must agree';
                msg = strcat(msg1, msg2);
                error(msg)
            end
            
            obj.chordLeadingEdges = chordLeadingEdges;
            obj.chordLengths = chordLengths;
            
            % for following calculations
            iChords = chordLengths(1:end-1);
            iPlusOneChords = chordLengths(2:end);
            iSpans = chordLeadingEdges(1:end-1, 2).';
            iPlusOneSpans = chordLeadingEdges(2:end, 2).';
            
            obj.spanSectionCount = length(chordLeadingEdges) - 1;
            obj.spanSectionLengths = iPlusOneSpans - iSpans;
            
            obj.span = 2*sum(obj.spanSectionLengths);
            obj.area = 2*sum(obj.spanSectionLengths/2 .* ...
                             (iPlusOneChords + iChords));
            obj.aspectRatio = obj.span^2/obj.area;
            
%             obj.mac = obj.span/obj.area * ...
%                       sum(iChords + (obj.span/4 - iSpans) .* ...
%                           (iPlusOneChords - iChords) ./ ...
%                           (iPlusOneSpans - iSpans));
            li = iChords;
            lip1 = iPlusOneChords;
            yi = iSpans;
            yip1 = iPlusOneSpans;
            m = (lip1 - li)./(yip1 - yi);
            s = obj.area;
            obj.mac = 2/s * ...
                      sum(li.^2.*yip1 + ...
                          2*m.*li.*(yip1.^2/2 - yi.*yip1) + ...
                          m.^2.*(yip1.^3/3 - yi.*yip1.^2 + yi.^2.*yip1) - ...
                          (li.^2.*yi + ...
                           2*m.*li.*(yi.^2/2 - yi.*yi) + ...
                           m.^2.*(yi.^3/3 - yi.*yi.^2 + yi.^2.*yi)));
                          
            obj.macY = interp1(chordLengths, ...
                               chordLeadingEdges(:,2), obj.mac);
            obj.macLeadingEdgeX = interp1(chordLengths, ...
                                          chordLeadingEdges(:,1), obj.mac);
            
            obj.N = N;
            
            obj.verifyDivisionToTiles();
            obj.createTiles();
            
        end
        
        
        function [cl, cmApex] = AeroCoeffsAtY(obj, y, uoo, varargin)
            % handle input
            clTimesC = false;
            cmTimesCSquared = false;
            returnSeparateXValues = false;
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'cl times c')
                    clTimesC = true;
                elseif strcmp(varargin{i}, 'cm times c squared')
                    cmTimesCSquared = true;
                elseif strcmp(varargin{i}, 'return separate x values')
                    returnSeparateXValues = true;
                end
            end
            
            xTileCount = obj.N(1);
            yTileCount = 2*obj.N(2);
            
            tilesAlongY = cell(xTileCount, length(y));
            for k = 1:length(y)
                for j = 1:yTileCount
                    if obj.tiles{1,j}.bottomLeft(2) <= y(k) && ...
                       y(k) <= obj.tiles{1,j}.topLeft(2)
                        tilesAlongY(:,k) = {obj.tiles{:,j}};
                        break
                    end
                end
            end
            
            if ~returnSeparateXValues
                [gammaValues, gammaXoValues] = deal(zeros(1, length(y)));
            else
                [gammaValues, gammaXoValues] = deal(zeros(xTileCount, length(y)));
            end
            for k = 1:length(y)
                for i = 1:xTileCount
                    currentHorseshoe = tilesAlongY{i,k}.horseshoe;
                    
                    gamma = currentHorseshoe.gamma;
                    horseshoeMiddleX = currentHorseshoe.middleX;
                    
                    if ~returnSeparateXValues
                        gammaValues(k) = gammaValues(k) + gamma;
                        gammaXoValues(k) = gammaXoValues(k) + ...
                                        gamma*horseshoeMiddleX;
                    else
                        gammaValues(i,k) = gamma;
                        gammaXoValues(i,k) = gamma*horseshoeMiddleX;
                    end
                end
            end
            
            if clTimesC
                cl = 2./uoo.*gammaValues;
            else
                c = obj.GetChord(y);
                cl = 2./(uoo.*c).*gammaValues;
            end
            if cmTimesCSquared
                cmApex = - 2./uoo.*gammaXoValues;
            else
                cSquared = obj.GetChord(y).^2;
                cmApex = - 2./(uoo.*cSquared).*gammaXoValues;
            end
        end
        
        
        function [cL, cMApex] = AeroCoeffs(obj, y, uoo)
            xTileCount = obj.N(1);
            yTileCount = 2*obj.N(2);
            
            tilesAlongY = cell(xTileCount, length(y));
            for k = 1:length(y)
                for j = 1:yTileCount
                    if obj.tiles{1,j}.bottomLeft(2) <= y(k) && ...
                       y(k) <= obj.tiles{1,j}.topLeft(2)
                        tilesAlongY(:,k) = {obj.tiles{:,j}};
                        break
                    end
                end
            end
            
            sumGammaDeltaY = 0;
            sumGammaDeltaYXo = 0;
            for k = 1:length(y)
                for i = 1:xTileCount
                    currentHorseshoe = tilesAlongY{i,k}.horseshoe;
                    
                    gamma = currentHorseshoe.gamma;
                    deltaY = currentHorseshoe.right(2) - ...
                             currentHorseshoe.left(2);
                    horseshoeMiddleX = currentHorseshoe.middleX;
                    
                    sumGammaDeltaY = sumGammaDeltaY + gamma*deltaY;
                    sumGammaDeltaYXo = sumGammaDeltaYXo + ...
                                       gamma*deltaY*horseshoeMiddleX;
                end
            end
            
            cL = 2/(obj.area * uoo) * sumGammaDeltaY;
            cMApex = - 2/(obj.area * obj.mac * uoo) * sumGammaDeltaYXo;
        end
        
        
        function handles = PlotSelf(obj, varargin)
            hold on
            
            % handle input
            rotate = false;
            normalize = false;
            origin = [0, 0];
            plotTiles = false;
            plotHorseshoes = false;
            plotControlPoints = false;
            plotMac = false;
            noFill = false;
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'rotate')
                    rotate = true;
                elseif strcmp(varargin{i}, 'normalize')
                    normalize = obj.span/2;
                elseif strcmp(varargin{i}, 'origin')
                    origin = varargin{i+1};
                elseif strcmp(varargin{i}, 'tiles')
                    plotTiles = true;
                elseif strcmp(varargin{i}, 'horseshoes')
                    plotHorseshoes = true;
                elseif strcmp(varargin{i}, 'control points')
                    plotControlPoints = true;
                elseif strcmp(varargin{i}, 'MAC')
                    plotMac = true;
                elseif strcmp(varargin{i}, 'no fill')
                    noFill = true;
                end
            end
            
            % plot parameters
            lw = 1.2;    % line width
            ms = 8;      % marker size
            load('colors.mat', 'colors')
            
            xTileCount = obj.N(1);
            yTileCount = 2*obj.N(2);
            
            outlineX = [obj.chordLeadingEdges(:,1).', ...
                        fliplr(obj.chordLeadingEdges(:,1).' + obj.chordLengths), ...
                        obj.chordLeadingEdges(:,1).' + obj.chordLengths, ...
                        fliplr(obj.chordLeadingEdges(:,1).')];
            outlineY = [obj.chordLeadingEdges(:,2).', ...
                        fliplr(obj.chordLeadingEdges(:,2).'), ...
                       -obj.chordLeadingEdges(:,2).', ...
                       -fliplr(obj.chordLeadingEdges(:,2).')];
             if rotate
                 temp = outlineX;
                 outlineX = outlineY;
                 outlineY = -temp;
                 clear temp
             end
             if normalize
                 outlineX = outlineX./(obj.span/2);
                 outlineY = outlineY./(obj.span/2);
             end
                 
             if ~noFill
                 patch(outlineX + origin(1), outlineY + origin(2), ...
                      colors.blue, 'LineStyle', 'None', 'FaceAlpha', 0.2)
             end
            
            if plotTiles
                for i = 1:xTileCount
                    for j = 1:yTileCount
                        currentTile = obj.tiles{i,j};
                        tileHandles = currentTile.PlotSelf( ...
                            'rotate', rotate, ...
                            'normalize', normalize, ...
                            'horseshoe', plotHorseshoes, ...
                            'control point', plotControlPoints);
                    end
                end
                
                handles.tileHandles = tileHandles;
            end
            
            handles.outline = plot(outlineX + origin(1), ...
                                   outlineY + origin(2), ...
                                   'LineWidth', lw, 'Color', colors.blue);
            
            if plotMac
                macTrailingEdgeX = obj.macLeadingEdgeX + obj.mac;
                macQuarterX = 3/4*obj.macLeadingEdgeX + ...
                              1/4*macTrailingEdgeX;
                handles.mac = plot([obj.macLeadingEdgeX, macTrailingEdgeX], ...
                                   obj.macY*[1, 1], ...
                                   'LineWidth', lw, 'Color', colors.red);
                handles.quarterMac = plot(macQuarterX, obj.macY, ...
                                          'x', 'LineWidth', lw, ...
                                          'Color', colors.red);
            end
        end
        
        
        function [c, leadingEdgeX] = GetChord(obj, y)
            y = abs(y);
            c = interp1(obj.chordLeadingEdges(:,2), ...
                        obj.chordLengths, y);
            leadingEdgeX = interp1(obj.chordLeadingEdges(:,2), ...
                                  obj.chordLeadingEdges(:,1), y);
        end
    end
    
    
    methods (Access = private)
        function verifyDivisionToTiles(obj)
            if length(obj.N) == 1
                i = 1;
                divisor1 = round(sqrt(obj.N));
                divisorsFound = false;

                while ~divisorsFound
                    if mod(obj.N, divisor1) == 0
                        divisor2 = obj.N/divisor1;
                        division(1) = min([divisor1, divisor2]);
                        division(2) = max([divisor1, divisor2]);

                        if division(2) >= obj.spanSectionCount
                            divisorsFound = true;
                        end
                    end

                    if divisor1 > obj.N
                        msg1 = 'Could not find a proper division of';
                        msg2 = sprintf(' wing into %d tile(s); ', obj.N);
                        msg3 = ' specify a different amount';
                        msg = strcat(msg1, msg2, msg3);
                        error(msg)
                    end

                    divisor1 = divisor1 + i;
                    i = - sign(i) - i;
                end

                obj.N = division;
            end
            
            % divide along span such that each section gets a proper amount
            % of tiles, according to its percentage of the total span
            obj.tileDivision = {obj.N(1), zeros(size(obj.spanSectionLengths))};
            halfSpan = sum(obj.spanSectionLengths);
            for i = 1:length(obj.spanSectionLengths)
                obj.tileDivision{2}(i) = ...
                    round(obj.N(2)*obj.spanSectionLengths(i)/halfSpan);
            end
            if sum(obj.tileDivision{2}) > obj.N(2)
                addedTileCount = (sum(obj.tileDivision{2}) - obj.N(2))*obj.N(1);
                msg1 = sprintf('%d tile(s) added to', addedTileCount);
                msg2 = ' avoid odd division of spanwise wing sections';
                msg = strcat(msg1, msg2);
                warning(msg)
                
                obj.N(2) = sum(obj.tileDivision{2});
            end
        end
        
        
        function createTiles(obj)
            obj.controlPointsX = zeros(obj.N(1), obj.N(2));
            obj.controlPointsY = zeros(obj.N(1), obj.N(2));
            chordTrailingEdges = obj.chordLeadingEdges + ...
                                 [obj.chordLengths.', zeros(size(obj.chordLengths.'))];
            
            % section along chord
            chordCount = length(obj.chordLengths);
            sectionedChords = cell(1, chordCount);
            for i = 1:chordCount
                sectionedChords{i}(:,1) = linspace(obj.chordLeadingEdges(i,1), ...
                                                   chordTrailingEdges(i,1), ...
                                                   obj.tileDivision{1} + 1);
                sectionedChords{i}(:,2) = obj.chordLeadingEdges(i,2);
            end
            
            % section along span
            for i = 1:chordCount-1
                % tile vertices x and y values
                xValuesTemp = zeros(obj.tileDivision{1}+1, ...
                                    obj.tileDivision{2}(i)+1);
                yValuesTemp = zeros(size(xValuesTemp));
                for j = 1:obj.tileDivision{1}+1
                    xValuesTemp(j,:) = linspace(sectionedChords{i}(j,1), ...
                                                sectionedChords{i+1}(j,1), ...
                                                obj.tileDivision{2}(i)+1);
                    yValuesTemp(j,:) = linspace(sectionedChords{i}(j,2), ...
                                                sectionedChords{i+1}(j,2), ...
                                                obj.tileDivision{2}(i)+1);
                end
                
                if i == 1
                    xValues = xValuesTemp;
                    yValues = yValuesTemp;
                elseif i > 1
                    % avoid overlapping vertices
                    xValuesTemp = xValuesTemp(:,2:end);
                    yValuesTemp = yValuesTemp(:,2:end);
                    
                    xValues = [xValues, xValuesTemp];
                    yValues = [yValues, yValuesTemp];
                end
            end
            
            % generate tiles
            yTilesPerWing  = obj.N(2);
            obj.tiles = cell(obj.N(1), 2*obj.N(2));
            for i = 1:obj.N(1)
                for j = 1:obj.N(2)
                    % positive y half of wing
                    bottomLeftPos = [xValues(i,j), yValues(i,j)];
                    bottomRightPos = [xValues(i+1,j), yValues(i+1,j)];
                    topRightPos = [xValues(i+1,j+1), yValues(i+1,j+1)];
                    topLeftPos = [xValues(i,j+1), yValues(i,j+1)];
                    
                    % negative y half of wing
                    bottomLeftNeg = [topLeftPos(1), -topLeftPos(2)];
                    bottomRightNeg = [topRightPos(1), -topRightPos(2)];
                    topRightNeg = [bottomRightPos(1), -bottomRightPos(2)];
                    topLeftNeg = [bottomLeftPos(1), -bottomLeftPos(2)];
                    
                    % generate tiles
                    obj.tiles{i,j} = ...
                        VLMTile(bottomLeftPos, bottomRightPos, ...
                                topRightPos, topLeftPos);
                    obj.tiles{i,j+yTilesPerWing} = ...
                        VLMTile(bottomLeftNeg, bottomRightNeg, ...
                                topRightNeg, topLeftNeg);
                    
                    % save control points of positive y half of wing
                    obj.controlPointsX(i,j) = obj.tiles{i,j}.controlPoint.x;
                    obj.controlPointsY(i,j) = obj.tiles{i,j}.controlPoint.y;
                    
%                     k = k + 1;
                end
            end
        end
    end    
end