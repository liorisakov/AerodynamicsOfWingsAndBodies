classdef VLMTile
    properties
        bottomLeft
        bottomRight
        topRight
        topLeft
        
        horseshoe
        controlPoint
    end
    
    methods
        function obj = VLMTile(bottomLeft, bottomRight, topRight, topLeft)
            % all vertices named according to top view when wing faces left
            obj.bottomLeft = bottomLeft;
            obj.bottomRight = bottomRight;
            obj.topRight = topRight;
            obj.topLeft = topLeft;
            
            obj.horseshoe.left = 3/4*bottomLeft + 1/4*bottomRight;
            obj.horseshoe.right = 3/4*topLeft + 1/4*topRight;
            obj.horseshoe.middleX = (obj.horseshoe.left(1) + ...
                                     obj.horseshoe.right(1))/2;
            obj.horseshoe.gamma = NaN;    % magnitude not yet calculated
            
            bottom3_4thsChord = 1/4*bottomLeft + 3/4*bottomRight;
            top3_4thsChord = 1/4*topLeft + 3/4*topRight;
            obj.controlPoint.x = (bottom3_4thsChord(1) + top3_4thsChord(1))/2;
            obj.controlPoint.y = (bottom3_4thsChord(2) + top3_4thsChord(2))/2;
        end
        
        
        function w = Downwash(obj, pointsX, pointsY, varargin)
            % handle input
            solveForGamma = false;
            for i = 1:2:length(varargin)
                if strcmp(varargin{i}, 'solve for gamma')
                    solveForGamma = varargin{i+1};
                end
            end
            
            if ~solveForGamma && isnan(obj.horseshoe.gamma)
                error("Vortex line magnitude hasn't been defined yet");
            end
            
%             P = points;
            A = obj.horseshoe.left;
            B = obj.horseshoe.right;
            
%             xP = P(:,1).';
%             yP = P(:,2).';
            xP = pointsX;
            yP = pointsY;
            xA = A(1);
            yA = A(2);
            xB = B(1);
            yB = B(2);
            
            w = 1/(4*pi) * ...
                (1./((xP-xA).*(yP-yB) - (xP-xB).*(yP-yA)) .* ...
                 (((xB-xA)*(xP-xA) + (yB-yA)*(yP-yA))./sqrt((xP-xA).^2 + (yP-yA).^2) - ... % vecnorm((P-A).') - ...
                  ((xB-xA)*(xP-xB) + (yB-yA)*(yP-yB))./sqrt((xP-xB).^2 + (yP-yB).^2)) + ... % vecnorm((P-B).') - ...
                 1./(yA-yP).*(1 + (xP-xA)./sqrt((xP-xA).^2 + (yP-yA).^2)) - ... % vecnorm((P-A).') - ...
                 1./(yB-yP).*(1 + (xP-xB)./sqrt((xP-xB).^2 + (yP-yB).^2))); % vecnorm((P-B).') - ...
             
            if ~solveForGamma
                w = obj.horseshoe.gamma * w;
            end
        end
        
        
        function handles = PlotSelf(obj, varargin)
            % handle input
            rotate = false;
            normalize = false;
            origin = [0, 0];
            plotHorseshoe = false;
            plotControlPoint = false;
            for i = 1:2:length(varargin)
                if strcmp(varargin{i}, 'rotate')
                    rotate = varargin{i+1};
                elseif strcmp(varargin{i}, 'normalize')
                    normalize = varargin{i+1};
                elseif strcmp(varargin{i}, 'origin')
                    origin = varargin{i+1};
                elseif strcmp(varargin{i}, 'horseshoe')
                    plotHorseshoe = varargin{i+1};
                elseif strcmp(varargin{i}, 'control point')
                    plotControlPoint = varargin{i+1};
                end
            end
            
            % plot parameters
            lw = 0.9;    % line width
%             ms = 7;      % marker size
            load('colors.mat', 'colors')
            
            % tile
            x = [obj.bottomLeft(1), obj.bottomRight(1), ...
                 obj.topRight(1), obj.topLeft(1), obj.bottomLeft(1)];
            y = [obj.bottomLeft(2), obj.bottomRight(2), ...
                 obj.topRight(2), obj.topLeft(2), obj.bottomLeft(2)];
            if normalize
                x = x/normalize;
                y = y/normalize;
            end
            if rotate
                handles.tile = ...
                    plot(y + origin(1), -x + origin(2), ...
                         'LineWidth', lw, 'Color', colors.blue);
            else
                handles.tile = ...
                    plot(x + origin(1), y + origin(2), ...
                         'LineWidth', lw, 'Color', colors.blue);
            end
            
            % horseshoe
            if plotHorseshoe
                middleX = [obj.horseshoe.left(1), obj.horseshoe.right(1)];
                middleY = [obj.horseshoe.left(2), obj.horseshoe.right(2)];
                leftX = [obj.horseshoe.left(1), 15];
                leftY = [obj.horseshoe.left(2), obj.horseshoe.left(2)];
                rightX = [obj.horseshoe.right(1), 15];
                rightY = [obj.horseshoe.right(2), obj.horseshoe.right(2)];
                
                if normalize
                    middleX = middleX/normalize;
                    middleY = middleY/normalize;
                    leftX = leftX/normalize;
                    leftY = leftY/normalize;
                    rightX = rightX/normalize;
                    rightY = rightY/normalize;
                end
                if rotate
                    handles.horseshoe = ...
                        plot(middleY + origin(1), - middleX + origin(2), ...
                             '--', 'Color', colors.red);
                    plot(leftY + origin(1), - leftX + origin(2), ...
                         '--', 'Color', colors.red)
                    plot(rightY + origin(1), - rightX + origin(2), ...
                         '--', 'Color', colors.red)
                else
                    handles.horseshoe = ...
                        plot(middleX + origin(1), middleY + origin(2), ...
                             '--', 'Color', colors.red);
                    plot(leftX + origin(1), leftY + origin(2), ...
                         '--', 'Color', colors.red)
                    plot(rightX + origin(1), rightY + origin(2), ...
                         '--', 'Color', colors.red)
                end
            end
            
            % control point
            if plotControlPoint
                controlPointX = obj.controlPoint.x;
                controlPointY = obj.controlPoint.y;
                
                if normalize
                    controlPointX = controlPointX/normalize;
                    controlPointY = controlPointY/normalize;
                end
                if rotate
                    handles.controlPoint = ...
                        plot(controlPointY + origin(1), ...
                             - controlPointX + origin(2), ...
                            '.', 'Color', colors.purple);
                else
                    handles.controlPoint = ...
                        plot(controlPointX + origin(1), ...
                             controlPointY + origin(2), ...
                            '.', 'Color', colors.purple);
                end
            end
        end    
    end
end