classdef VLMTile < handle
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
            obj.horseshoe.middle = (obj.horseshoe.left + ...
                                    obj.horseshoe.right)/2;
            obj.horseshoe.gamma = NaN;    % magnitude not yet calculated
            
            bottom3_4thsChord = 1/4*bottomLeft + 3/4*bottomRight;
            top3_4thsChord = 1/4*topLeft + 3/4*topRight;
            obj.controlPoint.x = (bottom3_4thsChord(1) + top3_4thsChord(1))/2;
            obj.controlPoint.y = (bottom3_4thsChord(2) + top3_4thsChord(2))/2;
        end
        
        
        function w = Downwash(obj, pointsX, pointsY, varargin)
            % handle input
            solveForGamma = false;
            trailingOnly = false;
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'solve for gamma')
                    solveForGamma = true;
                elseif strcmp(varargin{i}, 'trailing only')
                    trailingOnly = true;
                end
            end
            
            if ~solveForGamma && isnan(obj.horseshoe.gamma)
                error("Vortex line magnitude hasn't been defined yet");
            end
            
            A = obj.horseshoe.left;
            B = obj.horseshoe.right;
            
            xP = pointsX;
            yP = pointsY;
            xA = A(1);
            yA = A(2);
            xB = B(1);
            yB = B(2);
            
            if ~trailingOnly
                w = 1/(4*pi) * ...
                    (1./((xP-xA).*(yP-yB) - (xP-xB).*(yP-yA)) .* ...
                     (((xB-xA)*(xP-xA) + (yB-yA)*(yP-yA))./sqrt((xP-xA).^2 + (yP-yA).^2) - ...
                      ((xB-xA)*(xP-xB) + (yB-yA)*(yP-yB))./sqrt((xP-xB).^2 + (yP-yB).^2)) + ...
                     1./(yA-yP).*(1 + (xP-xA)./sqrt((xP-xA).^2 + (yP-yA).^2)) - ...
                     1./(yB-yP).*(1 + (xP-xB)./sqrt((xP-xB).^2 + (yP-yB).^2)));
            else
                w = 1/(4*pi) * ...
                    (1./(yA-yP).*(1 + (xP-xA)./sqrt((xP-xA).^2 + (yP-yA).^2)) - ...
                     1./(yB-yP).*(1 + (xP-xB)./sqrt((xP-xB).^2 + (yP-yB).^2)));
            end
             
            if ~solveForGamma
                w = obj.horseshoe.gamma * w;
            end
        end
        
        
        function handles = PlotSelf(obj, varargin)
            % plot parameters
            lw = 0.9;    % line width
            load('colors.mat', 'colors')
            
            % handle input
            rotate = false;
            normalize = false;
            origin = [0, 0];
            plotHorseshoe = false;
            plotControlPoint = false;
            color = colors.blue;
            zOffset = 0;
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
                elseif strcmp(varargin{i}, 'color')
                    color = varargin{i+1};
                elseif strcmp(varargin{i}, 'z offset')
                    zOffset = varargin{i+1};
                end
            end
            
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
                    plot3(y + origin(1), ...
                          - x + origin(2), ...
                          zoffset*ones(size(x)), ...
                          'LineWidth', lw, 'Color', color);
            else
                handles.tile = ...
                    plot3(x + origin(1), ...
                          y + origin(2), ...
                          zOffset*ones(size(x)), ...
                          'LineWidth', lw, 'Color', color);
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
                        plot3(middleY + origin(1), ...
                              - middleX + origin(2), ...
                              zOffset*ones(size(middleX)), ...
                              '--', 'Color', colors.red);
                    plot3(leftY + origin(1), ...
                          - leftX + origin(2), ...
                          zOffset*ones(size(leftX)), ...
                          '--', 'Color', colors.red)
                    plot3(rightY + origin(1), ...
                          - rightX + origin(2), ...
                          zOffset*ones(size(rightX)), ...
                          '--', 'Color', colors.red)
                else
                    handles.horseshoe = ...
                        plot3(middleX + origin(1), ...
                              middleY + origin(2), ...
                              zOffset*ones(size(middleX)), ...
                              '--', 'Color', colors.red);
                    plot3(leftX + origin(1), ...
                          leftY + origin(2), ...
                          zOffset*ones(size(leftX)), ...
                          '--', 'Color', colors.red)
                    plot3(rightX + origin(1), ...
                          rightY + origin(2), ...
                          zOffset*ones(size(rightX)), ...
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
                        plot3(controlPointY + origin(1), ...
                              - controlPointX + origin(2), ...
                              zOffset, ...
                              '.', 'Color', colors.yellow);
                else
                    handles.controlPoint = ...
                        plot3(controlPointX + origin(1), ...
                              controlPointY + origin(2), ...
                              zOffset, ...
                              '.', 'Color', colors.yellow);
                end
            end
        end    
    end
end