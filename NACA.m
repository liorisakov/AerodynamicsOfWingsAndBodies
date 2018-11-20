classdef NACA < handle
    
    properties
        digits
        m
        p
        t
        
        xU
        xL
        yU
        yL
        
        yU_func
        yL_func
    end
    
    methods
        function obj = NACA(digits)
        if ~(ischar(digits) || isstring(digits))
            msg1 = 'NACA profile must be inititalized with ';
            msg2 = 'digit string or character array';
            msg = strcat(msg1, msg2);
            error(msg)
        end
        
        obj.digits = digits;
        if length(digits) == 4
            first_digit = str2num(digits(1));
            second_digit = str2num(digits(2));
            last_digits = str2num(digits(3:4));
        end
        obj.m = first_digit/100;    % max camber
        obj.p = second_digit/10;    % max camber location
        obj.t = last_digits/100;    % max thickness as fraction of cord
        
        [obj.yU, obj.yL, obj.xU, obj.xL] = obj.generate_geometry();
        end
        
        
        function yt = thickness(obj, x)
            yt = 5 * obj.t * (0.2969 * sqrt(x) - 0.1260 * x - ...
                              0.3516 * x.^2 + 0.2843 * x.^3 - ...
                              0.1015 * x.^4);
        end
        
        
        function dyc_dx = camber_slope(obj, x)
            x_LE = x(x <= obj.p);    % section closer to leading edge
            x_TE = x(x > obj.p);     % section closer to trailing edge
            
            dyc_dx_LE = 2*obj.m/obj.p^2 * (obj.p - x_LE);
            dyc_dx_TE = 2*obj.m/(1 - obj.p)^2 * (obj.p - x_TE);
            if size(x, 1) > size(x, 2)    % concatenate correctly
                dyc_dx = [dyc_dx_LE; dyc_dx_TE];
            else
                dyc_dx = [dyc_dx_LE, dyc_dx_TE];
            end
        end
        
        
        function yc = camber_line(obj, x)
            x_LE = x(x <= obj.p);    % section closer to leading edge
            x_TE = x(x > obj.p);     % section closer to trailing edge
            
            yc_LE = obj.m/obj.p^2 * (2*obj.p*x_LE - x_LE.^2);
            yc_TE = obj.m/(1 - obj.p)^2 * ...
                    ((1 - 2*obj.p) + 2*obj.p*x_TE - x_TE.^2);
            if size(x, 1) > size(x, 2)    % concatenate correctly
                yc = [yc_LE; yc_TE];
            else
                yc = [yc_LE, yc_TE];
            end
        end
        
        
        function [yU, yL, xU, xL] = generate_geometry(obj, x)
            if nargin == 1
                x = linspace(0, 1, 200);
            end
            
            % geometry properties
            yt = obj.thickness(x);
            dyc_dx = obj.camber_slope(x);
            yc = obj.camber_line(x);
            
            % geometry
            theta = atan(dyc_dx);
            xU = x - yt.*sin(theta);
            xL = x + yt.*sin(theta);
            yU = yc + yt.*cos(theta);
            yL = yc - yt.*cos(theta);
        end
        
        
        function plot_self(obj, transparent)
            if nargin == 1
                transparent = false;
            end
            load('colors.mat', 'colors')

            hold on
            if ~transparent
                patch([obj.xU, fliplr(obj.xL)], [obj.yU, fliplr(obj.yL)], ...
                  white(1))
            end
            patch([obj.xU, fliplr(obj.xL)], [obj.yU, fliplr(obj.yL)], ...
                  colors.blue, 'FaceAlpha', 0.2)
            plot([obj.xU, fliplr(obj.xL)], [obj.yU, fliplr(obj.yL)], ...
                 'LineWidth', 1.2, 'Color', colors.blue); 
            
            axis equal
        end

    end
end