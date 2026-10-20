% =================================================================================================================
% STRIKE POINT FINDER:

% This code finds the strike points in a simplified tokamak geometry.

% The strike points are defined as the intersection of the separatrix with the divertor plates.
% The strike points are represented by a normalised parameter along the divertor plates, 
% in such a way that they are the same in a symmetric configuration.
% The separatrix is defined as that magnetic field line which goes through the X-point.
% The X-point is defined as that point where the magnetic field is zero.
% The divertor plates are modelled as lines between a start and end point.
% The coils and plasma current are modelled as infinitely thin wires.

% The magnetic field is determined analytically from the set current configuration.
% The equation for the magnetic field lines is determined analytically as well.

% The position of the X-point is determined via a quasi Monte-Carlo simulation (random walk).
% The position of the strike points is determined via the same method.
% NOTE: strike point determination might be better implemented using Newton's method. This may be implemented later.

% Current implementation is focussed towards strike-point calculation and visualisation.
% For an integrated implementation, the visualisation may be removed.
% (which parts may be removed is indicated in the code)

% This code was translated from my own, original Python code into MatLab by Claude AI
% - Raf van der Veeken
% ===================================================================================================================

% Define a random seed
seed = 30122005;
rng(seed);

% The following inputs are all used for testing the function
% They are not important for calling the function

% Define the static currents
a0 = [0, 0];
I0 = 2;

a1 = [1, 1];
I1 = -0.5;

a2 = [-1, 1];
I2 = -0.5;

a3 = [0, 4];
I3 = 4;

Currents = [I0, I1, I2, I3];
Locations = [a0; a1; a2; a3];

% Define the inner and outer divertor plates
InStart = [-1, 1];
InEnd = [0, 0];
InPlate = [InStart; InEnd];

OutStart = [0, 0];
OutEnd = [1, 1];
OutPlate = [OutStart; OutEnd];

% Test the function by calling it with the above inputs
[sIn, sOut] = StrikePoints(Currents, Locations, InPlate, OutPlate);
fprintf('(%.17g, %.17g)\n', sIn, sOut);

% FUNCTION DEFINITION
function [sIn, sOut] = StrikePoints(Currents, Locations, InPlate, OutPlate)
    % Current values used in the field-line definition (same values as I0..I3 above)
    I0 = Currents(1);
    I1 = Currents(2);
    I2 = Currents(3);
    I3 = Currents(4);

    % Initialize the X-point and strike point first guesses
    Xpoint = [0, 1];
    % NOTE: in Python these are views into Locations; the flags below emulate that aliasing
    InStrikePoint = Locations(3, :);
    inStrikeIsView = true;
    OutStrikePoint = Locations(2, :);
    outStrikeIsView = true;

    % Define the slope and intersect of the divertor plates
    InSlope = (InPlate(2, 2) - InPlate(1, 2)) / (InPlate(2, 1) - InPlate(1, 1));
    InIntersect = InPlate(1, 2) - InSlope * InPlate(1, 1);

    OutSlope = (OutPlate(2, 2) - OutPlate(1, 2)) / (OutPlate(2, 1) - OutPlate(1, 1));
    OutIntersect = OutPlate(1, 2) - OutSlope * OutPlate(1, 1);

    % ==================================================================================================================
    % DETERMINE THE X-POINT VIA QUASI MONTE-CARLO SIMULATION

    accuracy = 0.0001;
    d = 0.001;
    while true
        oldB = sqrt( ...
            (Currents(1) * ((Xpoint(2) - Locations(1, 2)) / ((Xpoint(1) - Locations(1, 1))^2 + (Xpoint(2) - Locations(1, 2))^2)) ...
                + Currents(2) * ((Xpoint(2) - Locations(2, 2)) / ((Xpoint(1) - Locations(2, 1))^2 + (Xpoint(2) - Locations(2, 2))^2)) ...
                + Currents(3) * ((Xpoint(2) - Locations(3, 2)) / ((Xpoint(1) - Locations(3, 1))^2 + (Xpoint(2) - Locations(3, 2))^2)) ...
                + Currents(4) * ((Xpoint(2) - Locations(4, 2)) / ((Xpoint(1) - Locations(4, 1))^2 + (Xpoint(2) - Locations(4, 2))^2)))^2 ...
                + (Currents(1) * ((Xpoint(1) - Locations(1, 1)) / ((Xpoint(1) - Locations(1, 1))^2 + (Xpoint(2) - Locations(1, 2))^2)) ...
                + Currents(2) * ((Xpoint(1) - Locations(2, 1)) / ((Xpoint(1) - Locations(2, 1))^2 + (Xpoint(2) - Locations(2, 2))^2)) ...
                + Currents(3) * ((Xpoint(1) - Locations(3, 1)) / ((Xpoint(1) - Locations(3, 1))^2 + (Xpoint(2) - Locations(3, 2))^2)) ...
                + Currents(4) * ((Xpoint(1) - Locations(4, 1)) / ((Xpoint(1) - Locations(4, 1))^2 + (Xpoint(2) - Locations(4, 2))^2)))^2 ...
            );

        dXpoint = d * (2 * rand(1, 2) - 1);
        newXpoint = Xpoint + dXpoint;

        newB = sqrt( ...
            (Currents(1) * ((newXpoint(2) - Locations(1, 2)) / ((newXpoint(1) - Locations(1, 1))^2 + (newXpoint(2) - Locations(1, 2))^2)) ...
                + Currents(2) * ((newXpoint(2) - Locations(2, 2)) / ((newXpoint(1) - Locations(2, 1))^2 + (newXpoint(2) - Locations(2, 2))^2)) ...
                + Currents(3) * ((newXpoint(2) - Locations(3, 2)) / ((newXpoint(1) - Locations(3, 1))^2 + (newXpoint(2) - Locations(3, 2))^2)) ...
                + Currents(4) * ((newXpoint(2) - Locations(4, 2)) / ((newXpoint(1) - Locations(4, 1))^2 + (newXpoint(2) - Locations(4, 2))^2)))^2 ...
                + (Currents(1) * ((newXpoint(1) - Locations(1, 1)) / ((newXpoint(1) - Locations(1, 1))^2 + (newXpoint(2) - Locations(1, 2))^2)) ...
                + Currents(2) * ((newXpoint(1) - Locations(2, 1)) / ((newXpoint(1) - Locations(2, 1))^2 + (newXpoint(2) - Locations(2, 2))^2)) ...
                + Currents(3) * ((newXpoint(1) - Locations(3, 1)) / ((newXpoint(1) - Locations(3, 1))^2 + (newXpoint(2) - Locations(3, 2))^2)) ...
                + Currents(4) * ((newXpoint(1) - Locations(4, 1)) / ((newXpoint(1) - Locations(4, 1))^2 + (newXpoint(2) - Locations(4, 2))^2)))^2 ...
            );

        if newB < oldB
            Xpoint = newXpoint;
        end

        if newB >= 0 && newB < accuracy
            break;
        end
    end

    % Determine the constant for the separatrix
    C = (Currents(1) * log((Xpoint(1) - Locations(1, 1))^2 + (Xpoint(2) - Locations(1, 2))^2) ...
        + Currents(2) * log((Xpoint(1) - Locations(2, 1))^2 + (Xpoint(2) - Locations(2, 2))^2) ...
        + Currents(3) * log((Xpoint(1) - Locations(3, 1))^2 + (Xpoint(2) - Locations(3, 2))^2) ...
        + Currents(4) * log((Xpoint(1) - Locations(4, 1))^2 + (Xpoint(2) - Locations(4, 2))^2));

    % ==================================================================================================================
    % DETERMINE THE STRIKE POINTS VIA QUASI MONTE-CARLO SIMULATION

    % Determine the inner strike point
    while true
        oldStrike = abs( ...
            Currents(1) * log((InStrikePoint(1) - Locations(1, 1))^2 + (InSlope * InStrikePoint(1) + InIntersect - Locations(1, 2))^2) ...
            + Currents(2) * log((InStrikePoint(1) - Locations(2, 1))^2 + (InSlope * InStrikePoint(1) + InIntersect - Locations(2, 2))^2) ...
            + Currents(3) * log((InStrikePoint(1) - Locations(3, 1))^2 + (InSlope * InStrikePoint(1) + InIntersect - Locations(3, 2))^2) ...
            + Currents(4) * log((InStrikePoint(1) - Locations(4, 1))^2 + (InSlope * InStrikePoint(1) + InIntersect - Locations(4, 2))^2) ...
            - C);

        % Determine the new strikepoint
        dInStrikePoint = d * (2 * rand(1, 2) - 1);
        newInStrikePoint = InStrikePoint + dInStrikePoint;

        % Make sure this strike point is within the desired range
        if newInStrikePoint(1) < InPlate(1, 1)
            newInStrikePoint(1) = InPlate(1, 1) - (newInStrikePoint(1) - InPlate(1, 1));
        end
        if newInStrikePoint(1) > InPlate(2, 1)
            newInStrikePoint(1) = InPlate(2, 1) + (InPlate(2, 1) - newInStrikePoint(1));
        end

        newStrike = abs( ...
            Currents(1) * log((newInStrikePoint(1) - Locations(1, 1))^2 + (InSlope * newInStrikePoint(1) + InIntersect - Locations(1, 2))^2) ...
            + Currents(2) * log((newInStrikePoint(1) - Locations(2, 1))^2 + (InSlope * newInStrikePoint(1) + InIntersect - Locations(2, 2))^2) ...
            + Currents(3) * log((newInStrikePoint(1) - Locations(3, 1))^2 + (InSlope * newInStrikePoint(1) + InIntersect - Locations(3, 2))^2) ...
            + Currents(4) * log((newInStrikePoint(1) - Locations(4, 1))^2 + (InSlope * newInStrikePoint(1) + InIntersect - Locations(4, 2))^2) ...
            - C ...
        );

        if newStrike < oldStrike
            InStrikePoint = newInStrikePoint;
            inStrikeIsView = false;
        end

        if newStrike >= 0 && newStrike < accuracy
            InStrikePoint(2) = InSlope * InStrikePoint(1) + InIntersect;
            if inStrikeIsView
                % Python's InStrikePoint is a view of Locations[2, :], so this write also changes Locations
                Locations(3, 2) = InStrikePoint(2);
            end
            break;
        end
    end

    % Determine the outer strike point
    while true
        oldStrike = abs( ...
            Currents(1) * log((OutStrikePoint(1) - Locations(1, 1))^2 + (OutSlope * OutStrikePoint(1) + OutIntersect - Locations(1, 2))^2) ...
            + Currents(2) * log((OutStrikePoint(1) - Locations(2, 1))^2 + (OutSlope * OutStrikePoint(1) + OutIntersect - Locations(2, 2))^2) ...
            + Currents(3) * log((OutStrikePoint(1) - Locations(3, 1))^2 + (OutSlope * OutStrikePoint(1) + OutIntersect - Locations(3, 2))^2) ...
            + Currents(4) * log((OutStrikePoint(1) - Locations(4, 1))^2 + (OutSlope * OutStrikePoint(1) + OutIntersect - Locations(4, 2))^2) ...
            - C);

        dOutStrikePoint = d * (2 * rand(1, 2) - 1);
        newOutStrikePoint = OutStrikePoint + dOutStrikePoint;

        % Make sure this strike point is within the desired range
        if newOutStrikePoint(1) < OutPlate(1, 1)
            newOutStrikePoint(1) = OutPlate(1, 1) - (newOutStrikePoint(1) - OutPlate(1, 1));
        end
        if newOutStrikePoint(1) > OutPlate(2, 1)
            newOutStrikePoint(1) = OutPlate(2, 1) + (OutPlate(2, 1) - newOutStrikePoint(1));
        end

        newStrike = abs( ...
            Currents(1) * log((newOutStrikePoint(1) - Locations(1, 1))^2 + (OutSlope * newOutStrikePoint(1) + OutIntersect - Locations(1, 2))^2) ...
            + Currents(2) * log((newOutStrikePoint(1) - Locations(2, 1))^2 + (OutSlope * newOutStrikePoint(1) + OutIntersect - Locations(2, 2))^2) ...
            + Currents(3) * log((newOutStrikePoint(1) - Locations(3, 1))^2 + (OutSlope * newOutStrikePoint(1) + OutIntersect - Locations(3, 2))^2) ...
            + Currents(4) * log((newOutStrikePoint(1) - Locations(4, 1))^2 + (OutSlope * newOutStrikePoint(1) + OutIntersect - Locations(4, 2))^2) ...
            - C ...
        );

        if newStrike < oldStrike
            OutStrikePoint = newOutStrikePoint;
            outStrikeIsView = false;
        end

        if newStrike >= 0 && newStrike < accuracy
            OutStrikePoint(2) = OutSlope * OutStrikePoint(1) + OutIntersect;
            if outStrikeIsView
                % Python's OutStrikePoint is a view of Locations[1, :], so this write also changes Locations
                Locations(2, 2) = OutStrikePoint(2);
            end
            break;
        end
    end

    % ==================================================================================================================

    % Determine the paramaterised positions of the strike points along the divertor plates, as measured from the y-axis
    sIn = sqrt( ...
        ((InStrikePoint(1) - InPlate(2, 1))^2 + (InStrikePoint(2) - InPlate(2, 2))^2) ...
        / ((InPlate(1, 1) - InPlate(2, 1))^2 + (InPlate(1, 2) - InPlate(2, 2))^2) ...
        );

    sOut = sqrt( ...
        ((OutStrikePoint(1) - OutPlate(1, 1))^2 + (OutStrikePoint(2) - OutPlate(1, 2))^2) ...
        / ((OutPlate(1, 1) - OutPlate(2, 1))^2 + (OutPlate(1, 2) - OutPlate(2, 2))^2) ...
        );

    % ==================================================================================================================
    % CREATE A PLOT
    % This part of the function can be deleted when integrating with
    % simulink. Plot creation is only used for diagnostic purposes

    % Define a grid of points
    delta = 0.0625;
    xmin = -2.5;
    xmax = 2.5;
    ymin = -2.5;
    ymax = 7.5;
    x = xmin + delta * (0:ceil((xmax - xmin) / delta) - 1);   % arange(xmin, xmax, delta)
    y = ymin + delta * (0:ceil((ymax - ymin) / delta) - 1);   % arange(ymin, ymax, delta)
    [X, Y] = meshgrid(x, y);

    % Define the field lines
    Z = (I0 * log((X - Locations(1, 1)).^2 + (Y - Locations(1, 2)).^2) ...
         + I1 * log((X - Locations(2, 1)).^2 + (Y - Locations(2, 2)).^2) ...
         + I2 * log((X - Locations(3, 1)).^2 + (Y - Locations(3, 2)).^2) ...
         + I3 * log((X - Locations(4, 1)).^2 + (Y - Locations(4, 2)).^2));

    % Create a plot of all the field lines
    figure;
    ax = gca;
    hold(ax, 'on');
    density = 0.25;
    nLevels = ceil(((C + 50 * density) - (C - 30 * density)) / density);
    levels = (C - 30 * density) + density * (0:nLevels - 1);   % arange(C - 30*density, C + 50*density, density)
    plotContourLevels(ax, x, y, Z, levels);

    % Highlight the separatrix by plotting it a second time
    levels = C;
    plotContourLevels(ax, x, y, Z, levels);

    % Draw the X-point and the positions of the currents
    plot(ax, Xpoint(1), Xpoint(2), 'o', 'Color', [1 0 0], 'MarkerFaceColor', [1 0 0], 'MarkerSize', 2);
    plot(ax, Locations(1, 1), Locations(1, 2), 'o', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0]);
    plot(ax, Locations(2, 1), Locations(2, 2), 'o', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0]);
    plot(ax, Locations(3, 1), Locations(3, 2), 'o', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0]);
    plot(ax, Locations(4, 1), Locations(4, 2), 'o', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0]);

    % Draw the divertor plates
    purple = [128, 0, 128] / 255;   % matplotlib 'purple'
    plot(ax, [InPlate(1, 1), InPlate(2, 1)], [InPlate(1, 2), InPlate(2, 2)], 'Color', purple);
    plot(ax, [OutPlate(1, 1), OutPlate(2, 1)], [OutPlate(1, 2), OutPlate(2, 2)], 'Color', purple);

    % Draw the stike points
    plot(ax, InStrikePoint(1), InStrikePoint(2), 'o', 'Color', [0 0 1], 'MarkerFaceColor', [0 0 1]);
    plot(ax, OutStrikePoint(1), OutStrikePoint(2), 'o', 'Color', [0 0 1], 'MarkerFaceColor', [0 0 1]);

    % Configure the axes
    axis(ax, 'equal');
    xlim(ax, [xmin, xmax]);
    ylim(ax, [ymin, ymax]);

end

function plotContourLevels(ax, x, y, Z, levels)
    % This function can be removed when implementing with simulink, as it
    % is only used for plotting

    % Emulates matplotlib ax.contour(Z, levels, origin='lower', colors='green', linewidths=1, extent=...)
    % followed by setp(CS.collections[-1], linewidth=2): the highest level in 'levels' is drawn with linewidth 2
    green = [0, 128, 0] / 255;   % matplotlib 'green'
    M = contourc(x, y, Z, levels);
    k = 1;
    while k <= size(M, 2)
        L = M(1, k);
        n = M(2, k);
        pts = M(:, k + 1:k + n);
        if L == levels(end)
            lw = 2;
        else
            lw = 1;
        end
        plot(ax, pts(1, :), pts(2, :), 'Color', green, 'LineWidth', lw);
        k = k + n + 1;
    end
end
