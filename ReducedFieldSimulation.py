# Import modules
import numpy as np
import matplotlib.pyplot as plt

# Define a random seed
seed = 30122005
np.random.seed(seed)

# Define the static currents
a0 = [0, 0]
I0 = 2

a1 = [1, 1]
I1 = -0.5

a2 = [-1, 1]
I2 = -0.5

a3 = [0, 4]
I3 = 4

Currents = np.array([I0, I1, I2, I3])
Locations = np.array([a0, a1, a2, a3])

# Define the inner and outer divertor plates
InStart = [-1, 1]
InEnd = [0, 0]
InPlate = np.array([InStart, InEnd])

OutStart = [0, 0]
OutEnd = [1, 1]
OutPlate = np.array([OutStart, OutEnd])

def StrikePoints(Currents, Locations, InPlate, OutPlate):
    # Define a grid of points
    delta = 0.0625
    xmin, xmax = -2.5, 2.5
    ymin, ymax = -2.5, 7.5
    x = np.arange(xmin, xmax, delta)
    y = np.arange(ymin, ymax, delta)
    X, Y = np.meshgrid(x, y)

    # Initialize the X-point and strike point first guesses
    Xpoint = np.array([0, 1])
    InStrikePoint = Locations[2, :]
    OutStrikePoint = Locations[1, :]

    # Define the slope and intersect of the divertor plates
    InSlope = (InPlate[1, 1] - InPlate[0, 1])/(InPlate[1, 0] - InPlate[0, 0])
    InIntersect = InPlate[0, 1] - InSlope * InPlate[0, 0]

    OutSlope = (OutPlate[1, 1] - OutPlate[0, 1]) / (OutPlate[1, 0] - OutPlate[0, 0])
    OutIntersect = OutPlate[0, 1] - OutSlope * OutPlate[0, 0]

    # ==================================================================================================================
    # DETERMINE THE X-POINT VIA QUASI MONTE-CARLO SIMULATION

    accuracy = 0.0001
    d = 0.001
    while True:
        oldB = np.sqrt(
            (Currents[0] * ((Xpoint[1] - Locations[0, 1])/((Xpoint[0] - Locations[0, 0])**2 + (Xpoint[1] - Locations[0, 1])**2))
                + Currents[1] * ((Xpoint[1] - Locations[1, 1])/((Xpoint[0] - Locations[1, 0])**2 + (Xpoint[1] - Locations[1, 1])**2))
                + Currents[2] * ((Xpoint[1] - Locations[2, 1])/((Xpoint[0] - Locations[2, 0])**2 + (Xpoint[1] - Locations[2, 1])**2))
                + Currents[3] * ((Xpoint[1] - Locations[3, 1])/((Xpoint[0] - Locations[3, 0])**2 + (Xpoint[1] - Locations[3, 1])**2)))**2
                + (Currents[0] * ((Xpoint[0] - Locations[0, 0])/((Xpoint[0] - Locations[0, 0])**2 + (Xpoint[1] - Locations[0, 1])**2))
                + Currents[1] * ((Xpoint[0] - Locations[1, 0])/((Xpoint[0] - Locations[1, 0])**2 + (Xpoint[1] - Locations[1, 1])**2))
                + Currents[2] * ((Xpoint[0] - Locations[2, 0])/((Xpoint[0] - Locations[2, 0])**2 + (Xpoint[1] - Locations[2, 1])**2))
                + Currents[3] * ((Xpoint[0] - Locations[3, 0])/((Xpoint[0] - Locations[3, 0])**2 + (Xpoint[1] - Locations[3, 1])**2)))**2
        )

        dXpoint = d * (2 * np.random.rand(2) - 1)
        newXpoint = Xpoint + dXpoint

        newB = np.sqrt(
            (Currents[0] * ((newXpoint[1] - Locations[0, 1])/((newXpoint[0] - Locations[0, 0])**2 + (newXpoint[1] - Locations[0, 1])**2))
                + Currents[1] * ((newXpoint[1] - Locations[1, 1])/((newXpoint[0] - Locations[1, 0])**2 + (newXpoint[1] - Locations[1, 1])**2))
                + Currents[2] * ((newXpoint[1] - Locations[2, 1])/((newXpoint[0] - Locations[2, 0])**2 + (newXpoint[1] - Locations[2, 1])**2))
                + Currents[3] * ((newXpoint[1] - Locations[3, 1])/((newXpoint[0] - Locations[3, 0])**2 + (newXpoint[1] - Locations[3, 1])**2)))**2
                + (Currents[0] * ((newXpoint[0] - Locations[0, 0])/((newXpoint[0] - Locations[0, 0])**2 + (newXpoint[1] - Locations[0, 1])**2))
                + Currents[1] * ((newXpoint[0] - Locations[1, 0])/((newXpoint[0] - Locations[1, 0])**2 + (newXpoint[1] - Locations[1, 1])**2))
                + Currents[2] * ((newXpoint[0] - Locations[2, 0])/((newXpoint[0] - Locations[2, 0])**2 + (newXpoint[1] - Locations[2, 1])**2))
                + Currents[3] * ((newXpoint[0] - Locations[3, 0])/((newXpoint[0] - Locations[3, 0])**2 + (newXpoint[1] - Locations[3, 1])**2)))**2
        )

        if newB < oldB:
            Xpoint = newXpoint

        if 0 <= newB < accuracy:
            break

    # Determine the constant for the separatrix
    C = (Currents[0] * np.log((Xpoint[0] - Locations[0, 0]) ** 2 + (Xpoint[1] - Locations[0, 1]) ** 2)
            + Currents[1] * np.log((Xpoint[0] - Locations[1, 0]) ** 2 + (Xpoint[1] - Locations[1, 1]) ** 2)
            + Currents[2] * np.log((Xpoint[0] - Locations[2, 0]) ** 2 + (Xpoint[1] - Locations[2, 1]) ** 2)
            + Currents[3] * np.log((Xpoint[0] - Locations[3, 0]) ** 2 + (Xpoint[1] - Locations[3, 1]) ** 2))

    # ==================================================================================================================
    # DETERMINE THE STRIKE POINTS VIA QUASI MONTE-CARLO SIMULATION

    # Determine the inner strike point
    while True:
        oldStrike = np.abs(
            Currents[0] * np.log((InStrikePoint[0] - Locations[0, 0])**2 + (InSlope * InStrikePoint[0] + InIntersect - Locations[0, 1])**2)
            + Currents[1] * np.log((InStrikePoint[0] - Locations[1, 0])**2 + (InSlope * InStrikePoint[0] + InIntersect - Locations[1, 1])**2)
            + Currents[2] * np.log((InStrikePoint[0] - Locations[2, 0])**2 + (InSlope * InStrikePoint[0] + InIntersect - Locations[2, 1])**2)
            + Currents[3] * np.log((InStrikePoint[0] - Locations[3, 0])**2 + (InSlope * InStrikePoint[0] + InIntersect - Locations[3, 1])**2)
            - C)

        # Determine the new strikepoint
        dInStrikePoint = d * (2 * np.random.rand(2) - 1)
        newInStrikePoint = InStrikePoint + dInStrikePoint

        # Make sure this strike point is within the desired range
        if newInStrikePoint[0] < InPlate[0, 0]:
            newInStrikePoint[0] = InPlate[0, 0] - (newInStrikePoint[0] - InPlate[0, 0])
        if newInStrikePoint[0] > InPlate[1, 0]:
            newInStrikePoint[0] = InPlate[1, 0] + (InPlate[1, 0] - newInStrikePoint[0])

        newStrike = np.abs(
            Currents[0] * np.log((newInStrikePoint[0] - Locations[0, 0]) ** 2 + (InSlope * newInStrikePoint[0] + InIntersect -Locations[0, 1]) ** 2)
            + Currents[1] * np.log((newInStrikePoint[0] - Locations[1, 0]) ** 2 + (InSlope * newInStrikePoint[0] + InIntersect - Locations[1, 1]) ** 2)
            + Currents[2] * np.log((newInStrikePoint[0] - Locations[2, 0]) ** 2 + (InSlope * newInStrikePoint[0] + InIntersect - Locations[2, 1]) ** 2)
            + Currents[3] * np.log((newInStrikePoint[0] - Locations[3, 0]) ** 2 + (InSlope * newInStrikePoint[0] + InIntersect - Locations[3, 1]) ** 2)
            - C
        )

        if newStrike < oldStrike:
            InStrikePoint = newInStrikePoint

        if 0 <= newStrike < accuracy:
            InStrikePoint[1] = InSlope * InStrikePoint[0] + InIntersect
            break

    # Determine the outer strike point
    while True:
        oldStrike = np.abs(
            Currents[0] * np.log((OutStrikePoint[0] - Locations[0, 0]) ** 2 + (OutSlope * OutStrikePoint[0] + OutIntersect - Locations[0, 1]) ** 2)
            + Currents[1] * np.log((OutStrikePoint[0] - Locations[1, 0]) ** 2 + (OutSlope * OutStrikePoint[0] + OutIntersect - Locations[1, 1]) ** 2)
            + Currents[2] * np.log((OutStrikePoint[0] - Locations[2, 0]) ** 2 + (OutSlope * OutStrikePoint[0] + OutIntersect - Locations[2, 1]) ** 2)
            + Currents[3] * np.log((OutStrikePoint[0] - Locations[3, 0]) ** 2 + (OutSlope * OutStrikePoint[0] + OutIntersect - Locations[3, 1]) ** 2)
            - C)

        dOutStrikePoint = d * (2 * np.random.rand(2) - 1)
        newOutStrikePoint = OutStrikePoint + dOutStrikePoint

        # Make sure this strike point is within the desired range
        if newOutStrikePoint[0] < OutPlate[0, 0]:
            newOutStrikePoint[0] = OutPlate[0, 0] - (newOutStrikePoint[0] - OutPlate[0, 0])
        if newOutStrikePoint[0] > OutPlate[1, 0]:
            newOutStrikePoint[0] = OutPlate[1, 0] + (OutPlate[1, 0] - newOutStrikePoint[0])

        newStrike = np.abs(
            Currents[0] * np.log(
                (newOutStrikePoint[0] - Locations[0, 0]) ** 2 + (OutSlope * newOutStrikePoint[0] + OutIntersect - Locations[0, 1]) ** 2)
            + Currents[1] * np.log(
                (newOutStrikePoint[0] - Locations[1, 0]) ** 2 + (OutSlope * newOutStrikePoint[0] + OutIntersect - Locations[1, 1]) ** 2)
            + Currents[2] * np.log(
                (newOutStrikePoint[0] - Locations[2, 0]) ** 2 + (OutSlope * newOutStrikePoint[0] + OutIntersect - Locations[2, 1]) ** 2)
            + Currents[3] * np.log(
                (newOutStrikePoint[0] - Locations[3, 0]) ** 2 + (OutSlope * newOutStrikePoint[0] + OutIntersect - Locations[3, 1]) ** 2)
            - C
        )

        if newStrike < oldStrike:
            OutStrikePoint = newOutStrikePoint

        if 0 <= newStrike < accuracy:
            OutStrikePoint[1] = OutSlope * OutStrikePoint[0] + OutIntersect
            break

    # ==================================================================================================================

    # Determine the paramaterised positions of the strike points along the divertor plates, as measured from the y-axis
    sIn = np.sqrt(
        ((InStrikePoint[0] - InPlate[1, 0])**2 + (InStrikePoint[1] - InPlate[1, 1])**2)
        / ((InPlate[0, 0] - InPlate[1, 0])**2 + (InPlate[0, 1] - InPlate[1, 1])**2)
    )

    sOut = np.sqrt(
        ((OutStrikePoint[0] - OutPlate[0, 0])**2 + (OutStrikePoint[1] - OutPlate[0, 1])**2)
        / ((OutPlate[0, 0] - OutPlate[1, 0])**2 + (OutPlate[0, 1] - OutPlate[1, 1])**2)
    )

    # ==================================================================================================================
    # Create a plot

    # Define the field lines
    Z = (I0 * np.log((X - Locations[0, 0]) ** 2 + (Y - Locations[0, 1]) ** 2)
         + I1 * np.log((X - Locations[1, 0]) ** 2 + (Y - Locations[1, 1]) ** 2)
         + I2 * np.log((X - Locations[2, 0]) ** 2 + (Y - Locations[2, 1]) ** 2)
         + I3 * np.log((X - Locations[3, 0]) ** 2 + (Y - Locations[3, 1]) ** 2))

    # Create a plot of all the field lines
    fig, ax = plt.subplots()
    density = 1
    levels = np.arange(C - 30 * density, C + 50 * density, density)
    CS = ax.contour(Z, levels,
                    origin='lower',
                    colors='green',
                    linewidths=1,
                    extent=(xmin, xmax, ymin, ymax))

    zc = CS.collections[-1]
    plt.setp(zc, linewidth=2)

    # Highlight the separatrix by plotting it a second time
    levels = np.array([C])
    CS = ax.contour(Z, levels,
                    origin='lower',
                    colors='green',
                    linewidths=1,
                    extent=(xmin, xmax, ymin, ymax))

    zc = CS.collections[-1]
    plt.setp(zc, linewidth=2)

    # Draw the X-point and the positions of the currents
    ax.plot(Xpoint[0], Xpoint[1], 'o', color='red', markersize=2)
    ax.plot(Locations[0, 0], Locations[0, 1], 'o', color='black')
    ax.plot(Locations[1, 0], Locations[1, 1], 'o', color='black')
    ax.plot(Locations[2, 0], Locations[2, 1], 'o', color='black')
    ax.plot(Locations[3, 0], Locations[3, 1], 'o', color='black')

    # Draw the divertor plates
    plt.plot([InPlate[0, 0], InPlate[1, 0]], [InPlate[0, 1], InPlate[1, 1]], color='purple')
    plt.plot([OutPlate[0, 0], OutPlate[1, 0]], [OutPlate[0, 1], OutPlate[1, 1]], color='purple')

    # Draw the stike points
    ax.plot(InStrikePoint[0], InStrikePoint[1], 'o', color='blue')
    ax.plot(OutStrikePoint[0], OutStrikePoint[1], 'o', color='blue')

    # Configure the axes
    ax.set_aspect('equal')
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.title(f"The strike points")

    plt.show()

    return (sIn, sOut)

print(StrikePoints(Currents, Locations, InPlate, OutPlate))