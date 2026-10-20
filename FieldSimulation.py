"""
STRIKE POINT FINDER:

This code finds the strike points in a simplified tokamak geometry.

The strike points are defined as the intersection of the separatrix with the divertor plates.
The separatrix is defined as that magnetic field line which goes through the X-point.
The X-point is defined as that point where the magnetic field is zero.
The divertor plates are modelled as lines: y = plate[0]*x + plate[1].
The coils and plasma current are modelled as infinitely thin wires.

The magnetic field is determined analytically from the set current configuration.
The equation for the magnetic field lines is determined analytically as well.

The position of the X-point is determined via a quasi Monte-Carlo simulation (random walk).
The position of the strike points is determined via the same method.
NOTE: strike point determination may be better implemented using Newton's method. This will be implemented later.

Current implementation is focussed towards strike-point calculation and visualisation.
For an integrated implementation, the visualisation for several frames may be removed.
This will require some adaptation of the code, but should in principle not be difficult.

This code is my own, original work.
- Raf van der Veeken
"""

# Import modules
import numpy as np
import matplotlib.pyplot as plt

# Define a random seed
seed = 30122005
np.random.seed(seed)

# Define a grid of points
delta = 0.0625
xmin, xmax = -2.5, 2.5
ymin, ymax = -2.5, 7.5
x = np.arange(xmin, xmax, delta)
y = np.arange(ymin, ymax, delta)
X, Y = np.meshgrid(x, y)

# Define the static currents
a1 = [1, 1]
I1 = -0.5

a2 = [-1, 1]
I2 = -0.5

a3 = [0, 4]
I3 = 4

# Define the slope and intersection with the y-axis for the divertor plates
InStart = [-1, 1]
InEnd = [0, 0]
InPlate = np.array([InStart, InEnd])

OutStart = [0, 0]
OutEnd = [1, 1]
OutPlate = np.array([OutStart, OutEnd])

InSlope = (InPlate[1, 1] - InPlate[0, 1])/(InPlate[1, 0] - InPlate[0, 0])
InIntersect = InPlate[0, 1] - InSlope * InPlate[0, 0]

OutSlope = (OutPlate[1, 1] - OutPlate[0, 1]) / (OutPlate[1, 0] - OutPlate[0, 0])
OutIntersect = OutPlate[0, 1] - OutSlope * OutPlate[0, 0]


InPlate = np.array([InSlope, InIntersect])
OutPlate = np.array([OutSlope, OutIntersect])

# Define the number of frames to capture
frames = 1

# Initialize the X-point and strike point first guesses
Xpoint = np.array([0, 1])
InStrikePoint = a2
OutStrikePoint = a1

# Find the strike points for every frame of the movie
for i in range(frames):

    # Define the varying current
    a0 = [0, 0]
    I0Amp = 1.0
    I0Avg = 2.0
    freq = 3/frames

    I0 = I0Avg - I0Amp*np.cos(2*np.pi*freq*i)

    """I0min = 1
    I0max = 3
    if i < frames / 2:
        I0 = I0min + ((I0max - I0min) / (frames / 2)) * i
    else:
        I0 = 2 * I0max - I0min - ((I0max - I0min) / (frames / 2)) * i"""

    # ===========================================================================================
    # DETERMINE X-POINT VIA QUASI MONTE-CARLO SIMULATION
    accuracy = 0.0001
    d = 0.001
    while True:
        oldB = np.sqrt(
            (I0 * (Xpoint[1]/(Xpoint[0]**2 + Xpoint[1]**2))
             + I1 * ((Xpoint[1] - a1[1])/((Xpoint[0] - a1[0])**2 + (Xpoint[1] - a1[1])**2))
             + I2 * ((Xpoint[1] - a2[1])/((Xpoint[0] - a2[0])**2 + (Xpoint[1] - a2[1])**2))
             + I3 * ((Xpoint[1] - a3[1])/((Xpoint[0] - a3[0])**2 + (Xpoint[1] - a3[1])**2)))**2
            + (I0 * (Xpoint[0]/(Xpoint[0]**2 + Xpoint[1]**2))
             + I1 * ((Xpoint[0] - a1[0])/((Xpoint[0] - a1[0])**2 + (Xpoint[1] - a1[1])**2))
             + I2 * ((Xpoint[0] - a2[0])/((Xpoint[0] - a2[0])**2 + (Xpoint[1] - a2[1])**2))
             + I3 * ((Xpoint[0] - a3[0])/((Xpoint[0] - a3[0])**2 + (Xpoint[1] - a3[1])**2)))**2
        )

        dXpoint = d * (2 * np.random.rand(2) - 1)
        newXpoint = Xpoint + dXpoint

        newB = np.sqrt(
            (I0 * (newXpoint[1]/(newXpoint[0]**2 + newXpoint[1]**2))
             + I1 * ((newXpoint[1] - a1[1])/((newXpoint[0] - a1[0])**2 + (newXpoint[1] - a1[1])**2))
             + I2 * ((newXpoint[1] - a2[1])/((newXpoint[0] - a2[0])**2 + (newXpoint[1] - a2[1])**2))
             + I3 * ((newXpoint[1] - a3[1])/((newXpoint[0] - a3[0])**2 + (newXpoint[1] - a3[1])**2)))**2
            + (I0 * (newXpoint[0]/(newXpoint[0]**2 + newXpoint[1]**2))
             + I1 * ((newXpoint[0] - a1[0])/((newXpoint[0] - a1[0])**2 + (newXpoint[1] - a1[1])**2))
             + I2 * ((newXpoint[0] - a2[0])/((newXpoint[0] - a2[0])**2 + (newXpoint[1] - a2[1])**2))
             + I3 * ((newXpoint[0] - a3[0])/((newXpoint[0] - a3[0])**2 + (newXpoint[1] - a3[1])**2)))**2
        )

        if newB < oldB:
            Xpoint = newXpoint

        if 0 <= newB < accuracy:
            #print("X-point found")
            #print()
            break

    # ===========================================================================================

    # Determine the constant for the separatrix
    C = (I0 * np.log((Xpoint[0] - a0[0]) ** 2 + (Xpoint[1] - a0[1]) ** 2)
         + I1 * np.log((Xpoint[0] - a1[0]) ** 2 + (Xpoint[1] - a1[1]) ** 2)
         + I2 * np.log((Xpoint[0] - a2[0]) ** 2 + (Xpoint[1] - a2[1]) ** 2)
         + I3 * np.log((Xpoint[0] - a3[0]) ** 2 + (Xpoint[1] - a3[1]) ** 2))

    # ===========================================================================================
    # DETERMINE STRIKE POINTS VIA QUASI MONTE-CARLO SIMULATION

    # Determine the inner strike point
    while True:
        oldStrike = np.abs(
            I0 * np.log((InStrikePoint[0] - a0[0])**2 + (InPlate[0] * InStrikePoint[0] + InPlate[1] - a0[1])**2)
            + I1 * np.log((InStrikePoint[0] - a1[0])**2 + (InPlate[0] * InStrikePoint[0] + InPlate[1] - a1[1])**2)
            + I2 * np.log((InStrikePoint[0] - a2[0])**2 + (InPlate[0] * InStrikePoint[0] + InPlate[1] - a2[1])**2)
            + I3 * np.log((InStrikePoint[0] - a3[0])**2 + (InPlate[0] * InStrikePoint[0] + InPlate[1] - a3[1])**2)
            - C)

        # Determine the new strikepoint
        dInStrikePoint = d * (2 * np.random.rand(2) - 1)
        newInStrikePoint = InStrikePoint + dInStrikePoint

        # Make sure this strike point is within the desired range
        if newInStrikePoint[0] < a2[0]:
            newInStrikePoint[0] = a2[0] - (newInStrikePoint[0] - a2[0])
        if newInStrikePoint[0] > 0:
            newInStrikePoint[0] = -newInStrikePoint[0]

        newStrike = np.abs(
            I0 * np.log((newInStrikePoint[0] - a0[0]) ** 2 + (InPlate[0] * newInStrikePoint[0] + InPlate[1] - a0[1]) ** 2)
            + I1 * np.log((newInStrikePoint[0] - a1[0]) ** 2 + (InPlate[0] * newInStrikePoint[0] + InPlate[1] - a1[1]) ** 2)
            + I2 * np.log((newInStrikePoint[0] - a2[0]) ** 2 + (InPlate[0] * newInStrikePoint[0] + InPlate[1] - a2[1]) ** 2)
            + I3 * np.log((newInStrikePoint[0] - a3[0]) ** 2 + (InPlate[0] * newInStrikePoint[0] + InPlate[1] - a3[1]) ** 2)
            - C
        )

        if newStrike < oldStrike:
            #print(newStrike)
            #print()
            InStrikePoint = newInStrikePoint

        if 0 <= newStrike < accuracy:
            InStrikePoint[1] = InPlate[0] * InStrikePoint[0] + InPlate[1]
            #print(f"Inner strike point found: {InStrikePoint}")
            #print()
            break

    # Determine the outer strike point
    while True:
        oldStrike = np.abs(
            I0 * np.log((OutStrikePoint[0] - a0[0]) ** 2 + (OutPlate[0] * OutStrikePoint[0] + OutPlate[1] - a0[1]) ** 2)
            + I1 * np.log((OutStrikePoint[0] - a1[0]) ** 2 + (OutPlate[0] * OutStrikePoint[0] + OutPlate[1] - a1[1]) ** 2)
            + I2 * np.log((OutStrikePoint[0] - a2[0]) ** 2 + (OutPlate[0] * OutStrikePoint[0] + OutPlate[1] - a2[1]) ** 2)
            + I3 * np.log((OutStrikePoint[0] - a3[0]) ** 2 + (OutPlate[0] * OutStrikePoint[0] + OutPlate[1] - a3[1]) ** 2)
            - C)

        dOutStrikePoint = d * (2 * np.random.rand(2) - 1)
        newOutStrikePoint = OutStrikePoint + dOutStrikePoint

        # Make sure this strike point is within the desired range
        if newOutStrikePoint[0] > a1[0]:
            newOutStrikePoint[0] = a1[0] - (newOutStrikePoint[0] - a1[0])
        if newOutStrikePoint[0] < 0:
            newOutStrikePoint[0] = -newOutStrikePoint[0]

        newStrike = np.abs(
            I0 * np.log(
                (newOutStrikePoint[0] - a0[0]) ** 2 + (OutPlate[0] * newOutStrikePoint[0] + OutPlate[1] - a0[1]) ** 2)
            + I1 * np.log(
                (newOutStrikePoint[0] - a1[0]) ** 2 + (OutPlate[0] * newOutStrikePoint[0] + OutPlate[1] - a1[1]) ** 2)
            + I2 * np.log(
                (newOutStrikePoint[0] - a2[0]) ** 2 + (OutPlate[0] * newOutStrikePoint[0] + OutPlate[1] - a2[1]) ** 2)
            + I3 * np.log(
                (newOutStrikePoint[0] - a3[0]) ** 2 + (OutPlate[0] * newOutStrikePoint[0] + OutPlate[1] - a3[1]) ** 2)
            - C
        )

        if newStrike < oldStrike:
            #print(newStrike)
            #print()
            OutStrikePoint = newOutStrikePoint

        if 0 <= newStrike < accuracy:
            OutStrikePoint[1] = OutPlate[0] * OutStrikePoint[0] + OutPlate[1]
            #print(f"Inner strike point found: {InStrikePoint}")
            #print()
            break

    # ===========================================================================================

    # Define the field lines
    Z = (I0 * np.log((X - a0[0])**2 + (Y - a0[1])**2)
         + I1 * np.log((X - a1[0])**2 + (Y - a1[1])**2)
         + I2 * np.log((X - a2[0])**2 + (Y - a2[1])**2)
         + I3 * np.log((X - a3[0])**2 + (Y - a3[1])**2))

    # Create a plot of all the field lines
    fig, ax = plt.subplots()
    density = 0.25
    levels = np.arange(C - 30*density, C + 80*density, density)
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
    ax.plot(a0[0], a0[1], 'o', color='black')
    ax.plot(a1[0], a1[1], 'o', color='black')
    ax.plot(a2[0], a2[1], 'o', color='black')
    ax.plot(a3[0], a3[1], 'o', color='black')

    # Draw the divertor plates
    plt.plot([InStart[0], InEnd[0]], [InStart[1], InEnd[1]], color='purple')
    plt.plot([OutStart[0], OutEnd[0]], [OutStart[1], OutEnd[1]], color='purple')

    # Draw the stike points
    ax.plot(InStrikePoint[0], InStrikePoint[1], 'o', color='blue')
    ax.plot(OutStrikePoint[0], OutStrikePoint[1], 'o', color='blue')

    # Configure the axes
    ax.set_aspect('equal')
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    #plt.title(f"X-point: ({Xpoint[0]}, {Xpoint[1]}), accuracy < {accuracy}")

    # Save the plots to a folder in order to splice them together later using ffmpeg in cmd
    if i < 10:
        location = f"C:\\Users\\Rafva\\Bureaublad\\school\\Universiteit\\Msc NF\\Year 1\\36FBA\\FieldLineMovie2\\Frame_00{i}"
    elif i < 100:
        location = f"C:\\Users\\Rafva\\Bureaublad\\school\\Universiteit\\Msc NF\\Year 1\\36FBA\\FieldLineMovie2\\Frame_0{i}"
    else:
        location = f"C:\\Users\\Rafva\\Bureaublad\\school\\Universiteit\\Msc NF\\Year 1\\36FBA\\FieldLineMovie2\\Frame_{i}"
    fig.savefig(location)
    plt.show()
    plt.close()

    # Show frame generation progress
    print(f"Progress: {(i + 1) / frames * 100}%")

