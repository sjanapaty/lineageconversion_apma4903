import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import convolve2d
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# run parameters
dy = 0.01
# spatial extent
Size = 50  # 10x10 grid
states = 6  # Number of species
Tmax = 10  # Maximum time
# numerical parameters
run_method = "RK45"  # You can use 'BDF' if the system is stiff
run_method = "BDF"  # You can use 'BDF' if the system is stiff
run_rtol = 1e-6
run_atol = 1e-8

# derived run parameters:
# Grid dimensions
dim_x, dim_y = Size, Size
# Time span for the simulation
t_span = (0, Tmax)  # Short time span for demonstration; adjust as needed

# chemical Parameters
params = {
    "gZ": 0.1,
    "gA": 0.01,
    "gB": 0.01,
    "gAA": 0.001,  # Assuming gA1 in your parameters corresponds to gAA in equations
    "dZ": 0.00001,
    "dA": 0.00001,
    "dB": 0.00001,
    "dAA": 0.00001,  # Assuming dA1 in your parameters corresponds to dAA in equations
    "pZA": 0.2,
    "pZB": 0.05,
    "pAAA": 0.2,  # Assuming pAA1 in your parameters corresponds to pAAA in equations
    "uAA": 0.1,  # Assuming uA1 in your parameters corresponds to uAA in equations
    "uB": 0.1,
    "kAA": 0.1,  # Assuming kA1 in your parameters corresponds to kAA in equations
    "kB": 0.1,
    "qAA": 0.3,
    "qB": 0.3,
    "DZ": 0.3,  # Assuming separate diffusion coefficients for each species
    "DA": 0.3,
    "DB": 0.3,
    "DAA": 0.3,  # Assuming DtAA in your parameters corresponds to DAA in equations
    "DtB": 0.3,  # Assuming DtB in your parameters corresponds to DAA in equations
}


# Laplacian kernel
laplacian_kernel = np.array([[0, 1, 0], [1, -4, 1], [0, 1, 0]])


# Function to compute the Laplacian using convolution
def compute_laplacian(u):
    return convolve2d(u, laplacian_kernel, mode="same", boundary="wrap")


# Reaction-diffusion system
def reaction_diffusion_system(t, y, params):
    # Reshape the 1D array back into (states, dim_x, dim_y)
    concentrations = y.reshape((states, dim_x, dim_y))

    # Compute the Laplacian for each state
    laplacians = np.array([compute_laplacian(concentrations[i]) for i in range(states)])

    # Initialize the derivative array
    dydt = np.zeros_like(concentrations)

    Z, A, AA, B, tAA, tB = concentrations

    # Reaction terms, corrected
    dydt[0] = (
        params["gZ"] * 1 - params["dZ"] * 1 - (params["pZA"] + params["pZB"]) * Z
    ) + params["DZ"] * laplacians[0]
    dydt[1] = (
        params["gA"] * 1 - params["dA"] * 1 + params["pZA"] * Z - params["pAAA"] * A
    ) + params["DA"] * laplacians[1]
    dydt[2] = (
        params["gAA"] * 1
        - params["dAA"] * 1
        + params["pAAA"] * A
        + params["uAA"] * tAA / (params["kAA"] + tAA) * B
        - params["uB"] * tB / (params["kB"] + tB) * AA
    ) + params["DAA"] * laplacians[2]
    dydt[3] = (
        params["gB"] * 1
        - params["dB"] * 1
        + params["pZB"] * Z
        + params["uB"] * tB / (params["kB"] + tB) * AA
        - params["uAA"] * tAA / (params["kAA"] + tAA) * B
    ) + params["DB"] * laplacians[3]
    dydt[4] = -params["qAA"] * tAA + params["DAA"] * laplacians[4]
    dydt[5] = -params["qB"] * tB + params["DtB"] * laplacians[5]

    # Flatten the dydt array to 1D for the ODE solver
    return dydt.ravel()


# Initial conditions: small random concentrations to avoid instabilities
y0 = np.random.rand(states * dim_x * dim_y) * dy



# Run the simulation
solution = solve_ivp(
    fun=reaction_diffusion_system,
    t_span=t_span,
    y0=y0,
    args=(params,),
    method=run_method,
    rtol=run_rtol,
    atol=run_atol,
)

# Reshape the solution to (states, dim_x, dim_y, time_points)
concentration_history = solution.y.reshape((states, dim_x, dim_y, -1))
# Calculate mean concentrations over time for each species
mean_concentrations = np.mean(concentration_history, axis=(1, 2))
# Get the time points from the solution
time_points = solution.t

# Assuming we want to plot the concentration of the first species over time
species_index = 0  # Index of the species to plot


# The code should now execute the equations including both the reaction terms and the diffusion terms as per the Laplacian computed using `convolve2d`. Each `dydt[i]` line corresponds to one of the equations from your image, with the appropriate parameters and Laplacian terms added.


# Species labels as per the image
species_labels = ["Z", "A", "B", "AA", "t_AA", "t_B"]

# Plot mean concentrations over time with correct labels
plt.figure(figsize=(14, 7))
for i, label in enumerate(species_labels):
    plt.plot(time_points, mean_concentrations[i], label=f"Species {label}")
plt.xlabel("Time")
plt.ylabel("Mean Concentration")
plt.legend()
plt.title("Mean Concentrations of Species Over Time")
plt.show()

for t_idx, time_point in enumerate(time_points):
    plt.imshow(concentration_history[species_index, :, :, t_idx], cmap="viridis")
    plt.title(f"Species {species_index} Concentration at time {time_point:.2f}")
    plt.colorbar()
    plt.show()


# Function to create CMYK overlay
def create_cmyk_frame(concentration_history, t_idx):
    cmy_image = np.zeros((*concentration_history[0, :, :, t_idx].shape, 3))
    for i in range(3):
        cmy_image[..., i] = concentration_history[i, :, :, t_idx] / np.max(
            concentration_history[i, :, :, t_idx]
        )

    # Convert CMY to RGB (ignoring black for overlay)
    rgb_image = 1 - cmy_image
    return rgb_image.clip(0, 1)

# ... (rest of your code for the CMYK frame creation)

# Initialize the figure for animation
fig, ax = plt.subplots()
print("creating cmyk frame...")
cmyk_frame = create_cmyk_frame(concentration_history, 0)
im = ax.imshow(cmyk_frame)
ax.set_title(f"CMYK Overlay at time {time_points[0]:.2f}")



# Function to update the frame for the animation
def update(frame):
    cmyk_frame = create_cmyk_frame(concentration_history, frame)
    # im.set_array(cmyk_frame)
    im.set_data(cmyk_frame)
    ax.set_title(f"CMYK Overlay at time {time_points[frame]:.2f}")
    return [im]


# Create the animation
print("Animating...")
ani = FuncAnimation(fig, update, frames=len(time_points), blit=True)

# To display the animation in a Jupyter notebook
from IPython.display import HTML

HTML(ani.to_html5_video())

# To save the animation to a file
# filename is a string with the name of the file to save to
# but the filename includes the variable "method"
filename = "cmyk_animation_%s.mp4" % run_method

ani.save(filename, writer="ffmpeg")

# Make sure to have ffmpeg installed for video saving and displaying
