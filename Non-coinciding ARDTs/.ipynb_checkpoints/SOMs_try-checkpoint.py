import numpy as np
import xarray as xr
from sklearn.preprocessing import MinMaxScaler
from sklearn.cluster import MiniBatchKMeans
from minisom import MiniSom
import matplotlib.pyplot as plt

# Load the data
potential_temperature = xr.open_mfdataset('west_europe_pot/cascade_bard_v1_pot/*.nc',parallel=True).sel(time='2000')
specific_humidity = xr.open_mfdataset('west_europe_q/cascade_bard_v1_q/*.nc',parallel=True).sel(time='2000')

# Combine the two variables into one feature vector
data = np.stack([potential_temperature['PT'].values, specific_humidity['Q'].values], axis=-1)

# Reshape the data to be 2D (time x longitude x latitude, features)
data = data.reshape(-1, data.shape[-1])

# Normalize the data (optional but recommended for SOM)
data = (data - data.mean(axis=0)) / data.std(axis=0)

# Define the SOM parameters (you can adjust these)
map_size = (10, 10)  # Size of the SOM grid
sigma = 1.0  # Initial neighborhood radius
learning_rate = 0.5  # Initial learning rate
num_epochs = 1000  # Number of training epochs

# Create the SOM
som = MiniSom(map_size[0], map_size[1], data.shape[1], sigma=sigma, learning_rate=learning_rate)

# Initialize the SOM weights
som.random_weights_init(data)

# Train the SOM
som.train(data, num_epochs)

# Get the mapped indices for each data point
mapped_data = np.array([som.winner(d) for d in data])

# Reshape the mapped data to match the original dataset shape
mapped_data = mapped_data.reshape(potential_temperature.shape[0], potential_temperature.shape[1], potential_temperature.shape[2], -1)

# Create a new xarray dataset to store the SOM results
som_dataset = xr.Dataset(
    {
        'mapped_indices': (('time', 'longitude', 'latitude', 'map'), mapped_data),
        'potential_temperature': (('time', 'longitude', 'latitude'), potential_temperature.values),
        'specific_humidity': (('time', 'longitude', 'latitude'), specific_humidity.values),
    },
    coords={
        'time': dataset['time'],
        'longitude': dataset['longitude'],
        'latitude': dataset['latitude'],
        'map': np.arange(map_size[0] * map_size[1]),
    },
)

# Save or visualize the SOM results as needed
print(som_dataset)