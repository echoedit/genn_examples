import csv
import matplotlib.pyplot as plt
import numpy as np

with open("stim_spikes.csv", "rb") as stim_csv_file, open("neuron_spikes.csv", "rb") as neuron_csv_file:
    stim_csv_reader = csv.reader(stim_csv_file, delimiter = ",")
    neuron_csv_reader = csv.reader(neuron_csv_file, delimiter = ",")

    # Skip headers
    stim_csv_reader.next()
    neuron_csv_reader.next()

    # Read data and zip into columns
    stim_data_columns = zip(*stim_csv_reader)
    neuron_data_columns = zip(*neuron_csv_reader)

    # Convert CSV columns to numpy
    stim_times = np.asarray(stim_data_columns[0], dtype=float)
    stim_neuron_id = np.asarray(stim_data_columns[1], dtype=int)
    neuron_times = np.asarray(neuron_data_columns[0], dtype=float)
    neuron_neuron_id = np.asarray(neuron_data_columns[1], dtype=int)

    # Create plot
    figure, axes = plt.subplots(2, sharex=True)

    # Plot spikes
    axes[0].scatter(stim_times, stim_neuron_id, s=2, edgecolors="none", label="Stim")
    axes[0].scatter(neuron_times, neuron_neuron_id, s=2, edgecolors="none", label="Neurons")
    axes[0].legend()

    # Plot rates
    bins = np.arange(0, 5000 + 1, 100)
    stim_rate = np.histogram(stim_times, bins=bins)[0] * (1000.0 / 100.0) * (1.0 / 30000.0)
    neuron_rate = np.histogram(neuron_times, bins=bins)[0] * (1000.0 / 100.0) * (1.0 / 30000.0)
    axes[1].plot(bins[0:-1], stim_rate, label="Stim")
    axes[1].plot(bins[0:-1], neuron_rate, label="Neurons")
    axes[1].legend()

    axes[0].set_title("Spikes")
    axes[1].set_title("Firing rates")

    axes[0].set_xlim((0, 5000))
    axes[0].set_ylim((0, 30000))

    axes[0].set_ylabel("Neuron number")
    axes[1].set_ylabel("Mean firing rate [Hz]")

    axes[1].set_xlabel("Time [ms]")

    # Show plot
    plt.show()

