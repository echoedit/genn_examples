import csv
import matplotlib.pyplot as plt
import numpy as np

with open("spikes.csv", "rb") as spikes_csv_file:
    spikes_csv_reader = csv.reader(spikes_csv_file, delimiter = ",")

    # Skip headers
    spikes_csv_reader.next()

    # Read data and zip into columns
    spikes_data_columns = zip(*spikes_csv_reader)

    # Convert CSV columns to numpy
    spike_times = np.asarray(spikes_data_columns[0], dtype=float)
    spike_neuron_id = np.asarray(spikes_data_columns[1], dtype=int)

    # Create plot
    figure, axes = plt.subplots(2, sharex=True)

    # Plot spikes
    axes[0].scatter(spike_times, spike_neuron_id, s=2, edgecolors="none")

    # Plot rates
    bins = np.arange(0, 10000 + 1, 10)
    rate = np.histogram(spike_times, bins=bins)[0] *  (1000.0 / 10.0) * (1.0 / 3200.0)
    axes[1].plot(bins[0:-1], rate)

    axes[0].set_title("Spikes")
    axes[1].set_title("Firing rates")

    axes[0].set_xlim((0, 10000))
    axes[0].set_ylim((0, 3200))

    axes[0].set_ylabel("Neuron number")
    axes[1].set_ylabel("Mean firing rate [Hz]")

    axes[1].set_xlabel("Time [ms]")

    # Show plot
    plt.show()

