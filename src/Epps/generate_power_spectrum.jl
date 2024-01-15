using FFTW, Plots

function generate_power_spectrum(average_epps_mean)
  Fs = 1000  # Sampling frequency
  t = 0:1/Fs:1  # Time vector
  # Compute FFT
  Y = fft(average_epps_mean)
  frequencies = fftfreq(length(average_epps_mean), 1 / Fs)

  # Compute power spectrum
  power_spectrum = abs.(Y) .^ 2 / length(average_epps_mean)
  return (power_spectrum, frequencies)
end
