using System;
using System.Numerics; // For complex numbers

// Class demonstrating various signal processing examples.
class SignalProcessingExamples
{
    // Entry point for the console application.
    static void Main()
    {
        // Test various signal processing functions by calling specific methods.
        TestConvolution();
        TestFourierTransform();
        TestLowPassFilter();
        TestSamplingAndReconstruction();
        TestWienerFilter();
    }

    // Tests convolution with a basic example.
    static void TestConvolution()
    {
        double[] signal = { 1, 2, 3, 4 }; // Define a signal array.
        double[] kernel = { 1, 0.5 }; // Define a kernel array.
        double[] result = Convolve(signal, kernel); // Perform convolution.

        Console.WriteLine("Convolution Result:");
        foreach (var value in result)
        {
            Console.Write(value + " "); // Print each value in the result.
        }
        Console.WriteLine(); // Move to a new line.
    }

    // Method to perform convolution.
    static double[] Convolve(double[] signal, double[] kernel)
    {
        int signalLength = signal.Length;
        int kernelLength = kernel.Length;
        double[] result = new double[signalLength + kernelLength - 1]; // Result array size.

        for (int i = 0; i < result.Length; i++)
        {
            result[i] = 0;
            for (int j = 0; j < kernelLength; j++)
            {
                if (i - j >= 0 && i - j < signalLength)
                {
                    result[i] += signal[i - j] * kernel[j]; // Convolution sum.
                }
            }
        }
        return result;
    }

    // Tests Fourier transform with a basic example.
    static void TestFourierTransform()
    {
        double[] signal = { 0, 1, 0, -1, 0, 1, 0, -1 }; // Define a simple sinusoidal signal.
        Complex[] dftResult = DFT(signal); // Compute the discrete Fourier transform.

        Console.WriteLine("DFT Result:");
        foreach (var value in dftResult)
        {
            Console.WriteLine(value); // Print each complex number in the result.
        }
    }

    // Method to perform a discrete Fourier transform.
    static Complex[] DFT(double[] signal)
    {
        int N = signal.Length;
        Complex[] result = new Complex[N];

        for (int k = 0; k < N; k++)
        {
            result[k] = Complex.Zero;
            for (int n = 0; n < N; n++)
            {
                double angle = -2 * Math.PI * k * n / N;
                result[k] += signal[n] * Complex.FromPolarCoordinates(1, angle); // Accumulate DFT sum.
            }
        }
        return result;
    }

    // Tests low-pass filtering with a basic example.
    static void TestLowPassFilter()
    {
        double[] signal = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }; // Define a linearly increasing signal.
        double[] kernel = CreateLowPassKernel(3, 0.33); // Create a low-pass filter kernel.
        double[] filteredSignal = Convolve(signal, kernel); // Apply convolution for filtering.

        Console.WriteLine("Filtered Signal:");
        foreach (var value in filteredSignal)
        {
            Console.Write(value + " "); // Print each value of the filtered signal.
        }
        Console.WriteLine(); // Move to a new line.
    }

    // Method to create a low-pass filter kernel.
    static double[] CreateLowPassKernel(int size, double cutoffFrequencyFactor)
    {
        double[] kernel = new double[size];
        double sum = 0;

        for (int i = 0; i < size; i++)
        {
            // Simple low-pass filter kernel using a sinc function approximation.
            kernel[i] = Math.Sin(2 * Math.PI * cutoffFrequencyFactor * (i - size / 2)) / (i - size / 2);
            if (i == size / 2) kernel[i] = 2 * Math.PI * cutoffFrequencyFactor; // Handle division by zero.
            sum += kernel[i];
        }

        // Normalize the kernel by dividing each element by the sum.
        for (int i = 0; i < size; i++)
        {
            kernel[i] /= sum;
        }

        return kernel;
    }

    // Tests sampling and reconstruction with a basic example.
    static void TestSamplingAndReconstruction()
    {
        int points = 32;
        double frequency = 1; // Frequency of the sine wave.
        double[] sampledSignal = new double[points];
        for (int i = 0; i < points; i++)
        {
            sampledSignal[i] = Math.Sin(2 * Math.PI * frequency * i / points); // Sample the sine wave.
        }

        Complex[] frequencyDomain = DFT(sampledSignal); // Perform DFT.
        double[] reconstructedSignal = IDFT(frequencyDomain); // Perform inverse DFT.

        Console.WriteLine("Original Signal:");
        foreach (var value in sampledSignal)
            Console.Write($"{value:F2} ");

        Console.WriteLine("\nReconstructed Signal:");
        foreach (var value in reconstructedSignal)
            Console.Write($"{value:F2} ");
        Console.WriteLine(); // Move to a new line.
    }

    // Method to perform an inverse discrete Fourier transform.
    static double[] IDFT(Complex[] frequencyDomain)
    {
        int N = frequencyDomain.Length;
        double[] result = new double[N];
        for (int n = 0; n < N; n++)
        {
            for (int k = 0; k < N; k++)
            {
                double angle = 2 * Math.PI * k * n / N;
                result[n] += (frequencyDomain[k] * Complex.FromPolarCoordinates(1, angle)).Real / N; // Accumulate IDFT sum.
            }
        }
        return result;
    }

    // Tests Wiener filtering with a basic example.
    static void TestWienerFilter()
    {
        int sampleCount = 64;
        double[] signal = new double[sampleCount];
        Random rand = new Random();
        for (int i = 0; i < sampleCount; i++)
        {
            signal[i] = Math.Sin(2 * Math.PI * i / 16) + 0.5 * (rand.NextDouble() - 0.5); // Signal with noise.
        }

        Complex[] signalDFT = DFT(signal); // Compute the DFT of the signal.
        double[] powerSpectrum = new double[sampleCount];
        for (int i = 0; i < sampleCount; i++)
        {
            powerSpectrum[i] = signalDFT[i].Magnitude; // Compute power spectrum.
        }

        double noisePower = 0.05; // Assume noise power is known.
        double[] filteredSignal = WienerFilter(signal, powerSpectrum, noisePower); // Apply Wiener filter.

        Console.WriteLine("Original Signal:");
        foreach (var value in signal)
            Console.Write($"{value:F2} ");

        Console.WriteLine("\nFiltered Signal:");
        foreach (var value in filteredSignal)
            Console.Write($"{value:F2} ");
        Console.WriteLine(); // Move to a new line.
    }

    // Method to apply Wiener filtering.
    static double[] WienerFilter(double[] inputSignal, double[] powerSpectrum, double noisePower)
    {
        int N = inputSignal.Length;
        Complex[] inputDFT = DFT(inputSignal);
        Complex[] outputDFT = new Complex[N];

        for (int i = 0; i < N; i++)
        {
            double signalPower = powerSpectrum[i];
            double filterFactor = signalPower / (signalPower + noisePower); // Compute filter factor.
            outputDFT[i] = inputDFT[i] * filterFactor; // Apply filter.
        }

        return IDFT(outputDFT); // Convert back to time domain.
    }
}