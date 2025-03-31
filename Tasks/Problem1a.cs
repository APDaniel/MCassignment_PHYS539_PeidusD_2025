using MCassignment_PHYS539_PeidusD_2025.Containers;
using MCassignment_PHYS539_PeidusD_2025.Engine;
using MCassignment_PHYS539_PeidusD_2025.Helpers;
using ScottPlot;
using static Tensorflow.Binding;
using Tensorflow;
using System.Numerics;


namespace MCassignment_PHYS539_PeidusD_2025.Tasks
{
    /************************************************************************************
    Class is a controller fot the Problem Ia. Convolution-Superposition calculations.
    I was limited with the deadline, so to simplify the workflow, 
    I perform MC for 10MeV 10e6 particles and then use it as a Kernel for each particle
    position calculated for the FS10x10, 10e7 aprticles. The MC elgorithm is very efficient,
    and does not take long to be calculated. Even on the 10+ years old laptop.
    ************************************************************************************/
    public static class Problem1a
    {
        public static void InterpolationSuperpositionTask(StatisticsToExport statistics)
        {


            double[] photonEnergy = [10.0]; //Leftover. Do not have time to fix. 

            //Create container for Convolution
            double voxelSizeCm = 0.05;
            var nx = 250;
            var ny = 250;
            var nz = 500;
            LargeWaterTank largeWaterTank = new(nx, ny, nz, voxelSizeCm);

            //Create container for Kernel MC calculations.
            WaterTankPhantom_UniformResolution voxels = new(250,250,250, voxelSizeCm); 

            Console.WriteLine("\nCalculating the amount of photons interacting per voxel...");

            //Calculate the amount of ptotons interacting in the large container.
            ConvolutionSuperpositionEngine engineConvolution = new();
            ConvolutionSuperpositionEngine.CalculateNumberOfInteractingPhotonsPerVoxel(largeWaterTank,100000,statistics);

            Console.WriteLine("Number of interactions calculated.\nInitiaizing convolution-superposition...");

            //Calculate maximim KN probability for the Rejection Method.
            //This is needed for determination of beta- start angle.
            var maxKlNichProb = maxKleinNichinaHelper.SampleMaxKleinNichina(1000);

            //Calculate a lookup table
            var lookupElScat = ElectronPositronScatteringAngleLookupTable.LookupElectronScatAngle(1000, voxelSizeCm);

            //Calculate primary Kernel. 10e6 particles, 10MeV.
            MonteCarloEngine engineMC = new();
            engineMC.RunPhotonKernelCalculation(voxels, 1000000, 10, maxKlNichProb, statistics, 0.05, lookupElScat);

            //Export peak slice for verification of the result.
            //Check that the kernel is non-zero.
            var test=largeWaterTank.FindPeakYSlice();
            largeWaterTank.ExportSliceY(test, "TEST1", 100000, 10);

            //Perform convolution superposition for 10e7 particles.
            var kernel = voxels.GetDoseGrid();
            var photonMap = largeWaterTank.GetDoseGrid();
            ManualConvolver3D.Convolve(largeWaterTank, voxels);

            //Check functionality manually.
            statistics.AveragePhotonTravel();

            //Expoort .CSV files for the peak slices for reporting: profiles and PDDs.
            var peakY=largeWaterTank.FindPeakYSlice();
            largeWaterTank.ExportSliceX(peakY, "CS_Xslice.CSV", 1, 10);
            largeWaterTank.ExportSliceY(peakY, "CS_Yslice.CSV", 1, 10);

            Console.WriteLine("Convolution-superposition complete.");
            
        }

    }
}
