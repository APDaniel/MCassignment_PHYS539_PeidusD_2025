using MCassignment_PHYS539_PeidusD_2025.Containers;
using MCassignment_PHYS539_PeidusD_2025.Engine;
using MCassignment_PHYS539_PeidusD_2025.Helpers;
using ScottPlot;
using static Tensorflow.Binding;
using Tensorflow;
using System.Numerics;


namespace MCassignment_PHYS539_PeidusD_2025.Tasks
{
    public static class Problem1a
    {
        public static void InterpolationSuperpositionTask(StatisticsToExport statistics)
        {
            //Calculate kernel for 10MeV first. Did not have time to write a method to store kernel from the problem#1, so...
            //...I just recalculate it. The MC code is sufficiently efficient. Does not take long to calculate. ~2-3 min.
            //Future development.

            double[] photonEnergy = [10.0]; //Leftover. Do not have time to fix. 
            //Create phantom
            double voxelSizeCm = 0.05;
            var nx = 250;
            var ny = 250;
            var nz = 500;
            WaterTankPhantom_UniformResolution voxels = new(250,250,250, voxelSizeCm); //Leftover. Do not have time to fix. 
            LargeWaterTank largeWaterTank = new(nx, ny, nz, voxelSizeCm);
            //Calculate the amount of ptotons interacting in a large water tank.
            //10e7 photons. FS=10x10cm2

            Console.WriteLine("\nCalculating the amount of photons interacting per voxel...");
            ConvolutionSuperpositionEngine engineConvolution = new();

            ConvolutionSuperpositionEngine.CalculateNumberOfInteractingPhotonsPerVoxel(largeWaterTank,100000,statistics);
            Console.WriteLine("Number of interactions calculated.\nInitiaizing convolution-superposition...");
            var maxKlNichProb = maxKleinNichinaHelper.SampleMaxKleinNichina(10);
            var lookupElScat = ElectronPositronScatteringAngleLookupTable.LookupElectronScatAngle(10, voxelSizeCm);
            MonteCarloEngine engineMC = new();

            engineMC.RunPhotonKernelCalculation(voxels, 100, 10, maxKlNichProb, statistics, 0.05, lookupElScat);
            var test=largeWaterTank.FindPeakYSlice();
            largeWaterTank.ExportSliceY(test, "TEST1", 100000, 10);


            var kernel = voxels.GetDoseGrid();
            var photonMap = largeWaterTank.GetDoseGrid();
            //Perform convolution superposition
            ManualConvolver3D.Convolve(largeWaterTank, voxels);
            
            /*for (int ix = 0; ix < voxels.Nx; ix++)
            {
                for (int iy = 0; iy < voxels.Ny; iy++)
                {
                    for (int iz = 0; iz < voxels.Nz; iz++)
                    {
                        var test = result[ix, iy, iz];
                        photonMap[ix, iy, iz] = result[ix, iy, iz];
                    }
                }
            }*/
            //largeWaterTank.ExportFullKernel("FullKernel", 10000000, 10);
            var peakY=largeWaterTank.FindPeakYSlice();
            largeWaterTank.ExportSliceX(peakY, "CS_Xslice.CSV", 1, 10);
            largeWaterTank.ExportSliceY(peakY, "CS_Yslice.CSV", 1, 10);
            Console.WriteLine("Convolution-superposition complete.");
            statistics.AveragePhotonTravel();
        }

    }
}
