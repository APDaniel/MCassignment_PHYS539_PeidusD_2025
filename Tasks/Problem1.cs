using MCassignment_PHYS539_PeidusD_2025.Containers;
using MCassignment_PHYS539_PeidusD_2025.Engine;
using MCassignment_PHYS539_PeidusD_2025.Helpers;
using ScottPlot;
using ScottPlot.WinForms;
using System.Drawing;
namespace MCassignment_PHYS539_PeidusD_2025.Tasks
{
    public static class Problem1
    {
        public static void KernelGenerationTask(StatisticsToExport statistics)
        {
            double[] photonEnergy = [2.0, 6.0, 10.0];
            int numPhotons = 1000000;
            //Create 8x8x8 cm phantom
            double voxelSizeCm = 0.05;
            var nx = 250;
            var ny = 250;
            var nz = 250;
            foreach (var _photonEnergy in photonEnergy)
            {

                Console.WriteLine($"\n{_photonEnergy}MeV photon simulation started...");
                var maxKlNichProb = maxKleinNichinaHelper.SampleMaxKleinNichina(_photonEnergy);
                var voxels = new WaterTankPhantom_UniformResolution(nx, ny, nz, voxelSizeCm);
                Dictionary<double, double> lookupElScat = new();
                lookupElScat= ElectronPositronScatteringAngleLookupTable.LookupElectronScatAngle(_photonEnergy, voxelSizeCm);
                var engine = new MonteCarloEngine();
                engine.RunPhotonKernelCalculation(voxels, numPhotons, _photonEnergy, maxKlNichProb, statistics,voxelSizeCm, lookupElScat);

                string baseName = $"Kernel_{(double)_photonEnergy}MeV.CSV";

                //voxels.ExportFullKernel(baseName, numPhotons, _photonEnergy);
                Console.WriteLine($"Kernel generation complete for {(double)_photonEnergy}MeV");

                //Plot histogram for reporting. Hard-coded, non-flexible, specified iteratively for each task due to the time constraint.
                var plt = new Plot();
                var data = statistics.GetPhotoelectricScatteringAngles();
                //plt.Add.Histogram(ScottPlot.Statistics.Histogram.WithBinCount(10,data), ScottPlot.Color.FromHex("#9e0afa"),true);
                //plt.SavePng("TestGraph.png", 400, 300);

                //Metrics to check functionality
                double meanComptonAngle = statistics.AverageAngleCompton();
                double meanPairAngle_e = statistics.AverageAnglePair_e();
                double meanPairAngle_p = statistics.AverageAnglePair_p();
                double meanHardCollisionAngle = statistics.AverageHardCollisionAngle();

                // Export statistics
                var filePath = "Statistics.csv"; // Add .csv extension
                using (var writer = new StreamWriter(filePath))
                {
                    // Optional: Write header
                    writer.WriteLine("Statistics");

                    // Write each value in a new line
                    int iterator = 0;
                    foreach (var angle in StatisticsToExport.scatteringAnglesComptonPhotons)
                    {
                        iterator = iterator + 1;
                        writer.WriteLine(angle.ToString());
                        if (iterator > 200000) break;
                    }
                }

                Console.WriteLine("CSV export completed: " + filePath);
                statistics.ClearComptonEnergy();
                statistics.ClearCompton();
                statistics.ClearPair_e();
                statistics.ClearPair_p();
                statistics.ClearHardCollisionAngle();
                


            }

                
            

        }
    }
}
