using MCassignment_PHYS539_PeidusD_2025.Containers;
using MCassignment_PHYS539_PeidusD_2025.Engine;
using MCassignment_PHYS539_PeidusD_2025.Helpers;

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
            var nx = 200;
            var ny = 200;
            var nz = 200;
            foreach (var _photonEnergy in photonEnergy)
            {

                Console.WriteLine($"\n{_photonEnergy}MeV photon simulation started...");
                var maxKlNichProb = maxKleinNichinaHelper.SampleMaxKleinNichina(_photonEnergy);
                var voxels = new WaterTankPhantom(nx, ny, nz, voxelSizeCm);
                var engine = new MonteCarloEngine();
                engine.RunKernelCalculation(voxels, numPhotons, _photonEnergy, maxKlNichProb, statistics);
                string baseName = $"Kernel_{(double)_photonEnergy}MeV.CSV";
                voxels.ExportFullKernel(baseName, numPhotons, _photonEnergy);
                Console.WriteLine($"Kernel generation complete for {(double)_photonEnergy}MeV");
                double meanComptonAngle = statistics.AverageAngleCompton();
                double meanPairAngle = statistics.AverageAnglePair();
                statistics.ClearCompton();
                statistics.ClearPair();
            }
            
        }
    }
}
