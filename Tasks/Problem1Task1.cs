using MCassignment_PHYS539_PeidusD_2025.Containers;
using MCassignment_PHYS539_PeidusD_2025.Engine;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MCassignment_PHYS539_PeidusD_2025.Tasks
{
    public static class Problem1Task1
    {
        public static void KernelGenerationTask()
        {
            double[] photonEnergy = [2.0, 6.0, 10.0];
            int numPhotons = 1000000;
            //Create 5x cm phantom
            double voxelSizeCm = 0.05;
            var nx = 150;
            var ny = 150;
            var nz = 150;
            foreach (var _photonEnergy in photonEnergy)
            {
                Console.WriteLine($"\n{_photonEnergy}MeV photon simulation started...");
                var voxels = new WaterTankPhantom(nx, ny, nz, voxelSizeCm);
                var engine = new MonteCarloEngine();
                engine.Run(voxels, numPhotons, _photonEnergy);
                string baseName = $"Kernel_{(double)_photonEnergy}MeV.CSV";
                voxels.ExportFullKernel(baseName, numPhotons, _photonEnergy);
                Console.WriteLine($"Kernel generation complete for {(double)_photonEnergy}MeV");
            }
            
        }
    }
}
