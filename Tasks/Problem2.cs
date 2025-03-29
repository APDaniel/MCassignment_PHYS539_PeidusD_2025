using MCassignment_PHYS539_PeidusD_2025.Containers;
using MCassignment_PHYS539_PeidusD_2025.Engine;

namespace MCassignment_PHYS539_PeidusD_2025.Tasks
{
    public static class Problem2
    {
        public static void ElectromBroadBeamTask(StatisticsToExport statistics)
        {
            double electronEnergy = 10.0;
            var numElectrons = 1000000;
            //Create 8x8x8 cm phantom
            double voxelSizeCm_XY = 0.5;
            double voxelSizeCm_Z = 0.2;
            var nx = 30; 
            var ny = 30; 
            var nz = 50; //resolution 0.2cm

            Console.WriteLine($"\n10MeV calculation for FS4x4 started...");
            var voxels = new WaterTankPhantom_ReducedResolution(nx, ny, nz, voxelSizeCm_XY, voxelSizeCm_Z);
            var electronsEnteredPerVoxel = numElectrons / ((nx* voxelSizeCm_XY) *(ny * voxelSizeCm_XY)); //Assuming FS of 4x4cm and 0.5cm grid XY

            MonteCarloEngine engine = new();
            engine.RunElectronBroadBeamSimulation(voxels,electronsEnteredPerVoxel,electronEnergy,statistics,numElectrons);
                
            string baseName = $"Electrons_DosePrimaryDepositedCAX_{(double)10}MeV.CSV";
            var yMaxSlice=voxels.FindPeakYSlice();
            //voxels.ExportFullKernel(baseName, numPhotons, _photonEnergy);
            Console.WriteLine($"Kernel generation complete for {(double)10}MeV");
   
            

        }
    }
}
