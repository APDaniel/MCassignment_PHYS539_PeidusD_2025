using MCassignment_PHYS539_PeidusD_2025.Containers;
using MCassignment_PHYS539_PeidusD_2025.Particles;
using System.Numerics;

namespace MCassignment_PHYS539_PeidusD_2025.Engine
{
    public class ConvolutionSuperpositionEngine()
    {
        public static void CalculateNumberOfInteractingPhotonsPerVoxel(LargeWaterTank voxels, double numberOfPhotons, StatisticsToExport statistics)
        {
            //Convert FS to mm
            var voxelSize=voxels.VoxelSizeCm;
            var amountOfPhotonsEnteringPerVoxel = 0.0;
            var fieldSize_cm = 10;
            amountOfPhotonsEnteringPerVoxel = numberOfPhotons/(Math.Pow(fieldSize_cm,2) / voxelSize);
            amountOfPhotonsEnteringPerVoxel = Math.Round(amountOfPhotonsEnteringPerVoxel); //This is required for debugging. Non-relevant for the submission.
            //Method calculates amount of photons per voxel for primary interactions
            SamplePhotonTransport(voxels, numberOfPhotons, amountOfPhotonsEnteringPerVoxel,statistics);
            
        }

        public static void SamplePhotonTransport(LargeWaterTank voxels, double numPhotons, double photPerVxl, StatisticsToExport statistics) 
        {
            var mu = 0.022; //Linear attenuation coefficient.
            

            for (int ix = 0; ix < voxels.Nx; ix++)
            {
                for (int iy=0;iy<voxels.Ny; iy++)
                {
                    for (int iz = 0; iz < voxels.Nz; iz++)
                    {
                            Random random = new();
                            var R1 = random.NextDouble();
                            var oneOverMu = 1 / mu;
                            var log = Math.Log(R1);
                            var distanceToInteraction = -oneOverMu * log;
                            Vector3 direction = new(0f, 0f, (float)distanceToInteraction);
                            Vector3 position = new(0, 0, 0);
                            position = position + direction;
                            statistics.SavePhotonTravelDistance(position.Z);
                            var doseGrid = voxels.GetDoseGrid();
                            var positionZinVoxels = position.Z / voxels.VoxelSizeCm;
                            if (positionZinVoxels < voxels.Nz)
                                voxels.DepositEnergy(position.X, (int)position.Y, (int)positionZinVoxels,1);
                                
                    }
                    
                }
            } 
        }
    }
}
